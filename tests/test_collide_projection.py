"""Per-leg collision fields, coupled geometry and Maya evaluation."""

import hashlib
from decimal import Decimal, localcontext
import itertools
import json
from pathlib import Path
import math
import struct
import unittest

import collide_reference as reference

try:
    import maya.api.OpenMaya as om
    import maya.cmds as cmds
except ImportError:
    om = cmds = None

from collide_reference import (
    Cylinder,
    add,
    finite,
    length,
    mul,
    profile_stations,
    sub,
)

if cmds is not None:
    from helpers import bind_rest_inputs, coordinates as geometry_coordinates, projection_rig


def coordinates(node, attribute, kind="mesh"):
    if kind == "surface":
        snapshot = surface_snapshot(node + "." + attribute)
        return [list(p) + [w] for p, w in zip(snapshot["cvs"], snapshot["weights"])]
    return geometry_coordinates(node, attribute, kind)


def points(rig):
    return [tuple(point[:3]) for point in coordinates(rig["node"], "outputGeometry[0]")]


def rest_points(rig):
    shape = cmds.listRelatives(rig["rest"], shapes=True, noIntermediate=True, fullPath=True)[0]
    return [tuple(point[:3]) for point in coordinates(shape, "outMesh")]


def point_hash(values):
    return hashlib.sha256(b"".join(struct.pack("<3d", *point) for point in values)).hexdigest()


def assert_displacements_close(test, original, actual, expected):
    test.assertEqual(len(actual), len(expected))
    for index, (p, a, b) in enumerate(zip(original, actual, expected)):
        actual_u, expected_u = sub(a, p), sub(b, p)
        tolerance = 1e-6 + 1e-9 * max(length(actual_u), length(expected_u))
        test.assertLessEqual(length(sub(actual_u, expected_u)), tolerance, "point {}".format(index))


def assert_topology_observations(test, observations, comparator):
    saved = {}
    for key, original, actual, expected in observations:
        with test.subTest(topology_check="movement", key=key):
            test.assertGreater(max(length(sub(a, p)) for a, p in zip(actual, original)), 1e-5)
        if key in saved:
            with test.subTest(topology_check="repeatability", key=key):
                test.assertEqual(actual, saved[key])
        else:
            saved[key] = actual
    with test.subTest(topology_check="closed_open"):
        test.assertGreater(max(length(sub(a, b)) for a, b in zip(saved[0, False], saved[0, True])), 1e-5)
    for key, original, actual, expected in observations:
        with test.subTest(topology_check="reference", key=key):
            comparator(test, original, actual, expected)


def node_cylinders(node):
    scale = cmds.getAttr(node + ".ringScale")[0]
    names = [station + "Radius" + axis for station in ("thigh", "knee", "calf", "ankle") for axis in ("X", "Z")]
    profile = {name: cmds.getAttr(node + "." + name) for name in names + ["thighPosition", "calfPosition"]}
    current, rest = {}, {}
    for leg_id, side in enumerate(("left", "right")):
        matrices = [
            [cmds.getAttr(node + "." + prefix + joint + "Matrix") for joint in ("Hip", "Knee", "Heel")]
            for prefix in (side, "rest" + side.capitalize())
        ]
        axis_index = cmds.getAttr(node + "." + side + "RingAxis")
        frames = []
        for group in matrices:
            if any(not finite(matrix[12:15]) for matrix in group):
                frames.append({})
                continue
            positions = [tuple(matrix[12:15]) for matrix in group]
            thigh, calf = [length(sub(b, a)) for a, b in zip(positions, positions[1:])]
            segments = {}
            if thigh + calf >= 1e-6:
                knee = thigh / (thigh + calf)
                knots = profile_stations(thigh, calf, **profile)
                for segment in range(1 + cmds.getAttr(node + ".skirtType")):
                    joint, target = group[segment : segment + 2]
                    if any(
                        not finite(matrix)
                        or not math.isfinite(om.MMatrix(matrix).det4x4())
                        or abs(om.MMatrix(matrix).det4x4()) < 1e-8
                        for matrix in (joint, target)
                    ):
                        continue
                    start, finish = (0, knee) if segment == 0 else (knee, 1)
                    extent = (thigh if segment == 0 else calf) * scale[1]
                    stations = [
                        ((t - start) * (thigh + calf) * scale[1], a * scale[0], b * scale[2])
                        for t, a, b in knots
                        if start <= t <= finish
                    ]
                    if (
                        not math.isfinite(extent)
                        or extent < 1e-6
                        or len(stations) < 2
                        or any(not finite(station) or min(station[1:]) <= 0 for station in stations)
                    ):
                        continue
                    rows = [tuple(joint[i : i + 3]) for i in (0, 4, 8)]
                    axis = om.MVector(*sub(positions[segment + 1], positions[segment])).normal()
                    y = om.MVector(*rows[axis_index % 3]).normal() * (-1 if axis_index >= 3 else 1)
                    x = om.MVector(*rows[(axis_index + 1) % 3]).normal()
                    rotation = y.rotateTo(axis)
                    y, x = y.rotateBy(rotation), x.rotateBy(rotation)
                    x = (x - y * (x * y)).normal()
                    z = (y ^ x).normal()
                    if not finite(tuple(x)) or not finite(tuple(z)) or abs(x * (y ^ z)) < 1e-8:
                        continue
                    segments[segment] = Cylinder(
                        positions[segment], tuple(axis), tuple(x), tuple(z), extent, stations, rows
                    )
            frames.append(segments)
        for segment, cylinder in frames[0].items():
            current[(leg_id, segment)] = cylinder
        for segment, cylinder in frames[1].items():
            rest[(leg_id, segment)] = cylinder
    return current, rest


def surface_snapshot(plug):
    transform = cmds.createNode("transform")
    try:
        shape = cmds.createNode("nurbsSurface", parent=transform)
        cmds.connectAttr(plug, shape + ".create", force=True)
        fn = om.MFnNurbsSurface(om.MSelectionList().add(shape).getDagPath(0))
        cvs = fn.cvPositions(om.MSpace.kObject)
        return dict(
            nu=fn.numCVsInU,
            nv=fn.numCVsInV,
            du=fn.degreeInU,
            dv=fn.degreeInV,
            ku=list(fn.knotsInU()),
            kv=list(fn.knotsInV()),
            formU=fn.formInU,
            formV=fn.formInV,
            cvs=[(p.x, p.y, p.z) for p in cvs],
            weights=[p.w for p in cvs],
        )
    finally:
        cmds.delete(transform)


def surface_node_reference(node, iterations=None, coupling=True, members=None):
    source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
    surface = surface_snapshot(source or node + ".input[0].inputGeometry")
    surface["closedU"] = cmds.getAttr(node + ".closedU")
    rest_source = cmds.connectionInfo(node + ".restGeometry", sourceFromDestination=True)
    rest_surface = surface_snapshot(rest_source or node + ".restGeometry")
    original, rps = surface["cvs"], rest_surface["cvs"]
    current, rest = node_cylinders(node)
    bell = cmds.getAttr(node + ".restBellMatrix")
    waist, axis = tuple(bell[12:15]), tuple(bell[4:7])
    falloff, envelope = [cmds.getAttr(node + "." + attr) for attr in ("falloff", "envelope")]
    if not finite(waist) or not finite(axis) or length(axis) < 1e-8 or not math.isfinite(falloff):
        return original, rps, original[:], [[] for _ in original], {}
    q = [cmds.getAttr(node + ".weightList[0].weights[{}]".format(i)) for i in range(len(original))]
    u = reference.surface_displacements(
        surface,
        original,
        rps,
        reference.make_segs(current, rest, bool(cmds.getAttr(node + ".skirtType"))),
        [rest[k] for k in sorted(rest)] if coupling else [],
        waist,
        axis,
        falloff,
        q=q,
        members=members,
        iters=iterations,
    )
    uu = surface["nu"] - surface["du"] if surface["formU"] == 3 else surface["nu"]
    vv = surface["nv"] - surface["dv"] if surface["formV"] == 3 else surface["nv"]
    expected = original[:]
    for i, p in enumerate(original):
        canonical = (i // surface["nv"] % uu) * surface["nv"] + i % surface["nv"] % vv
        if math.isfinite(envelope) and envelope != 0:
            candidate = add(p, mul(u[canonical], envelope))
            if finite(candidate):
                expected[i] = candidate
    legs = reference.make_segs(current, rest, bool(cmds.getAttr(node + ".skirtType")))
    fields = [reference.all_cons(legs, p, r, waist, axis, falloff) for p, r in zip(original, rps)]
    return original, rps, expected, fields, dict(u=u, indices=sorted(u))


def node_reference(node, kind="mesh", iterations=None, reverse=False, coupling=True, members=None):
    if kind == "surface":
        return surface_node_reference(node, iterations, coupling, members)
    iterations = 30 if iterations is None else iterations
    original = [tuple(p[:3]) for p in coordinates(node, "input[0].inputGeometry", kind)]
    rest_points = [tuple(p[:3]) for p in coordinates(node, "restGeometry", kind)]
    current, rest = node_cylinders(node)
    legs = reference.make_segs(current, rest, bool(cmds.getAttr(node + ".skirtType")))
    geometry = om.MSelectionList().add(node + ".input[0].inputGeometry").getPlug(0).asMObject()
    adjacency = {}
    vertex = om.MItMeshVertex(geometry)
    while not vertex.isDone():
        adjacency[vertex.index()] = sorted(vertex.getConnectedVertices())
        vertex.next()
    indices = sorted(adjacency)
    if reverse:
        indices.reverse()
    lookup = {index: i for i, index in enumerate(indices)}
    adj = [[lookup[j] for j in adjacency[index] if j in lookup] for index in indices]
    bell = cmds.getAttr(node + ".restBellMatrix")
    waist, waist_axis = tuple(bell[12:15]), tuple(bell[4:7])
    falloff, envelope = [cmds.getAttr(node + "." + attr) for attr in ("falloff", "envelope")]
    pts, rps = [original[i] for i in indices], [rest_points[i] for i in indices]
    if not finite(waist) or not finite(waist_axis) or length(waist_axis) < 1e-8 or not math.isfinite(falloff):
        return original, rest_points, original[:], [[] for _ in original], {}
    q = reference.point_weights(
        pts,
        rps,
        [cmds.getAttr(node + ".weightList[0].weights[{}]".format(i)) for i in indices],
        waist,
        waist_axis,
    )
    cons = [reference.all_cons(legs, p, rp, waist, waist_axis, falloff) if w else [] for p, rp, w in zip(pts, rps, q)]
    u = reference.couple(
        pts, rps, cons, adj, falloff if coupling else 0.0, [rest[key] for key in sorted(rest)], q, iterations
    )
    expected = original[:]
    all_cons = [[] for _ in original]
    for i, index in enumerate(indices):
        if q[i] and math.isfinite(envelope) and envelope != 0:
            candidate = add(original[index], mul(u[i], envelope))
            if finite(candidate):
                expected[index] = candidate
        all_cons[index] = cons[i]
    residual = max(
        [0.0] + [max(0.0, d - reference.dot(n, value)) for planes, value in zip(cons, u) for n, d, _ in planes]
    )
    return original, rest_points, expected, all_cons, dict(residual=residual, indices=indices, u=u)


@unittest.skipIf(cmds is None, "Maya is required")
class CollideProjectionTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")

    def assert_points_close(self, actual, expected, tolerance=1e-9):
        self.assertEqual(len(actual), len(expected))
        self.assertTrue(all(math.isfinite(v) for point in actual for v in point))
        for index, (a, b) in enumerate(zip(actual, expected)):
            self.assertLessEqual(length(sub(a, b)), tolerance, "vertex {}".format(index))

    def check_projection(self, rig):
        original, _, expected, fields, solves = node_reference(rig["node"])
        result = points(rig)
        assert_displacements_close(self, original, result, expected)
        if cmds.getAttr(rig["node"] + ".falloff") == 0:
            for p, result_point, planes in zip(original, result, fields):
                if all(d <= 0 for _, d, _ in planes):
                    self.assertEqual(result_point, p)
        print("collide reference max_error={}".format(max(length(sub(a, b)) for a, b in zip(result, expected))))

    def test_exterior_rest_identity_is_exact(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type)
            original, _, expected, fields, _ = node_reference(rig["node"])
            self.assertTrue(all(d <= 0 for entry in fields for _, d, _ in entry))
            self.assertEqual(points(rig), original)
            self.assertEqual(expected, original)

    def test_interior_rest_uses_reference_displacement(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type, angle=90, rest_angle=90)
            original, _, expected, fields, _ = node_reference(rig["node"])
            interior = [i for i, entry in enumerate(fields) if any(d > 0 for _, d, _ in entry)]
            self.assertTrue(interior)
            self.assertTrue(any(length(sub(expected[i], original[i])) > 1e-6 for i in interior))
            self.check_projection(rig)

    def test_horizontal_short_and_long_match_constraints(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type, angle=90)
            self.check_projection(rig)

    def test_elliptical_station_profile_matches_reference(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type, angle=65, scale=(2, 1, 1))
            for attr, value in (
                ("thighRadiusX", 1.4),
                ("kneeRadiusZ", 0.8),
                ("calfRadiusX", 1.8),
                ("ankleRadiusZ", 0.6),
            ):
                cmds.setAttr(rig["node"] + "." + attr, value)
            self.check_projection(rig)

    def test_unequal_leg_lengths(self):
        lengths = ((4, 6), (6, 3))
        for skirt_type in (0, 1):
            with self.subTest(skirt_type=skirt_type):
                cmds.file(new=True, force=True)
                rig = projection_rig(skirt_type, angle=90, lengths=lengths)
                self.check_projection(rig)

    def test_sweep_and_reverse(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type)
            cmds.setAttr(rig["node"] + ".falloff", 0.2)
            original = rest_points(rig)
            shape = cmds.listRelatives(rig["mesh"], shapes=True, noIntermediate=True, fullPath=True)[0]
            mesh = om.MFnMesh(om.MSelectionList().add(shape).getDagPath(0))
            edges = [mesh.getEdgeVertices(index) for index in range(mesh.numEdges)]
            hashes = {}
            previous = None
            maximum = 0.0
            neighbor_maximum = 0.0
            for angle in range(91):
                cmds.setAttr(rig["joints"]["leftHip"] + ".rotateX", -angle)
                current = points(rig)
                displacement = [sub(q, p) for p, q in zip(original, current)]
                neighbor_maximum = max(
                    neighbor_maximum, max(length(sub(displacement[a], displacement[b])) for a, b in edges)
                )
                if previous is not None:
                    maximum = max(maximum, max(length(sub(a, b)) for a, b in zip(current, previous)))
                hashes[angle] = point_hash(current)
                previous = current
            print(
                "collide sweep: skirt_type={}, max_displacement_per_degree={}, max_neighbor_displacement_difference={}".format(
                    skirt_type, maximum, neighbor_maximum
                )
            )
            with self.subTest(skirt_type=skirt_type, check="round_trip"):
                for angle in reversed(range(91)):
                    cmds.setAttr(rig["joints"]["leftHip"] + ".rotateX", -angle)
                    self.assertEqual(point_hash(points(rig)), hashes[angle], "angle {}".format(angle))

            previous = None
            local_maximum = 0.0
            for step in range(61):
                angle = 60 + step / 10.0
                cmds.setAttr(rig["joints"]["leftHip"] + ".rotateX", -angle)
                current = points(rig)
                if previous is not None:
                    change = max(length(sub(a, b)) for a, b in zip(current, previous))
                    local_maximum = max(local_maximum, change)
                previous = current
            print(
                "collide local sweep: skirt_type={}, max_displacement_per_tenth_degree={}".format(
                    skirt_type, local_maximum
                )
            )

    def test_bilateral_symmetry_and_constraint_order(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type, angle=60, both=True)
            original, result = rest_points(rig), points(rig)
            inputs, _, _, fields, _ = node_reference(rig["node"])
            for index, p in enumerate(original):
                reflected = (-p[0], p[1], p[2])
                other = min(range(len(original)), key=lambda j: length(sub(original[j], reflected)))
                # Mesh points are stored as floats, so mirrored coordinates match to about 1e-5.
                self.assertLessEqual(length(sub(original[other], reflected)), 1e-5)
                symmetry = length(sub(result[other], (-result[index][0], result[index][1], result[index][2])))
                self.assertTrue(math.isfinite(symmetry))
            sixty = node_reference(rig["node"], iterations=60)
            reversed_indices = node_reference(rig["node"], reverse=True)
            print(
                "coupling symmetry={}, K30-K60={}, index_reverse={}, residual={}".format(
                    max(
                        length(
                            sub(
                                result[
                                    min(
                                        range(len(original)),
                                        key=lambda j: length(sub(original[j], (-p[0], p[1], p[2]))),
                                    )
                                ],
                                (-result[i][0], result[i][1], result[i][2]),
                            )
                        )
                        for i, p in enumerate(original)
                    ),
                    max(length(sub(a, b)) for a, b in zip(result, sixty[2])),
                    max(length(sub(a, b)) for a, b in zip(result, reversed_indices[2])),
                    sixty[4]["residual"],
                )
            )
            node = rig["node"]
            values = {
                prefix + side + joint + "Matrix": cmds.getAttr(node + "." + prefix + side + joint + "Matrix")
                for prefix, sides in (("", ("left", "right")), ("rest", ("Left", "Right")))
                for side in sides
                for joint in ("Hip", "Knee", "Heel")
            }
            for prefix, sides in (("", ("left", "right")), ("rest", ("Left", "Right"))):
                for side, other in (sides, sides[::-1]):
                    for joint in ("Hip", "Knee", "Heel"):
                        name = prefix + side + joint + "Matrix"
                        plug = node + "." + name
                        for source in cmds.listConnections(plug, source=True, destination=False, plugs=True) or []:
                            cmds.disconnectAttr(source, plug)
                        cmds.setAttr(plug, *values[prefix + other + joint + "Matrix"], type="matrix")
            swapped = points(rig)
            for index, (before, after) in enumerate(zip(result, swapped)):
                self.assertLessEqual(length(sub(before, after)), 1e-6)
            self.check_projection(rig)

    def test_evaluation_modes_and_time_order(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type)
            for time, value in ((0, 0), (12, -90)):
                cmds.setKeyframe(
                    rig["joints"]["leftHip"] + ".rotateX",
                    time=time,
                    value=value,
                    inTangentType="linear",
                    outTangentType="linear",
                )
            references = {}
            for mode in ("off", "serial", "parallel"):
                cmds.evaluationManager(mode=mode)
                for time in (0, 4, 12, 7, 0):
                    with self.subTest(skirt_type=skirt_type, mode=mode, time=time):
                        cmds.currentTime(time)
                        result = points(rig)
                        self.assertTrue(all(math.isfinite(v) for p in result for v in p))
                        digest = point_hash(result)
                        if time in references:
                            self.assertEqual(digest, references[time])
                        else:
                            references[time] = digest
                if mode != "off":
                    self.assertFalse(cmds.evaluationManager(query=True, empty=True))
                    self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))

    def test_missing_and_mismatched_rest_geometry_recover(self):
        warnings = []

        def capture(message, kind, _data):
            if kind == om.MCommandMessage.kWarning and "restGeometry" in message:
                warnings.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)
        # A disconnected generic plug keeps its last data, so the missing case starts unbound.
        rig = projection_rig(angle=90, bind=False)
        original = rest_points(rig)
        plug = rig["node"] + ".restGeometry"
        for angle in (80, 90, 80):
            cmds.setAttr(rig["joints"]["leftHip"] + ".rotateX", -angle)
            self.assert_points_close(points(rig), original)
        other = cmds.polyCube(constructionHistory=False)[0]
        shape = cmds.listRelatives(other, shapes=True)[0]
        cmds.connectAttr(shape + ".outMesh", plug)
        self.assert_points_close(points(rig), original)
        self.assertEqual(len(warnings), 1)
        bind_rest_inputs(rig["node"], rig["rest"])
        cmds.setAttr(rig["joints"]["leftHip"] + ".rotateX", -90)
        active = points(rig)
        self.assertGreater(max(length(sub(p, q)) for p, q in zip(original, active)), 1e-6)

    def test_half_envelope_scales_displacement(self):
        for skirt_type in (0, 1):
            cmds.file(new=True, force=True)
            rig = projection_rig(skirt_type, angle=90)
            original, full = rest_points(rig), points(rig)
            cmds.setAttr(rig["node"] + ".envelope", 0.5)
            self.assert_points_close(
                points(rig), [add(p, mul(sub(q, p), 0.5)) for p, q in zip(original, full)], tolerance=1e-6
            )

    def test_nurbs_rest_geometry(self):
        surface = cmds.nurbsPlane(axis=(0, 1, 0), width=5.2, patchesU=4, patchesV=4, constructionHistory=False)[0]
        cmds.move(0, 10.5, 0, surface + ".cv[*][*]", relative=True, objectSpace=True)
        rig = projection_rig(angle=90, mesh=surface)
        node = rig["node"]
        rest_shape = cmds.listRelatives(rig["rest"], shapes=True, noIntermediate=True, fullPath=True)[0]
        before = coordinates(rest_shape, "local", "surface")
        after = coordinates(node, "outputGeometry[0]", "surface")
        self.assert_points_close([tuple(p[:3]) for p in after], node_reference(node, "surface")[2], 1e-6)
        self.assertTrue(all(math.isfinite(v) for p in after for v in p))
        self.assertGreater(max(length(sub(p[:3], q[:3])) for p, q in zip(before, after)), 1e-6)
        cmds.setAttr(node + ".envelope", 0)
        self.assert_points_close(coordinates(node, "outputGeometry[0]", "surface"), before)

    def test_periodic_nurbs_uses_geometry_indices(self):
        surface = cmds.cylinder(
            axis=(0, 1, 0),
            radius=2.6,
            heightRatio=4.2,
            sections=8,
            spans=4,
            degree=3,
            constructionHistory=False,
        )[0]
        cmds.move(0, 5.5, 0, surface + ".cv[*][*]", relative=True, objectSpace=True)
        shape = cmds.listRelatives(surface, shapes=True, noIntermediate=True, fullPath=True)[0]
        fn = om.MFnNurbsSurface(om.MSelectionList().add(shape).getDagPath(0))
        self.assertIn(om.MFnNurbsSurface.kPeriodic, (fn.formInU, fn.formInV))
        for skirt_type in (0, 1):
            geometry = cmds.duplicate(surface, returnRootsOnly=True)[0]
            rig = projection_rig(skirt_type, angle=90, mesh=geometry)
            original, _, expected, _, _ = node_reference(rig["node"], "surface")
            actual = [tuple(p[:3]) for p in coordinates(rig["node"], "outputGeometry[0]", "surface")]
            assert_displacements_close(self, original, actual, expected)
            self.assertTrue(any(length(sub(a, b)) > 1e-6 for a, b in zip(actual, original)))
            cmds.setAttr(rig["node"] + ".envelope", 0)
            self.assertEqual([tuple(p[:3]) for p in coordinates(rig["node"], "outputGeometry[0]", "surface")], original)


class CollideReferenceTests(unittest.TestCase):
    def cylinder(self, **overrides):
        values = dict(
            origin=(0, 0, 0),
            axis=(0, 0, 1),
            x_axis=(1, 0, 0),
            z_axis=(0, 1, 0),
            extent=4,
            stations=[(0, 1, 1), (4, 1, 1)],
        )
        values.update(overrides)
        return Cylinder(**values)

    def test_sample_single_plane_and_unit_row(self):
        for row in ({0: 1.0}, {0: 0.2, 1: 0.8}):
            for strength in (0.0, 1e-8, 0.2, 1.0):
                n = len(row)
                u = reference.sample_couple(
                    [reference.ZERO] * n,
                    [reference.ZERO] * n,
                    [row],
                    [[i] for i in range(n)],
                    [[] for _ in row],
                    [1.0] * n,
                    lambda p, r: [((1, 0, 0), 2.0 - p[0], strength)],
                    iters=200,
                )
                y = sum(b * u[i][0] for i, b in row.items())
                self.assertAlmostEqual(y, 5 * strength / (1 + 4 * strength) * 2, delta=1e-8)
                if n == 1:
                    self.assertLess(
                        length(sub(u[0], reference.correction([((1, 0, 0), 2, strength)], reference.ZERO))),
                        1e-10,
                    )

    def test_sample_local_minimum_independent_active_enumeration(self):
        terms = [
            (reference.unit(n), t, w, b)
            for n, t, w, b in (
                ((1, 0.2, -0.3), 1, 8, 0.6),
                ((-1, 1, 0), 0.4, 200, 0.4),
                ((0.2, -1, 1), 0.9, 200, 0.8),
                ((0, 0.4, -1), -0.2, 8, 0.7),
                ((1, -0.2, 1), 0.3, 40, 1),
            )
        ]
        neighbours = [(0.7, (0.4, -0.2, 0.1))]
        mass = 2

        def gaussian(matrix, rhs):
            a = [list(row) + [value] for row, value in zip(matrix, rhs)]
            for k in range(3):
                pivot = max(range(k, 3), key=lambda i: abs(a[i][k]))
                a[k], a[pivot] = a[pivot], a[k]
                scale = a[k][k]
                a[k] = [v / scale for v in a[k]]
                for i in range(3):
                    if i != k:
                        scale = a[i][k]
                        a[i] = [x - scale * y for x, y in zip(a[i], a[k])]
            return tuple(row[3] for row in a)

        candidates = []
        for active in itertools.product((False, True), repeat=len(terms)):
            a = [[2.7 if r == c else 0.0 for c in range(3)] for r in range(3)]
            rhs = [0.7 * v for v in neighbours[0][1]]
            for enabled, (n, t, w, b) in zip(active, terms):
                if enabled:
                    for r in range(3):
                        rhs[r] += w * b * t * n[r]
                        for c in range(3):
                            a[r][c] += w * b * b * n[r] * n[c]
            z = gaussian(a, rhs)
            hinges = [t - b * reference.dot(n, z) for n, t, w, b in terms]
            if all(h >= -1e-10 if enabled else h <= 1e-10 for enabled, h in zip(active, hinges)):
                candidates.append(z)
        self.assertTrue(candidates)
        trace = []
        actual = reference.solve_group(mass, neighbours, terms, trace=trace)
        self.assertLess(length(sub(actual, candidates[0])), 1e-8)
        self.assertTrue(all(b[1] <= a[1] for a, b in zip(trace, trace[1:])))
        gradient = add(mul(actual, 2.7), mul(neighbours[0][1], -0.7))
        for n, t, w, b in terms:
            gradient = sub(gradient, mul(n, w * b * max(0, t - b * reference.dot(n, actual))))
        self.assertLess(length(gradient), 1e-7)

    def test_sample_two_cv_active_switch_is_continuous(self):
        def run(b):
            v = [reference.ZERO, reference.ZERO]
            for _ in range(2):
                for i, target in enumerate((1.0, b)):
                    v[i] = reference.solve_group(1, [(1, v[1 - i])], [((1, 0, 0), target, 4, 1)], v[i])
            return v

        a, b = run(1 / 3 - 1e-8), run(1 / 3 + 1e-8)
        self.assertLess(max(length(sub(x, y)) for x, y in zip(a, b)), 1e-7)

    def sample_surface(self, **overrides):
        surface = dict(nu=4, nv=3, du=2, dv=1, ku=[0, 0, 0.3, 1, 1], kv=[0, 0.6, 1], formU=1, formV=1)
        surface.update(overrides)
        return surface

    def test_sample_basis_support_endpoints_rational_and_periodic(self):
        surface = self.sample_surface()
        parameters, rows = reference.surface_basis(surface)
        self.assertEqual(len(rows), 36)
        self.assertTrue(all(abs(sum(row.values()) - 1) < 1e-12 for row in rows))
        self.assertTrue(all(i < 3 for row in rows[-12:-6] for i in row))
        self.assertTrue(all(i >= 9 for row in rows[-6:] for i in row))
        self.assertEqual(reference.spline_basis([0, 0, 0.5, 0.5, 1, 1], 2, 5, 0.5), {2: 1.0})
        surface["weights"] = [1 + i / 10 for i in range(12)]
        _, rational = reference.surface_basis(surface)
        for row, rat in zip(rows, rational):
            total = sum(b * surface["weights"][i] for i, b in row.items())
            for i, b in row.items():
                self.assertAlmostEqual(rat[i], b * surface["weights"][i] / total)
        periodic = self.sample_surface(nu=5, du=2, ku=[-1, 0, 1, 2, 3, 4], formU=3)
        sites, aggregated = reference.surface_basis(periodic)
        self.assertEqual(len(sites), 36)
        for (u, v), row in zip(sites, aggregated):
            raw = {}
            for a, x in reference.spline_basis(periodic["ku"], 2, 5, u).items():
                for b, y in reference.spline_basis(periodic["kv"], 1, 3, v).items():
                    i = (a % 3) * 3 + b
                    raw[i] = raw.get(i, 0) + x * y
            self.assertEqual(row.keys(), raw.keys())
            for i in row:
                self.assertAlmostEqual(row[i], raw[i])
        for weight in (0, -1, float("nan")):
            surface["weights"][0] = weight
            with self.assertRaises(ValueError):
                reference.surface_basis(surface)

    def test_sample_basis_gauss3_nonuniform_span_weights_and_open_endpoints(self):
        surface = self.sample_surface(kv=[0, 0.2, 1.7])
        sites, rows, weights = reference.surface_basis(surface, with_quadrature=True)
        self.assertEqual(len(sites), len(rows))
        self.assertEqual(len(sites), len(weights))
        self.assertEqual(len(sites), 36)
        self.assertEqual([u for u, _ in sites[-12:-6]], [0] * 6)
        self.assertEqual([u for u, _ in sites[-6:]], [1] * 6)
        for start, end, offset in ((0, 0.2, 0), (0.2, 1.7, 3)):
            values = [sites[offset + k][1] for k in range(3)]
            quadrature = weights[offset:offset + 3]
            self.assertEqual(quadrature, [10 / 9, 16 / 9, 10 / 9])
            self.assertAlmostEqual(sum(quadrature), 4)
            self.assertAlmostEqual(sum(w * v for w, v in zip(quadrature, values)), 2 * (start + end))
            self.assertAlmostEqual(
                sum(w * v * v for w, v in zip(quadrature, values)),
                4 * (start * start + start * end + end * end) / 3,
            )
        self.assertEqual(weights[-12:-6], weights[:6])
        self.assertEqual(weights[-6:], weights[:6])

    def test_sample_fixed_support_and_regeneration_schedule(self):
        row = {0: 0.5, 1: 0.5}
        points = [(2, 0, 0), (0, 0, 0)]
        for sweeps, regenerations, expected in (
            (10, 1, [0, 5]),
            (11, 1, [0, 5]),
            (10, 0, [0]),
            (40, 3, [0, 10, 20, 30]),
            (30, 3, [0, 7, 15, 22]),
        ):
            events = []
            u = reference.sample_couple(
                points,
                points,
                [row],
                [[0], [1]],
                [[], []],
                [0, 1],
                lambda p, r: [((1, 0, 0), 3 - p[0], 1)],
                iters=sweeps,
                regenerations=regenerations,
                events=events,
            )
            self.assertEqual(u[0], reference.ZERO)
            self.assertAlmostEqual(u[1][0], 10 / 3)
            self.assertEqual([event[0] for event in events], expected)
            self.assertEqual(events[0][2], (1, 0, 0))
            for _, _, p in events[1:]:
                self.assertAlmostEqual((3 - p[0]) + (p[0] - 1), 2)
                self.assertAlmostEqual(p[0], 1 + u[1][0] / 2)

    def test_sample_closed_seam_paint_membership_and_open_ends(self):
        for form_u, form_v, closed in itertools.product((1, 2, 3), (1, 2, 3), (False, True)):
            with self.subTest(formU=form_u, formV=form_v, closedU=closed):
                surface = self.sample_surface(formU=form_u, formV=form_v, closedU=closed)
                groups, adjacency = reference.surface_groups(surface)
                mapping = {i: g for g, group in enumerate(groups) for i in group}
                if form_u == 3:
                    self.assertTrue(all(len(group) == 1 for group in groups))
                    self.assertEqual(set(mapping), set(adjacency))
                else:
                    vv = surface["nv"] - surface["dv"] if form_v == 3 else surface["nv"]
                    for v in range(vv):
                        self.assertEqual(mapping[v] == mapping[9 + v], closed)
                self.assertNotEqual(mapping[0], mapping[1])
                if form_v != 3:
                    self.assertNotEqual(mapping[0], mapping[2])
        groups, _ = reference.surface_groups(self.sample_surface(formU=2, formV=2))
        self.assertTrue(all(len(group) == 1 for group in groups))
        points = [(0, 0, 0)] * 2

        def run(groups, q):
            return reference.sample_couple(
                points,
                points,
                [{0: 1}, {1: 1}],
                groups,
                [[] for _ in groups],
                q,
                lambda p, r: [((1, 0, 0), 1 - p[0], 1)],
            )

        closed = run([[0, 1]], [0.2, 0.8])
        self.assertEqual(closed[0], closed[1])
        self.assertEqual(closed, run([[0, 1]], [0.2, 0.2]))
        self.assertEqual(run([[0, 1]], [0, 1]), {0: reference.ZERO, 1: reference.ZERO})
        opened = run([[0], [1]], [0.2, 0.8])
        self.assertNotEqual(opened[0], opened[1])
        surface = self.sample_surface(closedU=True)
        pts = [(0.5, 0, 2)] * 12
        c = self.cylinder()
        u = reference.surface_displacements(
            surface,
            pts,
            [(2, 0, 2)] * 12,
            reference.make_segs({(0, 0): c}, {(0, 0): c}, False),
            [c],
            reference.ZERO,
            (0, 0, 1),
            0.2,
            members=set(range(11)),
        )
        self.assertEqual(u[2], reference.ZERO)
        self.assertEqual(u[11], reference.ZERO)
        self.assertEqual(u[0], u[9])
        self.assertGreater(length(u[0]), 0)

    def test_sample_paint_applied_after_each_group_solve(self):
        q = [0.2, 0.6]
        expected = [0.0, 0.0]
        for sweeps in range(1, 11):
            # One half/half sample gives hinge weight 8 and local diagonal 3.
            expected[0] = q[0] * (5 - 2 * expected[1]) / 3
            expected[1] = q[1] * (5 - 2 * expected[0]) / 3
            actual = reference.sample_couple(
                [reference.ZERO] * 2,
                [reference.ZERO] * 2,
                [{0: 0.5, 1: 0.5}],
                [[0], [1]],
                [[], []],
                q,
                lambda p, r: [((1, 0, 0), 1, 1)],
                iters=sweeps,
                regenerations=0,
            )
            for i in range(2):
                self.assertAlmostEqual(actual[i][0], expected[i], places=14)
                self.assertEqual(actual[i][1:], (0.0, 0.0))

    def test_sample_periodic_paint_and_membership_are_canonical(self):
        c = self.cylinder()
        for pu, pv in ((True, False), (False, True), (True, True)):
            surface = self.sample_surface(formU=3 if pu else 1, formV=3 if pv else 1)
            nu, nv = surface["nu"], surface["nv"]
            uu, vv = nu - surface["du"] if pu else nu, nv - surface["dv"] if pv else nv

            def canonical(i):
                return (i // nv % uu) * nv + i % nv % vv

            duplicate = next(i for i in range(nu * nv) if canonical(i) != i)
            original = canonical(duplicate)

            def solve(paint, members=None):
                return reference.surface_displacements(
                    surface,
                    [(0.5, 0, 2)] * (nu * nv),
                    [(2, 0, 2)] * (nu * nv),
                    [[reference.Seg(c, c, 0, True, True)]],
                    [],
                    reference.ZERO,
                    (0, 0, 1),
                    0.2,
                    q=paint,
                    members=members,
                )

            paint = [1.0] * (nu * nv)
            full = solve(paint)
            self.assertEqual(len(full), uu * vv)
            self.assertGreater(length(full[original]), 0.1)
            paint[duplicate] = 0.2
            self.assertEqual(full, solve(paint))
            paint[original] = 0.2
            reduced = solve(paint)
            self.assertLess(length(reduced[original]), length(full[original]))
            for duplicate_paint in (0.0, 1.0):
                paint[duplicate] = duplicate_paint
                self.assertEqual(reduced, solve(paint))
            paint[original] = 0.0
            self.assertEqual(solve(paint)[original], reference.ZERO)
            real_members = {i for i in range(nu * nv) if canonical(i) == i}
            self.assertEqual(full, solve([1.0] * (nu * nv), real_members))
            fixed = solve([1.0] * (nu * nv), real_members - {original})
            self.assertEqual(fixed[original], reference.ZERO)

    def test_sample_all_outside_and_zero_paint(self):
        for q, distance in (([1, 1], -1), ([0, 0], 1)):
            u = reference.sample_couple(
                [(0, 0, 0)] * 2,
                [(0, 0, 0)] * 2,
                [{0: 0.3, 1: 0.7}],
                [[0], [1]],
                [[(1, 0.5)], [(0, 0.5)]],
                q,
                lambda p, r: [((1, 0, 0), distance - p[0], 1)],
            )
            self.assertEqual(u, {0: reference.ZERO, 1: reference.ZERO})

    def test_default_collision_schedule_and_wrapper_propagation(self):
        points = [(0, 0, 0)]
        rows, groups, edges, q = [{0: 1}], [[0]], [[]], [1]
        original_solve = reference.solve_group
        for iters, regenerations, expected_sweeps in ((None, None, 3), (10, 1, 10), (8, 1, 8)):
            solves, events = [], []

            def generate(position, rest):
                return []

            def solve(*args):
                solves.append(args)
                return reference.ZERO

            try:
                reference.solve_group = solve
                reference.sample_couple(
                    points, points, rows, groups, edges, q, generate,
                    iters=iters, regenerations=regenerations, events=events,
                )
            finally:
                reference.solve_group = original_solve
            self.assertEqual(len(solves), expected_sweeps)
            expected_events = [0] if iters is None else ([0, 5] if iters == 10 else [0, 4])
            self.assertEqual([event[0] for event in events], expected_events)

        calls = []
        def sample_spy(*args, **kwargs):
            calls.append((args, kwargs))
            return original_sample(*args, **kwargs)

        original_sample = reference.sample_couple
        reference.sample_couple = sample_spy
        try:
            data_path = Path(__file__).resolve().parent / "fixtures/long_skirt_poses/frame_1.json"
            data = json.loads(data_path.read_text())
            reference.fixture_output(data)
            self.assertIsNone(calls[-1][0][7])
            self.assertIsNone(calls[-1][0][9])
            surface = dict(data["input"], closedU=True)
            current, rest = reference.fixture_cylinders(data)
            reference.surface_displacements(
                surface, [tuple(p[:3]) for p in surface["cvs"]],
                [tuple(p[:3]) for p in data["rest"]["cvs"]],
                reference.make_segs(current, rest, bool(data["skirtType"])),
                [rest[k] for k in sorted(rest)], tuple(data["matrices"]["restBellMatrix"][12:15]),
                tuple(data["matrices"]["restBellMatrix"][4:7]), data["falloff"],
            )
            self.assertIsNone(calls[-1][0][7])
            self.assertIsNone(calls[-1][0][9])
        finally:
            reference.sample_couple = original_sample

    def test_topology_assertion_helper_runs_all_checks_after_reference_failure(self):
        observations = [
            ((0, False), [(0, 0, 0)], [(9, 0, 0)], [(1, 0, 0)]),
            ((0, True), [(0, 0, 0)], [(2, 0, 0)], [(2, 0, 0)]),
            ((0, False), [(0, 0, 0)], [(1, 0, 0)], [(1, 0, 0)]),
        ]
        calls = []
        class Probe(unittest.TestCase):
            def __init__(self):
                super().__init__()
                self.greater_calls = 0
                self.equal_calls = 0

            def runTest(self):
                def compare(test, original, actual, expected):
                    calls.append(actual)
                    if len(calls) == 1:
                        test.fail("reference mismatch")
                assert_topology_observations(self, observations, compare)

            def assertGreater(self, *args, **kwargs):
                self.greater_calls += 1
                return super().assertGreater(*args, **kwargs)

            def assertEqual(self, *args, **kwargs):
                self.equal_calls += 1
                return super().assertEqual(*args, **kwargs)

        result = unittest.TestResult()
        probe = Probe()
        probe.run(result)
        self.assertEqual(len(calls), 3)
        self.assertEqual(len(result.failures), 2)
        self.assertEqual(probe.greater_calls, 4)
        self.assertEqual(probe.equal_calls, 1)
        self.assertTrue(any("reference mismatch" in failure[1] for failure in result.failures))
        self.assertTrue(all("KeyError" not in failure[1] for failure in result.failures))
        self.assertTrue(any("(9, 0, 0)" in failure[1] for failure in result.failures))

    def test_sample_same_index_continuity(self):
        selected = {"taper-side", "corner-side", "corner-cap", "corner-inside", "station-kink"}
        cases = [c for c in continuity_cases() if c["name"] in selected]
        uniform = self.cylinder()
        cases.append(
            dict(
                name="far-side",
                center=0,
                path=lambda t: (-0.5, t, 2),
                rest=(2, 0, 2),
                segments=reference.make_segs({(0, 0): uniform}, {(0, 0): uniform}, False)[0],
            )
        )
        for case in cases:

            def run(t):
                pts = [case["path"](case["center"] + (i - 2) * 0.00025 + t) for i in range(5)]
                return reference.sample_couple(
                    pts,
                    [case["rest"]] * 5,
                    [{i: 0.25, i + 1: 0.75} for i in range(4)],
                    [[i] for i in range(5)],
                    [[(j, 0.2) for j in (i - 1, i + 1) if 0 <= j < 5] for i in range(5)],
                    [1] * 5,
                    lambda p, r: reference.leg_constraints(case["segments"], p, r, reference.ZERO, (0, 0, 1), 0.2)
                    + [(n, d - reference.dot(n, p), w) for n, d, w in case.get("push", [])],
                )

            quotients = []
            for epsilon in (1e-5, 5e-6, 2.5e-6):
                a, b = run(-epsilon), run(epsilon)
                difference = max(length(sub(a[i], b[i])) for i in a)
                quotients.append(difference / (2 * epsilon))
            print("sample continuity", case["name"], quotients)
            for coarse, fine in zip(quotients, quotients[1:]):
                self.assertLessEqual(fine, max(1e-5, 1.5 * coarse))

        def paint_run(q):
            return reference.sample_couple(
                [(0, 0, 0)] * 2,
                [(0, 0, 0)] * 2,
                [{0: 0.5, 1: 0.5}],
                [[0], [1]],
                [[(1, 1)], [(0, 1)]],
                [q, 1],
                lambda p, r: [((1, 0, 0), 1 - p[0], 1)],
            )

        fixed = paint_run(0)
        ratios = []
        for epsilon in (1e-5, 5e-6, 2.5e-6):
            actual = paint_run(epsilon)
            ratios.append(max(length(sub(actual[i], fixed[i])) for i in fixed) / epsilon)
        for coarse, fine in zip(ratios, ratios[1:]):
            self.assertLessEqual(fine, max(1e-5, 1.5 * coarse))

    def test_sample_regenerated_active_switch_continuity(self):
        def run(t):
            events = []
            original = [(t, 0, 0), (t, 0.1, 0)]
            u = reference.sample_couple(
                original,
                original,
                [{0: 1}, {1: 1}],
                [[0], [1]],
                [[(1, 1)], [(0, 1)]],
                [1, 1],
                lambda p, r: [((1, 0, 0), 1 - p[0], 1), ((1, 0, 0), 0.8 - p[0], 1)],
                events=events,
                regenerations=1,
            )
            return u, events

        ratios = []
        for epsilon in (1e-5, 5e-6, 2.5e-6):
            a, events = run(-epsilon)
            b, _ = run(epsilon)
            self.assertLess(events[0][2][0], 0.8)
            self.assertGreater(events[2][2][0], 0.8)
            ratios.append(max(length(sub(a[i], b[i])) for i in a) / (2 * epsilon))
        for coarse, fine in zip(ratios, ratios[1:]):
            self.assertLessEqual(fine, max(1e-5, 1.5 * coarse))

    def test_sample_synced_fixtures(self):
        directory = Path(__file__).resolve().parent / "fixtures/long_skirt_poses"
        golden = json.loads(
            (Path(__file__).resolve().parent / "fixtures/collide_golden_maya2026.json").read_text()
        )
        self.assertEqual(golden["schema"], "yddColliders.collideGolden.gauss3.v1")
        cases = {(case["frame"], case["closedU"]): case for case in golden["cases"]}
        self.assertEqual(set(cases), {(frame, closed) for frame in (0, 1, 2, 10) for closed in (False, True)})
        for frame in (0, 1, 2, 10):
            with self.subTest(frame=frame):
                input_path = directory / ("frame_%d.json" % frame)
                data = json.loads(input_path.read_text())
                result = reference.fixture_output(data, closed=True)
                nv, nu = data["input"]["nv"], data["input"]["nu"]
                for v in range(nv):
                    self.assertLessEqual(length(sub(result[v], result[(nu - 1) * nv + v])), 1e-12)
                if frame == 0:
                    self.assertEqual(result, [tuple(p) for p in data["input"]["cvs"]])
                opened = reference.fixture_output(data, closed=False)
                if frame == 0:
                    self.assertEqual(opened, result)
                else:
                    self.assertGreater(max(length(sub(opened[v], opened[(nu - 1) * nv + v])) for v in range(nv)), 1e-5)
                for closed, expected in ((True, result), (False, opened)):
                    case = cases[frame, closed]
                    self.assertEqual(case["inputSha256"], hashlib.sha256(input_path.read_bytes()).hexdigest())
                    self.assertEqual(case["cvCount"], len(expected))
                    difference = max(length(sub(a, b)) for a, b in zip(expected, case["cvs"]))
                    self.assertLessEqual(difference, 1e-6, "golden scalar reference drift")

    def test_transport_true_inverse_and_fallback_thresholds(self):
        shear = ((1, 0.4, 0), (0, 1, 0.2), (0.1, 0, 1))
        current = self.cylinder()
        rest = self.cylinder(joint=shear)
        side, mode = reference.transported_side((1, 0, 0), current, rest)
        self.assertEqual(mode, "joint")
        rotation = reference.joint_rotation(shear)
        transported = reference.row_vector((1, 0, 0), reference.inverse(rotation))
        self.assertLessEqual(length(sub(reference.row_vector(transported, rotation), (1, 0, 0))), 1e-12)
        expected = reference.unit(reference.radial(transported, current.axis))
        self.assertLessEqual(length(sub(side, expected)), 1e-12)
        transpose = tuple(zip(*rotation))
        wrong = reference.unit(reference.radial(reference.row_vector((1, 0, 0), transpose), current.axis))
        self.assertGreater(length(sub(side, wrong)), 0.1)
        for det, valid in ((0.499999, False), (0.500001, True)):
            matrix = ((1, 0, 0), (math.sqrt(1 - det * det), det, 0), (0, 0, 1))
            self.assertEqual(reference.joint_rotation(matrix) is not None, valid)
        for norm, valid in ((0.999e-8, False), (1e-8, True)):
            self.assertEqual(reference.joint_rotation(((norm, 0, 0), (0, 1, 0), (0, 0, 1))) is not None, valid)
        for epsilon, expected_mode in ((0.999e-6, "quaternion"), (1.001e-6, "joint")):
            z = math.sqrt(1 - epsilon * epsilon)
            current = self.cylinder(joint=((0, epsilon, z), (1, 0, 0), (0, z, -epsilon)))
            side, mode = reference.transported_side((1, 0, 0), current, self.cylinder())
            self.assertEqual(mode, expected_mode)
            self.assertLessEqual(length(sub(side, (1, 0, 0) if mode == "quaternion" else (0, 1, 0))), 1e-12)

    def test_station_interpolation_and_ellipse_radius(self):
        cylinder = self.cylinder(stations=[(0, 1, 2), (2, 2, 1), (4, 3, 2)])
        for z, expected in (
            (-1, (1, 2)),
            (0, (1, 2)),
            (1, (1.5, 1.5)),
            (2, (2, 1)),
            (3, (2.5, 1.5)),
            (5, (3, 2)),
        ):
            self.assertEqual(cylinder.section_radii(z), expected)
        for name, priority in (("ankle", 2.4), ("knee", 1.6)):
            knots = reference.profile_stations(
                thighPosition=1,
                calfPosition=1,
                thighRadiusX=1.2,
                kneeRadiusX=1.6,
                calfRadiusX=2,
                ankleRadiusX=2.4,
            )
            self.assertEqual(dict((s, a) for s, a, _ in knots)[1 if name == "ankle" else 0.5], priority)
        for radius in (0, -1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                self.cylinder(stations=[(0, radius, 1), (4, 1, 1)])

    def test_weighted_active_sets_and_shift(self):
        for weight in (0.0, 1e-20, 1e-13, 0.01, 0.5, 1.0):
            planes = [((1.0, 0.0, 0.0), 2.0, weight)]
            delta = reference.correction(planes, (0.5, 0.0, 0.0))
            self.assertAlmostEqual(
                delta[0], 5 * weight / (1 + 4 * weight) * 1.5, delta=1e-30 if weight < 1e-12 else 1e-12
            )
        for d in (0.999e-9, 1.001e-9):
            self.assertAlmostEqual(reference.project_w([((1, 0, 0), d, 1)])[0], d, delta=1e-20)
        cases = [
            ([], reference.ZERO),
            ([((1, 0, 0), -1, 1)], reference.ZERO),
            ([((1, 0, 0), 1, 1), ((0, 1, 0), 2, 1)], (1, 2, 0)),
            ([((1, 0, 0), 1, 1), ((-1, 0, 0), 1, 1)], reference.ZERO),
            ([((1, 0, 0), 1, 1)] * 4, (20 / 17, 0, 0)),
        ]
        for planes, expected in cases:
            for order in (planes, list(reversed(planes))):
                self.assertLessEqual(length(sub(reference.project_w(order), expected)), 1e-12)
        self.assertEqual(reference.correction([((1, 0, 0), 1, 1)], (2, 0, 0)), reference.ZERO)

    def test_local_geometry_and_ellipse_normal(self):
        c = self.cylinder()
        for r in (0, 0.999999e-8, 1e-8, 1.000001e-8):
            self.assertEqual(reference.local(c, (0, r, 2))["e"], (1, 0, 0) if r < 1e-8 else (0, 1, 0))
        segment = reference.Seg(c, c, 0, True, True)
        for point in ((0.2, 0, 2), (0.95, 0, 2)):
            d, who = reference.union_exit([segment], point, (1, 0, 0))
            self.assertIs(who, segment)
            self.assertAlmostEqual(d, 1 - point[0], delta=2 * math.hypot(4, 1) * 2**-40)
        ell = self.cylinder(stations=[(0, 2, 1), (4, 2, 1)])
        radius = math.sqrt(1.6)
        point = (radius / math.sqrt(2), radius / math.sqrt(2), 2)
        n = reference.side_normal(ell, reference.local(ell, point), 0.2)
        expected = reference.unit((point[0] / 4, point[1], 0))
        self.assertLess(length(sub(n, expected)), 1e-12)
        self.assertGreater(length(sub(n, reference.unit((1, 1, 0)))), 0.1)

    def assert_intervals_close(self, actual, expected, tolerance=1e-12):
        merged = []
        for start, end in sorted(interval[:2] for interval in actual):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        self.assertEqual(len(merged), len(expected), (merged, expected))
        for interval, target in zip(merged, expected):
            for value, wanted in zip(interval, target):
                self.assertAlmostEqual(value, wanted, delta=tolerance)

    def test_polynomial_roots_degrees_reentry_and_tangencies(self):
        cases = [
            ([0, 0, 0, 0, 0], [(0, 5)]),
            ([-1, 0, 0, 0, 0], [(0, 5)]),
            ([1, 0, 0, 0, 0], []),
            ([-2, 1, 0, 0, 0], [(0, 2)]),
            ([4, -4, 1, 0, 0], [(2, 2)]),
            ([-4, 4, -1, 0, 0], [(0, 5)]),
            ([-24, 50, -35, 10, -1], [(0, 1), (2, 3), (4, 5)]),
            ([4, 0, -4, 0, 1], [(math.sqrt(2), math.sqrt(2))]),
            ([1, 0, -6, 0, 9], [(math.sqrt(1 / 3), math.sqrt(1 / 3))]),
            ([4, 0, -12, 0, 9], [(math.sqrt(2 / 3), math.sqrt(2 / 3))]),
            ([1 + 2**-48, 0, -6, 0, 9], []),
            ([1 - 2**-48, 0, -6, 0, 9], [(math.sqrt((1 - 2**-24) / 3), math.sqrt((1 + 2**-24) / 3))]),
        ]
        for coefficients, expected in cases:
            with self.subTest(coefficients=coefficients):
                self.assert_intervals_close(reference.polynomial_intervals(coefficients, 0, 5), expected)
        for scale in (1e-20, 1, 1e20):
            roots = reference.polynomial_roots([-scale, 0, 0, 0, scale * 1e-16], 0, 20000)
            self.assertEqual(len(roots), 1)
            self.assertAlmostEqual(roots[0], 10000, delta=1e-11)

    def test_ray_short_intervals_at_arbitrary_positions(self):
        for position in (0.0, 0.00037, 0.017, 0.31, 0.735, 0.99):
            with self.subTest(position=position):
                c = self.cylinder(origin=(0, 0, position), extent=1e-5, stations=[(0, 1, 1), (1e-5, 1, 1)])
                self.assert_intervals_close(
                    reference.ray_intervals(c, (0.5, 0, 0), (0, 0, 1), 1), [(position, position + 1e-5)]
                )
        c = self.cylinder()
        segment = reference.Seg(c, c, 0, True, True)
        for depth in (1e-8, 1e-6, 0.03, 0.12, 0.14):
            distance, who = reference.union_exit([segment], (1 - depth, 0, 2), (1, 0, 0))
            self.assertIs(who, segment)
            self.assertAlmostEqual(distance, depth, delta=1e-13)

    def test_ray_non_similar_ellipses_quartic(self):
        c = self.cylinder(stations=[(0, 1, 2), (4, 5, 4)])
        with localcontext() as context:
            context.prec = 70
            lo, hi = Decimal(0), Decimal(4)
            for _ in range(240):
                t = (lo + hi) / 2
                value = 1 / (1 + t) ** 2 + 1 / (2 + t / 2) ** 2 - 1
                if value > 0:
                    lo = t
                else:
                    hi = t
            entry = float((lo + hi) / 2)
        self.assert_intervals_close(reference.ray_intervals(c, (1, 1, 0), (0, 0, 1), 5), [(entry, 4)], 5 / 2**46)
        self.assert_intervals_close(reference.ray_intervals(c, (1, 1, 4), (0, 0, -1), 5), [(0, 4 - entry)], 5 / 2**46)

    def test_ray_tangent_axial_degeneracy_and_caps(self):
        c = self.cylinder(stations=[(0, 1, 2), (4, 3, 1)])
        self.assert_intervals_close(reference.ray_intervals(c, (-2, 1.5, 2), (1, 0, 0), 5), [(2, 2)])
        segment = reference.Seg(c, c, 0, True, True)
        self.assertEqual(reference.union_exit([segment], (0, 1.5, 2), (1, 0, 0)), (0, segment))
        for z, radius in ((0, 1), (2, 2), (4, 3)):
            self.assert_intervals_close(reference.ray_intervals(c, (0, 0, z), (1, 0, 0), 5), [(0, radius)])
        for z in (-1e-12, 4 + 1e-12):
            self.assertEqual(reference.ray_intervals(c, (0, 0, z), (1, 0, 0), 5), [])
        cylinder = self.cylinder()
        self.assert_intervals_close(reference.ray_intervals(cylinder, (1, 0, 0), (0, 0, 1), 5), [(0, 4)])
        self.assert_intervals_close(reference.ray_intervals(cylinder, (0, 0, 2), (0, 0, -1e-12), 5), [(0, 5)])

    def test_ray_boundary_keeps_zero_distance_constraint(self):
        c = self.cylinder()
        segment = reference.Seg(c, c, 0, True, True)
        point = (1, 0, 2)
        self.assertEqual(reference.union_exit([segment], point, (1, 0, 0)), (0, segment))
        self.assertEqual(reference.union_exit([segment], (2, 0, 2), (1, 0, 0)), (0, None))
        constraints = reference.leg_constraints([segment], point, (2, 0, 2), reference.ZERO, (0, 0, 1), 0.2)
        self.assertEqual(constraints, [((1, 0, 0), 0, 1)])
        for a, b in ((1, 1), (2, 1), (1.7, 0.9), (16.8, 14.1)):
            c = self.cylinder(stations=[(0, a, b), (4, a, b)])
            segment = reference.Seg(c, c, 0, True, True)
            for j in range(1, 100):
                angle = j * math.pi / 50
                point = (a * math.cos(angle), b * math.sin(angle), 2)
                if reference.inside_closed(c, point):
                    normal = reference.unit((point[0] / (a * a), point[1] / (b * b), 0))
                    distance, who = reference.union_exit([segment], point, normal)
                    self.assertIs(who, segment)
                    self.assertAlmostEqual(distance, 0, delta=1e-13)
                    constraints = reference.leg_constraints(
                        [segment], point, (2 * a, 0, 2), reference.ZERO, (0, 0, 1), 0.2
                    )
                    self.assertEqual(len(constraints), 1)
                    self.assertAlmostEqual(constraints[0][1], 0, delta=1e-13)

    def test_ray_reentry_bridged_through_knee(self):
        thigh = self.cylinder(stations=[(0, 1, 1), (2, 0.25, 0.25), (4, 1, 1)])
        calf = self.cylinder(origin=(0, 0, 1), extent=2, stations=[(0, 1, 1), (2, 1, 1)])
        segments = [reference.Seg(thigh, thigh, 0, True, False), reference.Seg(calf, calf, 1, False, True)]
        self.assert_intervals_close(reference.ray_intervals(thigh, (0.5, 0, 0), (0, 0, 1), 6), [(0, 4 / 3), (8 / 3, 4)])
        distance, who = reference.union_exit(segments[:1], (0.5, 0, 0), (0, 0, 1))
        self.assertAlmostEqual(distance, 4 / 3)
        distance, who = reference.union_exit(segments, (0.5, 0, 0), (0, 0, 1))
        self.assertEqual(distance, 4)
        self.assertEqual(who.id, 0)

    def test_ray_merge_tolerance_both_sides(self):
        a = self.cylinder(extent=1, stations=[(0, 1, 1), (1, 1, 1)])
        for factor in (0.99, 1.01):
            gap = factor * 1e-7 * (2 * math.sqrt(2) + 2)
            b = self.cylinder(origin=(0, 0, 1 + gap), extent=1, stations=[(0, 1, 1), (1, 1, 1)])
            segments = [reference.Seg(a, a, 0, True, False), reference.Seg(b, b, 1, False, True)]
            expected = 2 + gap if factor < 1 else 1
            for order in (segments, list(reversed(segments))):
                distance, who = reference.union_exit(order, (0.5, 0, 0), (0, 0, 1))
                self.assertAlmostEqual(distance, expected, delta=1e-12)
                self.assertEqual(who.id, 1 if factor < 1 else 0)

    def test_straight_long_leg_knee_is_not_an_exit(self):
        a = self.cylinder(origin=(0, 0, -4))
        b = self.cylinder()
        segments = [reference.Seg(a, a, 0, True, False), reference.Seg(b, b, 1, False, True)]
        d, who = reference.union_exit(segments, (0.5, 0, -0.1), (0, 0, 1))
        self.assertEqual(who.id, 1)
        self.assertAlmostEqual(d, 4.1, delta=1e-10)
        for z in (-1e-6, 0, 1e-6):
            cons = reference.leg_constraints(segments, (0.5, 0, z), (2, 0, -2), reference.ZERO, (0, 0, 1), 0.2)
            self.assertEqual(len(cons), 1)
            self.assertAlmostEqual(cons[0][0][2], 0, delta=1e-10)
            self.assertAlmostEqual(cons[0][1], 0.5, delta=1e-10)

    def test_outside_omission_and_small_side_gradient(self):
        thigh = self.cylinder(origin=(0, 0, -4), extent=3, stations=[(0, 1, 1), (3, 1, 1)])
        calf = self.cylinder(extent=3, stations=[(0, 1, 1), (3, 1, 1)])
        segments = [reference.Seg(thigh, thigh, 0, True, False), reference.Seg(calf, calf, 1, False, True)]
        for kappa in (1, 0, -0.2):
            self.assertEqual(
                reference.leg_constraints(segments, (0, 0, -0.5), (2, 0, -2), reference.ZERO, (0, 0, 1), kappa),
                [],
            )
        large = self.cylinder(stations=[(0, 1e10, 2e10), (4, 1e10, 2e10)])
        normal = reference.side_normal(large, reference.local(large, (1, 1, 2)), 0.2)
        self.assertLess(length(sub(normal, reference.unit((4, 1, 0)))), 1e-12)

    def test_periodic_and_open_adjacency(self):
        open_grid = reference.surface_adjacency(3, 4, 1, 1, False, False)
        self.assertEqual(open_grid[0], [1, 4])
        for pu, pv in ((True, False), (False, True), (True, True)):
            adjacency = reference.surface_adjacency(6, 7, 3, 3, pu, pv)
            nu, nv = 3 if pu else 6, 4 if pv else 7
            self.assertEqual(len(adjacency), nu * nv)
            for i, neighbours in adjacency.items():
                self.assertEqual(neighbours, sorted(set(neighbours)))
                self.assertLess(i // 7, nu)
                self.assertLess(i % 7, nv)
                for j in neighbours:
                    self.assertIn(i, adjacency[j])
            if pu:
                self.assertIn((nu - 1) * 7, adjacency[0])
            if pv:
                self.assertIn(nv - 1, adjacency[0])

    def test_fixed_ids_and_partial_legs(self):
        c = self.cylinder()
        legs = reference.make_segs({(0, 1): c, (1, 0): c}, {(0, 1): c, (1, 0): c})
        self.assertEqual([[s.id for s in leg] for leg in legs], [[1], [0]])
        for leg in legs:
            values = [reference.local(s.c, (0.5, 0, 2)) for s in leg]
            self.assertEqual(reference.exposures(leg, values, 0.2), ([1.0], [1.0]))

    def test_staircase_equalities(self):
        c = self.cylinder(stations=[(0, 1, 1), (2, 2, 2), (4, 1, 1)])
        for kappa in (0, -0.2, 2):
            self.assertEqual(reference.slopes(c, -1e6, kappa), reference.slopes(c, 0, kappa))
            self.assertEqual(reference.slopes(c, 1e6, kappa), reference.slopes(c, 4, kappa))
        for kappa in (0, -0.2):
            self.assertEqual(reference.slopes(c, 2, kappa), (-0.5, -0.5))
            straight = self.cylinder()
            seg = reference.Seg(straight, straight, 0, True, True)
            normal = reference.seg_normal(seg, (0.5, 0, 0.5), kappa, 1, 1)
            self.assertEqual(normal, (1, 0, 0))
            self.assertEqual(reference.seg_normal(seg, (0, 0, 2), kappa, 1, 1), reference.ZERO)
            values = [dict(phi=0), dict(phi=0)]
            self.assertEqual(reference.exposures([seg, seg], values, kappa), ([1, 0], [0, 1]))
            nu, _, sn = reference.leg_normal_and_anchor([seg, seg], (1, 0, 2), kappa)
            self.assertEqual(nu, (1, 0, 0))
            self.assertEqual(sn, 1)

    def test_coupling_dirichlet_validation_and_order_measurements(self):
        points = [(1 + u, v, 2) for u in range(3) for v in range(3)]
        adj = reference.surface_adjacency(3, 3, 1, 1, False, False)
        cons = [[] for _ in points]
        cons[4] = [((1, 0, 0), 1, 1)]
        c = self.cylinder()
        weights = [1.0] * 9
        weights[1] = 0
        u = reference.couple(points, points, cons, adj, 0.8, [c], weights)
        self.assertEqual(u[1], reference.ZERO)
        self.assertGreater(length(u[3]), 0)
        for q in (1e-3, 5e-4, 2.5e-4):
            weights[1] = q
            other = reference.couple(points, points, cons, adj, 0.8, [c], weights)
            self.assertLess(max(length(sub(a, b)) for a, b in zip(u, other)), q)
        weights[1] = 0
        sixty = reference.couple(points, points, cons, adj, 0.8, [c], weights, 60)
        rev_adj = {8 - i: [8 - j for j in neighbours] for i, neighbours in adj.items()}
        reverse = list(
            reversed(reference.couple(points[::-1], points[::-1], cons[::-1], rev_adj, 0.8, [c], weights[::-1]))
        )
        residual = max(max(0, d - reference.dot(n, u[i])) for i in range(9) for n, d, _ in cons[i])
        print(
            "grid K30-K60={}, index_reverse={}, residual={}".format(
                max(length(sub(a, b)) for a, b in zip(u, sixty)),
                max(length(sub(a, b)) for a, b in zip(u, reverse)),
                residual,
            )
        )
        invalid = points[:]
        invalid[1] = (math.nan, 0, 0)
        result = reference.couple(invalid, invalid, cons, adj, 0.8, [c], weights)
        self.assertTrue(all(reference.finite(v) for v in result))
        self.assertEqual(result[1], reference.ZERO)
        self.assertEqual(reference.point_weights(points[:3], points[:3], [math.nan, -1, 2]), [0, 0, 1])
        self.assertEqual(reference.couple(points, points, cons, adj, 0, [], weights)[4], (1, 0, 0))

    def test_continuously_moving_rest_and_reverse_evaluation(self):
        for case in (case for case in continuity_cases() if case["name"] in ("far-side-push", "ellipse-axis")):
            maxima = []
            for h in (1e-3, 5e-4, 2.5e-4):
                times = [case["center"] + i * h for i in range(-8, 9)]

                def evaluate(t):
                    point = case["path"](t)
                    rest = add(case["rest"], (0, 0.25 * (t - case["center"]), 0))
                    cons = (
                        reference.all_cons([case["segments"]], point, rest, reference.ZERO, (0, 0, 1), 1) + case["push"]
                    )
                    delta = reference.correction(cons, reference.ZERO)
                    return add(delta, reference.correction(cons, delta))

                values = [evaluate(t) for t in times]
                self.assertEqual(values, list(reversed([evaluate(t) for t in reversed(times)])))
                maxima.append(max(length(sub(a, b)) / h for a, b in zip(values, values[1:])))
            self.assertLess(max(maxima), 64)
            for a, b in zip(maxima, maxima[1:]):
                self.assertLessEqual(b, max(1e-5, 1.5 * a))
            print("moving rest {} final quotients={}".format(case["name"], maxima))

    def test_continuity_sweeps(self):
        """These finite fixtures peak at 22.603 locally and 39.417 after final correction.

        The retained bound 64 and h-halving checks exclude the recorded reentrant exit jump.
        """
        for case in continuity_cases():
            measurements = continuity_measurements(case)
            print("continuity {} {}".format(case["name"], measurements))
            if case.get("exit_jump"):
                continue
            for series in measurements.values():
                self.assertTrue(all(math.isfinite(value) for value in series))
                self.assertLess(max(series), 64)
                for coarse, fine in zip(series, series[1:]):
                    self.assertLessEqual(fine, max(1e-5, 1.5 * coarse))


def continuity_cases():
    def cyl(origin=(0, 0, 0), axis=(0, 0, 1), x=(1, 0, 0), z=(0, 1, 0), extent=4, stations=((0, 1, 1), (4, 1, 1))):
        return Cylinder(origin, axis, x, z, extent, stations)

    def leg(*cylinders):
        return [reference.Seg(c, c, i, i == 0, i == len(cylinders) - 1) for i, c in enumerate(cylinders)]

    result = []

    def add_case(name, segments, path, lo, hi, count, rest=(2, 0, 2), push=(), center=0, exit_jump=False):
        result.append(
            dict(
                name=name,
                segments=segments,
                path=path,
                lo=lo,
                hi=hi,
                count=count,
                rest=rest,
                push=list(push),
                center=center,
                exit_jump=exit_jump,
            )
        )

    taper = leg(cyl(stations=((0, 2, 2), (4, 1, 1))))
    uniform = leg(cyl())
    a = cyl(origin=(0, 0, -4))
    straight = leg(a, cyl())
    axis = (0.5, 0, math.sqrt(3) / 2)
    x = reference.unit(reference.cross((0, 1, 0), axis))
    bent = leg(a, cyl(axis=axis, x=x, z=reference.cross(axis, x)))
    bent[1].c0 = straight[1].c
    negative_z = [((0, 0, -1), 1, 1)]
    diagonal = math.sqrt(0.5)
    add_case("taper-side-push", taper, lambda r: (r, 0, 0.2), 1.90, 2.0, 200, (3, 0, 2), negative_z, 1.95)
    add_case("taper-side", taper, lambda r: (r, 0, 2), 1.4, 1.6, 200, (3, 0, 2), center=1.5)
    add_case(
        "far-side-push",
        uniform,
        lambda x: (x, 0, 2),
        -0.9,
        0.9,
        360,
        push=[((diagonal, 0, diagonal), 1, 1)],
        center=-0.5,
    )
    push = [((-diagonal, 0, diagonal), 1, 1)]
    add_case("corner-side", uniform, lambda z: (1.02, 0, z), -0.5, 0.5, 200, push=push)
    add_case("corner-cap", uniform, lambda r: (r, 0, -0.02), 0.5, 1.5, 200, push=push, center=1)
    add_case("corner-inside", uniform, lambda z: (0.98, 0, z), -0.5, 0.5, 200, push=push)
    add_case("corner-diagonal", uniform, lambda t: (1 + t, 0, -t), -0.3, 0.3, 200, push=push)
    add_case("straight-switch", straight, lambda r: (r, 0, -1), 0.9, 1.1, 200, (2, 0, -2), [((0, 0, 1), 2, 1)], 1)
    add_case("straight-knee", straight, lambda z: (0.5, 0, z), -1, 1, 200, (2, 0, -2))
    add_case("bent-outer", bent, lambda z: (-0.5, 0, z), -0.5, 0.5, 200, (2, 0, -2))
    add_case("bent-inner", bent, lambda z: (0.5, 0, z), -0.5, 0.5, 200, (2, 0, -2))
    add_case("bent-extended", bent, lambda z: (-0.5, 0, z), -0.2, 0.6, 400, (2, 0, -2))
    ellipse = leg(cyl(stations=((0, 2, 1), (4, 2, 1))))
    add_case("ellipse-axis", ellipse, lambda y: (0, y, 1.5), -0.5, 0.5, 200, (3, 0, 1.5))
    kink = leg(cyl(stations=((0, 1.5, 1.5), (2, 2, 2), (4, 1.5, 1.5))))
    add_case("station-kink", kink, lambda z: (0.5, 0, z), 1.5, 2.5, 200, (3, 0, 2), center=2)
    add_case("exposed-knee-p1", bent, lambda z: (-0.5, 0, z), -0.2, 0.2, 400, (2, 0, -2), negative_z)
    zc = 1 / (2 * math.sqrt(3))
    add_case("exposed-knee-p2", bent, lambda t: (-0.5, 0, zc + t), -0.2, 0.2, 400, (2, 0, -2), negative_z)
    short = leg(cyl(extent=0.1, stations=((0, 1, 1), (0.1, 1, 1))))
    add_case("short-cap-axis", short, lambda z: (0, 0, z), -0.05, 0.05, 400, (2, 0, 0.05), [((0, 0, 1), 1, 1)])
    gap = leg(
        cyl(origin=(0, 0, -4), extent=3, stations=((0, 1, 1), (3, 1, 1))),
        cyl(extent=3, stations=((0, 1, 1), (3, 1, 1))),
    )
    add_case("gap-axis", gap, lambda y: (0, y, -0.5), -0.1, 0.1, 400, (2, 0, -2), [((-1, 0, 0), 1, 1)])
    right_angle = leg(a, cyl(axis=(1, 0, 0), x=(0, 1, 0), z=(0, 0, 1)))
    right_angle[1].c0 = straight[1].c
    add_case("right-angle-boundary", right_angle, lambda z: (-0.1, 0, z), -8e-6, 8e-6, 16, (2, 0, -2), negative_z)
    add_case(
        "right-angle-exit-jump",
        right_angle,
        lambda z: (-0.1, 0, z),
        -0.2,
        0.2,
        400,
        (2, 0, -2),
        negative_z,
        exit_jump=True,
    )
    return result


def continuity_measurements(case):
    measurements = {name: [] for name in ("local", "final", "boundary_local", "boundary_final")}
    for divisor in (1, 2, 4):
        count = case["count"] * divisor
        h = (case["hi"] - case["lo"]) / count
        for boundary in (False, True):
            ts = (
                [case["center"] + i * h for i in range(-8, 9)]
                if boundary
                else [case["lo"] + (case["hi"] - case["lo"]) * i / count for i in range(count + 1)]
            )
            samples = []
            for t in ts:
                p = case["path"](t)
                cons = (
                    reference.all_cons([case["segments"]], p, case["rest"], reference.ZERO, (0, 0, 1), 1.0)
                    + case["push"]
                )
                local = reference.correction(cons, reference.ZERO)
                final = add(local, reference.correction(cons, local))
                samples.append((local, final))
            for component, name in enumerate(("local", "final")):
                key = ("boundary_" if boundary else "") + name
                measurements[key].append(
                    max(length(sub(a[component], b[component])) / h for a, b in zip(samples, samples[1:]))
                )
    return measurements
