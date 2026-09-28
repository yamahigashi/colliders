import math
import unittest
import copy
import hashlib
import json
from pathlib import Path
import struct

import collide_reference as cr

import maya.cmds as cmds
import maya.api.OpenMaya as om

from collide_reference import (
    Cylinder,
    add,
    dot,
    joint_rotation,
    length,
    mul,
    rest_contact,
    sub,
    transported_side,
)
from test_collide_projection import (
    assert_topology_observations,
    assert_displacements_close,
    node_cylinders,
    node_reference,
    surface_snapshot,
)
from helpers import bind_rest_inputs, duplicate_rest_geometry, coordinates as geometry_coordinates


def coordinates(node, attribute, kind="mesh"):
    if kind == "surface":
        snapshot = surface_snapshot(node + "." + attribute)
        return [list(p) + [w] for p, w in zip(snapshot["cvs"], snapshot["weights"])]
    return geometry_coordinates(node, attribute, kind)


def matrix(x=0, y=0, z=0):
    return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]


def distance(point, reference):
    return math.sqrt(sum((value - target) ** 2 for value, target in zip(point, reference)))


class DeformerTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def mesh(self, points):
        fn = om.MFnMesh()
        obj = fn.create(
            om.MPointArray([om.MPoint(*point) for point in points]), [len(points)], list(range(len(points)))
        )
        return om.MFnDagNode(obj).fullPathName()

    def points(self, mesh):
        shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
        path = om.MSelectionList().add(shape).getDagPath(0)
        return [tuple(p)[:3] for p in om.MFnMesh(path).getPoints()]

    def world_points(self, mesh):
        shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
        path = om.MSelectionList().add(shape).getDagPath(0)
        return [tuple(p)[:3] for p in om.MFnMesh(path).getPoints(om.MSpace.kWorld)]

    def transform_matrix(self, translation, rotation, scale):
        transform = om.MTransformationMatrix()
        transform.setTranslation(om.MVector(*translation), om.MSpace.kTransform)
        transform.setRotation(om.MEulerRotation(*(math.radians(value) for value in rotation)))
        transform.setScale((scale, scale, scale), om.MSpace.kTransform)
        return list(transform.asMatrix())

    def assert_world_is_one_transform(self, mesh, object_points):
        world_matrix = om.MMatrix(cmds.getAttr(mesh + ".worldMatrix[0]"))
        expected = [tuple(om.MPoint(*point) * world_matrix)[:3] for point in object_points]
        self.assert_points_close(self.world_points(mesh), expected, tolerance=1e-9)

    def collision(self, mesh, long=False):
        geometry = mesh[0].split(".", 1)[0] if isinstance(mesh, (list, tuple)) else mesh
        rest = duplicate_rest_geometry(geometry)
        node = cmds.deformer(mesh, type="yddSkirtCollideDeformer")[0]
        for side, x in [("left", 0), ("right", 10)]:
            for joint, y in [("Hip", 0), ("Knee", 1), ("Heel", 2)]:
                cmds.setAttr(
                    node + "." + side + joint + "Matrix",
                    *matrix(x, y if side == "left" else 0),
                    type="matrix",
                )
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        cmds.setAttr(node + ".skirtType", int(long))
        bind_rest_inputs(node, rest)
        return node

    def collision_reference(self, node, original):
        actual_input, _, expected, _, _ = node_reference(node)
        self.assertEqual(actual_input, original)
        return expected

    def assert_collision_reference(self, mesh, node, original):
        expected = self.collision_reference(node, original)
        actual = self.points(mesh)
        assert_displacements_close(self, original, actual, expected)
        return actual

    def wave(self, mesh):
        node = cmds.deformer(mesh, type="yddSkirtWaveDeformer")[0]
        for attr, value in [
            ("amplitude", 0.7),
            ("idleAmplitude", 0.3),
            ("impulseX", 0),
            ("impulseZ", 0),
            ("noiseAmplitude", 0),
            ("wavePhaseV", 0.23),
            ("wavePhaseU", 0.17),
            ("waveCountV", 1.3),
            ("waveCountU", 2),
            ("skew", 0.3),
            ("sharpness", 0.3),
        ]:
            cmds.setAttr(node + "." + attr, value)
        for i, position in enumerate([0, 1]):
            prefix = node + ".amplitudeRamp[%d]." % i
            cmds.setAttr(prefix + "amplitudeRamp_Position", position)
            cmds.setAttr(prefix + "amplitudeRamp_FloatValue", 1)
            cmds.setAttr(prefix + "amplitudeRamp_Interp", 1)
        return node

    def test_collision_rest_caps_and_strength(self):
        for long in (False, True):
            with self.subTest(long=long):
                end = 2 if long else 1
                mesh = self.mesh([(0.25, y, 0) for y in (-0.001, 0, 0.001, 0.5, end - 0.001, end, end + 0.001)])
                original = self.points(mesh)
                node = self.collision(mesh, long)
                result = self.assert_collision_reference(mesh, node, original)
                cmds.setAttr(node + ".falloff", 0)
                hard = self.assert_collision_reference(mesh, node, original)
                for index in (0, 6):
                    self.assertEqual(hard[index], original[index])
                for index in (1, 5):
                    self.assertGreater(distance(hard[index], original[index]), 0.1)
                cmds.setAttr(node + ".falloff", 0.2)
                self.assertGreater(distance(result[3], original[3]), 0.1)
                full = self.points(mesh)
                cmds.setAttr(node + ".envelope", 0.5)
                self.assert_points_close(
                    self.points(mesh), [add(p, mul(sub(q, p), 0.5)) for p, q in zip(original, full)]
                )
                cmds.setAttr(node + ".envelope", 0)
                self.assertEqual(self.points(mesh), original)
                cmds.setAttr(node + ".envelope", 1)
                cmds.setAttr(node + ".weightList[0].weights[3]", 0)
                self.assertEqual(self.points(mesh)[3], original[3])

    def set_profile(self, node, **values):
        for name, value in values.items():
            cmds.setAttr(node + "." + name, value)

    def test_leg_profile_attribute_contract(self):
        radii = [station + "Radius" + axis for station in ("thigh", "knee", "calf", "ankle") for axis in ("X", "Z")]
        for node_type in ("yddSkirtBellCollider", "yddSkirtCollideDeformer"):
            node = cmds.createNode(node_type)
            for name in radii + ["thighPosition", "calfPosition"]:
                with self.subTest(node_type=node_type, attribute=name):
                    radius = name in radii
                    self.assertEqual(cmds.getAttr(node + "." + name, type=True), "double")
                    self.assertEqual(cmds.getAttr(node + "." + name), 1 if radius else 0.5)
                    self.assertTrue(cmds.getAttr(node + "." + name, keyable=True))
                    self.assertEqual(cmds.attributeQuery(name, node=node, minimum=True), [0.001 if radius else 0])
                    self.assertEqual(cmds.attributeQuery(name, node=node, maxExists=True), not radius)
                    if not radius:
                        self.assertEqual(cmds.attributeQuery(name, node=node, maximum=True), [1])

    def test_leg_profile_station_radii_increase_support_and_displacement(self):
        for station, long, y in (
            ("thigh", False, 0.5),
            ("knee", False, 0.5),
            ("calf", True, 1.5),
            ("ankle", True, 1.5),
        ):
            with self.subTest(station=station):
                mesh = self.mesh([(0.2, y, 0), (0, y, 0.2), (0.2, y, 0.1)])
                original = self.points(mesh)
                node = self.collision(mesh, long)
                baseline = self.points(mesh)
                self.set_profile(node, **{station + "RadiusX": 1.4, station + "RadiusZ": 1.8})
                result = self.assert_collision_reference(mesh, node, original)
                self.assertGreater(max(distance(a, b) for a, b in zip(baseline, result)), 1e-5)

    def test_leg_profile_station_positions_change_axial_support(self):
        for station, long, start in (("thigh", False, 0), ("calf", True, 1)):
            mesh = self.mesh([(0.2, start + y, 0) for y in (0.2, 0.5, 0.8)])
            original = self.points(mesh)
            node = self.collision(mesh, long)
            self.set_profile(node, **{station + "RadiusX": 2, station + "RadiusZ": 2, station + "Position": 0.25})
            early = self.assert_collision_reference(mesh, node, original)
            cmds.setAttr(node + "." + station + "Position", 0.75)
            late = self.assert_collision_reference(mesh, node, original)
            self.assertGreater(distance(early[0], late[0]), 1e-5)
            self.assertGreater(distance(early[2], late[2]), 1e-5)
            self.assertTrue(all(math.isfinite(v) for point in (early[1], late[1]) for v in point))

    def test_leg_profile_coincident_thigh_and_calf_stations(self):
        for thigh_position in (0, 1e-9, 1):
            for calf_position in (0, 1 - 1e-9, 1):
                with self.subTest(thigh_position=thigh_position, calf_position=calf_position):
                    mesh = self.mesh([(0.2, y, 0) for y in (0.2, 0.5, 0.8, 1.2, 1.5, 1.8)])
                    original = self.points(mesh)
                    node = self.collision(mesh, long=True)
                    self.set_profile(
                        node,
                        thighPosition=thigh_position,
                        calfPosition=calf_position,
                        thighRadiusX=1.2,
                        thighRadiusZ=1.4,
                        kneeRadiusX=1.6,
                        kneeRadiusZ=1.8,
                        calfRadiusX=2.0,
                        calfRadiusZ=2.2,
                        ankleRadiusX=2.4,
                        ankleRadiusZ=2.6,
                    )
                    self.assert_collision_reference(mesh, node, original)

    def test_leg_profile_long_adds_only_the_lower_segment(self):
        mesh = self.mesh([(0.2, y, 0) for y in (1.2, 1.5, 1.8)])
        original = self.points(mesh)
        node = self.collision(mesh)
        self.assertEqual(self.points(mesh), original)
        cmds.setAttr(node + ".skirtType", 1)
        self.set_profile(node, calfRadiusX=1.8, calfRadiusZ=1.8)
        result = self.assert_collision_reference(mesh, node, original)
        self.assertNotEqual(result, original)

    def test_leg_profile_degenerate_lengths_and_effective_length_guard(self):
        for knee_y, heel_y, scale_y in ((0, 0, 1), (0, 2, 1), (1, 1, 1), (1, 2, 1e-7)):
            with self.subTest(knee=knee_y, heel=heel_y, scale_y=scale_y):
                mesh = self.mesh([(0.2, y, 0) for y in (0.2, 0.5, 1.5)])
                original = self.points(mesh)
                node = self.collision(mesh, long=True)
                for prefix in ("left", "restLeft"):
                    cmds.setAttr(node + "." + prefix + "KneeMatrix", *matrix(0, knee_y), type="matrix")
                    cmds.setAttr(node + "." + prefix + "HeelMatrix", *matrix(0, heel_y), type="matrix")
                cmds.setAttr(node + ".ringScale1", scale_y)
                self.set_profile(node, kneeRadiusX=1.7, kneeRadiusZ=1.7, ankleRadiusX=2.3, ankleRadiusZ=2.3)
                result = self.assert_collision_reference(mesh, node, original)
                if heel_y == 0 or scale_y < 1e-6:
                    self.assertEqual(result, original)

    def test_leg_profile_length_scale_moves_station_support(self):
        mesh = self.mesh([(0.2, y, 0) for y in (0.5, 1.0, 1.5)])
        original = self.points(mesh)
        node = self.collision(mesh)
        self.set_profile(node, thighRadiusX=1.5, thighRadiusZ=1.5)
        self.assert_collision_reference(mesh, node, original)
        cmds.setAttr(node + ".ringScale1", 2)
        after = self.assert_collision_reference(mesh, node, original)
        cmds.setAttr(node + ".ringScale1", 1)
        cmds.setAttr(node + ".falloff", 0)
        self.assertEqual(self.points(mesh)[2], original[2])
        cmds.setAttr(node + ".ringScale1", 2)
        cmds.setAttr(node + ".falloff", 0.2)
        self.assertGreater(distance(after[2], original[2]), 1e-5)

    def test_kappa_depth_ramp_and_rest_attribute_contract(self):
        node = cmds.createNode("yddSkirtCollideDeformer")
        self.assertEqual(cmds.getAttr(node + ".falloff", type=True), "float")
        self.assertAlmostEqual(cmds.getAttr(node + ".falloff"), 0.2)
        self.assertEqual(cmds.attributeQuery("falloff", node=node, minimum=True), [0])
        self.assertEqual(cmds.attributeQuery("falloff", node=node, maximum=True), [1])
        self.assertTrue(cmds.attributeQuery("restGeometry", node=node, hidden=True))
        self.assertFalse(cmds.attributeQuery("restGeometry", node=node, storable=True))
        self.assertFalse(cmds.attributeQuery("restGeometry", node=node, readable=True))
        for name in ("Bell", "LeftHip", "LeftKnee", "LeftHeel", "RightHip", "RightKnee", "RightHeel"):
            attr = "rest" + name + "Matrix"
            self.assertTrue(cmds.attributeQuery(attr, node=node, hidden=True))
            self.assertEqual(cmds.getAttr(node + "." + attr, type=True), "matrix")
            self.assertEqual(cmds.getAttr(node + "." + attr), matrix())

        mesh = self.mesh([(0.95, 0.5, 0), (0.9, 0.5, 0), (0.8, 0.5, 0)])
        original = self.points(mesh)
        node = self.collision(mesh)
        for kappa in (0.05, 0.2, 0.8):
            cmds.setAttr(node + ".falloff", kappa)
            self.assert_collision_reference(mesh, node, original)
            current, rest = node_cylinders(node)
            leg = cr.make_segs(current, rest, False)[0]
            for point in original:
                values = [cr.local(seg.c, point) for seg in leg]
                exp0, exp1 = cr.exposures(leg, values, kappa)
                direction = cr.direction(leg, values, exp0, exp1, point, point, (0, 0, 0), (0, 1, 0), kappa, True)
                d, who = cr.union_exit(leg, point, cr.unit(direction))
                self.assertIsNotNone(who)
                self.assertGreater(d, 0)

    def collision_from_opposite_rest_side(self, radii):
        mesh = self.mesh([(-radius, 0.5, 0) for radius in radii])
        original = self.points(mesh)
        node = self.collision(mesh)
        rest_shape = cmds.listConnections(node + ".restGeometry", source=True, destination=False)[0]
        for index in range(len(original)):
            cmds.xform(rest_shape + ".vtx[{}]".format(index), objectSpace=True, translation=(2, 0.5, 0))
        return mesh, node, original

    def test_collision_shallow_deep_and_hard_mode(self):
        mesh, node, original = self.collision_from_opposite_rest_side((0.95, 0.1, 1.1))
        for kappa in (0.2, 0):
            cmds.setAttr(node + ".falloff", kappa)
            expected = self.assert_collision_reference(mesh, node, original)
            self.assertLessEqual(expected[0][0], -1 + 1e-6)
            self.assertGreaterEqual(expected[1][0], 1 - 1e-6)
            if kappa == 0:
                self.assertEqual(self.points(mesh)[2], original[2])
        mesh = self.mesh([(0.5, y, 0) for y in (-1e-5, 0, 1e-5)])
        original = self.points(mesh)
        node = self.collision(mesh)
        cmds.setAttr(node + ".falloff", 0)
        result = self.assert_collision_reference(mesh, node, original)
        self.assertEqual(result[0], original[0])
        self.assertGreater(distance(result[1], original[1]), 0.4)
        self.assertGreater(distance(result[2], original[2]), 0.4)

    def test_collision_exterior_barriers_remain(self):
        mesh, node, original = self.collision_from_opposite_rest_side((1.1, 2.1, 1.5))
        self.assertEqual(self.points(mesh), original)
        _, _, _, fields, _ = node_reference(node)
        for entry in fields:
            self.assertLess(entry[0][1], 0)
        cylinder = Cylinder((0, 0, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0), 4, [(0, 1, 1), (4, 1, 1)])
        other = Cylinder((-1.8, 0, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0), 4, [(0, 1, 1), (4, 1, 1)])
        legs = [[cr.Seg(c, c, 0, True, True)] for c in (cylinder, other)]
        planes = cr.all_cons(legs, (-1.1, 0, 2), (2, 2, 2), (0, 0, 0), (0, 0, 1), 0.2)
        self.assertEqual(len(planes), 2)
        self.assertLess(planes[0][1], 0)
        delta = cr.project_w(planes)
        alone = cr.project_w(planes[1:])
        self.assertGreater(dot(planes[0][0], delta), dot(planes[0][0], alone))
        self.assertGreater(length(sub(delta, alone)), 1e-3)

    def test_collision_invalid_cylinders_are_excluded_and_recover(self):
        mesh = self.mesh([(0.2, y, 0) for y in (0.2, 0.5, 0.8)])
        original = self.points(mesh)
        node = self.collision(mesh)
        active = self.points(mesh)
        for attribute in ("leftHipMatrix", "restLeftHipMatrix", "leftKneeMatrix", "restLeftKneeMatrix"):
            valid = cmds.getAttr(node + "." + attribute)
            for bad in ([0.0] * 16, matrix(0, 1e-7) if "Knee" in attribute else matrix(0, 1 - 1e-7)):
                cmds.setAttr(node + "." + attribute, *bad, type="matrix")
                self.assertEqual(self.points(mesh), original)
                self.assert_collision_reference(mesh, node, original)
            cmds.setAttr(node + "." + attribute, *valid, type="matrix")
            self.assertEqual(self.points(mesh), active)
        for scale in ((0, 1, 1), (-1, 1, 1), (1, 1, 0), (1, 1, -1)):
            cmds.setAttr(node + ".ringScale", *scale, type="double3")
            self.assertEqual(self.points(mesh), original)
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        self.assertEqual(self.points(mesh), active)

    def sweep_node(self, fixture, samples):
        mesh = self.mesh([sample["point"] for sample in samples])
        node = self.collision(mesh, long=len(fixture["current"]) == 2)
        for prefix in ("left", "restLeft"):
            first = fixture["current"][0]
            knee = add(first.origin, mul(first.axis, first.extent))
            last = fixture["current"][-1]
            heel = (
                add(last.origin, mul(last.axis, last.extent))
                if len(fixture["current"]) == 2
                else add(knee, mul(first.axis, 4))
            )
            for joint, position in zip(("Hip", "Knee", "Heel"), (first.origin, knee, heel)):
                value = matrix(*position)
                if fixture["name"] == "ellipse_45":
                    value[:12] = [0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0]
                cmds.setAttr(node + "." + prefix + joint + "Matrix", *value, type="matrix")
        cmds.setAttr(node + ".leftRingAxis", 1 if fixture["name"] == "ellipse_45" else 0)
        if fixture["name"] == "ellipse_45":
            cmds.setAttr(node + ".ringScale", 2, 1, 1, type="double3")
        bell = matrix()
        bell[:12] = [1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0]
        cmds.setAttr(node + ".restBellMatrix", *bell, type="matrix")
        rest_shape = cmds.listConnections(node + ".restGeometry", source=True, destination=False)[0]
        for index, sample in enumerate(samples):
            cmds.xform(rest_shape + ".vtx[{}]".format(index), objectSpace=True, translation=sample["rest_point"])
        return mesh, node

    def test_collision_continuity_sweeps(self):
        fixture = self.simple_sweep_fixture()
        paths = [
            lambda t: (1 + t, 0, 2),
            lambda t: (0.5, 0, t),
            lambda t: (0.5, 0, 4 + t),
            lambda t: (t, 0, 2),
            lambda t: (-0.5, t, 2),
        ]
        for path_index, path in enumerate(paths):
            maxima = []
            for h in (1e-3, 5e-4, 2.5e-4):
                samples = [dict(point=path(i * h), rest_point=(2, 0, 2)) for i in range(-8, 9)]
                mesh, node = self.sweep_node(fixture, samples)
                original, _, expected, _, _ = node_reference(node)
                actual = self.assert_collision_reference(mesh, node, original)
                deltas = [sub(a, p) for a, p in zip(actual, original)]
                maximum = max(
                    length(sub(a, b)) / length(sub(p, q))
                    for a, b, p, q in zip(deltas, deltas[1:], original, original[1:])
                )
                maxima.append(maximum)
                self.assertLess(maximum, 64)
                cmds.delete(mesh, node)
            print("node continuity path={} quotients={}".format(path_index, maxima))
            for coarse, fine in zip(maxima, maxima[1:]):
                self.assertLessEqual(fine, max(1e-5, 1.5 * coarse))

    def simple_sweep_fixture(self):
        c = Cylinder((0, 0, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0), 4, [(0, 1, 1), (4, 1, 1)])
        return dict(name="uniform", current=[c], rest=[c], rest_point=(2, 0, 2))

    def test_collision_threshold_samples_match_reference(self):
        fixture = self.simple_sweep_fixture()
        positions = [(0, r, 2) for r in (0, 0.999999e-8, 1e-8, 1.000001e-8)]
        positions += [(1, 0, 2), (0.5, 0, 0), (0.5, 0, 4)]
        g = (3 - math.sqrt(5)) / 4
        positions.append((-(1 - g), 0, -0.05))
        samples = [dict(point=point, rest_point=fixture["rest_point"]) for point in positions]
        mesh, node = self.sweep_node(fixture, samples)
        original, _, _, _, _ = node_reference(node)
        self.assert_collision_reference(mesh, node, original)

    def test_collision_short_ray_exits_match_analytic_cylinder(self):
        fixture = self.simple_sweep_fixture()
        depths = (1e-5, 0.00037, 0.017, 0.03, 0.12, 0.14)
        samples = [dict(point=(1 - depth, 0, 2), rest_point=(2, 0, 2)) for depth in depths]
        mesh, node = self.sweep_node(fixture, samples)
        cmds.setAttr(node + ".falloff", 0)
        original, _, expected, _, _ = node_reference(node)
        actual = self.points(mesh)
        for p in actual:
            self.assertAlmostEqual(p[0], 1.0, delta=1e-6)
            self.assertAlmostEqual(p[1], 0.0, delta=1e-6)
            self.assertAlmostEqual(p[2], 2.0, delta=1e-6)
        assert_displacements_close(self, original, actual, expected)

    def test_collision_nonfinite_inputs_do_not_deform_and_recover(self):
        mesh = self.mesh([(0.2, y, 0) for y in (0.2, 0.5, 0.8)])
        original = self.points(mesh)
        node = self.collision(mesh)
        active = self.points(mesh)
        for name in ("falloff", "envelope", "weightList[0].weights[0]"):
            plug = om.MSelectionList().add(node + "." + name).getPlug(0)
            previous = plug.asFloat()
            for value in (math.nan, math.inf, -math.inf):
                plug.setFloat(value)
                self.assertFalse(math.isfinite(plug.asFloat()))
                result = self.points(mesh)
                self.assertEqual(result[0], original[0])
                if name.startswith("weight"):
                    self.assert_collision_reference(mesh, node, original)
                else:
                    self.assertEqual(result[1:], original[1:])
            plug.setFloat(previous)
            self.assertEqual(self.points(mesh), active)
        for name, component in (("restBellMatrix", 5), ("leftHipMatrix", 0), ("restLeftHipMatrix", 0)):
            previous = cmds.getAttr(node + "." + name)
            for bad in (math.nan, math.inf):
                value = list(previous)
                value[component] = bad
                data = om.MFnMatrixData().create(om.MMatrix(value))
                om.MSelectionList().add(node + "." + name).getPlug(0).setMObject(data)
                self.assertEqual(self.points(mesh), original)
            cmds.setAttr(node + "." + name, *previous, type="matrix")
            self.assertEqual(self.points(mesh), active)
        for name in ("ringScale0", "ringScale2"):
            plug = om.MSelectionList().add(node + "." + name).getPlug(0)
            for bad in (math.nan, math.inf):
                plug.setDouble(bad)
                self.assertEqual(self.points(mesh), original)
            plug.setDouble(1)
            self.assertEqual(self.points(mesh), active)
        for attribute in ("input[0].inputGeometry", "restGeometry"):
            for bad in (math.nan, math.inf, -math.inf):
                mesh = self.mesh([(0.2, y, 0) for y in (0.2, 0.5, 0.8)])
                original = self.points(mesh)
                node = self.collision(mesh)
                active = self.points(mesh)
                shape = cmds.listConnections(node + "." + attribute, source=True, destination=False, shapes=True)[0]
                fn = om.MFnMesh(om.MSelectionList().add(shape).getDagPath(0))
                saved = fn.getPoints()
                modified = om.MPointArray(saved)
                modified[0] = om.MPoint(bad, saved[0].y, saved[0].z)
                fn.setPoints(modified)
                result = self.points(mesh)
                if attribute == "restGeometry":
                    self.assertEqual(result[0], original[0])
                else:
                    self.assertFalse(math.isfinite(result[0][0]))
                    self.assertEqual(result[0][1:], original[0][1:])
                expected = node_reference(node)[2]
                self.assert_points_close(result[1:], expected[1:])

    def test_collision_joint_transport_and_quaternion_fallback(self):
        rest_point = (-1.73204, -2.50071, 1.73204)
        mesh = self.mesh([rest_point, rest_point, rest_point])
        node = self.collision(mesh, long=True)
        poses = {
            "left": ((1, -1.5, 0), (-1.64, -1.73, 2.29), (-2.66, -3.22, -1.28)),
            "restLeft": ((1, -1.5, 0), (0.9, -5, 0), (0.8, -9, 0)),
        }
        for prefix, positions in poses.items():
            for joint, position in zip(("Hip", "Knee", "Heel"), positions):
                cmds.setAttr(node + "." + prefix + joint + "Matrix", *matrix(*position), type="matrix")
        cmds.setAttr(node + ".leftRingAxis", 0)
        bell = matrix()
        bell[5] = bell[10] = -1
        cmds.setAttr(node + ".restBellMatrix", *bell, type="matrix")
        current, rest = node_cylinders(node)
        b0 = rest_contact(rest_point, rest[(0, 1)], (0, 0, 0), (0, -1, 0))
        side, mode = transported_side(b0, current[(0, 1)], rest[(0, 1)])
        quaternion_side, _ = transported_side(b0, current[(0, 1)], rest[(0, 1)], True)
        self.assertEqual(mode, "joint")
        self.assertGreater(length(sub(side, quaternion_side)), 0.1)
        point = add(add(current[(0, 1)].origin, mul(current[(0, 1)].axis, current[(0, 1)].extent / 2)), mul(side, 0.1))
        input_shape = cmds.listConnections(
            node + ".input[0].inputGeometry", source=True, destination=False, shapes=True
        )[0]
        fn = om.MFnMesh(om.MSelectionList().add(input_shape).getDagPath(0))
        fn.setPoints(om.MPointArray([om.MPoint(*point) for _ in range(3)]))
        self.points(mesh)
        original, rp, expected, fields, _ = node_reference(node)
        self.assertLessEqual(length(sub(original[0], point)), 1e-6)
        self.assert_collision_reference(mesh, node, original)
        joint_delta = sub(expected[0], original[0])
        quaternion_current = copy.deepcopy(current)
        for c in quaternion_current.values():
            c.joint = ((0, 0, 0),) * 3
        legs = cr.make_segs(quaternion_current, rest)
        planes = cr.all_cons(legs, original[0], rp[0], (0, 0, 0), (0, -1, 0), cmds.getAttr(node + ".falloff"))
        quaternion_delta = cr.couple(
            original,
            rp,
            [planes] * len(original),
            [[j for j in range(len(original)) if j != i] for i in range(len(original))],
            cmds.getAttr(node + ".falloff"),
            [rest[key] for key in sorted(rest)],
            [1] * len(original),
        )[0]
        self.assertGreater(length(sub(joint_delta, quaternion_delta)), 1e-3)
        rotation = om.MVector(*b0).rotateTo(om.MVector(*current[(0, 1)].axis)).asMatrix()
        value = list(rotation)
        value[12:15] = poses["left"][1]
        cmds.setAttr(node + ".leftKneeMatrix", *value, type="matrix")
        _, _, _, fields, _ = node_reference(node)
        current, rest = node_cylinders(node)
        self.assertEqual(transported_side(b0, current[(0, 1)], rest[(0, 1)])[1], "quaternion")
        self.assert_collision_reference(mesh, node, original)
        rotation = list(om.MEulerRotation(0, 0, math.pi / 4).asMatrix())
        for rotating in (("left",), ("restLeft",), ("left", "restLeft")):
            rotations = {prefix: rotation if prefix in rotating else matrix() for prefix in ("left", "restLeft")}
            for prefix, unscaled in rotations.items():
                value = unscaled[:]
                value[12:15] = poses[prefix][1]
                cmds.setAttr(node + "." + prefix + "KneeMatrix", *value, type="matrix")
            baseline = self.points(mesh)
            self.assert_collision_reference(mesh, node, original)
            for prefix, unscaled in rotations.items():
                value = unscaled[:]
                for row, scale in enumerate((2, 1, 0.5)):
                    for column in range(3):
                        value[4 * row + column] *= scale
                value[12:15] = poses[prefix][1]
                cmds.setAttr(node + "." + prefix + "KneeMatrix", *value, type="matrix")
                normalized = joint_rotation([value[i : i + 3] for i in (0, 4, 8)])
                self.assert_points_close(normalized, [unscaled[i : i + 3] for i in (0, 4, 8)], 1e-12)
                self.assert_points_close(self.points(mesh), baseline, 1e-6)
                self.assert_collision_reference(mesh, node, original)

    def test_collision_degenerate_rest_side_keeps_constraints(self):
        fixture = self.simple_sweep_fixture()
        samples = [
            dict(point=point, rest_point=(0, 2, 0)) for point in ((0, 0, 2), (0.5, 0, 2), (2, 0, 2), (0, 0, -0.05))
        ]
        mesh, node = self.sweep_node(fixture, samples)
        for joint, position in zip(("Hip", "Knee", "Heel"), ((0, 2, 0), (0, 6, 0), (0, 10, 0))):
            cmds.setAttr(node + ".restLeft" + joint + "Matrix", *matrix(*position), type="matrix")
        original, _, _, fields, _ = node_reference(node)
        current, rest = node_cylinders(node)
        self.assertIsNone(rest_contact((0, 2, 0), rest[(0, 0)], (0, 0, 0), (0, 0, 1)))
        for entry in fields:
            self.assertTrue(entry)
        self.assert_collision_reference(mesh, node, original)

    def test_collision_partial_long_leg_and_independent_rest_radii(self):
        mesh = self.mesh([(0.2, 0.5, 0), (2, 0.5, 0), (2, 1.5, 0)])
        original = self.points(mesh)
        node = self.collision(mesh, long=True)
        cmds.setAttr(node + ".falloff", 0.8)
        self.set_profile(node, calfRadiusX=3, calfRadiusZ=3, ankleRadiusX=3, ankleRadiusZ=3)
        bad_heel = matrix(0, 2)
        bad_heel[0] = 0
        cmds.setAttr(node + ".leftHeelMatrix", *bad_heel, type="matrix")
        current, rest = node_cylinders(node)
        self.assertNotIn((0, 1), current)
        self.assertIn((0, 1), rest)
        retained = self.assert_collision_reference(mesh, node, original)
        cmds.setAttr(node + ".restLeftHeelMatrix", *bad_heel, type="matrix")
        removed = self.assert_collision_reference(mesh, node, original)
        self.assertGreater(max(distance(a, b) for a, b in zip(retained, removed)), 1e-6)
        cmds.setAttr(node + ".restLeftHeelMatrix", *matrix(0, 2), type="matrix")
        self.assert_points_close(self.points(mesh), retained)
        cmds.setAttr(node + ".leftHeelMatrix", *matrix(0, 2), type="matrix")
        bad_hip = matrix()
        bad_hip[0] = 0
        cmds.setAttr(node + ".leftHipMatrix", *bad_hip, type="matrix")
        current, rest = node_cylinders(node)
        self.assertEqual([segment.id for segment in cr.make_segs(current, rest)[0]], [1])
        self.assert_collision_reference(mesh, node, original)

    def test_collision_couples_neighbours_on_mesh(self):
        mesh = self.mesh([(0.2, 0.5, 0), (2, 0.5, 0), (2, 0.5, 2), (0.2, 0.5, 2)])
        original = self.points(mesh)
        node = self.collision(mesh)
        cmds.setAttr(node + ".falloff", 0.8)
        result = self.assert_collision_reference(mesh, node, original)
        uncoupled = node_reference(node, coupling=False)[2]
        self.assertEqual(uncoupled[1:], original[1:])
        self.assertGreater(distance(result[1], original[1]), 1e-5)
        self.assertGreater(distance(result[3], original[3]), 1e-5)

    def test_collision_couples_neighbours_on_nurbs_grid(self):
        for periodic in (False, True):
            with self.subTest(periodic=periodic):
                if periodic:
                    surface = cmds.cylinder(
                        axis=(0, 1, 0),
                        radius=1.3,
                        heightRatio=2,
                        sections=8,
                        spans=4,
                        degree=3,
                        constructionHistory=False,
                    )[0]
                    cmds.move(0, 0.5, 0, surface + ".cv[*][*]", relative=True, objectSpace=True)
                else:
                    surface = cmds.nurbsPlane(
                        axis=(0, 1, 0), width=4, degree=1, patchesU=2, patchesV=2, constructionHistory=False
                    )[0]
                    cmds.move(0.2, 0.5, 0, surface + ".cv[*][*]", relative=True, objectSpace=True)
                node = self.collision(surface, long=periodic)
                cmds.setAttr(node + ".falloff", 0.8)
                if periodic:
                    for joint, y in (("Hip", 0), ("Knee", 1), ("Heel", 2)):
                        cmds.setAttr(node + ".left" + joint + "Matrix", *matrix(1.1, y), type="matrix")
                original, _, expected, planes, _ = node_reference(node, "surface")
                result = [tuple(p[:3]) for p in coordinates(node, "outputGeometry[0]", "surface")]
                assert_displacements_close(self, original, result, expected)
                uncoupled = node_reference(node, "surface", coupling=False)[2]
                self.assertGreater(max(distance(a, b) for a, b in zip(result, uncoupled)), 1e-5)
                self.assertTrue(
                    any(
                        all(d <= 0 for _, d, _ in cons) and distance(p, q) > 1e-5
                        for cons, p, q in zip(planes, original, result)
                    )
                )
                if periodic:
                    snapshot = surface_snapshot(node + ".input[0].inputGeometry")
                    nu, nv = snapshot["nu"], snapshot["nv"]
                    uu = nu - snapshot["du"] if snapshot["formU"] == 3 else nu
                    uv = nv - snapshot["dv"] if snapshot["formV"] == 3 else nv
                    self.assertTrue(uu < nu or uv < nv)
                    adjacency = cr.surface_adjacency(
                        nu, nv, snapshot["du"], snapshot["dv"], snapshot["formU"] == 3, snapshot["formV"] == 3
                    )
                    wrap = (uu - 1) * nv if uu < nu else uv - 1
                    self.assertIn(wrap, adjacency[0])
                    self.assertIn(0, adjacency[wrap])
                    for u in range(nu):
                        for v in range(nv):
                            if u >= uu or v >= uv:
                                self.assert_points_close([result[u * nv + v]], [result[(u % uu) * nv + v % uv]], 1e-12)

    def test_collision_periodic_canonical_displacement_and_paint(self):
        surface = cmds.cylinder(
            axis=(0, 1, 0),
            radius=1.3,
            heightRatio=2,
            sections=8,
            spans=2,
            degree=3,
            constructionHistory=False,
        )[0]
        cmds.move(0, 0.5, 0, surface + ".cv[*][*]", relative=True, objectSpace=True)
        node = self.collision(surface, long=True)
        cmds.setAttr(node + ".falloff", 0.8)
        for joint, y in (("Hip", 0), ("Knee", 1), ("Heel", 2)):
            cmds.setAttr(node + ".left" + joint + "Matrix", *matrix(1.1, y), type="matrix")
        snapshot = surface_snapshot(node + ".input[0].inputGeometry")
        self.assertIn(om.MFnNurbsSurface.kPeriodic, (snapshot["formU"], snapshot["formV"]))
        nu, nv = snapshot["nu"], snapshot["nv"]
        uu = nu - snapshot["du"] if snapshot["formU"] == 3 else nu
        vv = nv - snapshot["dv"] if snapshot["formV"] == 3 else nv
        duplicates = {i: (i // nv % uu) * nv + i % nv % vv for i in range(nu * nv) if i // nv >= uu or i % nv >= vv}
        self.assertTrue(duplicates)
        full = None
        for paint in (1.0, 0.2, 0.0, 1.0):
            for i in set(duplicates.values()) | set(duplicates):
                cmds.setAttr(node + ".weightList[0].weights[%d]" % i, paint)
            original, _, expected, _, _ = node_reference(node, "surface")
            actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
            for i, canonical in duplicates.items():
                self.assertLessEqual(distance(actual[i], actual[canonical]), 1e-12)
                if paint == 0:
                    self.assertEqual(actual[canonical], original[canonical])
            if paint == 1:
                self.assertGreater(max(distance(actual[i], original[i]) for i in duplicates.values()), 1e-5)
                if full is not None:
                    self.assertEqual(actual, full)
                full = actual
            self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry"), snapshot)
            assert_displacements_close(self, original, actual, expected)
            for duplicate_paint in (0.0, 1.0):
                for i in duplicates:
                    cmds.setAttr(node + ".weightList[0].weights[%d]" % i, duplicate_paint)
                changed = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                self.assertEqual(changed, actual)
                expected = node_reference(node, "surface")[2]
                assert_displacements_close(self, original, changed, expected)

        for pattern in ("uniform", "varying"):
            with self.subTest(paint=pattern):
                for i in range(nu * nv):
                    paint = 0.2 if pattern == "uniform" else 0.1 + 0.1 * (i % 9)
                    cmds.setAttr(node + ".weightList[0].weights[%d]" % i, paint)
                original, _, expected, _, _ = node_reference(node, "surface")
                actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                assert_displacements_close(self, original, actual, expected)
                for duplicate_paint in (0.0, 1.0):
                    for i in duplicates:
                        cmds.setAttr(node + ".weightList[0].weights[%d]" % i, duplicate_paint)
                    changed = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                    self.assertEqual(changed, actual)
                    for i, canonical in duplicates.items():
                        self.assertLessEqual(distance(changed[i], changed[canonical]), 1e-12)
                self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry"), snapshot)

    def test_collision_rational_surface_cartesian_cvs_match_reference(self):
        directory = Path(__file__).resolve().parent / "fixtures/long_skirt_poses"
        data = json.loads((directory / "frame_1.json").read_text())
        for key in ("input", "rest"):
            surface_data = data[key]
            nu, nv = surface_data["nu"], surface_data["nv"]
            surface_data["weights"] = [0.5 + 0.5 * (((i // nv) % (nu - 1) + i % nv) % 4) for i in range(nu * nv)]
        surface = self.fixture_surface(data["input"])
        rest = self.fixture_surface(data["rest"])
        node = cmds.deformer(surface, type="yddSkirtCollideDeformer")[0]
        rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
        cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
        self.configure_fixture(node, data)
        cmds.setAttr(node + ".closedU", True)
        snapshot = surface_snapshot(node + ".input[0].inputGeometry")
        self.assertEqual(snapshot["weights"], data["input"]["weights"])
        self.assertGreater(len(set(snapshot["weights"])), 1)
        self.assert_points_close(snapshot["cvs"], data["input"]["cvs"], 1e-12)
        source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True).split(".", 1)[0]
        iterator = om.MItGeometry(om.MSelectionList().add(source).getDagPath(0))
        while not iterator.isDone():
            self.assert_points_close([tuple(iterator.position())[:3]], [snapshot["cvs"][iterator.index()]], 1e-12)
            iterator.next()
        original, _, expected, _, _ = node_reference(node, "surface")
        output = surface_snapshot(node + ".outputGeometry[0]")
        actual = output["cvs"]
        self.assertEqual(output["weights"], snapshot["weights"])
        self.assertGreater(max(distance(a, p) for a, p in zip(actual, original)), 1e-3)
        nu, nv = snapshot["nu"], snapshot["nv"]
        for v in range(nv):
            self.assertLessEqual(distance(actual[v], actual[(nu - 1) * nv + v]), 1e-12)
        assert_displacements_close(self, original, actual, expected)
        source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
        source_shape = source.split(".", 1)[0]
        source_fn = om.MFnNurbsSurface(om.MSelectionList().add(source_shape).getDagPath(0))
        first_output = actual
        for weights in (list(reversed(snapshot["weights"])), snapshot["weights"]):
            source_fn.setCVPositions(om.MPointArray([om.MPoint(*p, w) for p, w in zip(original, weights)]))
            source_fn.updateSurface()
            original, _, expected, _, _ = node_reference(node, "surface")
            output = surface_snapshot(node + ".outputGeometry[0]")
            self.assertEqual(output["weights"], weights)
            assert_displacements_close(self, original, output["cvs"], expected)
            if weights == snapshot["weights"]:
                self.assertEqual(output["cvs"], first_output)
            else:
                self.assertGreater(max(distance(a, b) for a, b in zip(output["cvs"], first_output)), 1e-5)
            self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry")["cvs"], original)
            geometry = om.MSelectionList().add(node + ".outputGeometry[0]").getPlug(0).asMObject()
            output_fn = om.MFnNurbsSurface(geometry)
            sites, rows = cr.surface_basis(output)
            for (u, v), row in zip(sites, rows):
                evaluated = tuple(output_fn.getPointAtParam(u, v, om.MSpace.kObject))[:3]
                expected_point = tuple(sum(b * output["cvs"][i][k] for i, b in row.items()) for k in range(3))
                self.assertLessEqual(distance(evaluated, expected_point), 1e-8)
        for envelope in (0.5, 0):
            cmds.setAttr(node + ".envelope", envelope)
            output = surface_snapshot(node + ".outputGeometry[0]")
            self.assertEqual(output["weights"], snapshot["weights"])
            scaled = [add(p, mul(sub(a, p), envelope)) for p, a in zip(original, first_output)]
            assert_displacements_close(self, original, output["cvs"], scaled)

    def test_collision_unit_weight_full_cv_iterator_output(self):
        fixture = Path(__file__).resolve().parent / "fixtures/long_skirt_poses/frame_1.json"
        data = json.loads(fixture.read_text())
        surface = self.fixture_surface(data["input"])
        rest = self.fixture_surface(data["rest"])
        node = cmds.deformer(surface, type="yddSkirtCollideDeformer")[0]
        rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
        cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
        self.configure_fixture(node, data)

        input_before = surface_snapshot(node + ".input[0].inputGeometry")
        self.assertEqual(input_before["weights"], [1.0] * len(input_before["cvs"]))
        source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
        source_shape = source.split(".", 1)[0]
        source_iter = om.MItGeometry(om.MSelectionList().add(source_shape).getDagPath(0))
        indices = []
        while not source_iter.isDone():
            indices.append(source_iter.index())
            source_iter.next()
        self.assertEqual(indices, list(range(len(input_before["cvs"]))))

        cmds.dgdirty(node)
        first = surface_snapshot(node + ".outputGeometry[0]")
        self.assertEqual(first["weights"], input_before["weights"])
        self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry"), input_before)
        self.assertGreater(max(sum((a - b) ** 2 for a, b in zip(p, q))
                               for p, q in zip(first["cvs"], input_before["cvs"])), 1e-10)
        cmds.dgdirty(node)
        self.assertEqual(surface_snapshot(node + ".outputGeometry[0]"), first)
        self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry"), input_before)

        previous_mode = cmds.evaluationManager(query=True, mode=True)[0]
        try:
            cmds.evaluationManager(mode="parallel")
            cmds.dgdirty(node)
            parallel = surface_snapshot(node + ".outputGeometry[0]")
            self.assertEqual(parallel["weights"], first["weights"])
            self.assert_points_close(parallel["cvs"], first["cvs"], tolerance=1e-9)
            self.assertEqual(surface_snapshot(node + ".input[0].inputGeometry"), input_before)
        finally:
            cmds.evaluationManager(mode=previous_mode)

    def test_collision_surface_identical_input_output_reuse(self):
        data = dict(
            nu=4,
            nv=3,
            du=2,
            dv=1,
            ku=[0, 0, 0.3, 1, 1],
            kv=[0, 0.6, 1],
            formU=1,
            formV=1,
            cvs=[(0.2 + 0.12 * u, 0.2 + 0.2 * v, 0.15) for u in range(4) for v in range(3)],
            weights=[0.5 + 0.5 * (i % 4) for i in range(12)],
        )
        surface = self.fixture_surface(data)
        components = [surface + ".cv[%d][%d]" % divmod(i, data["nv"]) for i in range(len(data["cvs"]))]
        node = self.collision(components, long=True)
        cmds.setAttr(node + ".falloff", 0.8)
        for i in range(len(components)):
            cmds.setAttr(node + ".weightList[0].weights[%d]" % i, 1)
        tag_expression = cmds.getAttr(node + ".input[0].componentTagExpression")
        reduced_tag = cmds.componentTag([c for i, c in enumerate(components) if i != 4], newTagName="reduced", create=True)
        source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
        source = source.split(".", 1)[0]
        source_fn = om.MFnNurbsSurface(om.MSelectionList().add(source).getDagPath(0))
        original = surface_snapshot(node + ".input[0].inputGeometry")

        def evaluate():
            cmds.dgdirty(node)
            return coordinates(node, "outputGeometry[0]", "surface")

        baseline = evaluate()
        self.assertEqual(evaluate(), baseline)
        self.assertEqual([p[3] for p in baseline], original["weights"])
        self.assertGreater(max(distance(a[:3], b) for a, b in zip(baseline, original["cvs"])), 1e-5)
        self.assertGreater(distance(baseline[4][:3], original["cvs"][4]), 1e-5)

        def check_change(change, restore):
            try:
                change()
                changed = evaluate()
                self.assertGreater(max(distance(a[:3], b[:3]) for a, b in zip(changed, baseline)), 1e-6)
                self.assertEqual(evaluate(), changed)
            finally:
                restore()
            self.assertEqual(evaluate(), baseline)
            self.assertEqual(evaluate(), baseline)

        for attribute, value in (
            ("weightList[0].weights[4]", 0),
            ("falloff", 0.2),
            ("closedU", True),
            ("envelope", 0.5),
            ("thighRadiusX", 1.6),
            ("calfRadiusZ", 1.4),
            ("skirtType", 0),
        ):
            with self.subTest(input=attribute):
                plug = node + "." + attribute
                saved = cmds.getAttr(plug)
                check_change(lambda: cmds.setAttr(plug, value), lambda: cmds.setAttr(plug, saved))

        for attribute, value, kind in (
            ("restLeftHipMatrix", matrix(0.35, 0, 0.1), "matrix"),
            ("leftKneeMatrix", matrix(0.3, 1, 0.2), "matrix"),
            ("restBellMatrix", matrix(0.1, 0.4, 0.2), "matrix"),
            ("ringScale", (1.3, 1, 0.8), "double3"),
        ):
            with self.subTest(input=attribute):
                plug = node + "." + attribute
                saved = cmds.getAttr(plug)
                if kind == "double3":
                    saved = saved[0]
                check_change(
                    lambda: cmds.setAttr(plug, *value, type=kind),
                    lambda: cmds.setAttr(plug, *saved, type=kind),
                )

        def set_input(cvs, weights):
            source_fn.setCVPositions(om.MPointArray([om.MPoint(*p, w) for p, w in zip(cvs, weights)]))
            source_fn.updateSurface()

        changed_cvs = list(original["cvs"])
        changed_cvs[4] = add(changed_cvs[4], (0.12, 0, 0.1))
        changed_weights = list(original["weights"])
        changed_weights[4] = 1.75
        for label, cvs, weights in (
            ("input CV", changed_cvs, original["weights"]),
            ("rational weight", original["cvs"], changed_weights),
        ):
            with self.subTest(input=label):
                check_change(
                    lambda: set_input(cvs, weights),
                    lambda: set_input(original["cvs"], original["weights"]),
                )

        with self.subTest(input="membership"):
            expression = node + ".input[0].componentTagExpression"
            check_change(
                lambda: cmds.setAttr(expression, reduced_tag, type="string"),
                lambda: cmds.setAttr(expression, tag_expression, type="string"),
            )

    def test_collision_surface_topology_cache_invalidation(self):
        base = dict(nu=4, nv=3, du=2, dv=1, ku=[0, 0, 0.3, 1, 1], kv=[0, 0.6, 1], formU=1, formV=1)
        variants = [
            base,
            dict(base, ku=[0, 0, 0.7, 1, 1]),
            dict(base, kv=[0, 0.25, 1]),
            dict(base, du=1, ku=[0, 0.2, 0.6, 1]),
            dict(base, dv=2, kv=[0, 0, 1, 1]),
            dict(base, nu=5, ku=[0, 0, 0.3, 0.7, 1, 1]),
            dict(base, nv=4, kv=[0, 0.2, 0.6, 1]),
            dict(base, nu=6, ku=[-1, 0, 1, 2, 3, 4, 5], formU=3),
            dict(base, nv=4, kv=[-1, 0, 1, 2], formV=3),
        ]
        shapes = []
        for data in variants:
            nu, nv = data["nu"], data["nv"]
            uu = nu - data["du"] if data["formU"] == 3 else nu
            vv = nv - data["dv"] if data["formV"] == 3 else nv
            data["cvs"] = [(0.2 + 0.12 * (u % uu), 0.2 + 0.2 * (v % vv), 0.15) for u in range(nu) for v in range(nv)]
            surface = self.fixture_surface(data)
            shapes.append(cmds.listRelatives(surface, shapes=True, noIntermediate=True, fullPath=True)[0])
        driven = self.fixture_surface(base)
        node = self.collision(driven)
        cmds.setAttr(node + ".falloff", 0.8)
        source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
        source_shape = source.split(".", 1)[0]
        observations = []
        for index, closed in [(i, False) for i in range(len(shapes))] + [(0, True), (0, False), (0, True)]:
            with self.subTest(topology=index, closedU=closed):
                cmds.connectAttr(shapes[index] + ".local", source_shape + ".create", force=True)
                cmds.connectAttr(shapes[index] + ".local", node + ".restGeometry", force=True)
                cmds.setAttr(node + ".closedU", closed)
                original, _, expected, _, _ = node_reference(node, "surface")
                actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                key = index, closed
                observations.append((key, original, actual, expected))
        assert_topology_observations(self, observations, assert_displacements_close)

    def fixture_surface(self, data):
        fn = om.MFnNurbsSurface()
        cvs = [
            om.MPoint(*(list(p[:3]) + [data.get("weights", [1] * len(data["cvs"]))[i]]))
            for i, p in enumerate(data["cvs"])
        ]
        obj = fn.create(
            om.MPointArray(cvs),
            data["ku"],
            data["kv"],
            data["du"],
            data["dv"],
            data["formU"],
            data["formV"],
            any(cv.w != 1.0 for cv in cvs),
        )
        path = om.MFnDagNode(obj).fullPathName()
        return path if cmds.nodeType(path) == "transform" else cmds.listRelatives(path, parent=True, fullPath=True)[0]

    def configure_fixture(self, node, data):
        for name, value in data["matrices"].items():
            cmds.setAttr(node + "." + name, *value, type="matrix")
        cmds.setAttr(node + ".ringScale", *data["ringScale"], type="double3")
        for name in (
            "thighPosition",
            "calfPosition",
            "leftRingAxis",
            "rightRingAxis",
            "skirtType",
            "falloff",
            "envelope",
        ):
            cmds.setAttr(node + "." + name, data[name])
        names = [station + "Radius" + axis for station in ("thigh", "knee", "calf", "ankle") for axis in ("X", "Z")]
        for name, value in zip(names, data["radii"]):
            cmds.setAttr(node + "." + name, value)

    def test_collision_sample_synced_fixtures_and_open_seams(self):
        directory = Path(__file__).resolve().parent / "fixtures/long_skirt_poses"
        fixtures = {f: json.loads((directory / ("frame_%d.json" % f)).read_text()) for f in (0, 1, 2, 10)}
        for closed in (True, False):
            data = fixtures[0]
            surface = self.fixture_surface(data["input"])
            rest = self.fixture_surface(data["rest"])
            node = cmds.deformer(surface, type="yddSkirtCollideDeformer")[0]
            self.assertFalse(cmds.getAttr(node + ".closedU"))
            if closed:
                cmds.setAttr(node + ".closedU", True)
            rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
            cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
            source = cmds.connectionInfo(node + ".input[0].inputGeometry", sourceFromDestination=True)
            source_shape = source.split(".", 1)[0]
            fn = om.MFnNurbsSurface(om.MSelectionList().add(source_shape).getDagPath(0))
            saved = {}
            for frame in (0, 1, 10, 2, 1, 0, 10):
                cmds.currentTime(frame)
                data = fixtures[frame]
                fn.setCVPositions(om.MPointArray([om.MPoint(*p) for p in data["input"]["cvs"]]))
                fn.updateSurface()
                self.configure_fixture(node, data)
                original, _, expected, _, _ = node_reference(node, "surface")
                full = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                nv, nu = data["input"]["nv"], data["input"]["nu"]
                displacement = [sub(a, p) for a, p in zip(full, original)]
                seam = max(length(sub(displacement[v], displacement[(nu - 1) * nv + v])) for v in range(nv))
                if closed:
                    self.assertLessEqual(seam, 1e-12)
                elif frame != 0:
                    self.assertGreater(seam, 1e-5)
                if frame == 0:
                    self.assertEqual(full, original)
                assert_displacements_close(self, original, full, expected)
                if frame in saved:
                    assert_displacements_close(self, original, full, saved[frame])
                saved[frame] = full
                print(
                    "sample node",
                    frame,
                    "closed",
                    closed,
                    "waist",
                    [max(length(displacement[u * nv + v]) for u in range(nu)) for v in (0, 1)],
                    "max",
                    max(map(length, displacement)),
                    "cv[12][1]",
                    displacement[12 * nv + 1],
                )
                for envelope in (0, 0.5, 1):
                    cmds.setAttr(node + ".envelope", envelope)
                    actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                    scaled = [add(p, mul(d, envelope)) for p, d in zip(original, displacement)]
                    assert_displacements_close(self, original, actual, scaled)

    def test_collide_golden(self):
        directory = Path(__file__).resolve().parent / "fixtures"
        golden_path = directory / "collide_golden_maya2026.json"
        golden = json.loads(golden_path.read_text())
        self.assertEqual(golden["schema"], "yddColliders.collideGolden.gauss3.v1")
        input_dir = Path(__file__).resolve().parent / "fixtures/long_skirt_poses"
        fixtures = {
            frame: json.loads((input_dir / "frame_{}.json".format(frame)).read_text())
            for frame in (0, 1, 2, 10)
        }
        expected_keys = {(frame, closed) for frame in (0, 1, 2, 10) for closed in (False, True)}
        cases_by_key = {(case["frame"], case["closedU"]): case for case in golden["cases"]}
        self.assertEqual(len(golden["cases"]), 8)
        self.assertEqual(set(cases_by_key), expected_keys)
        for case in golden["cases"]:
            frame, closed = case["frame"], case["closedU"]
            data = fixtures[frame]
            stored = case["cvs"]
            self.assertEqual(case["cvCount"], len(stored))
            self.assertEqual(case["cvCount"], len(data["input"]["cvs"]))
            self.assertTrue(all(math.isfinite(value) for point in stored for value in point))
            stored_hash = hashlib.sha256(
                b"".join(struct.pack("<3d", *point) for point in stored)
            ).hexdigest()
            self.assertEqual(stored_hash, case["sha256LittleEndian3d"])
            input_path = input_dir / "frame_{}.json".format(frame)
            self.assertEqual(hashlib.sha256(input_path.read_bytes()).hexdigest(), case["inputSha256"])
            with self.subTest(frame=frame, closedU=closed):
                cmds.file(new=True, force=True)
                surface = self.fixture_surface(data["input"])
                rest = self.fixture_surface(data["rest"])
                node = cmds.deformer(surface, type="yddSkirtCollideDeformer")[0]
                cmds.setAttr(node + ".closedU", closed)
                rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
                cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
                self.configure_fixture(node, data)
                cmds.currentTime(frame)
                original, _, expected, _, _ = node_reference(node, "surface")
                actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                self.assertEqual(len(actual), case["cvCount"])
                self.assertTrue(all(math.isfinite(value) for point in actual for value in point))
                assert_displacements_close(self, original, actual, expected)
                if cmds.about(version=True).startswith("2026"):
                    actual_hash = hashlib.sha256(
                        b"".join(struct.pack("<3d", *point) for point in actual)
                    ).hexdigest()
                    self.assertEqual(actual_hash, case["sha256LittleEndian3d"])
                nv, nu = data["input"]["nv"], data["input"]["nu"]
                displacement = [sub(a, p) for a, p in zip(actual, original)]
                seam = max(length(sub(displacement[v], displacement[(nu - 1) * nv + v])) for v in range(nv))
                if closed:
                    self.assertLessEqual(seam, 1e-12)
                elif frame != 0:
                    self.assertGreater(seam, 1e-5)

    def test_collision_sample_partial_membership_and_paint(self):
        directory = Path(__file__).resolve().parent / "fixtures/long_skirt_poses"
        data = json.loads((directory / "frame_1.json").read_text())
        for partial in (False, True):
            surface = self.fixture_surface(data["input"])
            rest = self.fixture_surface(data["rest"])
            nv, nu = data["input"]["nv"], data["input"]["nu"]
            members = set(range((nu - 1) * nv)) if partial else set(range(nu * nv))
            components = [surface + ".cv[%d][%d]" % divmod(i, nv) for i in sorted(members)]
            node = self.collision(components)
            cmds.setAttr(node + ".closedU", True)
            rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
            cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
            self.configure_fixture(node, data)
            cmds.setAttr(node + ".weightList[0].weights[1]", 0.2)
            cmds.setAttr(node + ".weightList[0].weights[%d]" % ((nu - 1) * nv + 1), 0.8)
            original, _, expected, _, _ = node_reference(node, "surface", members=members)
            actual = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
            assert_displacements_close(self, original, actual, expected)
            if partial:
                for v in range(nv):
                    self.assertEqual(actual[v], original[v])
                    self.assertEqual(actual[(nu - 1) * nv + v], original[(nu - 1) * nv + v])
            else:
                displacement = [sub(a, p) for a, p in zip(actual, original)]
                self.assertLess(length(sub(displacement[1], displacement[(nu - 1) * nv + 1])), 1e-12)
                cmds.setAttr(node + ".weightList[0].weights[%d]" % ((nu - 1) * nv + 1), 0.2)
                assert_displacements_close(self, original, surface_snapshot(node + ".outputGeometry[0]")["cvs"], actual)

    def test_collision_zero_weight_points_stay_fixed(self):
        mesh = self.mesh([(0.2, 0.5, 0), (0.4, 0.5, 0), (2, 0.5, 0), (2, 0.5, 2)])
        original = self.points(mesh)
        node = self.collision(mesh)
        cmds.setAttr(node + ".falloff", 0.8)
        cmds.setAttr(node + ".weightList[0].weights[0]", 0)
        fixed = self.assert_collision_reference(mesh, node, original)
        self.assertEqual(self.points(mesh)[0], original[0])
        self.assertGreater(distance(fixed[1], original[1]), 0.1)
        quotients = []
        for q in (1e-3, 5e-4, 2.5e-4):
            cmds.setAttr(node + ".weightList[0].weights[0]", q)
            result = self.assert_collision_reference(mesh, node, original)
            difference = max(distance(a, b) for a, b in zip(result, fixed))
            quotients.append(difference / q)
            self.assertLess(difference, 2 * q)
        print("paint approaching zero quotients={}".format(quotients))

    def test_collision_long_leg_knee_is_not_an_exit(self):
        mesh = self.mesh([(0.5, y, 0) for y in (1 - 1e-5, 1, 1 + 1e-5)])
        original = self.points(mesh)
        node = self.collision(mesh, long=True)
        result = self.assert_collision_reference(mesh, node, original)
        for p, q in zip(original, result):
            self.assertAlmostEqual(q[0], 1, delta=1e-6)
            self.assertAlmostEqual(q[1], p[1], delta=1e-6)
        current, rest = node_cylinders(node)
        leg = cr.make_segs(current, rest)[0]
        d, who = cr.union_exit(leg, (0.5, 0.9, 0), (0, 1, 0))
        self.assertEqual(who.id, 1)
        self.assertAlmostEqual(d, 1.1, delta=1e-10)

    def test_wave_zero_paths_and_default_impulse(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        cmds.setAttr(node + ".impulseX", -0.5)
        cmds.setAttr(node + ".impulsePosition", 0.5)
        for attr in ["amplitude", "envelope"]:
            previous = cmds.getAttr(node + "." + attr)
            cmds.setAttr(node + "." + attr, 0)
            self.assertEqual(self.points(mesh), original)
            cmds.setAttr(node + "." + attr, previous)
        cmds.setAttr(node + ".idleAmplitude", 0)
        self.assertNotEqual(self.points(mesh), original)
        cmds.setAttr(node + ".impulseX", 0)
        self.assertEqual(self.points(mesh), original)
        cmds.setAttr(node + ".noiseAmplitude", 0.000001)
        self.assertEqual(self.points(mesh), original)

    def test_wave_layers_and_zero_weight_hem_height(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 2, 1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        for complexity in [0, 1, 0.4]:
            with self.subTest(complexity=complexity):
                cmds.setAttr(node + ".idleComplexity", complexity)
                cmds.setAttr(node + ".weightList[0].weights[2]", 1)
                full = self.points(mesh)
                cmds.setAttr(node + ".weightList[0].weights[2]", 0)
                result = self.points(mesh)
                self.assertEqual(result[:2], full[:2])
                self.assertEqual(result[2], original[2])
                golden = (1 + math.sqrt(5)) * 0.5

                def shaped(phase):
                    sine = math.sin(phase + 0.3 * math.sin(phase))
                    return math.copysign(abs(sine) ** 1.9, sine)

                for p, actual in zip(original[:2], result[:2]):
                    x, y, z = p
                    v, theta = y / 2, math.atan2(z, x)
                    primary = 0.5 * (shaped(2 * math.pi * (0.23 - 1.3 * v)) + shaped(2 * (theta - 2 * math.pi * 0.17)))
                    secondary = 0.5 * (
                        shaped(2 * math.pi * (golden * 0.23 - 1.3 * v) + 0.5 * math.pi)
                        + shaped(2 * (theta - 2 * math.pi * golden * 0.17) + 0.5 * math.pi)
                    )
                    displacement = 0.7 * 0.3 * ((1 - complexity) * primary + complexity * secondary)
                    radius = math.hypot(x, z)
                    self.assertAlmostEqual(actual[0], x + x / radius * displacement, places=6)
                    self.assertAlmostEqual(actual[2], z + z / radius * displacement, places=6)

    def test_partial_membership_and_reverse_inputs(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1), (0.8, 0.7, 0.3)])
        original = self.points(mesh)
        node = self.wave([mesh + ".vtx[0]", mesh + ".vtx[2:3]"])
        results = {}
        for phase in [0.1, 0.7, 1.2, 0.7, 0.1]:
            cmds.setAttr(node + ".wavePhaseV", phase)
            result = self.points(mesh)
            self.assertEqual(result[1], original[1])
            self.assertNotEqual(result[0], original[0])
            self.assertNotEqual(result[2], original[2])
            self.assertNotEqual(result[3], original[3])
            self.assertTrue(all(math.isfinite(value) for p in result for value in p))
            if phase in results:
                self.assertEqual(result, results[phase])
            results[phase] = result
        cmds.setAttr(node + ".bellMatrix", *([0] * 16), type="matrix")
        self.assertEqual(self.points(mesh), original)

    def assert_points_close(self, actual, expected, tolerance=1e-6):
        self.assertEqual(len(actual), len(expected))
        for point, reference in zip(actual, expected):
            for value, target in zip(point, reference):
                self.assertAlmostEqual(value, target, delta=tolerance)

    def assert_values_close(self, actual, expected, tolerance=1e-6):
        self.assertEqual(len(actual), len(expected))
        for value, target in zip(actual, expected):
            self.assertAlmostEqual(value, target, delta=tolerance)

    def set_wave_values(self, node, **values):
        for name, value in values.items():
            cmds.setAttr(node + "." + name, value)

    def test_wave_expression_attribute_contract(self):
        node = cmds.createNode("yddSkirtWaveDeformer")
        defaults = {
            "idleAmplitudeV": 1,
            "idleAmplitudeU": 1,
            "idleDirectionality": 0,
            "idleDirectionX": -1,
            "idleDirectionZ": 0,
            "idleDirectionSpace": 0,
            "phaseSpread": 0,
            "impulseAmount": 1,
        }
        for name, value in defaults.items():
            with self.subTest(attribute=name):
                self.assertEqual(cmds.getAttr(node + "." + name), value)
                self.assertEqual(cmds.getAttr(node + "." + name, keyable=True), name != "idleDirectionSpace")
        self.assertTrue(cmds.getAttr(node + ".idleDirectionSpace", channelBox=True))
        self.assertEqual(cmds.attributeQuery("idleDirectionSpace", node=node, listEnum=True), ["World:Bell Local"])
        for name in ["idleDirectionX", "idleDirectionZ"]:
            self.assertFalse(cmds.attributeQuery(name, node=node, minExists=True))
            self.assertFalse(cmds.attributeQuery(name, node=node, maxExists=True))
            self.assertEqual(cmds.attributeQuery(name, node=node, softRange=True), [-1, 1])
        for name, maximum in [
            ("idleAmplitudeV", 2),
            ("idleAmplitudeU", 2),
            ("phaseSpread", 0.5),
            ("impulseAmount", 2),
        ]:
            self.assertEqual(cmds.attributeQuery(name, node=node, minimum=True), [0])
            self.assertFalse(cmds.attributeQuery(name, node=node, maxExists=True))
            self.assertEqual(cmds.attributeQuery(name, node=node, softMax=True), [maximum])

    def test_wave_direction_symmetry_normalization_and_zero(self):
        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        self.set_wave_values(node, idleAmplitudeU=0, idleDirectionality=1, idleComplexity=0, waveCountV=0)
        result = self.points(mesh)
        self.assertNotEqual(result[0], original[0])
        self.assertAlmostEqual(result[0][0] - 1, result[1][0] + 1, delta=1e-6)
        self.assertEqual(result[2:], original[2:])
        cmds.setAttr(node + ".idleDirectionX", -17)
        self.assertEqual(self.points(mesh), result)
        for direction in [0, -0.000001]:
            cmds.setAttr(node + ".idleDirectionX", direction)
            self.assertEqual(self.points(mesh), original)
        self.set_wave_values(node, idleDirectionX=-1, idleDirectionality=0)
        breathing = self.points(mesh)
        self.assertAlmostEqual(breathing[0][0] - 1, -(breathing[1][0] + 1), delta=1e-6)
        self.set_wave_values(node, idleDirectionality=0.5)
        self.assert_points_close(
            self.points(mesh), [tuple((a + b) / 2 for a, b in zip(p, q)) for p, q in zip(result, breathing)]
        )

    def test_wave_direction_spaces_rotation_projection_and_scale(self):
        local = [(1, 0.3, 0), (-1, 0.3, 0), (0, 1, 1), (0, 1, -1)]
        for angle in [0, math.pi / 4, math.pi / 2]:
            rotation = om.MEulerRotation(0, 0, angle).asMatrix()
            for scale in [0.5, 1, 2]:
                with self.subTest(angle=angle, scale=scale):
                    bell = om.MMatrix(
                        [entry * scale if index < 12 else entry for index, entry in enumerate(list(rotation))]
                    )
                    mesh = self.mesh([tuple(om.MPoint(*point) * bell)[:3] for point in local])
                    original = self.points(mesh)
                    node = self.wave(mesh)
                    self.set_wave_values(node, idleAmplitudeU=0, idleDirectionality=1, idleComplexity=0, waveCountV=0)
                    cmds.setAttr(node + ".bellMatrix", *list(bell), type="matrix")
                    cmds.setAttr(node + ".idleDirectionSpace", 1)
                    local_result = self.points(mesh)
                    cmds.setAttr(node + ".idleDirectionSpace", 0)
                    world_result = self.points(mesh)
                    local_delta = distance(local_result[0], original[0])
                    world_delta = distance(world_result[0], original[0])
                    self.assertGreater(local_delta, 0.01 * scale)
                    self.assertAlmostEqual(world_delta, local_delta * abs(math.cos(angle)), delta=1e-6)
                    if angle == 0:
                        self.assert_points_close(world_result, local_result)
                    if angle == math.pi / 2:
                        self.assert_points_close(world_result, original)
                    if scale == 0.5:
                        reference_delta = local_delta / scale
                    else:
                        self.assertAlmostEqual(local_delta / scale, reference_delta, delta=1e-6)

    def rotation_matrix(self, x=0, y=0, z=0):
        return list(om.MEulerRotation(x, y, z).asMatrix())

    def set_rotation(self, node, values):
        cmds.setAttr(node + ".evaluationToWorldRotation", *values, type="matrix")

    def directional_wave(self, mesh):
        node = self.wave(mesh)
        self.set_wave_values(node, idleAmplitudeU=0, idleDirectionality=1, idleComplexity=0, waveCountV=0)
        return node

    def test_wave_rotation_attribute_contract(self):
        node = cmds.createNode("yddSkirtWaveDeformer")
        name = "evaluationToWorldRotation"
        self.assertEqual(cmds.attributeQuery(name, node=node, shortName=True), "etwr")
        self.assertEqual(cmds.getAttr(node + "." + name), matrix())
        self.assertTrue(cmds.attributeQuery(name, node=node, hidden=True))
        self.assertTrue(cmds.attributeQuery(name, node=node, storable=True))
        self.assertTrue(cmds.attributeQuery(name, node=node, connectable=True))
        self.assertFalse(cmds.attributeQuery(name, node=node, keyable=True))
        source = cmds.createNode("transform")
        cmds.connectAttr(source + ".matrix", node + "." + name)
        self.assertEqual(cmds.listConnections(node + "." + name, source=True, destination=False), [source])
        cmds.setAttr(source + ".rotateY", 30)
        self.assertNotEqual(cmds.getAttr(node + "." + name), matrix())

    def test_wave_rotation_identity_default_and_rotated_direction_equivalence(self):
        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1), (0.6, 0.5, -0.8)])
        original = self.points(mesh)
        node = self.directional_wave(mesh)
        baseline = self.points(mesh)
        self.assertNotEqual(baseline, original)
        self.set_rotation(node, matrix())
        self.assertEqual(self.points(mesh), baseline)
        for angle in [math.radians(30), math.pi / 2, math.radians(137)]:
            with self.subTest(angle=angle):
                rotation = om.MEulerRotation(0, angle, 0).asMatrix()
                self.set_wave_values(node, idleDirectionX=-1, idleDirectionZ=0)
                self.set_rotation(node, list(rotation))
                rotated = self.points(mesh)
                self.assertNotEqual(rotated, baseline)
                expected_direction = om.MVector(-1, 0, 0) * rotation.inverse()
                self.set_rotation(node, matrix())
                self.set_wave_values(node, idleDirectionX=expected_direction.x, idleDirectionZ=expected_direction.z)
                self.assert_points_close(rotated, self.points(mesh), tolerance=1e-9)
        # impulse keeps its input strength through the rotation
        self.set_wave_values(
            node,
            idleAmplitude=0,
            impulseAmount=1,
            impulseX=-0.5,
            impulseZ=0.25,
            impulsePosition=0.5,
            directionality=0.7,
        )
        self.set_rotation(node, matrix())
        impulse_baseline = self.points(mesh)
        self.assertNotEqual(impulse_baseline, original)
        rotation = om.MEulerRotation(0, math.radians(75), 0).asMatrix()
        self.set_rotation(node, list(rotation))
        rotated = self.points(mesh)
        self.assertNotEqual(rotated, impulse_baseline)
        expected_impulse = om.MVector(-0.5, 0, 0.25) * rotation.inverse()
        self.set_rotation(node, matrix())
        self.set_wave_values(node, impulseX=expected_impulse.x, impulseZ=expected_impulse.z)
        self.assert_points_close(rotated, self.points(mesh), tolerance=1e-9)

    def test_wave_rotation_tilt_keeps_y_and_projects_without_renormalizing(self):
        mesh = self.mesh([(0, 1, 1), (0, 1, -1), (1, 1, 0), (-1, 1, 0)])
        original = self.points(mesh)
        node = self.directional_wave(mesh)
        self.set_wave_values(node, idleDirectionX=0, idleDirectionZ=1)
        reference = self.points(mesh)
        reference_delta = distance(reference[0], original[0])
        self.assertGreater(reference_delta, 0.01)
        for angle in [math.radians(45), math.radians(60), math.pi / 2]:
            with self.subTest(angle=angle):
                self.set_rotation(node, self.rotation_matrix(x=angle))
                tilted = self.points(mesh)
                self.assertAlmostEqual(
                    distance(tilted[0], original[0]), reference_delta * abs(math.cos(angle)), delta=1e-6
                )
                self.assertEqual(tilted[2:], original[2:])
        self.assert_points_close(self.points(mesh), original)

    def test_wave_rotation_bell_local_ignores_rotation(self):
        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1)])
        node = self.directional_wave(mesh)
        self.set_wave_values(
            node, idleDirectionSpace=1, impulseSpace=1, impulseAmount=1, impulseX=-0.5, impulsePosition=0.5
        )
        local = self.points(mesh)
        for values in [
            self.rotation_matrix(y=1.0),
            [2 if i in (0, 5, 10) else 0 for i in range(16)],
            [float("nan")] * 16,
        ]:
            self.set_rotation(node, values)
            self.assertEqual(self.points(mesh), local)
        self.set_wave_values(node, idleDirectionSpace=0)
        self.assertNotEqual(self.points(mesh), local)

    def test_wave_rotation_invalid_passthrough_and_recovery(self):
        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1)])
        original = self.points(mesh)
        node = self.directional_wave(mesh)
        baseline = self.points(mesh)
        self.assertNotEqual(baseline, original)
        uniform = [2 if i in (0, 5, 10) else 0 for i in range(16)]
        uniform[15] = 1
        shear = matrix()
        shear[4] = 0.01
        reflection = matrix()
        reflection[0] = -1
        translated = matrix(x=1e-3)
        perspective = matrix()
        perspective[3] = 1e-3
        homogeneous = matrix()
        homogeneous[15] = 1 + 1e-6
        nonfinite = matrix()
        nonfinite[5] = float("nan")
        infinite = matrix()
        infinite[12] = float("inf")
        slightly_scaled = matrix()
        slightly_scaled[0] = 1 + 5e-6
        invalid = {
            "uniform_scale": uniform,
            "shear": shear,
            "reflection": reflection,
            "translation": translated,
            "perspective": perspective,
            "homogeneous": homogeneous,
            "nan": nonfinite,
            "inf": infinite,
            "slightly_scaled": slightly_scaled,
        }
        for label, values in invalid.items():
            with self.subTest(matrix=label):
                self.set_rotation(node, values)
                self.assertEqual(self.points(mesh), original)
                self.set_rotation(node, matrix())
                self.assertEqual(self.points(mesh), baseline)
        within = self.rotation_matrix(y=0.3)
        within[0] += 2e-7
        within[12] = 5e-9
        self.set_rotation(node, within)
        self.assertNotEqual(self.points(mesh), original)
        self.assertNotEqual(self.points(mesh), baseline)
        self.set_rotation(node, matrix(x=1e-7))
        self.assertEqual(self.points(mesh), original)

    def test_wave_rotation_survives_save_and_reload(self):
        import os
        import tempfile

        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1)])
        node = self.directional_wave(mesh)
        rotation = self.rotation_matrix(y=math.radians(40))
        self.set_rotation(node, rotation)
        expected = self.points(mesh)
        source = cmds.createNode("transform", name="rotationSource")
        cmds.setAttr(source + ".rotateY", 40)
        connected = cmds.createNode("yddSkirtWaveDeformer", name="connectedWave")
        cmds.connectAttr(source + ".matrix", connected + ".evaluationToWorldRotation")
        directory = tempfile.mkdtemp()
        try:
            for extension, file_type in [("ma", "mayaAscii"), ("mb", "mayaBinary")]:
                with self.subTest(extension=extension):
                    path = os.path.join(directory, "rotation." + extension).replace(os.sep, "/")
                    cmds.file(rename=path)
                    cmds.file(save=True, force=True, type=file_type)
                    cmds.file(new=True, force=True)
                    cmds.file(path, open=True, force=True)
                    self.assert_values_close(
                        cmds.getAttr(node + ".evaluationToWorldRotation"), rotation, tolerance=1e-12
                    )
                    self.assert_points_close(self.points(mesh), expected, tolerance=1e-12)
                    self.assertEqual(
                        cmds.listConnections(connected + ".evaluationToWorldRotation", source=True, destination=False),
                        [source],
                    )
            cmds.file(new=True, force=True)
            fresh = cmds.createNode("yddSkirtWaveDeformer")
            self.assertEqual(cmds.getAttr(fresh + ".evaluationToWorldRotation"), matrix())
        finally:
            for name in os.listdir(directory):
                os.remove(os.path.join(directory, name))
            os.rmdir(directory)

    def test_wave_independent_gains_and_phases(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        self.set_wave_values(node, idleDirectionality=1, phaseSpread=0.3)
        combined = self.points(mesh)
        cmds.setAttr(node + ".idleAmplitudeU", 0)
        vertical = self.points(mesh)
        cmds.setAttr(node + ".wavePhaseU", 0.9)
        self.assertEqual(self.points(mesh), vertical)
        self.set_wave_values(node, idleAmplitudeU=1, idleAmplitudeV=0, wavePhaseU=0.17)
        around = self.points(mesh)
        self.set_wave_values(node, wavePhaseV=0.89, idleDirectionX=0.2, idleDirectionZ=0.9, phaseSpread=0.5)
        self.assertEqual(self.points(mesh), around)
        self.assert_points_close(
            combined,
            [tuple(v + u - p for v, u, p in zip(a, b, c)) for a, b, c in zip(vertical, around, original)],
        )
        cmds.setAttr(node + ".idleAmplitudeU", 0)
        self.assertEqual(self.points(mesh), original)
        self.set_wave_values(node, idleAmplitudeV=1, waveCountU=0)
        vertical_only = self.points(mesh)
        self.set_wave_values(node, idleAmplitudeU=9, wavePhaseU=-7.8)
        self.assertEqual(self.points(mesh), vertical_only)

    def test_wave_phase_spread_columns_seam_loop_and_noise_independence(self):
        epsilon = 1e-7
        coordinates = [
            (math.cos(theta), height, math.sin(theta))
            for theta in [0.2, 1.1, -math.pi + epsilon, math.pi - epsilon]
            for height in [0.25, 1]
        ]
        mesh = self.mesh(coordinates)
        original = self.points(mesh)
        node = self.wave(mesh)
        self.set_wave_values(node, idleAmplitudeU=0, waveCountV=0, phaseSpread=0, idleComplexity=0)
        unspread = self.points(mesh)
        self.set_wave_values(node, phaseSpread=0.4, noisePhase=0.37)
        spread = self.points(mesh)
        self.assertNotEqual(spread, unspread)
        for index in range(0, len(original), 2):
            self.assertAlmostEqual(
                distance(spread[index], original[index]),
                distance(spread[index + 1], original[index + 1]),
                delta=1e-6,
            )
        self.assertAlmostEqual(distance(spread[4], original[4]), distance(spread[6], original[6]), delta=1e-6)
        cmds.setAttr(node + ".noiseFrequencyV", 9)
        self.assertEqual(self.points(mesh), spread)
        cmds.setAttr(node + ".wavePhaseV", 1.23)
        self.assert_points_close(self.points(mesh), spread)
        self.set_wave_values(node, wavePhaseV=0.23, idleComplexity=1)
        golden_spread = self.points(mesh)
        cmds.setAttr(node + ".phaseSpread", 0)
        self.assertNotEqual(self.points(mesh), golden_spread)
        self.set_wave_values(node, phaseSpread=0.4, noiseFrequencyU=0, idleComplexity=0)
        uniform = self.points(mesh)
        for point, reference in zip(uniform, original):
            self.assertAlmostEqual(distance(point, reference), distance(uniform[0], original[0]), delta=1e-6)

    def test_wave_send_amount_and_expression_zero_paths(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        self.set_wave_values(node, idleAmplitude=0, impulseAmount=0, impulseX=-0.5, impulsePosition=0.5)
        self.assertEqual(self.points(mesh), original)
        for direction, position in [(15, 0.2), (-3, 1), (1, 0.5)]:
            self.set_wave_values(node, impulseX=direction, impulseZ=direction, impulsePosition=position)
            self.assertEqual(self.points(mesh), original)
        cmds.setAttr(node + ".impulseAmount", 1)
        single = self.points(mesh)
        self.assertNotEqual(single, original)
        cmds.setAttr(node + ".impulseAmount", 2)
        doubled = self.points(mesh)
        self.assert_points_close(
            doubled,
            [tuple(2 * a - b for a, b in zip(point, reference)) for point, reference in zip(single, original)],
        )
        self.set_wave_values(node, impulseAmount=1, impulseX=2, impulseZ=2)
        self.assert_points_close(self.points(mesh), doubled)
        self.set_wave_values(node, idleAmplitude=0.3, idleDirectionality=1, phaseSpread=0.4)
        for attribute in ["amplitude", "envelope"]:
            value = cmds.getAttr(node + "." + attribute)
            cmds.setAttr(node + "." + attribute, 0)
            self.assertEqual(self.points(mesh), original)
            cmds.setAttr(node + "." + attribute, value)
        cmds.setAttr(node + ".weightList[0].weights[1]", 0)
        self.assertEqual(self.points(mesh)[1], original[1])
        for index in [0, 1]:
            cmds.setAttr(node + ".amplitudeRamp[%d].amplitudeRamp_FloatValue" % index, 0)
        self.assertEqual(self.points(mesh), original)

    def test_wave_expression_same_frame_dirty(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1)])
        node = self.wave(mesh)
        bell = om.MEulerRotation(0, 0.4, 0.2).asMatrix()
        cmds.setAttr(node + ".bellMatrix", *list(bell), type="matrix")
        self.set_wave_values(node, idleDirectionality=0.7, impulseX=-0.5, impulsePosition=0.5)
        for name, value in [
            ("idleAmplitudeV", 0.2),
            ("idleAmplitudeU", 0.2),
            ("idleDirectionality", 0.1),
            ("idleDirectionX", 1),
            ("idleDirectionZ", 1),
            ("idleDirectionSpace", 1),
            ("phaseSpread", 0.4),
            ("impulseAmount", 0),
        ]:
            with self.subTest(attribute=name):
                previous = cmds.getAttr(node + "." + name)
                before = self.points(mesh)
                cmds.setAttr(node + "." + name, value)
                self.assertNotEqual(self.points(mesh), before)
                cmds.setAttr(node + "." + name, previous)
                self.assertEqual(self.points(mesh), before)

    def test_wave_expression_evaluation_order_and_retime(self):
        mesh = self.mesh([(1, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1, 1), (-1, 0.8, 0.2)])
        node = self.wave([mesh + ".vtx[0]", mesh + ".vtx[2:3]"])
        values = {
            "wavePhaseV": (0.1, 1.3),
            "wavePhaseU": (0.3, 0.8),
            "idleAmplitudeV": (0.3, 1.1),
            "idleAmplitudeU": (1, 0.2),
            "idleDirectionality": (0.2, 0.9),
            "idleDirectionX": (-1, 0.4),
            "idleDirectionZ": (0.3, 1),
            "phaseSpread": (0.1, 0.4),
            "impulseAmount": (0, 1.2),
            "impulsePosition": (0.1, 1),
            "noisePhase": (0.2, 0.9),
            "noiseFrequencyU": (1.2, 1.8),
            "amplitudeRamp[1].amplitudeRamp_FloatValue": (0.6, 1),
        }
        cmds.setAttr(node + ".impulseX", -0.5)
        for name, (start, end) in values.items():
            for time, value in [(1, start), (9, end)]:
                cmds.setKeyframe(
                    node + "." + name, time=time, value=value, inTangentType="linear", outTangentType="linear"
                )
        times = [1, 2.5, 5, 8.25, 9]
        references = {}
        try:
            for mode in ["off", "serial", "parallel"]:
                cmds.evaluationManager(mode=mode)
                for time in times + list(reversed(times)) + [8.25, 1, 5, 2.5]:
                    cmds.currentTime(time)
                    result = self.points(mesh)
                    if time in references:
                        self.assert_points_close(result, references[time])
                    else:
                        references[time] = result
            curves = cmds.listConnections(node, source=True, destination=False, type="animCurve") or []
            self.assertEqual(len(set(curves)), len(values))
            cmds.scaleKey(list(set(curves)), timeScale=2, timePivot=0)
            for time in times:
                cmds.currentTime(time * 2)
                self.assert_points_close(self.points(mesh), references[time])
        finally:
            cmds.evaluationManager(mode="off")

    def assert_geometry_transform_invariance(self, mesh, node, baseline):
        for label, values in (
            ("translation", ((4, -3, 2), (0, 0, 0), 1.0, None)),
            ("rotation", ((0, 0, 0), (21, -17, 9), 1.0, None)),
            ("scale", ((0, 0, 0), (0, 0, 0), 1.7, None)),
            ("opm", ((0, 0, 0), (0, 0, 0), 1.0, self.transform_matrix((1.2, -0.7, 2.1), (17, -23, 11), 1.3))),
        ):
            with self.subTest(transform=label):
                cmds.setAttr(mesh + ".translate", *values[0], type="double3")
                cmds.setAttr(mesh + ".rotate", *values[1], type="double3")
                cmds.setAttr(mesh + ".scale", values[2], values[2], values[2], type="double3")
                cmds.setAttr(mesh + ".offsetParentMatrix", *(values[3] or matrix()), type="matrix")
                cmds.dgdirty(node)
                self.assertEqual(self.points(mesh), baseline)
                self.assert_world_is_one_transform(mesh, baseline)

    def test_collision_object_space_is_invariant_under_geometry_transform(self):
        mesh = self.mesh([(0.25, 0.2, 0), (0.25, 1, 0), (0.1, 0.5, 0.4)])
        original = self.points(mesh)
        node = self.collision(mesh)
        self.assert_collision_reference(mesh, node, original)
        baseline = self.points(mesh)
        self.assertNotEqual(baseline, original)
        self.assert_geometry_transform_invariance(mesh, node, baseline)

    def test_wave_bell_local_object_space_is_invariant_under_geometry_transform(self):
        mesh = self.mesh([(1, 0.3, 0), (-1, 0.3, 0), (0, 1, 1), (0, 1, -1)])
        original = self.points(mesh)
        node = self.wave(mesh)
        self.set_wave_values(node, idleAmplitudeU=0, idleDirectionality=1, idleDirectionSpace=1, impulseSpace=1)
        cmds.setAttr(node + ".evaluationToWorldRotation", *matrix(), type="matrix")
        baseline = self.points(mesh)
        self.assertNotEqual(baseline, original)
        self.assert_geometry_transform_invariance(mesh, node, baseline)
