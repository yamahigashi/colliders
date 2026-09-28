"""Panel geometry, material coordinates, validation, and Maya evaluation."""

import itertools
import math
import os
import struct
import tempfile
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import output_object, transform
from test_waist_ring import point_hash, thigh_cylinder, outside_cylinder, waist_frame


def plug(node, attribute):
    return om.MSelectionList().add(node + "." + attribute).getPlug(0)


def doubles(values):
    return struct.pack("<{}d".format(len(values)), *values)


class SkirtPanelTests(unittest.TestCase):
    def setUp(self):
        previous = cmds.evaluationManager(query=True, mode=True)[0]
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode=previous)
        self.surface_data = []
        self.reader_shapes = {}
        self.shape_baselines = {}
        self.metadata_baselines = {}
        self.errors = []
        self.warnings = []

        def capture(message, kind, _data):
            if kind == om.MCommandMessage.kError:
                self.errors.append(message)
            elif kind == om.MCommandMessage.kWarning:
                self.warnings.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)

    def rig(self, skirt_type=1, angle=0, both=False, waist_height=11.5):
        node = cmds.createNode("yddSkirtBellCollider")
        waist = transform((0, waist_height, 0), (180, 0, 0))
        cmds.connectAttr(waist + ".worldMatrix[0]", node + ".bellMatrix")
        joints = {}
        for side, x in (("left", -1), ("right", 1)):
            hip = transform((x, 10, 0))
            knee = transform((x, 5, 0))
            heel = transform((x, 0, 0))
            cmds.parent(knee, heel, hip)
            for part, joint in (("Hip", hip), ("Knee", knee), ("Heel", heel)):
                joints[side + part] = joint
                cmds.connectAttr(joint + ".worldMatrix[0]", node + "." + side + part + "Matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 4)
            if side == "left" or both:
                cmds.setAttr(hip + ".rotateX", -angle)
        for name, value in (("skirtType", skirt_type), ("bellAxis", 1), ("bellSubdivision", 16)):
            cmds.setAttr(node + "." + name, value)
        cmds.setAttr(node + ".bellScale", 2, 1, 2, type="double3")
        cmds.setAttr(node + ".ringScale", 2, 1, 2, type="double3")
        return node, waist, joints

    def seam(self, node, index, u, height=0, enabled=True):
        base = node + ".seams[{}].".format(index)
        cmds.setAttr(base + "enabled", enabled)
        cmds.setAttr(base + "materialU", u)
        cmds.setAttr(base + "startHeight", height)

    def hem(self, node, index, us, heights):
        base = node + ".panelHems[{}].".format(index)
        for name, values in (("uSamples", us), ("heightSamples", heights)):
            cmds.setAttr(base + name, list(values), type="doubleArray")

    def surface(self, node, attribute="outputSurface"):
        # Read through a connected shape: a data object taken straight from an array
        # child plug with asMObject can be released before it is read.
        key = (node, attribute)
        reader = self.reader_shapes.get(key)
        if reader is None or not cmds.objExists(reader):
            reader = cmds.createNode("nurbsSurface", name="reader_" + attribute.replace("[", "_").replace("]", "").replace(".", "_"))
            cmds.connectAttr(node + "." + attribute, reader + ".create")
            self.reader_shapes[key] = reader
        path = om.MSelectionList().add(reader).getDagPath(0)
        fn = om.MFnNurbsSurface(path)
        fn.numCVsInU  # force the evaluation of the connection before the caller reads
        return fn

    def ids(self, node):
        array = plug(node, "outputPatches")
        array.evaluateNumElements()
        return tuple(sorted(array.getExistingArrayAttributeIndices()))

    def patch(self, node, index):
        return self.surface(node, "outputPatches[{}].surface".format(index))

    def array(self, node, attribute):
        data = plug(node, attribute).asMObject()
        return tuple(om.MFnDoubleArrayData(data).array())

    def metadata(self, node, index):
        base = "outputPatches[{}].".format(index)
        return (
            plug(node, base + "materialUStart").asDouble(),
            plug(node, base + "materialUEnd").asDouble(),
            self.array(node, base + "vBreaks"),
            self.array(node, base + "hemUSamples"),
            self.array(node, base + "hemHeightSamples"),
        )

    def snapshot(self, node):
        result = []
        for index in self.ids(node):
            surface = self.patch(node, index)
            result.append(
                (
                    index,
                    self.metadata(node, index),
                    tuple(surface.knotsInU()),
                    tuple(surface.knotsInV()),
                    surface.numCVsInU,
                    surface.numCVsInV,
                    point_hash(surface.cvPositions()),
                )
            )
        return tuple(result)

    def row(self, surface, v):
        # Copy each point while the MPointArray is alive: element access returns
        # references that do not extend the array's lifetime.
        points = surface.cvPositions()
        return [om.MPoint(points[u * surface.numCVsInV + v]) for u in range(surface.numCVsInU)]

    def shape(self, node, attribute):
        shape = cmds.createNode("nurbsSurface")
        cmds.connectAttr(node + "." + attribute, shape + ".create")
        self.shape_baselines[shape] = self.shape_points(shape)
        if attribute.startswith("outputPatches["):
            base = attribute.rsplit(".", 1)[0]
            values = {}
            for name in ("materialUStart", "materialUEnd", "vBreaks", "hemUSamples", "hemHeightSamples"):
                destination = "panel" + name[0].upper() + name[1:]
                if name in ("materialUStart", "materialUEnd"):
                    cmds.addAttr(shape, longName=destination, attributeType="double")
                else:
                    cmds.addAttr(shape, longName=destination, dataType="doubleArray")
                cmds.connectAttr(node + "." + base + "." + name, shape + "." + destination)
                values[destination] = cmds.getAttr(shape + "." + destination)
            self.metadata_baselines[shape] = values
        return shape

    def shape_points(self, shape):
        path = om.MSelectionList().add(shape).getDagPath(0)
        return tuple(tuple(p) for p in om.MFnNurbsSurface(path).cvPositions())

    def assert_position(self, a, b, tolerance=1e-6):
        self.assertLessEqual((om.MPoint(a) - om.MPoint(b)).length(), tolerance)

    def assert_open(self, surface):
        self.assertEqual((surface.degreeInU, surface.degreeInV), (3, 1))
        self.assertEqual((surface.formInU, surface.formInV), (om.MFnNurbsSurface.kOpen, om.MFnNurbsSurface.kOpen))
        self.assertEqual(surface.knotDomainInU, (0.0, 1.0))
        self.assertEqual(surface.knotDomainInV, (0.0, 1.0))
        knots = tuple(surface.knotsInU())
        self.assertEqual(knots[:3], (0.0,) * 3)
        self.assertEqual(knots[-3:], (1.0,) * 3)
        self.assertGreaterEqual(surface.numCVsInU, 4)
        for value in set(knots[3:-3]):
            self.assertLessEqual(knots.count(value), 3)
        self.assertTrue(all(math.isfinite(c) for p in surface.cvPositions() for c in p))
        self.assertTrue(all(p.w == 1 for p in surface.cvPositions()))

    def assert_reserved_patch(self, node):
        self.assertEqual(self.ids(node), (0,))
        surface = self.patch(node, 0)
        beta = surface.knotsInV()[-2]
        self.assertEqual(self.metadata(node, 0), (0.0, 1.0, (beta,), (), ()))

    def request(self, node, attribute):
        target = plug(node, attribute)
        if target.isArray:
            target.evaluateNumElements()
        elif target.isCompound:
            target.child(0).asMObject()
        elif attribute.endswith("surface") or attribute == "outputSurface":
            target.asMObject()
        else:
            target.asDouble()

    def assert_invalid_preserves(self, node, attribute, condition, shapes, reason=None):
        expected = [self.shape_baselines[shape] for shape in shapes]
        nodes = set(cmds.ls(long=True))
        connections = cmds.listConnections(node, connections=True, plugs=True) or []
        for attempt in range(3):
            with self.subTest(request=attribute, attempt=attempt):
                self.errors.clear()
                if attempt:
                    # Maya may keep an array child clean after a failed evaluation; a
                    # dirty propagation from any input must re-run the validation.
                    cmds.dgdirty(node)
                try:
                    self.request(node, attribute)
                except RuntimeError:
                    pass
                actual = [self.shape_points(shape) for shape in shapes]
                for shape in shapes:
                    for name, value in self.metadata_baselines.get(shape, {}).items():
                        self.assertEqual(cmds.getAttr(shape + "." + name), value)
                messages = [message for message in self.errors if "condition " in message]
                self.assertTrue(messages, self.errors)
                self.assertIn(node.split("|")[-1], messages[0])
                self.assertIn("condition {}:".format(condition), messages[0])
                if reason:
                    self.assertIn(reason, messages[0])
                self.assertEqual([point_hash(p) for p in actual], [point_hash(p) for p in expected])
                self.assertEqual(set(cmds.ls(long=True)), nodes)
                self.assertEqual(cmds.listConnections(node, connections=True, plugs=True) or [], connections)

    def test_seamless_new_path_regression(self):
        for skirt_type, smooth, follow, tightness, state in itertools.product((0, 1), (0, 0.5), (0, 1), (0, 0.5, 1), (0, 1)):
            with self.subTest(skirt_type=skirt_type, smooth=smooth, follow=follow, tightness=tightness, state=state):
                node, _, _ = self.rig(skirt_type, angle=60)
                for name, value in (("smoothness", smooth), ("follow", follow), ("tightness", tightness), ("nodeState", state)):
                    cmds.setAttr(node + "." + name, value)
                reference = self.surface(node)
                expected = point_hash(reference.cvPositions())
                identity_values = tuple(om.MMatrix())
                for configuration in ("empty", "disabled", "disabled_columns", "unit_columns", "default_hem"):
                    with self.subTest(configuration=configuration):
                        if configuration == "disabled":
                            self.seam(node, 9, float("nan"), float("inf"), False)
                        elif configuration == "disabled_columns":
                            cmds.setAttr(node + ".columnMaterialU", [float("nan")], type="doubleArray")
                            cmds.setAttr(node + ".columnOffsetMatrix[12]", *identity_values, type="matrix")
                            cmds.removeMultiInstance(node + ".columnOffsetMatrix[12]", b=True)
                        elif configuration == "unit_columns":
                            for index, u in ((1, 0.05), (8, 0.4), (20, 0.8)):
                                values = [0.5] * (index + 1)
                                values[index] = u
                                cmds.setAttr(node + ".columnMaterialU", values, type="doubleArray")
                                identity_values = tuple(om.MMatrix())
                                cmds.setAttr(
                                    node + ".columnOffsetMatrix[{}]".format(index),
                                    *identity_values,
                                    type="matrix",
                                )
                        elif configuration == "default_hem":
                            self.hem(node, 0, (0, 0.3, 1), (1, 1, 1))
                        surface = self.surface(node)
                        self.assertEqual((surface.degreeInU, surface.degreeInV), (3, 1))
                        self.assertEqual(
                            (surface.formInU, surface.formInV), (om.MFnNurbsSurface.kPeriodic, om.MFnNurbsSurface.kOpen)
                        )
                        self.assert_reserved_patch(node)
                        self.assertEqual(point_hash(surface.cvPositions()), expected)
                        for v in range(surface.numCVsInV):
                            row = self.row(surface, v)
                            self.assertEqual(point_hash(row[:3]), point_hash(row[-3:]))
                        lo, hi = surface.knotDomainInU
                        for v in surface.knotsInV():
                            self.assert_position(surface.getPointAtParam(lo, v), surface.getPointAtParam(hi, v), 1e-9)

    def test_seamless_reserved_patch(self):
        cases = itertools.product((0, 1), (0, 60), ((0, 0), (0.5, 1)), (0, 0.5, 1), (8, 16, 64), (0, 1))
        for skirt_type, angle, (smooth, follow), tightness, subdivision, state in cases:
            with self.subTest(
                skirt_type=skirt_type,
                angle=angle,
                smooth=smooth,
                follow=follow,
                tightness=tightness,
                subdivision=subdivision,
                state=state,
            ):
                node, _, _ = self.rig(skirt_type, angle=angle)
                cmds.setAttr(node + ".smoothness", smooth)
                cmds.setAttr(node + ".follow", follow)
                cmds.setAttr(node + ".tightness", tightness)
                cmds.setAttr(node + ".bellSubdivision", subdivision)
                cmds.setAttr(node + ".nodeState", state)
                periodic = self.surface(node)
                self.assert_reserved_patch(node)
                surface = self.patch(node, 0)
                self.assert_open(surface)
                expected_knots = (0.0,) * 3 + tuple(i / subdivision for i in range(1, subdivision)) + (1.0,) * 3
                self.assertEqual(len(surface.knotsInU()), len(expected_knots))
                for actual, expected in zip(surface.knotsInU(), expected_knots):
                    self.assertAlmostEqual(actual, expected, delta=1e-12)
                self.assertEqual(tuple(surface.knotsInV()), tuple(periodic.knotsInV()))
                for u, v in itertools.product((i / 10 for i in range(11)), (0, 0.3, 0.7, 1)):
                    self.assert_position(surface.getPointAtParam(u, v), periodic.getPointAtParam(u, v), 1e-9)
                if state or tightness or subdivision != 16:
                    # Keep the node: reader shapes are keyed by node name and a
                    # recreated node with the same name would read stale data.
                    continue
                breaks = self.metadata(node, 0)[2]
                self.assertTrue(breaks)
                self.seam(node, 0, 0.25, 0)
                seam_breaks = self.metadata(node, 1)[2]
                self.assertEqual(len(breaks), len(seam_breaks))
                for actual, expected in zip(breaks, seam_breaks):
                    self.assertAlmostEqual(actual, expected, delta=1e-9)

    def test_column_offset_frame_and_displacements(self):
        node, _, _ = self.rig(skirt_type=1)
        cmds.setAttr(node + ".nodeState", 1)
        baseline = self.surface(node)
        baseline_points = baseline.cvPositions()
        baseline_hash = point_hash(baseline.cvPositions())
        baseline_patch_points = self.patch(node, 0).cvPositions()
        shape = self.shape(node, "outputSurface")
        identity_values = tuple(om.MMatrix())
        for count in (0, 1, 16, 32):
            for index in cmds.getAttr(node + ".columnOffsetMatrix", multiIndices=True) or []:
                cmds.removeMultiInstance(node + ".columnOffsetMatrix[{}]".format(index), b=True)
            material_u = [0.5] * (count + 1)
            for index in range(count):
                material_u[index] = (index + 1.0) / (count + 2.0)
                cmds.setAttr(
                    node + ".columnOffsetMatrix[{}]".format(index), *identity_values, type="matrix"
                )
            cmds.setAttr(node + ".columnMaterialU", material_u, type="doubleArray")
            self.assertEqual(point_hash(self.surface(node).cvPositions()), baseline_hash)
        for index in cmds.getAttr(node + ".columnOffsetMatrix", multiIndices=True) or []:
            cmds.removeMultiInstance(node + ".columnOffsetMatrix[{}]".format(index), b=True)
        cmds.setAttr(node + ".columnMaterialU", [0.5, 0.5, 0.5, 0.5, 0.25], type="doubleArray")
        translated = om.MMatrix()
        translated.setElement(3, 0, 0.2)
        translated_values = tuple(translated)
        cmds.setAttr(node + ".columnOffsetMatrix[4]", *translated_values, type="matrix")
        result = self.surface(node)
        self.assertEqual(result.numCVsInU, baseline.numCVsInU)
        self.assertEqual(tuple(result.knotsInU()), tuple(baseline.knotsInU()))
        self.assertNotEqual(point_hash(result.cvPositions()), baseline_hash)

        def assert_uniform_translation(surface, distance):
            displacements = [surface.cvPosition(u, v) - baseline_points[u * baseline.numCVsInV + v]
                             for u in range(surface.numCVsInU) for v in range(surface.numCVsInV)]
            first = displacements[0]
            for delta in displacements[1:]:
                self.assertLessEqual((delta - first).length(), 1e-9)
            self.assertAlmostEqual(first.length(), distance, delta=1e-9)
            bell_matrix = om.MMatrix(cmds.getAttr(node + ".bellMatrix"))
            bell_axis = om.MVector(bell_matrix[4], bell_matrix[5], bell_matrix[6]).normal()
            self.assertAlmostEqual(first * bell_axis, 0.0, delta=1e-9)

        assert_uniform_translation(result, 0.2)
        for index in cmds.getAttr(node + ".columnOffsetMatrix", multiIndices=True) or []:
            cmds.removeMultiInstance(node + ".columnOffsetMatrix[{}]".format(index), b=True)
        cmds.setAttr(node + ".columnMaterialU", [0.25, 0.0, 0.75], type="doubleArray")
        cmds.setAttr(node + ".columnOffsetMatrix[0]", *translated_values, type="matrix")
        cmds.setAttr(node + ".columnOffsetMatrix[2]", *identity_values, type="matrix")
        # The base row control points lie on a circle at uniform material
        # spacing, so on those CVs the translation field, linear in material
        # U between the two columns, is linear in azimuth about the bell axis:
        # 0.2 where the translated column sits, 0 on the opposite side. CVs that
        # knot insertion adds at the clamped ends leave the circle and are
        # skipped.
        patch = self.patch(node, 0)
        subdivision = cmds.getAttr(node + ".bellSubdivision")
        bell_values = cmds.getAttr(node + ".bellMatrix")
        centre = om.MVector(bell_values[12], bell_values[13], bell_values[14])
        axis = om.MVector(bell_values[4], bell_values[5], bell_values[6]).normal()
        reference = om.MVector(bell_values[0], bell_values[1], bell_values[2])
        reference = (reference - axis * (reference * axis)).normal()
        binormal = axis ^ reference

        def azimuth(point):
            radial = om.MVector(point) - centre
            radial -= axis * (radial * axis)
            return math.atan2(radial * binormal, radial * reference), radial.length()

        checked = 0
        for v in range(patch.numCVsInV):
            rows = []
            for index in range(patch.numCVsInU):
                rest = baseline_patch_points[index * patch.numCVsInV + v]
                angle, radius = azimuth(rest)
                delta = patch.cvPosition(index, v) - rest
                rows.append((angle, radius, delta.length()))
            radius_max = max(radius for _, radius, _ in rows)
            # Base rows pass through float rounding, so the circle test and the
            # azimuth-derived expectation carry that tolerance.
            on_circle = [row for row in rows if abs(row[1] - radius_max) <= 1e-5 * radius_max]
            # Clamping the reserved patch at material 0 replaces the circle
            # point there with the curve point, so one circle CV is missing.
            self.assertGreaterEqual(len(on_circle), subdivision - 1)
            peak_angle, _, peak = max(on_circle, key=lambda row: row[2])
            self.assertAlmostEqual(peak, 0.2, delta=1e-9)
            for angle, _, length in on_circle:
                difference = abs(angle - peak_angle)
                difference = min(difference, 2.0 * math.pi - difference)
                self.assertAlmostEqual(length, 0.2 * (1.0 - difference / math.pi), delta=1e-6)
                checked += 1
        self.assertGreaterEqual(checked, (subdivision - 1) * patch.numCVsInV)

    def test_column_offset_material_interpolation(self):
        node, _, _ = self.rig(skirt_type=0)
        cmds.setAttr(node + ".nodeState", 1)
        baseline_hash = point_hash(self.surface(node).cvPositions())
        shape = self.shape(node, "outputSurface")
        cmds.setAttr(node + ".columnMaterialU", [float("nan")], type="doubleArray")
        identity_values = tuple(om.MMatrix())
        cmds.setAttr(node + ".columnOffsetMatrix[2]", *identity_values, type="matrix")
        cmds.removeMultiInstance(node + ".columnOffsetMatrix[2]", b=True)
        self.surface(node)
        self.assertEqual(cmds.getAttr(node + ".columnOffsetMatrix", multiIndices=True) or [], [])
        self.assertEqual(point_hash(self.surface(node).cvPositions()), baseline_hash)
        cmds.setAttr(node + ".columnOffsetMatrix[2]", *identity_values, type="matrix")
        for value in (-1e-12, 1.0, float("nan"), float("inf")):
            with self.subTest(material_u=value):
                cmds.setAttr(node + ".columnMaterialU", [0.0, 0.0, value], type="doubleArray")
                self.assert_invalid_preserves(node, "outputSurface", 17, [shape])
        cmds.setAttr(node + ".columnMaterialU", [0.0, 0.0], type="doubleArray")
        self.assert_invalid_preserves(node, "outputSurface", 17, [shape])
        cmds.setAttr(node + ".columnMaterialU", [0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25 + 0.5e-9], type="doubleArray")
        cmds.setAttr(node + ".columnOffsetMatrix[9]", *identity_values, type="matrix")
        self.assert_invalid_preserves(node, "outputSurface", 18, [shape])
        cmds.setAttr(node + ".columnMaterialU", [0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75], type="doubleArray")
        self.assertEqual(self.surface(node).numCVsInV, 3)

    def test_material_u_and_patch_ids(self):
        for configuration in (((0, 0.0),), ((19, 0.25), (0, 0.73)), ((19, 0.0), (0, 0.25), (7, 0.73))):
            with self.subTest(seams=configuration):
                node, _, _ = self.rig()
                cmds.setAttr(node + ".nodeState", 1)
                reference = self.surface(node)
                lo, hi = reference.knotDomainInU
                self.assertNotEqual(reference.knotsInU()[0], lo)
                for index, u in configuration:
                    self.seam(node, index, u, 1)
                ordered = sorted(configuration, key=lambda item: item[1])
                self.assertEqual(self.ids(node), tuple(sorted(i + 1 for i, _ in ordered)))
                for rank, (index, a) in enumerate(ordered):
                    b = ordered[(rank + 1) % len(ordered)][1] + (rank + 1 == len(ordered))
                    surface = self.patch(node, index + 1)
                    self.assert_open(surface)
                    self.assertEqual(self.metadata(node, index + 1)[:2], (a, b))
                    for u, v in itertools.product((0, 0.13, 0.5, 0.91, 1), (0, 0.2, 0.7, 1)):
                        s = (a + u * (b - a)) % 1.0
                        self.assert_position(surface.getPointAtParam(u, v), reference.getPointAtParam(lo + s * (hi - lo), v))
                    for u, s in ((0, a), (surface.numCVsInU - 1, b)):
                        for v, t in enumerate(surface.knotsInV()):
                            self.assert_position(
                                surface.cvPosition(u, v), reference.getPointAtParam(lo + (s % 1) * (hi - lo), t)
                            )

    def test_start_heights_and_state(self):
        for heights in ((0, 0), (4.0 / 11.5, 4.0 / 11.5), (0.45, 0.45), (0.25, 0.7), (1, 1), (0.5, 0.5 + 0.75e-9)):
            with self.subTest(heights=heights):
                node, _, _ = self.rig(angle=65)
                for index, u, height in zip((0, 8), (0.3, 0.8), heights):
                    self.seam(node, index, u, height)
                before = self.snapshot(node)
                for state in (0, 1, 0):
                    cmds.setAttr(node + ".nodeState", state)
                    snapshot = self.snapshot(node)
                    self.assertEqual(tuple(p[:-1] for p in snapshot), tuple(p[:-1] for p in before))
                    if state == 0:
                        self.assertEqual(snapshot, before)
                    for index, other, height in ((1, 9, heights[0]), (9, 1, heights[1])):
                        start, end = self.patch(node, index), self.patch(node, other)
                        for v in start.knotsInV():
                            if height > 0 and v <= height + 1e-9:
                                self.assert_position(start.getPointAtParam(0, v), end.getPointAtParam(1, v), 1e-9)
                                row = tuple(start.knotsInV()).index(v)
                                self.assertEqual(
                                    doubles(tuple(start.cvPosition(0, row))),
                                    doubles(tuple(end.cvPosition(end.numCVsInU - 1, row))),
                                )
                    internal = sorted(set(h for h in heights if 0 < h < 1))
                    if heights == (0.5, 0.5 + 0.75e-9):
                        internal = [0.5]
                    if internal:
                        self.assertEqual(self.metadata(node, 1)[2], tuple(internal))
                cmds.setAttr(node + ".nodeState", 1)
                waist, direction = waist_frame(node)
                for index in self.ids(node):
                    surface = self.patch(node, index)
                    for v, t in enumerate(surface.knotsInV()):
                        for u, p in enumerate(self.row(surface, v)):
                            self.assertAlmostEqual(
                                (p - waist) * direction,
                                11.5 * t,
                                delta=1e-6,
                                msg="patch {} row {} cv {} point {}".format(index, v, u, tuple(p)),
                            )

    def test_path_combinations(self):
        for smooth, follow in ((0.5, 0), (0, 1), (0.5, 1)):
            with self.subTest(smooth=smooth, follow=follow):
                node, _, joints = self.rig(angle=0, waist_height=10.25)
                cmds.setAttr(node + ".bellScale", 1.5, 1, 1.5, type="double3")
                cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
                cmds.setAttr(node + ".smoothness", smooth)
                cmds.setAttr(node + ".follow", follow)
                self.assert_reserved_patch(node)
                self.seam(node, 0, 5.0 / 16, 0)
                self.seam(node, 8, 13.0 / 16, 0)
                surfaces = {i: self.patch(node, i) for i in self.ids(node)}
                opposite = max(surfaces, key=lambda i: surfaces[i].getPointAtParam(0.5, 0).x)
                active = next(i for i in surfaces if i != opposite)
                before = {i: point_hash(s.cvPositions()) for i, s in surfaces.items()}
                near = []
                cylinder = thigh_cylinder(node, "right")
                for j, p in enumerate(surfaces[opposite].cvPositions()):
                    local = p * cylinder["inverse"]
                    if 0 < local.y < 1 and abs(math.hypot(local.x, local.z) - 1) < 0.6:
                        near.append(j)
                self.assertTrue(near, "The unchanged set must include CVs near the stationary ring")
                points_before = surfaces[opposite].cvPositions()
                near_before = point_hash([om.MPoint(points_before[j]) for j in near])
                cmds.setAttr(joints["leftHip"] + ".rotateX", -90)
                changed = {i: self.patch(node, i) for i in surfaces}
                self.assertNotEqual(point_hash(changed[active].cvPositions()), before[active])
                self.assertEqual(point_hash(changed[opposite].cvPositions()), before[opposite])
                points_after = changed[opposite].cvPositions()
                self.assertEqual(point_hash([om.MPoint(points_after[j]) for j in near]), near_before)
                for p in changed[opposite].cvPositions():
                    self.assertTrue(outside_cylinder(p, cylinder, 1e-6))
                cmds.setAttr(joints["rightHip"] + ".rotateX", -90)
                self.assertNotEqual(point_hash(self.patch(node, opposite).cvPositions()), before[opposite])
                self.assertNotEqual(point_hash(self.patch(node, active).cvPositions()), before[active])

    def test_has_no_effect_inserted_rows(self):
        for calf in (5.0, 5e-5):
            with self.subTest(calf=calf):
                node, _, joints = self.rig()
                cmds.setAttr(node + ".nodeState", 1)
                cmds.setAttr(node + ".bellScale1", 1.7)
                for side, x in (("left", -1), ("right", 1)):
                    cmds.xform(joints[side + "Heel"], worldSpace=True, translation=(x, 5 - calf, 0))
                height = 6.5 + calf
                t = (6.5 + calf * 0.5) / height
                self.seam(node, 0, 0.3, t)
                surface = self.patch(node, 1)
                row = tuple(surface.knotsInV()).index(t)
                waist, direction = waist_frame(node)
                expected = 6.5 + 0.5 * max(calf, 1e-4) * 1.7
                for p in self.row(surface, row):
                    self.assertAlmostEqual((p - waist) * direction, expected, delta=1e-6)

    def test_path_combinations_cut_end_injection(self):
        self.skipTest(
            "The node exposes no direct-field injection input; cut-end coefficients are tested through the C++ row API"
        )

    def test_waist_and_long_interactions(self):
        for tightness in (0, 0.5, 1):
            with self.subTest(tightness=tightness):
                node, _, _ = self.rig(angle=90)
                cmds.setAttr(node + ".tightness", tightness)
                self.seam(node, 0, 0.3, 0)
                cmds.setAttr(node + ".nodeState", 1)
                raw = point_hash(self.row(self.patch(node, 1), 0))
                cmds.setAttr(node + ".nodeState", 0)
                baseline = None
                for follow in (0, 1):
                    cmds.setAttr(node + ".follow", follow)
                    surface = self.patch(node, 1)
                    row = self.row(surface, 0)
                    self.assertNotEqual(point_hash(row), raw)
                    if baseline is None:
                        baseline = point_hash(row)
                    else:
                        self.assertEqual(point_hash(row), baseline)
                    for p in row:
                        for side in ("left", "right"):
                            self.assertTrue(outside_cylinder(p, thigh_cylinder(node, side), 2e-6))
                    samples = {(u, v): surface.getPointAtParam(u, v) for u in (0, 0.2, 0.5, 0.8, 1) for v in surface.knotsInV()}
                    self.seam(node, 8, 0.7, 0.5)
                    self.seam(node, 15, 0.9, 0.65)
                    for (u, v), expected in samples.items():
                        if v <= 0.5:
                            material = 0.3 + u
                            # A sample at the cut boundary belongs to the bank it was taken from:
                            # match the material position directly first and wrap only when no
                            # patch contains it, so the two banks of the seam are not conflated.
                            hits = []
                            for index in self.ids(node):
                                a, b = self.metadata(node, index)[:2]
                                if a <= material <= b:
                                    hits.append((index, a, b, material))
                            if not hits:
                                for index in self.ids(node):
                                    a, b = self.metadata(node, index)[:2]
                                    if a <= material + 1 <= b:
                                        hits.append((index, a, b, material + 1))
                            positions = [self.patch(node, index).getPointAtParam((s - a) / (b - a), v) for index, a, b, s in hits]
                            for position in positions:
                                self.assertLessEqual((om.MPoint(position) - om.MPoint(expected)).length(), 2e-6,
                                                     ((u, v), hits, tuple(expected)[:3], tuple(position)[:3]))
                    for index in (8, 15):
                        cmds.removeMultiInstance(node + ".seams[{}]".format(index), b=True)
                    self.seam(node, 8, 0.7, 0.5)
                    self.seam(node, 15, 0.9, 0.65)
                    physical = (0.0, 4.0 / 11.5, 6.5 / 11.5, 1.0)
                    original_rows = {}
                    for index in self.ids(node):
                        surface = self.patch(node, index)
                        for t in physical:
                            row = tuple(surface.knotsInV()).index(t)
                            original_rows[index, t] = point_hash(self.row(surface, row))
                    cmds.setAttr(node + ".seams[8].startHeight", 0.53)
                    cmds.setAttr(node + ".seams[15].startHeight", 0.7)
                    for (index, t), expected in original_rows.items():
                        surface = self.patch(node, index)
                        row = tuple(surface.knotsInV()).index(t)
                        self.assertEqual(point_hash(self.row(surface, row)), expected)
                    for index in (8, 15):
                        cmds.removeMultiInstance(node + ".seams[{}]".format(index), b=True)

    def test_standard_samplers(self):
        for waist_height in (11.5, 14):
            with self.subTest(waist_height=waist_height):
                node, _, _ = self.rig(angle=60, waist_height=waist_height)
                self.seam(node, 0, 0.3, 0.4)
                self.seam(node, 8, 0.8, 0.4)
                for index in self.ids(node):
                    attribute = ".outputPatches[{}].surface".format(index)
                    surface = self.patch(node, index)
                    posi = cmds.createNode("pointOnSurfaceInfo")
                    cmds.connectAttr(node + attribute, posi + ".inputSurface")
                    cmds.setAttr(posi + ".turnOnPercentage", False)
                    pin = cmds.createNode("uvPin")
                    for name in ("deformedGeometry", "originalGeometry"):
                        cmds.connectAttr(node + attribute, pin + "." + name)
                    cmds.setAttr(pin + ".normalizedIsoParms", False)
                    cmds.setAttr(pin + ".normalAxis", 1)
                    cmds.setAttr(pin + ".tangentAxis", 0)
                    for u, v in itertools.product((0, 0.23, 0.7, 1), (0.1, 0.4, 0.4000001, 0.8, 1)):
                        with self.subTest(panel=index, u=u, v=v):
                            cmds.setAttr(posi + ".parameterU", u)
                            cmds.setAttr(posi + ".parameterV", v)
                            cmds.setAttr(pin + ".coordinate[0].coordinateU", u)
                            cmds.setAttr(pin + ".coordinate[0].coordinateV", v)
                            expected = surface.getPointAtParam(u, v)
                            self.assert_position(cmds.getAttr(posi + ".position")[0], expected)
                            matrix = cmds.getAttr(pin + ".outputMatrix[0]")
                            self.assert_position(matrix[12:15], expected)
                            tangent_v = v + 1e-10 if v == 0.4 else v
                            tangent = surface.tangents(u, tangent_v)[0].normal()
                            normal = surface.normal(u, tangent_v).normal()
                            self.assertLessEqual(
                                (om.MVector(cmds.getAttr(posi + ".normalizedTangentU")[0]) - tangent).length(), 1e-6
                            )
                            self.assertLessEqual(
                                (om.MVector(cmds.getAttr(posi + ".normalizedNormal")[0]) - normal).length(), 1e-6
                            )
                            axes = [om.MVector(matrix[k : k + 3]) for k in (0, 4, 8)]
                            self.assertLessEqual((axes[0] - tangent).length(), 1e-6)
                            self.assertLessEqual((axes[1] - normal).length(), 1e-6)
                            for a in range(3):
                                for b in range(3):
                                    self.assertAlmostEqual(axes[a] * axes[b], float(a == b), delta=1e-6)

    def test_hem_mapping(self):
        for start in (None, 0.4, 1):
            for profile in ((0.9, 0.9, 0.9), (0.9, 0.8, 0.9)):
                with self.subTest(start=start, profile=profile):
                    node, _, _ = self.rig()
                    cmds.setAttr(node + ".nodeState", 1)
                    periodic = self.surface(node)
                    if start is not None:
                        self.seam(node, 0, 0.3, start)
                    index = 0 if start is None else 1
                    if start is None:
                        self.hem(node, 0, (0, 0.5, 1), (1, 1, 1 - 1e-12))
                    independent = self.patch(node, index)
                    beta = 0.4 if start == 0.4 else 6.5 / 11.5
                    waist, direction = waist_frame(node)
                    ku = list(independent.knotsInU())
                    closed = independent.formInU == om.MFnNurbsSurface.kPeriodic
                    full_knots = (
                        [ku[0] - (ku[1] - ku[0])] + ku + [ku[-1] + (ku[-1] - ku[-2])] if closed else [ku[0]] + ku + [ku[-1]]
                    )
                    lo, hi = independent.knotDomainInU
                    expected_points = om.MPointArray(independent.cvPositions())
                    for u in range(independent.numCVsInU):
                        g = (sum(full_knots[u + 1 : u + 4]) / 3.0 - lo) / (hi - lo)
                        x = g
                        h = (
                            profile[0] + (profile[1] - profile[0]) * x * 2
                            if x <= 0.5
                            else profile[1] + (profile[2] - profile[1]) * (x - 0.5) * 2
                        )
                        knots_v = list(independent.knotsInV())
                        physical = [k for k in knots_v if any(abs(k - q) <= 1e-9 for q in (0.0, 4.0 / 11.5, 6.5 / 11.5, 1.0))]
                        for v, t in enumerate(knots_v):
                            if t <= beta:
                                continue
                            mapped = beta + (t - beta) * (h - beta) / (1 - beta)
                            lower = max(k for k in physical if k <= mapped)
                            upper = min(k for k in physical if k >= mapped)
                            if upper <= lower:
                                continue
                            row_lower, row_upper = knots_v.index(lower), knots_v.index(upper)
                            base_lower = om.MPoint(independent.cvPosition(u, row_lower))
                            base_upper = om.MPoint(independent.cvPosition(u, row_upper))
                            lam = (mapped - lower) / (upper - lower)
                            cv = u * independent.numCVsInV + v
                            expected_points[cv] = om.MPoint(base_lower + (base_upper - base_lower) * lam)
                    data = om.MFnNurbsSurfaceData().create()
                    self.surface_data.append(data)
                    expected = om.MFnNurbsSurface()
                    expected.create(
                        expected_points,
                        independent.knotsInU(),
                        independent.knotsInV(),
                        3,
                        1,
                        independent.formInU,
                        om.MFnNurbsSurface.kOpen,
                        False,
                        data,
                    )
                    self.hem(node, index, (0, 0.5, 1), profile)
                    surface = self.patch(node, index)
                    self.assert_open(surface)
                    meta = self.metadata(node, index)
                    self.assertEqual(meta[3:], ((0.0, 0.5, 1.0), profile))
                    self.assertEqual(meta[:2], (0, 1) if start is None else (0.3, 1.3))
                    self.assertEqual(meta[2], (beta,))
                    last = -1.0
                    for v, t in enumerate(surface.knotsInV()):
                        self.assertGreater(t, last)
                        last = t
                        for u in (0, 0.13, 0.5, 0.77, 1):
                            material = (meta[0] + u * (meta[1] - meta[0])) % 1
                            parameter = lo + material * (hi - lo) if closed else u
                            self.assert_position(surface.getPointAtParam(u, t), expected.getPointAtParam(parameter, t))
                    if not closed:
                        self.assertEqual(tuple(surface.knotsInU()), tuple(expected.knotsInU()))
                        self.assertEqual(surface.numCVsInU, expected.numCVsInU)
                        for a, b in zip(surface.cvPositions(), expected.cvPositions()):
                            self.assert_position(a, b)
                    shape = self.shape(node, "outputPatches[{}].surface".format(index))
                    before = self.snapshot(node)
                    self.hem(node, index, (0, 1), (beta, 0.9))
                    self.assert_invalid_preserves(node, "outputPatches", 8, [shape])
                    self.hem(node, index, (0, 0.5, 1), profile)
                    self.assertEqual(self.snapshot(node), before)

    def test_uncut_hem_boundary_tolerance(self):
        base = 0.875
        for delta in (0.9999999e-9, 1e-9, 1.0000002e-9):
            with self.subTest(delta=delta):
                node, _, _ = self.rig()
                cmds.setAttr(node + ".nodeState", 1)
                self.seam(node, 8, 0.3, 1)
                self.seam(node, 0, 0.8, 0)
                self.hem(node, 1, (0, 1), (base, base))
                self.hem(node, 9, (0, 1), (base, base))
                shapes = [self.shape(node, "outputPatches[{}].surface".format(i)) for i in (1, 9)]
                expected = tuple(
                    self.patch(node, 1).cvPosition(self.patch(node, 1).numCVsInU - 1, self.patch(node, 1).numCVsInV - 1)
                )
                self.hem(node, 9, (0, 1), (base + delta, base))
                if (base + delta) - base > 1e-9:
                    self.assert_invalid_preserves(node, "outputPatches", 8, shapes)
                else:
                    left, right = self.patch(node, 1), self.patch(node, 9)
                    self.assertEqual(doubles(tuple(left.cvPosition(left.numCVsInU - 1, left.numCVsInV - 1))), doubles(expected))
                    self.assertEqual(doubles(tuple(right.cvPosition(0, right.numCVsInV - 1))), doubles(expected))

    def test_validation_and_recovery(self):
        invalids = (
            ("bellSubdivision", 2, 1),
            ("bellSubdivision", 4097, 1),
            ("skirtType", 2, 2),
            ("leftRingAxis", 6, 3),
            ("rightRingAxis", -1, 3),
            ("bellAxis", 6, 3),
            ("bellScale0", float("nan"), 9),
            ("followRange", -1, 13),
            ("followRange", float("nan"), 13),
            ("followRange", float("inf"), 13),
        )
        for panel, request in (
            (False, "outputSurface"),
            (False, "outputPatches"),
            (False, "outputPatches[0]"),
            (False, "outputPatches[0].surface"),
            (False, "outputPatches[0].materialUStart"),
            (True, "outputPatches"),
            (True, "outputPatches[1]"),
            (True, "outputPatches[1].surface"),
            (True, "outputPatches[1].materialUStart"),
        ):
            for name, value, condition in invalids:
                with self.subTest(panel=panel, request=request, attribute=name, value=value):
                    node, _, _ = self.rig()
                    if panel:
                        self.seam(node, 0, 0.3, 0.4)
                    output = "outputPatches[1].surface" if panel else "outputSurface"
                    saved = cmds.getAttr(node + "." + name)
                    driver = cmds.createNode("addDoubleLinear")
                    cmds.setAttr(driver + ".input1", saved)
                    cmds.connectAttr(driver + ".output", node + "." + name)
                    shape = self.shape(node, output)
                    shapes = [shape]
                    if not panel:
                        self.assert_reserved_patch(node)
                        shapes.append(self.shape(node, "outputPatches[0].surface"))
                    before = point_hash(self.shape_points(shape))
                    # A connected driver bypasses attribute range limits.
                    cmds.setAttr(driver + ".input1", value)
                    reason = "bellScale" if name == "bellScale0" else name
                    self.assert_invalid_preserves(node, request, condition, shapes, reason)
                    cmds.setAttr(driver + ".input1", saved)
                    self.assertEqual(point_hash(self.shape_points(shape)), before)
                    if not panel:
                        self.assert_reserved_patch(node)
                        self.assertEqual(self.shape_points(shapes[1]), self.shape_baselines[shapes[1]])
        for attribute, value in (
            ("materialU", -0.1),
            ("materialU", 1),
            ("materialU", float("nan")),
            ("startHeight", -0.1),
            ("startHeight", 1.1),
            ("startHeight", float("inf")),
        ):
            with self.subTest(attribute=attribute, value=value):
                node, _, _ = self.rig()
                self.seam(node, 0, 0.3, 0.4)
                shape = self.shape(node, "outputPatches[1].surface")
                before = self.snapshot(node)
                cmds.setAttr(node + ".seams[0]." + attribute, value)
                self.assert_invalid_preserves(node, "outputPatches", 4, [shape], "seams[0]")
                self.seam(node, 0, 0.3, 0.4)
                self.assertEqual(self.snapshot(node), before)
        for us, heights in (
            ((0,), (1,)),
            ((0, 1), (1,)),
            ((0.1, 1), (1, 1)),
            ((0, 0.5), (1, 1)),
            ((0, 0.5, 0.5, 1), (1, 1, 1, 1)),
            ((0, float("nan"), 1), (1, 1, 1)),
            ((0, 1), (0, 1)),
            ((0, 1), (1.1, 1)),
            ((0, 1), (float("inf"), 1)),
        ):
            with self.subTest(us=us, heights=heights):
                node, _, _ = self.rig()
                self.seam(node, 0, 0.3, 0.4)
                shape = self.shape(node, "outputPatches[1].surface")
                before = self.snapshot(node)
                self.hem(node, 1, us, heights)
                self.assert_invalid_preserves(node, "outputPatches", 6, [shape], "panelHems[1]")
                self.hem(node, 1, (), ())
                self.assertEqual(self.snapshot(node), before)

    def test_validation_per_plug_and_order(self):
        node, _, _ = self.rig()
        driver = cmds.createNode("addDoubleLinear")
        cmds.setAttr(driver + ".input1", 16)
        cmds.connectAttr(driver + ".output", node + ".bellSubdivision")
        shape = self.shape(node, "outputSurface")
        before = point_hash(self.shape_points(shape))
        for state in (0, 1):
            for start in (0, 1):
                with self.subTest(state=state, start=start):
                    cmds.setAttr(node + ".nodeState", state)
                    self.seam(node, 0, 0.3, start)
                    cmds.setAttr(driver + ".input1", 2)
                    self.assert_invalid_preserves(node, "outputSurface", 1, [shape])
                    cmds.setAttr(driver + ".input1", 16)
                    cmds.setAttr(node + ".seams[0].materialU", -1)
                    self.assert_invalid_preserves(node, "outputSurface", 4, [shape])
                    cmds.setAttr(node + ".seams[0].materialU", 0.3)
                    self.assert_invalid_preserves(node, "outputSurface", 11, [shape])
        cmds.removeMultiInstance(node + ".seams[0]", b=True)
        cmds.setAttr(node + ".nodeState", 0)
        self.assertEqual(point_hash(self.shape_points(shape)), before)
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(reference=value):
                cmds.setAttr(node + ".referenceMaterialHeight", value)
                self.errors.clear()
                self.surface(node)
                self.assert_reserved_patch(node)
                self.assertFalse(self.errors)
                self.assert_invalid_preserves(node, "outputReferenceHeight", 10, [shape])
        cmds.setAttr(node + ".referenceMaterialHeight", 1e-8)
        self.seam(node, 0, 0.3, 0.4)
        panel_shape = self.shape(node, "outputPatches[1].surface")
        panels = self.snapshot(node)
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(panel_reference=value):
                cmds.setAttr(node + ".referenceMaterialHeight", value)
                self.errors.clear()
                self.assertEqual(self.snapshot(node), panels)
                self.assertFalse(self.errors)
                self.assert_invalid_preserves(node, "outputReferenceHeight", 10, [panel_shape])
        cmds.setAttr(node + ".referenceMaterialHeight", 1e-8)
        self.seam(node, 0, float("nan"), 0)
        cmds.setAttr(driver + ".input1", 2)
        self.errors.clear()
        self.assertEqual(cmds.getAttr(node + ".outputReferenceHeight"), 1e-8)
        self.assertFalse(self.errors)
        self.assertTrue(cmds.isDirty(node + ".outputSurface", datablock=True))
        self.assertTrue(cmds.isDirty(node + ".outputPatches", datablock=True))

    def test_logical_index_limit(self):
        for kind in ("seams", "panelHems"):
            with self.subTest(kind=kind):
                node, _, _ = self.rig()
                self.seam(node, 0, 0.3, 0.4)
                shape = self.shape(node, "outputPatches[1].surface")
                before = self.snapshot(node)
                if kind == "seams":
                    self.seam(node, 2147483647, 0.7, 0.4)
                else:
                    self.hem(node, 2147483647, (), ())
                self.assert_invalid_preserves(node, "outputPatches", 15, [shape], kind)
                cmds.removeMultiInstance(node + ".{}[2147483647]".format(kind), b=True)
                self.assertEqual(self.snapshot(node), before)

    def test_seam_spacing_and_height_clusters(self):
        for spacing in (0.5e-9, 1e-9, 1.0001e-9):
            for wrap in (False, True):
                with self.subTest(spacing=spacing, wrap=wrap):
                    node, _, _ = self.rig()
                    self.seam(node, 0, 0, 0.5)
                    shape = self.shape(node, "outputPatches[1].surface")
                    u = 1 - spacing if wrap else spacing
                    self.seam(node, 8, u, 0.5 + 0.75e-9)
                    actual_spacing = 1 - u if wrap else u
                    if actual_spacing <= 1e-9:
                        self.assert_invalid_preserves(node, "outputPatches", 5, [shape])
                    else:
                        self.assertEqual(self.ids(node), (1, 9))
                        self.assertEqual(self.metadata(node, 1)[2], (0.5,))
                        for index in self.ids(node):
                            self.assert_open(self.patch(node, index))
        node, _, _ = self.rig()
        for index, u, height in ((0, 0, 0.5), (4, 0.3, 0.5 + 0.75e-9), (9, 0.7, 0.5 + 1.5e-9)):
            self.seam(node, index, u, height)
        self.assertEqual(self.metadata(node, 1)[2], (0.5, 0.5 + 1.5e-9))
        node, _, _ = self.rig()
        physical = 4.0 / 11.5
        self.seam(node, 0, 0.3, physical - 0.75e-9)
        self.seam(node, 8, 0.8, physical)
        self.assertEqual(self.metadata(node, 1)[2], (physical,))
        self.assertEqual(tuple(self.patch(node, 1).knotsInV()).count(physical), 1)

    def test_narrow_cubic_panel(self):
        node, _, _ = self.rig()
        cmds.setAttr(node + ".bellSubdivision", 3)
        self.seam(node, 0, 0.3, 0)
        self.seam(node, 8, 0.3 + 2e-9, 0)
        for index in self.ids(node):
            self.assert_open(self.patch(node, index))

    def test_zero_length_physical_stage(self):
        for request in ("outputSurface", "outputPatches"):
            with self.subTest(request=request):
                node, _, joints = self.rig()
                shape = self.shape(node, "outputSurface")
                before = point_hash(self.shape_points(shape))
                for side in ("left", "right"):
                    cmds.xform(joints[side + "Knee"], worldSpace=True, translation=(-1 if side == "left" else 1, 10, 0))
                self.assert_invalid_preserves(node, request, 14, [shape])
                for side in ("left", "right"):
                    cmds.xform(joints[side + "Knee"], worldSpace=True, translation=(-1 if side == "left" else 1, 5, 0))
                self.assertEqual(point_hash(self.shape_points(shape)), before)

    def test_sparse_arrays_and_persistence(self):
        node, _, _ = self.rig(angle=60)
        self.seam(node, 0, 0.1, 0.3)
        self.seam(node, 19, 0.7, 0.3)
        self.hem(node, 1, (0, 1), (0.9, 0.85))
        self.hem(node, 20, (0, 1), (0.8, 0.9))
        cmds.setAttr(node + ".followRange", 0.125)
        cmds.setAttr(node + ".referenceMaterialHeight", 13.5)
        shapes = [self.shape(node, "outputPatches[{}].surface".format(i)) for i in (1, 20)]
        expected = self.snapshot(node)
        for request in ("outputPatches", "outputPatches[20]", "outputPatches[1].surface", "outputPatches[20].materialUStart"):
            cmds.dgdirty(node)
            self.request(node, request)
            self.assertEqual(self.snapshot(node), expected)
        self.warnings.clear()
        for index in (19, 0, 19):
            cmds.setAttr(node + ".seams[{}].enabled".format(index), False)
            self.assertEqual(self.ids(node), (1,) if index == 19 else (20,))
            cmds.setAttr(node + ".seams[{}].materialU".format(index), float("nan"))
            self.ids(node)
            self.seam(node, index, 0.7 if index == 19 else 0.1, 0.3)
            self.assertEqual(self.snapshot(node), expected)
        warnings = [m for m in self.warnings if node in m and "panelHems" in m]
        self.assertEqual(len(warnings), 1, warnings)
        cmds.removeMultiInstance(node + ".seams[19]", b=True)
        self.assertEqual(self.ids(node), (1,))
        self.seam(node, 27, 0.7, 0.3)
        self.assertEqual(self.ids(node), (1, 28))
        cmds.removeMultiInstance(node + ".seams[27]", b=True)
        self.seam(node, 19, 0.7, 0.3)
        self.assertEqual(self.snapshot(node), expected)
        connections = [cmds.connectionInfo(s + ".create", sourceFromDestination=True) for s in shapes]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "panels.ma").replace(os.sep, "/")
            cmds.file(rename=path)
            cmds.file(save=True, force=True, type="mayaAscii")
            cmds.file(path, open=True, force=True)
            self.assertEqual(self.snapshot(node), expected)
            self.assertEqual(cmds.getAttr(node + ".followRange"), 0.125)
            self.assertEqual(cmds.getAttr(node + ".outputReferenceHeight"), 13.5)
            self.assertEqual([cmds.connectionInfo(s + ".create", sourceFromDestination=True) for s in shapes], connections)
            self.assertEqual(cmds.getAttr(node + ".seams", multiIndices=True), [0, 19])
            self.assertEqual(cmds.getAttr(node + ".panelHems", multiIndices=True), [1, 20])
            for index in (0, 19):
                cmds.setAttr(node + ".seams[{}].enabled".format(index), False)
            self.assert_reserved_patch(node)
            for index in (0, 19):
                cmds.setAttr(node + ".seams[{}].enabled".format(index), True)
            self.assertEqual(self.snapshot(node), expected)
            cmds.file(new=True, force=True)

    def test_evaluation_modes_and_time_order(self):
        node, _, joints = self.rig()
        self.seam(node, 0, 0.3, 0.4)
        self.seam(node, 8, 0.8, 0.4)
        for index in (1, 9):
            self.shape(node, "outputPatches[{}].surface".format(index))
        for time, value in ((0, 0), (12, -65)):
            cmds.setKeyframe(
                joints["leftHip"] + ".rotateX", time=time, value=value, inTangentType="linear", outTangentType="linear"
            )
        references = {}
        for mode in ("off", "serial", "parallel"):
            cmds.evaluationManager(mode=mode)
            for start in (0.4, 0.6, 0.4):
                cmds.setAttr(node + ".seams[8].startHeight", start)
                for time in (0, 4, 12, 7, 0):
                    with self.subTest(mode=mode, start=start, time=time):
                        cmds.currentTime(time)
                        self.request(node, "outputPatches[9].surface")
                        result = self.snapshot(node)
                        key = (start, time)
                        if key in references:
                            self.assertEqual(result, references[key])
                        else:
                            references[key] = result
                        for index in (1, 9):
                            for child in (
                                "surface",
                                "materialUStart",
                                "materialUEnd",
                                "vBreaks",
                                "hemUSamples",
                                "hemHeightSamples",
                            ):
                                self.assertFalse(
                                    cmds.isDirty(node + ".outputPatches[{}].{}".format(index, child), datablock=True)
                                )
                if mode != "off":
                    self.assertFalse(cmds.evaluationManager(query=True, empty=True))
                    self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))

    def test_seamless_evaluation_modes_and_request_order(self):
        node, _, joints = self.rig()
        attributes = ("outputSurface", "outputPatches[0].surface")
        shapes = [self.shape(node, attribute) for attribute in attributes]
        for time, value in ((0, 0), (12, -65)):
            cmds.setKeyframe(
                joints["leftHip"] + ".rotateX", time=time, value=value, inTangentType="linear", outTangentType="linear"
            )
        references = {}
        for mode in ("off", "serial", "parallel"):
            cmds.evaluationManager(mode=mode)
            for time in (0, 4, 12, 7, 0):
                cmds.currentTime(time)
                for order in (attributes, attributes[::-1]):
                    with self.subTest(mode=mode, time=time, order=order):
                        cmds.dgdirty(node)
                        self.request(node, order[0])
                        # Maya reports the array parent and element plugs as dirty in
                        # the datablock even after every child was set clean.
                        for attribute in (
                            "outputSurface",
                            "outputPatches[0].surface",
                            "outputPatches[0].materialUStart",
                            "outputPatches[0].materialUEnd",
                            "outputPatches[0].vBreaks",
                            "outputPatches[0].hemUSamples",
                            "outputPatches[0].hemHeightSamples",
                        ):
                            self.assertFalse(cmds.isDirty(node + "." + attribute, datablock=True), attribute)
                        self.request(node, order[1])
                        self.assert_reserved_patch(node)
                        result = tuple(point_hash(self.shape_points(shape)) for shape in shapes)
                        if time in references:
                            self.assertEqual(result, references[time])
                        else:
                            references[time] = result
            if mode != "off":
                self.assertFalse(cmds.evaluationManager(query=True, empty=True))
                self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))

    def test_parallel_warn_once(self):
        node, _, joints = self.rig()
        for index, u in ((0, 0.1), (8, 0.35), (19, 0.7)):
            self.seam(node, index, u, 0.4)
            self.shape(node, "outputPatches[{}].surface".format(index + 1))
        for time, value in ((0, 0), (12, -65)):
            cmds.setKeyframe(joints["leftHip"] + ".rotateX", time=time, value=value)
        self.hem(node, 99, (float("nan"),), ())
        self.warnings.clear()
        cmds.evaluationManager(mode="parallel")
        for time in (0, 4, 12, 7, 0):
            cmds.currentTime(time)
            self.assertEqual(self.ids(node), (1, 9, 20))
            self.assertEqual(len(self.snapshot(node)), 3)
            cmds.dgdirty(node)
        self.assertFalse(cmds.evaluationManager(query=True, empty=True))
        self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))
        warnings = [m for m in self.warnings if node in m and "panelHems" in m]
        self.assertEqual(len(warnings), 1, warnings)

    def test_evaluation_solve_count(self):
        self.skipTest("The node exposes no solve/build counter; clean plugs and complete arrays cannot prove exactly one solve")

    def test_prepare_for_draw_patch_data_path(self):
        references = {}
        for mode in ("off", "serial", "parallel"):
            cmds.evaluationManager(mode=mode)
            node, _, _ = self.rig()
            cmds.setAttr(node + ".nodeState", 1)
            old_shape = self.shape(node, "outputSurface")
            old_points = point_hash(self.shape_points(old_shape))
            cmds.disconnectAttr(node + ".outputSurface", old_shape + ".create")
            self.seam(node, 0, 0.1, 0)
            self.seam(node, 8, 0.3, 0)
            self.errors.clear()
            elements = []
            for index in self.ids(node):
                surface = self.patch(node, index)
                points = tuple(tuple(p) for p in surface.cvPositions())
                nu, nv = surface.numCVsInU, surface.numCVsInV
                lines = tuple((points[(u - 1) * nv + v], points[u * nv + v]) for v in range(nv) for u in range(1, nu))
                self.assertEqual(len(lines), (nu - 1) * nv)
                elements.append((index, nu, nv, point_hash(points), lines))
            self.assertFalse(self.errors)
            self.assertNotEqual(elements[0][1], elements[1][1])
            self.assertEqual(point_hash(self.shape_points(old_shape)), old_points)
            if references:
                self.assertEqual(elements, references["elements"])
            else:
                references["elements"] = elements
            self.assertTrue(cmds.isDirty(node + ".outputSurface", datablock=True))
            self.assert_invalid_preserves(node, "outputSurface", 11, [])
            shapes = [self.shape(node, "outputPatches[{}].surface".format(i)) for i in (1, 9)]
            saved = self.snapshot(node)
            cmds.setAttr(node + ".seams[0].materialU", -1)
            self.assert_invalid_preserves(node, "outputPatches", 4, shapes)
            cmds.setAttr(node + ".seams[0].materialU", 0.1)
            self.assertEqual(self.snapshot(node), saved)
            cmds.removeMultiInstance(node + ".seams[0]", b=True)
            self.seam(node, 19, 0.1, 0)
            renamed = self.snapshot(node)
            self.assertEqual(tuple(p[0] for p in renamed), (9, 20))
            self.assertEqual(next(p[1:] for p in renamed if p[0] == 20), saved[0][1:])
