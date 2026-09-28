"""Surface fitting topology, geometry, validation, and evaluation tests."""

import hashlib
import math
import os
import struct
import tempfile
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import coordinates, output_object, skirt


def distance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def point_hash(points):
    digest = hashlib.sha256()
    for point in points:
        digest.update(struct.pack("<4d", *point))
    return digest.hexdigest()


class SurfaceFitTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")
        self.errors = []
        self.surface_data = []

        def capture(message, message_type, _client_data):
            if message_type == om.MCommandMessage.kError:
                self.errors.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)

    def fitted_skirt(self, skirt_type=1, spans=4):
        cmds.file(new=True, force=True)
        collider, _, surface = skirt()
        cmds.setAttr(collider + ".skirtType", skirt_type)
        fit = cmds.createNode("yddSkirtSurfaceFit")
        cmds.setAttr(fit + ".spansV", spans)
        cmds.connectAttr(collider + ".outputSurface", fit + ".inputSurface")
        shape = cmds.listRelatives(surface, shapes=True, fullPath=True)[0]
        cmds.connectAttr(fit + ".outputSurface", shape + ".create", force=True)
        return collider, fit

    def points(self, node):
        return coordinates(node, "outputSurface", "surface")

    def shape_points(self, fit):
        shapes = cmds.listConnections(
            fit + ".outputSurface", source=False, destination=True, shapes=True, type="nurbsSurface"
        ) or []
        self.assertEqual(len(shapes), 1)
        dag_path = om.MSelectionList().add(shapes[0]).getDagPath(0)
        points = om.MFnNurbsSurface(dag_path).cvPositions(om.MSpace.kObject)
        return [[p.x, p.y, p.z, p.w] for p in points]

    def surface(self, node):
        data = output_object(node, "outputSurface")
        # The function set does not own the data object; keep it alive for the test.
        self.surface_data.append(data)
        return om.MFnNurbsSurface(data)

    def sampler(self, node):
        sampler = cmds.createNode("pointOnSurfaceInfo")
        cmds.connectAttr(node + ".outputSurface", sampler + ".inputSurface")
        cmds.setAttr(sampler + ".turnOnPercentage", False)
        return sampler

    def sample(self, sampler, u, v, attribute="position"):
        cmds.setAttr(sampler + ".parameterU", u)
        cmds.setAttr(sampler + ".parameterV", v)
        return cmds.getAttr(sampler + "." + attribute)[0]

    def assert_finite(self, points):
        self.assertTrue(points)
        self.assertTrue(all(math.isfinite(value) for point in points for value in point))

    def assert_invalid_preserves(self, fit, expected, reason, repeats=3):
        nodes = set(cmds.ls(long=True))
        connections = cmds.listConnections(fit, connections=True, plugs=True) or []
        spans = cmds.getAttr(fit + ".spansV")
        for attempt in range(repeats):
            with self.subTest(attempt=attempt):
                self.errors.clear()
                try:
                    output_object(fit, "outputSurface")
                except RuntimeError:
                    pass
                points = self.shape_points(fit)
                self.assertTrue(
                    any("yddSkirtSurfaceFit:" in message and reason in message for message in self.errors),
                    self.errors,
                )
                self.assertEqual(points, expected)
                self.assertEqual(cmds.getAttr(fit + ".spansV"), spans)
                self.assertEqual(set(cmds.ls(long=True)), nodes)
                self.assertEqual(cmds.listConnections(fit, connections=True, plugs=True) or [], connections)

    def test_output_topology_and_knots(self):
        for skirt_type in (1, 0):
            collider, fit = self.fitted_skirt(skirt_type)
            source = self.surface(collider)
            for spans in (1, 4, 9):
                with self.subTest(skirt_type=skirt_type, spans=spans):
                    cmds.setAttr(fit + ".spansV", spans)
                    result = self.surface(fit)
                    self.assertEqual(result.numCVsInV, spans + 3)
                    self.assertEqual(result.degreeInV, 3)
                    self.assertEqual(result.formInV, om.MFnNurbsSurface.kOpen)
                    self.assertEqual(result.numCVsInU, source.numCVsInU)
                    self.assertEqual(result.degreeInU, source.degreeInU)
                    self.assertEqual(result.formInU, source.formInU)
                    self.assertEqual(list(result.knotsInU()), list(source.knotsInU()))
                    self.assertEqual(result.knotDomainInU, source.knotDomainInU)
                    self.assertEqual(tuple(result.knotDomainInV), (0.0, 1.0))
                    self.assertEqual(
                        list(result.knotsInV()),
                        [0.0] * 3 + [i / float(spans) for i in range(1, spans)] + [1.0] * 3,
                    )
                    self.assert_finite(self.points(fit))

    def test_convex_bounds_and_endpoints(self):
        for skirt_type in (1, 0):
            collider, fit = self.fitted_skirt(skirt_type)
            source = self.surface(collider)
            self.assertEqual(source.degreeInV, 1)
            before = self.points(collider)
            source_sampler = self.sampler(collider)
            output_sampler = self.sampler(fit)
            umin, umax = source.knotDomainInU
            vmin, vmax = source.knotDomainInV
            for spans in (1, 4, 9):
                with self.subTest(skirt_type=skirt_type, spans=spans):
                    cmds.setAttr(fit + ".spansV", spans)
                    after = self.points(fit)
                    m = spans + 3
                    nv = source.numCVsInV
                    for u in range(source.numCVsInU):
                        for component in range(3):
                            values = [point[component] for point in before[u * nv:(u + 1) * nv]]
                            for point in after[u * m:(u + 1) * m]:
                                tolerance = 1e-9 * max(1.0, abs(min(values)), abs(max(values)))
                                self.assertGreaterEqual(point[component], min(values) - tolerance)
                                self.assertLessEqual(point[component], max(values) + tolerance)
                    for step in range(9):
                        u = umin + (umax - umin) * step / 8.0
                        for output_v, input_v in ((0, vmin), (1, vmax)):
                            actual = self.sample(output_sampler, u, output_v)
                            expected = self.sample(source_sampler, u, input_v)
                            self.assertLessEqual(distance(actual, expected), 1e-6)

    def test_shape_difference(self):
        collider, fit = self.fitted_skirt()
        source = self.surface(collider)
        umin, umax = source.knotDomainInU
        vmin, vmax = source.knotDomainInV
        source_sampler = self.sampler(collider)
        output_sampler = self.sampler(fit)
        leg_lengths = []
        for side in ("left", "right"):
            positions = []
            for part in ("Hip", "Knee", "Heel"):
                matrix = cmds.getAttr(collider + "." + side + part + "Matrix")
                positions.append(matrix[12:15])
            leg_lengths.append(distance(positions[0], positions[1]) + distance(positions[1], positions[2]))
        limit = 0.1 * sum(leg_lengths) / 2.0
        maximum = 0.0
        for ui in range(33):
            u = umin + (umax - umin) * ui / 32.0
            for vi in range(33):
                fraction = vi / 32.0
                actual = self.sample(output_sampler, u, fraction)
                expected = self.sample(source_sampler, u, vmin + fraction * (vmax - vmin))
                delta = distance(actual, expected)
                self.assertTrue(math.isfinite(delta))
                maximum = max(maximum, delta)
        print("yddSkirtSurfaceFit maximum shape difference: %.12g; limit: %.12g" % (maximum, limit))
        self.assertLessEqual(maximum, limit)

    def test_periodic_seam(self):
        for skirt_type in (1, 0):
            _, fit = self.fitted_skirt(skirt_type)
            sampler = self.sampler(fit)
            for spans in (1, 4, 9):
                cmds.setAttr(fit + ".spansV", spans)
                umin, umax = self.surface(fit).knotDomainInU
                for v in (0.0, 0.5, 1.0):
                    with self.subTest(skirt_type=skirt_type, spans=spans, v=v):
                        first = self.sample(sampler, umin, v)
                        last = self.sample(sampler, umax, v)
                        self.assertLessEqual(distance(first, last), 1e-9)
                        first = self.sample(sampler, umin, v, "tangentU")
                        last = self.sample(sampler, umax, v, "tangentU")
                        first_length = distance(first, (0, 0, 0))
                        last_length = distance(last, (0, 0, 0))
                        self.assertGreater(first_length, 0)
                        self.assertGreater(last_length, 0)
                        self.assertLessEqual(
                            distance([value / first_length for value in first], [value / last_length for value in last]),
                            1e-6,
                        )

    def test_invalid_input_and_recovery(self):
        collider, fit = self.fitted_skirt(spans=1)
        for invalid in (0, 257):
            with self.subTest(spans=invalid):
                cmds.setAttr(fit + ".spansV", 1)
                before = self.shape_points(fit)
                self.assert_finite(before)
                self.assertEqual(len(before), self.surface(collider).numCVsInU * 4)
                cmds.setAttr(fit + ".spansV", invalid)
                self.assert_invalid_preserves(fit, before, "spansV must be between 1 and 256")
                cmds.setAttr(fit + ".spansV", 4)
                after = self.shape_points(fit)
                self.assert_finite(after)
                self.assertNotEqual(after, before)
                self.assertEqual(len(after), self.surface(collider).numCVsInU * 7)
        before = self.shape_points(fit)
        cmds.disconnectAttr(collider + ".outputSurface", fit + ".inputSurface")
        self.assert_invalid_preserves(fit, before, "inputSurface must be a readable NURBS surface")
        cmds.connectAttr(collider + ".outputSurface", fit + ".inputSurface")
        after = self.shape_points(fit)
        self.assert_finite(after)
        self.assertEqual(after, before)
        # Two violated conditions: the lower-numbered one (spansV) is reported.
        cmds.setAttr(fit + ".spansV", 300)
        cmds.disconnectAttr(collider + ".outputSurface", fit + ".inputSurface")
        self.assert_invalid_preserves(fit, before, "spansV must be between 1 and 256")
        cmds.setAttr(fit + ".spansV", 4)
        cmds.connectAttr(collider + ".outputSurface", fit + ".inputSurface")
        self.assertEqual(self.shape_points(fit), before)

    def test_closed_v_input_is_rejected(self):
        collider, fit = self.fitted_skirt(spans=4)
        before = self.shape_points(fit)
        plane = cmds.nurbsPlane(axis=(0, 1, 0), width=4, lengthRatio=2, patchesU=3, patchesV=5, degree=3)[0]
        closed = cmds.closeSurface(plane, direction=1, constructionHistory=False, replaceOriginal=True)[0]
        closed_shape = cmds.listRelatives(closed, shapes=True, fullPath=True)[0]
        self.assertEqual(om.MFnNurbsSurface(om.MSelectionList().add(closed_shape).getDagPath(0)).formInV,
                         om.MFnNurbsSurface.kPeriodic)
        cmds.connectAttr(closed_shape + ".local", fit + ".inputSurface", force=True)
        self.assert_invalid_preserves(fit, before, "V form must be open")
        cmds.connectAttr(collider + ".outputSurface", fit + ".inputSurface", force=True)
        self.assertEqual(self.shape_points(fit), before)

    def test_minimum_height(self):
        for skirt_type in (1, 0):
            collider, fit = self.fitted_skirt(skirt_type)
            cmds.setAttr(collider + ".height", 0.01)
            for spans in (1, 4, 9):
                with self.subTest(skirt_type=skirt_type, spans=spans):
                    cmds.setAttr(fit + ".spansV", spans)
                    points = self.points(fit)
                    self.assert_finite(points)
                    self.assertEqual(len(points), self.surface(collider).numCVsInU * (spans + 3))

    def test_evaluation_modes_and_time_order(self):
        collider, fit = self.fitted_skirt()
        hip = cmds.listConnections(collider + ".leftHipMatrix", source=True, destination=False)[0]
        for time, value in ((0, 0), (12, 65)):
            cmds.setKeyframe(hip + ".rotateY", time=time, value=value, inTangentType="linear", outTangentType="linear")
        references = {}
        try:
            for mode in ("off", "serial", "parallel"):
                cmds.evaluationManager(mode=mode)
                for time in (0, 4, 12, 7, 0):
                    with self.subTest(mode=mode, time=time):
                        cmds.currentTime(time)
                        points = self.points(fit)
                        self.assert_finite(points)
                        result = point_hash(points)
                        if time in references:
                            self.assertEqual(result, references[time])
                        else:
                            references[time] = result
                if mode != "off":
                    self.assertFalse(cmds.evaluationManager(query=True, empty=True))
                    self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))
        finally:
            cmds.evaluationManager(mode="off")

    def test_save_and_reload(self):
        collider, fit = self.fitted_skirt(spans=9)
        expected = self.points(fit)
        source = cmds.connectionInfo(fit + ".inputSurface", sourceFromDestination=True)
        destinations = cmds.listConnections(fit + ".outputSurface", source=False, destination=True, plugs=True)
        directory = tempfile.mkdtemp()
        try:
            path = os.path.join(directory, "surface_fit.ma").replace(os.sep, "/")
            cmds.file(rename=path)
            cmds.file(save=True, force=True, type="mayaAscii")
            cmds.file(new=True, force=True)
            cmds.file(path, open=True, force=True)
            self.assertEqual(cmds.getAttr(fit + ".spansV"), 9)
            self.assertEqual(cmds.connectionInfo(fit + ".inputSurface", sourceFromDestination=True), source)
            self.assertTrue(cmds.isConnected(collider + ".outputSurface", fit + ".inputSurface"))
            self.assertEqual(
                cmds.listConnections(fit + ".outputSurface", source=False, destination=True, plugs=True), destinations
            )
            self.assertEqual(self.points(fit), expected)
        finally:
            cmds.file(new=True, force=True)
            for name in os.listdir(directory):
                os.remove(os.path.join(directory, name))
            os.rmdir(directory)
