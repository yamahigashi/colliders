"""Wave height normalization tests for yddSkirtWaveDeformer."""

import hashlib
import math
import struct
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import wave_settings


def point_hash(points):
    digest = hashlib.sha256()
    for point in points:
        digest.update(struct.pack("<3d", *point[:3]))
    return digest.hexdigest()


class WaveHeightTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")
        self.errors = []

        def capture(message, message_type, _client_data):
            if message_type == om.MCommandMessage.kError:
                self.errors.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)

    def mesh(self, points):
        fn = om.MFnMesh()
        obj = fn.create(om.MPointArray([om.MPoint(*point) for point in points]), [len(points)], list(range(len(points))))
        return om.MFnDagNode(obj).fullPathName()

    def points(self, mesh):
        shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
        path = om.MSelectionList().add(shape).getDagPath(0)
        return [tuple(p)[:3] for p in om.MFnMesh(path).getPoints()]

    def set_point(self, mesh, index, xyz):
        shapes = cmds.listRelatives(mesh, shapes=True, fullPath=True)
        intermediate = [s for s in shapes if cmds.getAttr(s + ".intermediateObject")]
        shape = (intermediate or shapes)[0]
        path = om.MSelectionList().add(shape).getDagPath(0)
        om.MFnMesh(path).setPoint(index, om.MPoint(*xyz), om.MSpace.kObject)

    def set_short(self, node, attr, value):
        om.MSelectionList().add(node + "." + attr).getPlug(0).setShort(value)

    def set_double(self, node, attr, value):
        om.MSelectionList().add(node + "." + attr).getPlug(0).setDouble(value)

    def wave(self, mesh, condition="active"):
        node = cmds.deformer(mesh, type="yddSkirtWaveDeformer")[0]
        wave_settings(node, condition)
        return node

    def assert_error_preserves(self, mesh, node, expected, reason, repeats=3):
        for attempt in range(repeats):
            with self.subTest(attempt=attempt):
                self.errors.clear()
                try:
                    cmds.dgeval(node + ".outputGeometry[0]")
                except RuntimeError:
                    pass
                self.assertTrue(
                    any("yddSkirtWaveDeformer" in message and reason in message for message in self.errors),
                    self.errors,
                )
                self.assertEqual(self.points(mesh), expected)

    def test_current_geometry_far_side_follows_max_height(self):
        mesh = self.mesh([(-1.0, 1.0, 0.3), (-1.0, 0.4, -0.3), (1.0, 1.0, 0.3), (1.0, 1.0, -0.3)])
        node = self.wave(mesh)
        cmds.setAttr(node + ".heightNormalization", 0)
        baseline = self.points(mesh)
        far_hash = point_hash(baseline[:2])
        self.set_point(mesh, 3, (1.0, 4.0, -0.3))
        raised = self.points(mesh)
        self.assertNotEqual(point_hash(raised[:2]), far_hash)

        self.set_point(mesh, 3, (1.0, 1.0, -0.3))
        cmds.setAttr(node + ".weightList[0].weights[3]", 0)
        weighted_base = self.points(mesh)
        far_hash = point_hash(weighted_base[:2])
        self.set_point(mesh, 3, (1.0, 4.0, -0.3))
        weighted_raised = self.points(mesh)
        self.assertNotEqual(point_hash(weighted_raised[:2]), far_hash)
        for actual, wanted in zip(weighted_raised[3], (1.0, 4.0, -0.3)):
            self.assertAlmostEqual(actual, wanted, places=6)

    def test_reference_height_keeps_far_side(self):
        mesh = self.mesh([(-1.0, 1.0, 0.3), (-1.0, 0.4, -0.3), (1.0, 1.0, 0.3), (1.0, 1.0, -0.3)])
        node = self.wave(mesh)
        cmds.setAttr(node + ".heightNormalization", 1)
        cmds.setAttr(node + ".referenceHeight", 1)
        baseline = self.points(mesh)
        far_hash = point_hash(baseline[:2])
        self.set_point(mesh, 3, (1.0, 4.0, -0.3))
        raised = self.points(mesh)
        self.assertEqual(point_hash(raised[:2]), far_hash)

    def test_invalid_height_inputs_preserve_and_recover(self):
        mesh = self.mesh([(1.0, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1.0, 1.0)])
        undeformed = self.points(mesh)
        node = self.wave(mesh)
        cmds.setAttr(node + ".heightNormalization", 1)
        cmds.setAttr(node + ".referenceHeight", 1)
        valid = self.points(mesh)
        self.assertNotEqual(valid, undeformed)

        self.set_short(node, "heightNormalization", 2)
        self.assert_error_preserves(mesh, node, undeformed, "heightNormalization must be 0 or 1")
        self.set_short(node, "heightNormalization", 1)

        for value in (0.0, -1.0, float("inf"), float("nan")):
            with self.subTest(referenceHeight=value):
                cmds.setAttr(node + ".referenceHeight", 1)
                self.assertEqual(self.points(mesh), valid)
                self.set_double(node, "referenceHeight", value)
                self.assert_error_preserves(mesh, node, undeformed, "referenceHeight must be finite and positive")

        cmds.setAttr(node + ".referenceHeight", 1)
        recovered = self.points(mesh)
        self.assertEqual(recovered, valid)

        cmds.setAttr(node + ".amplitude", 0)
        self.set_double(node, "referenceHeight", -2.0)
        self.assert_error_preserves(mesh, node, undeformed, "referenceHeight must be finite and positive")
        cmds.setAttr(node + ".amplitude", 1)

        cmds.setAttr(node + ".envelope", 0)
        self.set_double(node, "referenceHeight", 0.0)
        self.assert_error_preserves(mesh, node, undeformed, "referenceHeight must be finite and positive")
        cmds.setAttr(node + ".envelope", 1)
        cmds.setAttr(node + ".referenceHeight", 1)

        cmds.setAttr(node + ".bellMatrix", *[0.0] * 16, type="matrix")
        self.set_double(node, "referenceHeight", float("nan"))
        self.assert_error_preserves(mesh, node, undeformed, "referenceHeight must be finite and positive")
        cmds.setAttr(node + ".bellMatrix", *[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], type="matrix")
        cmds.setAttr(node + ".referenceHeight", 1)
        self.assertEqual(self.points(mesh), valid)

        cmds.setAttr(node + ".heightNormalization", 0)
        self.set_double(node, "referenceHeight", -8.0)
        self.errors.clear()
        current = self.points(mesh)
        self.assertTrue(all("referenceHeight must be finite and positive" not in message for message in self.errors))
        self.assertNotEqual(current, [(1.0, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1.0, 1.0)])

    def test_invalid_mode_precedes_early_returns(self):
        identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        for condition in ("matrix", "amplitude", "envelope"):
            with self.subTest(condition=condition):
                mesh = self.mesh([(1.0, 0.2, 0.1), (0.7, 0.5, 0.5), (0.1, 1.0, 1.0)])
                undeformed = self.points(mesh)
                node = self.wave(mesh)
                valid = self.points(mesh)
                self.assertNotEqual(valid, undeformed)
                if condition == "matrix":
                    cmds.setAttr(node + ".bellMatrix", *[0.0] * 16, type="matrix")
                else:
                    cmds.setAttr(node + "." + condition, 0)
                self.set_short(node, "heightNormalization", 2)
                self.set_double(node, "referenceHeight", -1.0)
                self.assert_error_preserves(
                    mesh, node, undeformed, "heightNormalization must be 0 or 1"
                )
                self.assertFalse(any("referenceHeight must" in message for message in self.errors))
                self.set_short(node, "heightNormalization", 0)
                self.errors.clear()
                self.assertEqual(self.points(mesh), undeformed)
                self.assertFalse(self.errors)
                if condition == "matrix":
                    cmds.setAttr(node + ".bellMatrix", *identity, type="matrix")
                else:
                    cmds.setAttr(node + "." + condition, 1)
                self.assertEqual(self.points(mesh), valid)

    def test_reference_height_below_minimum_is_not_clamped(self):
        mesh = self.mesh([(1.0, 5e-6, 0.0), (1.0, 5e-7, 0.0), (1.0, 5e-7, 0.5)])
        original = self.points(mesh)
        node = self.wave(mesh, "impulse_only")
        cmds.setAttr(node + ".heightNormalization", 1)
        cmds.setAttr(node + ".referenceHeight", 1e-6)
        cmds.setAttr(node + ".impulsePosition", 0.5)
        cmds.setAttr(node + ".impulseWidth", 0.25)
        result = self.points(mesh)
        self.assertEqual(result[0], original[0])
        self.assertNotEqual(result[1], original[1])
