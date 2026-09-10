"""Maya regression tests for surface attachments and plane output evaluation."""

import unittest

import maya.cmds as cmds

import yddColliders


class SurfaceAttachmentTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def assertVectorAlmostEqual(self, actual, expected):
        for value, reference in zip(actual, expected):
            self.assertAlmostEqual(value, reference, places=5)

    def create_skirt(self, skirt_type):
        joints = {}
        for side, x in (("left", -1), ("right", 1)):
            for part, y in (("Hip", 10), ("Knee", 5), ("Heel", 0)):
                obj = cmds.createNode("transform", name=side + part)
                cmds.setAttr(obj + ".translate", x, y, 0)
                joints[side + part + "Obj"] = obj
        node, _, surface = yddColliders.createSkirtBellCollider(skirtType=skirt_type, **joints)
        return node, surface

    def check_attachments(self, surface, joints, u_num, v_num):
        shape = cmds.listRelatives(surface, shapes=True, fullPath=True)[0]
        # Query the live domain independently of pointOnSurfaceInfo percentage mode.
        min_u, max_u = cmds.getAttr(shape + ".minMaxRangeU")[0]
        min_v, max_v = cmds.getAttr(shape + ".minMaxRangeV")[0]
        for i in range(u_num):
            for j in range(v_num):
                u = i / float(u_num)
                v = j / float(v_num - 1) if v_num > 1 else 0.5
                expected = cmds.pointOnSurface(
                    shape,
                    parameterU=min_u + u * (max_u - min_u),
                    parameterV=min_v + v * (max_v - min_v),
                    position=True,
                )
                actual = cmds.xform(
                    joints[i * v_num + j],
                    query=True,
                    worldSpace=True,
                    translation=True,
                )
                self.assertVectorAlmostEqual(actual, expected)
                self.assertGreater(sum(value * value for value in actual), 0.01)

    def test_grid_follows_initial_shortened_and_lengthened_surface(self):
        for skirt_type in (0, 1):
            with self.subTest(skirt_type=skirt_type):
                cmds.file(new=True, force=True)
                node, surface = self.create_skirt(skirt_type)
                joints = yddColliders.attachJointsToSurface(surface, 3, 3)
                self.assertEqual(len(joints), 9)
                for i in range(3):
                    for j in range(3):
                        posi = "skirt_{}_{}_posi".format(i, j)
                        self.assertTrue(cmds.getAttr(posi + ".turnOnPercentage"))
                        self.assertAlmostEqual(cmds.getAttr(posi + ".parameterU"), i / 3.0)
                        self.assertAlmostEqual(cmds.getAttr(posi + ".parameterV"), j / 2.0)
                for height in (1.0, 0.1, 1.0):
                    cmds.setAttr(node + ".height", height)
                    self.check_attachments(surface, joints, 3, 3)

    def test_single_v_stays_at_middle_after_height_changes(self):
        node, surface = self.create_skirt(1)
        joints = yddColliders.attachJointsToSurface(surface, 1, 1)
        self.assertEqual(cmds.getAttr("skirt_0_0_posi.parameterV"), 0.5)
        for height in (1.0, 0.1, 1.0):
            cmds.setAttr(node + ".height", height)
            self.check_attachments(surface, joints, 1, 1)


class PlaneOutputTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def test_parent_and_each_child_can_be_evaluated_first(self):
        for first in (None, 0, 1, 2):
            with self.subTest(first=first):
                node = cmds.createNode("yddPlaneCollider")
                cmds.setAttr(node + ".inputPosition", 1, -2, 3)
                expected = (1.0, 0.0, 3.0)
                if first is None:
                    self.assertEqual(cmds.getAttr(node + ".outputPosition")[0], expected)
                else:
                    self.assertEqual(cmds.getAttr(node + ".outputPosition" + str(first)), expected[first])
                self.assertEqual(cmds.getAttr(node + ".outputPosition")[0], expected)
                cmds.setAttr(node + ".inputPosition", 4, 5, 6)
                for index, value in enumerate((4.0, 5.0, 6.0)):
                    self.assertEqual(cmds.getAttr(node + ".outputPosition" + str(index)), value)

    def test_downstream_child_connections_evaluate_first_and_update(self):
        node = cmds.createNode("yddPlaneCollider")
        sink = cmds.createNode("multiplyDivide")
        for index, axis in enumerate("XYZ"):
            cmds.connectAttr(node + ".outputPosition" + str(index), sink + ".input1" + axis)
        for position, expected in (
            ((1, -2, 3), (1.0, 0.0, 3.0)),
            ((4, 5, 6), (4.0, 5.0, 6.0)),
            ((-3, -4, -5), (-3.0, 0.0, -5.0)),
        ):
            cmds.setAttr(node + ".inputPosition", *position)
            for axis, value in zip("XYZ", expected):
                self.assertEqual(cmds.getAttr(sink + ".output" + axis), value)
            self.assertEqual(cmds.getAttr(node + ".outputPosition")[0], expected)
