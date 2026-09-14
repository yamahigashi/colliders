"""Independent Maya registration contract for the internal ydd plugin identity."""

import unittest

import maya.cmds as cmds
import maya.api.OpenMaya as om

import yddColliders


class RegistrationTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def test_plugin_metadata_and_registered_nodes(self):
        self.assertEqual(cmds.pluginInfo("yddColliders", query=True, vendor=True), "yamahigashi")
        self.assertEqual(cmds.pluginInfo("yddColliders", query=True, version=True), "4.0.0")
        self.assertEqual(
            set(cmds.pluginInfo("yddColliders", query=True, dependNode=True)),
            {"yddBellCollider", "yddPlaneCollider", "yddSkirtBellCollider", "yddSkirtCollideDeformer", "yddSkirtWaveDeformer"},
        )

    def test_node_names_ids_kinds_and_draw_classification(self):
        for name, type_id, kind, classification in (
            ("yddBellCollider", 0x0007DD00, om.MFn.kLocator, "drawdb/geometry/yddBellCollider"),
            ("yddPlaneCollider", 0x0007DD01, om.MFn.kLocator, "drawdb/geometry/yddPlaneCollider"),
            ("yddSkirtBellCollider", 0x0007DD02, om.MFn.kLocator, "drawdb/geometry/yddSkirtBellCollider"),
            ("yddSkirtCollideDeformer", 0x0007DD03, om.MFn.kGeometryFilt, None),
            ("yddSkirtWaveDeformer", 0x0007DD04, om.MFn.kGeometryFilt, None),
        ):
            with self.subTest(node=name):
                node = cmds.createNode(name)
                obj = om.MSelectionList().add(node).getDependNode(0)
                fn = om.MFnDependencyNode(obj)
                self.assertEqual(fn.typeName, name)
                self.assertEqual(fn.typeId.id(), type_id)
                self.assertTrue(obj.hasFn(kind))
                if classification is not None:
                    self.assertIn(classification, om.MFnDependencyNode.classification(name).split(":"))

    def test_public_creator_uses_prefixed_custom_node_stems(self):
        self.assertEqual(yddColliders.WINDOW_ID, "yddBellColliderUI")
        node, _, _, _ = yddColliders.createBellCollider(prefix="hero_")
        self.assertEqual(cmds.nodeType(node), "yddBellCollider")
        self.assertEqual(node.rsplit("|", 1)[-1], "hero_yddBellCollider")
        joints = {}
        for side, x in (("left", -1), ("right", 1)):
            for part, y in (("Hip", 10), ("Knee", 5), ("Heel", 0)):
                joint = cmds.createNode("transform")
                cmds.setAttr(joint + ".translate", x, y, 0)
                joints[side + part + "Obj"] = joint
        node, _, _ = yddColliders.createSkirtBellCollider(prefix="hero_", **joints)
        self.assertEqual(cmds.nodeType(node), "yddSkirtBellCollider")
        self.assertEqual(node.rsplit("|", 1)[-1], "hero_yddSkirtBellCollider")
