"""DG dependencies for display-only controls (viewport buffers use a C++ test)."""

import unittest

import maya.cmds as cmds
import maya.api.OpenMaya as om


class DrawDependencyTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def test_display_controls_do_not_affect_geometry(self):
        for node_type, attributes in (
            ("bellCollider", ("ringSubdivision", "drawColor", "drawOpacity")),
            ("skirtBellCollider", ("ringSubdivision",)),
        ):
            node = cmds.createNode(node_type)
            selection = om.MSelectionList()
            selection.add(node)
            fn = om.MFnDependencyNode(selection.getDependNode(0))
            for attribute in attributes:
                with self.subTest(node_type=node_type, attribute=attribute):
                    self.assertEqual(len(fn.getAffectedAttributes(fn.attribute(attribute))), 0)
