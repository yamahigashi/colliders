"""Input-domain regression tests for collider and deformer compute paths.

Connected drivers bypass attribute limits so compute must validate its inputs.
"""

import math
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds


def matrix(x=0.0, y=0.0, z=0.0):
    return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]


def mesh_points(mesh):
    shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
    path = om.MSelectionList().add(shape).getDagPath(0)
    return [tuple(point)[:3] for point in om.MFnMesh(path).getPoints()]


def surface_points(node):
    plug = om.MSelectionList().add(node + ".outputSurface").getPlug(0)
    data = plug.asMObject()
    return [tuple(point) for point in om.MFnNurbsSurface(data).cvPositions()]


def bell_points(node, attribute):
    plug = om.MSelectionList().add(node + "." + attribute).getPlug(0)
    data = plug.asMObject()
    fn = om.MFnMesh(data) if attribute == "outputBellMesh" else om.MFnNurbsCurve(data)
    return [tuple(point)[:3] for point in (fn.getPoints() if attribute == "outputBellMesh" else fn.cvPositions())]


def finite_points(points):
    return points and all(math.isfinite(value) for point in points for value in point)


def evaluate(call):
    """Evaluate a plug, allowing Maya to report invalid input as a command error.

    The C++ contract is ``MS::kInvalidParameter``. Maya may surface that as a
    RuntimeError, or may return a cached/null data object, depending on the
    output plug and evaluation mode. Either is acceptable here; recovery after
    restoring a valid value is the assertion that matters.
    """
    try:
        return call(), None
    except RuntimeError as error:
        return None, error


class InputValidationTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        self.errors = []

        def capture(message, message_type, _client_data):
            if message_type == om.MCommandMessage.kError:
                self.errors.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)

    def bell(self):
        node = cmds.createNode("yddBellCollider")
        cmds.setAttr(node + ".bellMatrix", *matrix(), type="matrix")
        cmds.setAttr(node + ".ringMatrix[0]", *matrix(0, 2, 0), type="matrix")
        return node

    def skirt(self):
        node = cmds.createNode("yddSkirtBellCollider")
        for side, x in (("left", -1), ("right", 1)):
            for part, y in (("Hip", 4), ("Knee", 2), ("Heel", 0)):
                cmds.setAttr(node + "." + side + part + "Matrix", *matrix(x, y), type="matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".bellMatrix", *matrix(0, 4, 0), type="matrix")
        cmds.setAttr(node + ".bellScale", 2, 1, 2, type="double3")
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        cmds.setAttr(node + ".skirtType", 1)
        return node

    def collision_mesh(self):
        fn = om.MFnMesh()
        obj = fn.create(
            om.MPointArray([om.MPoint(0.2, 0.25, 0), om.MPoint(0.2, 1, 0), om.MPoint(0.2, 2.1, 0)]),
            [3],
            [0, 1, 2],
        )
        return om.MFnDagNode(obj).fullPathName()

    def collision(self, mesh):
        node = cmds.deformer(mesh, type="yddSkirtCollideDeformer")[0]
        for side, x in (("left", 0), ("right", 4)):
            for part, y in (("Hip", 0), ("Knee", 1), ("Heel", 2)):
                cmds.setAttr(node + "." + side + part + "Matrix", *matrix(x, y), type="matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        return node

    def assert_recovers(self, node, attribute, invalid_values, read_value):
        target = node + "." + attribute
        valid = cmds.getAttr(target)
        driver = cmds.createNode("addDoubleLinear")
        cmds.setAttr(driver + ".input1", valid)
        cmds.connectAttr(driver + ".output", target)
        baseline = read_value()
        self.assertTrue(finite_points(baseline))
        for invalid in invalid_values:
            with self.subTest(attribute=attribute, invalid=invalid):
                try:
                    cmds.setAttr(driver + ".input1", invalid)
                    self.assertEqual(cmds.getAttr(target), invalid)
                    self.errors.clear()
                    evaluate(read_value)
                    self.assertTrue(
                        any(attribute in message and "must be" in message for message in self.errors),
                        self.errors,
                    )
                finally:
                    cmds.setAttr(driver + ".input1", valid)
                recovered = read_value()
                self.assertTrue(finite_points(recovered))
                self.assertEqual(len(recovered), len(baseline))
                for actual, expected in zip(recovered, baseline):
                    for a, b in zip(actual, expected):
                        self.assertAlmostEqual(a, b, places=6)
        cmds.disconnectAttr(driver + ".output", target)
        cmds.setAttr(target, valid)
        cmds.delete(driver)

    def test_bell_subdivision_outputs_recover_from_invalid_values(self):
        node = self.bell()
        for attribute in ("outputBellMesh", "outputCurve"):
            self.assert_recovers(
                node, "bellSubdivision", (-1, 0, 2, 4097, 2**31 - 1),
                lambda: bell_points(node, attribute),
            )

    def test_skirt_domains_recover_from_invalid_values(self):
        node = self.skirt()
        for attribute, invalid_values in (
            ("bellSubdivision", (-1, 0, 2, 4097, 2**31 - 1)),
            ("skirtType", (-1, 2)),
            ("leftRingAxis", (-1, 6)),
            ("rightRingAxis", (-1, 6)),
            ("bellAxis", (-1, 6)),
        ):
            self.assert_recovers(node, attribute, invalid_values, lambda: surface_points(node))

    def test_skirt_deformer_domains_recover_from_invalid_values(self):
        mesh = self.collision_mesh()
        node = self.collision(mesh)
        for attribute in ("skirtType", "leftRingAxis", "rightRingAxis"):
            invalid_values = (-1, 2) if attribute == "skirtType" else (-1, 6)
            self.assert_recovers(node, attribute, invalid_values, lambda: mesh_points(mesh))

    def test_plane_normal_axis_recovers_from_invalid_values(self):
        node = cmds.createNode("yddPlaneCollider")
        cmds.setAttr(node + ".planeMatrix", *matrix(), type="matrix")
        cmds.setAttr(node + ".inputPosition", 1, -2, 3)
        self.assert_recovers(
            node, "normalAxis", (-1, 6),
            lambda: cmds.getAttr(node + ".outputPosition"),
        )

    def test_subdivision_limits_and_valid_endpoints(self):
        for node in (self.bell(), self.skirt()):
            for attribute in ("bellSubdivision", "ringSubdivision"):
                self.assertEqual(cmds.attributeQuery(attribute, node=node, minimum=True), [3])
                self.assertEqual(cmds.attributeQuery(attribute, node=node, maximum=True), [4096])
            for value in (3, 4096):
                with self.subTest(node=node, subdivision=value):
                    cmds.setAttr(node + ".bellSubdivision", value)
                    if cmds.nodeType(node) == "yddBellCollider":
                        points = bell_points(node, "outputBellMesh")
                        self.assertEqual(len(points), 2 * value + 1)
                    else:
                        points = surface_points(node)
                        self.assertEqual(len(points), (value + 3) * 4)
                    self.assertTrue(finite_points(points))

    def test_all_six_plane_axes_project_the_selected_coordinate(self):
        node = cmds.createNode("yddPlaneCollider")
        for axis in range(6):
            sign = -1 if axis < 3 else 1
            source = [sign * value for value in (1, 2, 3)]
            expected = list(source)
            expected[axis % 3] = 0
            cmds.setAttr(node + ".normalAxis", axis)
            cmds.setAttr(node + ".inputPosition", *source)
            self.assertEqual(list(cmds.getAttr(node + ".outputPosition")[0]), expected)
