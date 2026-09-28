"""Waist contact geometry and stateless Maya evaluation tests."""

import hashlib
import math
import struct
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

import yddColliders
from helpers import coordinates, output_object, skirt, transform


def matrix_row(matrix, row):
    return om.MVector(*(matrix[row * 4 + column] for column in range(3)))


def waist_frame(node):
    matrix = om.MMatrix(cmds.getAttr(node + ".bellMatrix"))
    axis = cmds.getAttr(node + ".bellAxis")
    direction = matrix_row(matrix, axis % 3)
    if axis >= 3:
        direction = -direction
    direction = direction.normal() if direction.length() >= 1e-4 else om.MVector(0, -1, 0)
    return om.MPoint(matrix_row(matrix, 3)), direction


def radius(point, waist, direction):
    offset = point - waist
    return (offset - direction * (offset * direction)).length()


def thigh_cylinder(node, side):
    joint = om.MMatrix(cmds.getAttr(node + "." + side + "HipMatrix"))
    hip = om.MPoint(matrix_row(joint, 3))
    knee = om.MPoint(matrix_row(om.MMatrix(cmds.getAttr(node + "." + side + "KneeMatrix")), 3))
    axis = cmds.getAttr(node + "." + side + "RingAxis")
    y = matrix_row(joint, axis % 3)
    if axis >= 3:
        y = -y
    y = y.normal()
    ux = matrix_row(joint, (axis + 1) % 3).normal()
    target = knee - hip
    if target.length() > 1e-6:
        rotation = y.rotateTo(target).asMatrix()
        y = y * rotation
        ux = ux * rotation
    x = (ux - y * (ux * y)).normal()
    z = (y ^ x).normal()
    lengths = []
    for leg in ("left", "right"):
        h = om.MPoint(matrix_row(om.MMatrix(cmds.getAttr(node + "." + leg + "HipMatrix")), 3))
        k = om.MPoint(matrix_row(om.MMatrix(cmds.getAttr(node + "." + leg + "KneeMatrix")), 3))
        lengths.append((k - h).length())
    scale = cmds.getAttr(node + ".ringScale")[0]
    x *= scale[0]
    y *= sum(lengths) * 0.5 * scale[1]
    z *= scale[2]
    matrix = om.MMatrix([
        x.x, x.y, x.z, 0, y.x, y.y, y.z, 0,
        z.x, z.y, z.z, 0, hip.x, hip.y, hip.z, 1,
    ])
    return {"inverse": matrix.inverse(), "hip": hip, "normal": y.normal()}


def plane_distance(point, cylinder):
    return (point - cylinder["hip"]) * cylinder["normal"]


def outside_cylinder(point, cylinder, tolerance):
    if plane_distance(point, cylinder) <= 0.0:
        return True
    local = point * cylinder["inverse"]
    return math.hypot(local.x, local.z) >= 1.0 - tolerance


def point_hash(points):
    digest = hashlib.sha256()
    for point in points:
        digest.update(struct.pack("<4d", *point))
    return digest.hexdigest()


class WaistRingTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")
        self.surface_data = []

    def surface(self, node):
        data = output_object(node, "outputSurface")
        self.surface_data.append(data)
        return om.MFnNurbsSurface(data)

    def rows(self, node):
        surface = self.surface(node)
        points = surface.cvPositions()
        nu, nv = surface.numCVsInU, surface.numCVsInV
        self.assertTrue(all(math.isfinite(value) for point in points for value in (point.x, point.y, point.z)))
        rows = [[om.MPoint(points[(nu - 1 - j) * nv + v]) for j in range(nu)] for v in range(nv)]
        count = cmds.getAttr(node + ".bellSubdivision")
        for row in rows:
            for j in range(3):
                self.assertEqual(tuple(row[count + j]), tuple(row[j]))
        return [row[:count] for row in rows]

    def rig(self, skirt_type=1, angle=90, both=False, waist_height=11.5):
        cmds.file(new=True, force=True)
        root = transform()
        joints = {}
        for side, x in (("left", -1), ("right", 1)):
            hip = transform((x, 10, 0))
            knee = transform((x, 5, 0))
            heel = transform((x, 0, 0))
            cmds.parent(knee, heel, hip)
            cmds.parent(hip, root)
            joints[side + "HipObj"] = hip
            joints[side + "KneeObj"] = knee
            joints[side + "HeelObj"] = heel
        node, locator, _ = yddColliders.createSkirtBellCollider(skirtType=skirt_type, **joints)
        cmds.xform(locator, worldSpace=True, translation=(0, waist_height, 0))
        cmds.parent(locator, root)
        cmds.setAttr(node + ".bellSubdivision", 16)
        cmds.setAttr(node + ".bellScale", 2, 1, 2, type="double3")
        cmds.setAttr(node + ".ringScale", 2, 1, 2, type="double3")
        # A negative X rotation takes the downward knee vector toward positive Z.
        cmds.setAttr(joints["leftHipObj"] + ".rotateX", -angle)
        if both:
            cmds.setAttr(joints["rightHipObj"] + ".rotateX", -angle)
        return node, locator, root, joints

    def baseline_and_solution(self, node):
        original = cmds.getAttr(node + ".nodeState")
        try:
            cmds.setAttr(node + ".nodeState", 1)
            baseline = self.rows(node)
            cmds.setAttr(node + ".nodeState", 0)
            return baseline, self.rows(node)
        finally:
            cmds.setAttr(node + ".nodeState", original)

    def tolerance(self, node):
        return 1e-6 * cmds.getAttr(node + ".ringScale")[0][0]

    def assert_row_unchanged(self, expected, actual, tolerance):
        for j, (before, after) in enumerate(zip(expected, actual)):
            with self.subTest(cv=j):
                self.assertLessEqual((after - before).length(), tolerance)

    def test_upper_contact(self):
        for skirt_type in (0, 1):
            for smooth in (0.0, 0.5):
                with self.subTest(skirt_type=skirt_type, smooth=smooth):
                    node, _, _, joints = self.rig(skirt_type)
                    cmds.setAttr(node + ".smoothness", smooth)
                    tolerance = self.tolerance(node)
                    for name, expected in (("leftHipObj", (-1, 10, 0)),
                                           ("rightHipObj", (1, 10, 0)),
                                           ("leftKneeObj", (-1, 10, 5)),
                                           ("rightKneeObj", (1, 5, 0)),
                                           ("leftHeelObj", (-1, 10, 10)),
                                           ("rightHeelObj", (1, 0, 0))):
                        position = cmds.xform(joints[name], query=True, worldSpace=True, translation=True)
                        self.assertLessEqual((om.MPoint(position) - om.MPoint(expected)).length(), tolerance)
                    before, after = self.baseline_and_solution(node)
                    cylinders = [thigh_cylinder(node, side) for side in ("left", "right")]
                    front = min(before[0], key=lambda p: (p - om.MPoint(0, 11.5, 2)).length())
                    self.assertLessEqual((front - om.MPoint(0, 11.5, 2)).length(), tolerance)
                    self.assertFalse(outside_cylinder(front, cylinders[0], tolerance))
                    self.assertGreater(max((q - p).length() for p, q in zip(before[0], after[0])), tolerance)
                    for side, cylinder in zip(("left", "right"), cylinders):
                        for j, point in enumerate(after[0]):
                            with self.subTest(side=side, cv=j):
                                self.assertTrue(outside_cylinder(point, cylinder, tolerance))
                    if smooth == 0.0:
                        waist, direction = waist_frame(node)
                        normal = cylinders[0]["normal"]
                        forward = (normal - direction * (normal * direction)).normal()
                        back_count = 0
                        for j, (p, q) in enumerate(zip(before[0], after[0])):
                            offset = p - waist
                            rho = (offset - direction * (offset * direction)).normal()
                            if rho * forward < 0.0:
                                back_count += 1
                                with self.subTest(back_cv=j):
                                    self.assertLessEqual((q - p).length(), tolerance)
                        self.assertGreater(back_count, 0)

    def test_hip_upstream_is_unchanged(self):
        for skirt_type in (0, 1):
            with self.subTest(skirt_type=skirt_type):
                node, _, _, _ = self.rig(skirt_type)
                cmds.setAttr(node + ".smoothness", 0)
                before, after = self.baseline_and_solution(node)
                cylinder = thigh_cylinder(node, "left")
                count = 0
                for j, (p, q) in enumerate(zip(before[0], after[0])):
                    if plane_distance(p, cylinder) <= 0.0:
                        count += 1
                        with self.subTest(cv=j):
                            self.assertLessEqual((q - p).length(), self.tolerance(node))
                self.assertGreater(count, 0)

    def test_higher_waist_upper_is_unchanged(self):
        for skirt_type in (0, 1):
            for smooth in (0.0, 0.5):
                with self.subTest(skirt_type=skirt_type, smooth=smooth):
                    node, _, _, _ = self.rig(skirt_type, waist_height=12)
                    cmds.setAttr(node + ".smoothness", smooth)
                    before, after = self.baseline_and_solution(node)
                    self.assert_row_unchanged(before[0], after[0], self.tolerance(node))

    def test_sweep_and_reverse(self):
        for skirt_type in (0, 1):
            node, _, _, joints = self.rig(skirt_type, angle=0)
            cmds.setAttr(node + ".smoothness", 0)
            reference = {}
            previous = None
            limit = 0.1 * cmds.getAttr(node + ".ringScale")[0][0]
            for angle in list(range(91)) + list(range(90, -1, -1)):
                cmds.setAttr(joints["leftHipObj"] + ".rotateX", -angle)
                upper = self.rows(node)[0]
                result = point_hash(coordinates(node, "outputSurface", "surface"))
                with self.subTest(skirt_type=skirt_type, angle=angle):
                    if previous is not None:
                        self.assertLessEqual(max((a - b).length() for a, b in zip(upper, previous)), limit)
                    if angle in reference:
                        self.assertEqual(result, reference[angle])
                    else:
                        reference[angle] = result
                previous = upper

    def test_evaluation_modes_and_time_order(self):
        for skirt_type in (0, 1):
            node, _, _, joints = self.rig(skirt_type, angle=0)
            cmds.setAttr(node + ".smoothness", 0.5)
            for time, value in ((0, 0), (12, -90)):
                cmds.setKeyframe(joints["leftHipObj"] + ".rotateX", time=time, value=value,
                                 inTangentType="linear", outTangentType="linear")
            references = {}
            for mode in ("off", "serial", "parallel"):
                cmds.evaluationManager(mode=mode)
                for time in (0, 4, 12, 7, 0):
                    with self.subTest(skirt_type=skirt_type, mode=mode, time=time):
                        cmds.currentTime(time)
                        points = coordinates(node, "outputSurface", "surface")
                        self.assertTrue(all(math.isfinite(value) for point in points for value in point))
                        result = point_hash(points)
                        if time in references:
                            self.assertEqual(result, references[time])
                        else:
                            references[time] = result
                if mode != "off":
                    self.assertFalse(cmds.evaluationManager(query=True, empty=True))
                    self.assertFalse(cmds.evaluationManager(query=True, fallbackTriggered=True))

    def test_symmetry(self):
        for skirt_type in (0, 1):
            with self.subTest(skirt_type=skirt_type):
                node, _, _, _ = self.rig(skirt_type, angle=60, both=True)
                _, after = self.baseline_and_solution(node)
                for j, point in enumerate(after[0]):
                    reflected = after[0][(16 // 2 - j) % 16]
                    self.assertLessEqual((reflected - om.MPoint(-point.x, point.y, point.z)).length(),
                                         self.tolerance(node))

    def test_root_tilt_displacement(self):
        for skirt_type in (0, 1):
            for angle, both in ((60, True), (90, False)):
                with self.subTest(skirt_type=skirt_type, angle=angle, both=both):
                    node, _, root, _ = self.rig(skirt_type, angle=angle, both=both)
                    before, after = self.baseline_and_solution(node)
                    inverse = om.MMatrix(cmds.getAttr(node + ".bellMatrix")).inverse()
                    displacements = [(q - p) * inverse for p, q in zip(before[0], after[0])]
                    cmds.setAttr(root + ".rotateZ", 30)
                    tilted_before, tilted_after = self.baseline_and_solution(node)
                    tilted_inverse = om.MMatrix(cmds.getAttr(node + ".bellMatrix")).inverse()
                    for j, (p, q) in enumerate(zip(tilted_before[0], tilted_after[0])):
                        self.assertLessEqual(((q - p) * tilted_inverse - displacements[j]).length(),
                                             self.tolerance(node))

    def test_static_noncontact(self):
        for count in (3, 16, 64):
            for height in (10, 12):
                for skirt_type in (0, 1):
                    with self.subTest(subdivision=count, waist_height=height, skirt_type=skirt_type):
                        cmds.file(new=True, force=True)
                        node, locator, _ = skirt()
                        cmds.xform(locator, worldSpace=True, translation=(0, height, 0))
                        cmds.setAttr(node + ".bellSubdivision", count)
                        cmds.setAttr(node + ".skirtType", skirt_type)
                        before, after = self.baseline_and_solution(node)
                        self.assert_row_unchanged(before[0], after[0], self.tolerance(node))

    def test_node_state_returns_raw_rings_and_restores_solution(self):
        for mode in ("off", "serial", "parallel"):
            for skirt_type in (0, 1):
                for tightness in (0.0, 0.5, 1.0):
                    with self.subTest(mode=mode, skirt_type=skirt_type, tightness=tightness):
                        node, _, _, _ = self.rig(skirt_type)
                        cmds.evaluationManager(mode=mode)
                        cmds.setAttr(node + ".smoothness", 0.5)
                        cmds.setAttr(node + ".follow", 1)
                        cmds.setAttr(node + ".tightness", tightness)
                        normal = self.rows(node)
                        surface = self.surface(node)
                        knots = (list(surface.knotsInU()), list(surface.knotsInV()))
                        waist, direction = waist_frame(node)
                        matrix = om.MMatrix(cmds.getAttr(node + ".bellMatrix"))
                        x = matrix_row(matrix, 0).normal()
                        x = (x - direction * (x * direction)).normal()
                        z = (x ^ direction).normal()
                        # This fixture has a 1.5 waist offset and two 5-unit leg segments.
                        levels = (0.0, 4.0, 6.5, 11.5) if skirt_type else (0.0, 4.0, 6.5)
                        count = cmds.getAttr(node + ".bellSubdivision")
                        expected = [
                            [waist + direction * level
                             + x * (2.0 * math.cos(2.0 * math.pi * j / count))
                             + z * (2.0 * math.sin(2.0 * math.pi * j / count))
                             for j in range(count)]
                            for level in levels
                        ]
                        cmds.setAttr(node + ".nodeState", 1)
                        raw = self.rows(node)
                        self.assertEqual(len(raw), len(expected))
                        for before, after in zip(expected, raw):
                            self.assert_row_unchanged(before, after, 2e-6)
                        surface = self.surface(node)
                        self.assertEqual(knots, (list(surface.knotsInU()), list(surface.knotsInV())))
                        self.assertGreater(
                            max((a - b).length() for row_a, row_b in zip(normal, raw)
                                for a, b in zip(row_a, row_b)), 1e-5
                        )
                        cmds.setAttr(node + ".nodeState", 0)
                        restored = self.rows(node)
                        self.assertEqual([point_hash(row) for row in normal],
                                         [point_hash(row) for row in restored])

    def test_upper_node_state_bypass_and_follow_isolation(self):
        for skirt_type in (0, 1):
            for smooth in (0.0, 0.5):
                with self.subTest(skirt_type=skirt_type, smooth=smooth):
                    node, _, _, _ = self.rig(skirt_type)
                    cmds.setAttr(node + ".smoothness", smooth)
                    reference = {}
                    for follow in (0.0, 1.0):
                        cmds.setAttr(node + ".follow", follow)
                        for state in (1, 0, 1):
                            cmds.setAttr(node + ".nodeState", state)
                            upper = self.rows(node)[0]
                            result = point_hash(upper)
                            if state in reference:
                                self.assertEqual(result, reference[state])
                            else:
                                reference[state] = result
                            if state == 1:
                                waist, direction = waist_frame(node)
                                for point in upper:
                                    self.assertAlmostEqual(radius(point, waist, direction), 2.0,
                                                           delta=self.tolerance(node))
                                    self.assertAlmostEqual((point - waist) * direction, 0.0,
                                                           delta=self.tolerance(node))
