"""Continuous first-pass lift on synthetic skirt surfaces."""

import json
import math
from pathlib import Path
import struct
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import transform


SUBDIVISION = 16
R0 = math.sqrt((2.0 * 1.001) ** 2 - 1.0)
REST_REFERENCE = Path(__file__).parent / "fixtures" / "continuous_lift_activation_rest.json"


def point_values(point):
    return (point.x, point.y, point.z, point.w)


def packed_points(points):
    return b"".join(struct.pack("<4d", *point_values(point)) for point in points)


def matrix_axis(matrix, row):
    return om.MVector(*(matrix[4 * row + column] for column in range(3)))


def matrix_origin(matrix):
    return om.MPoint(*(matrix[12 + column] for column in range(3)))


def frame(x, y, z, origin):
    return om.MMatrix(
        (
            x.x,
            x.y,
            x.z,
            0.0,
            y.x,
            y.y,
            y.z,
            0.0,
            z.x,
            z.y,
            z.z,
            0.0,
            origin.x,
            origin.y,
            origin.z,
            1.0,
        )
    )


def smoothstep(x):
    t = min(1.0, max(0.0, x))
    return (t * t) * (3.0 - 2.0 * t)


def sphere_line_hits(point, direction, centre, radius):
    direction = direction.normal()
    offset = point - centre
    b = 2.0 * (direction.x * offset.x + direction.y * offset.y + direction.z * offset.z)
    c = offset.x**2 + offset.y**2 + offset.z**2 - radius**2
    delta = b**2 - 4.0 * c
    if delta <= 0.0:
        return []
    return [point + direction * ((-b + sign * math.sqrt(delta)) / 2.0) for sign in (1.0, -1.0)]


def lift_geometry(bell, ring, base_row):
    result = dict(found=False, extended=False, a=0.0, byHeight=0.0, byMargin=0.0)
    axis = matrix_axis(bell, 1)
    if not base_row or axis.length() <= 1e-5:
        return result
    normal = axis.normal()
    origin = matrix_origin(bell)
    inverse = bell.inverse()
    translation = matrix_origin(ring)
    direction = matrix_axis(ring, 1)
    ring_normal = direction.normal()
    ring_inverse = ring.inverse()
    projected = direction - normal * (direction * normal)
    if projected.length() <= 1e-3:
        return result
    projected_origin = translation - normal * ((translation - origin) * normal)
    hits = sphere_line_hits(projected_origin * inverse, projected * inverse, om.MPoint(), 1.001)
    if not hits:
        return result
    bell_hit = hits[0] * bell + axis
    signed_projection = projected * (1.0 if ring_normal * normal > 0.0 else -1.0)
    ring_projection = signed_projection - ring_normal * (signed_projection * ring_normal)
    ring_projection_length = ring_projection.length()
    offset = om.MVector()
    if ring_projection_length > 1e-5:
        local_length = (ring_projection * ring_inverse).length()
        d = ring_projection_length / local_length if local_length > 1e-5 else 1.0
        offset = ring_projection.normal() * d
    line_point = translation + offset
    radius = (bell_hit - translation).length()
    ring_hit = None
    for hit in sphere_line_hits(line_point, direction, translation, radius):
        if (hit - translation) * direction > 0.0:
            ring_hit = hit
    if ring_hit is None and ring_projection_length > 1e-5:
        ring_hit = line_point
        result["extended"] = True
    if ring_hit is None:
        return result
    delta = (((ring_hit - origin) * normal) - ((bell_hit - origin) * normal)) / axis.length()
    result.update(found=True, collisionDelta=delta, R=radius, rOffset=(line_point - translation).length())
    if not delta < 0.0:
        return result
    matrix = ring_inverse.inverse()
    a_axis, b_axis = matrix_axis(matrix, 0), matrix_axis(matrix, 2)
    r = max(a_axis.length(), b_axis.length())
    length = direction.length()
    up = -normal
    height = math.hypot(a_axis * up, b_axis * up)
    contributions = []
    for raw in base_row:
        point = om.MPoint(raw)
        if point.w != 1.0:
            point.cartesianize()
        local = om.MVector(point - translation) * ring_inverse
        axial = local.y * length
        fade = smoothstep((axial + 0.25 * r) / (0.25 * r)) if r > 1e-12 else float(axial >= 0.0)
        penetration = 1.0 - math.hypot(local.x, local.z)
        centre = translation + direction * local.y
        h = height - (point - centre) * up
        g = smoothstep((penetration + 0.10) / 0.10)
        b = g * smoothstep((h - 0.25 * r) / (0.5 * r)) if r > 1e-12 else 0.0
        if fade <= 0.0 or g <= 0.0 or r <= 1e-12:
            continue
        contribution = b * fade
        contributions.append(contribution)
        if contribution > result["byHeight"]:
            result["byHeight"] = contribution
    if not result["extended"] and r > 1e-12:
        result["byMargin"] = smoothstep((radius - result["rOffset"]) / (0.5 * r))
    amplitude = max(result["byHeight"], result["byMargin"])
    direction_factor = smoothstep(projected.length() / (0.05 * length)) if length > 1e-12 else 1.0
    if length > 1e-12:
        amplitude *= direction_factor
    result.update(a=amplitude, r=r, H=height, directionFactor=direction_factor, contributions=contributions)
    return result


class LiftFixture:
    def __init__(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.node = cmds.createNode("yddSkirtBellCollider")
        self.waist = transform((0.0, 12.5, 0.0))
        cmds.connectAttr(self.waist + ".worldMatrix[0]", self.node + ".bellMatrix")
        self.hips = {}
        for side, x in (("left", 1.0), ("right", -1.0)):
            hip, knee, heel = (transform((x, y, 0.0)) for y in (10.0, 5.0, 0.0))
            cmds.parent(knee, heel, hip)
            self.hips[side] = hip
            for part, joint in (("Hip", hip), ("Knee", knee), ("Heel", heel)):
                cmds.connectAttr(joint + ".worldMatrix[0]", self.node + "." + side + part + "Matrix")
            cmds.setAttr(self.node + "." + side + "RingAxis", 4)
        for attribute, value in dict(
            bellAxis=4,
            skirtType=0,
            height=0.25,
            bellSubdivision=SUBDIVISION,
            falloff=-1.0,
            follow=0.0,
            smoothness=0.0,
            tightness=1.0,
        ).items():
            cmds.setAttr(self.node + "." + attribute, value)
        for attribute in ("bellScale", "ringScale"):
            cmds.setAttr(self.node + "." + attribute, 2.0, 1.0, 2.0, type="double3")
        for level in ("thigh", "knee", "calf", "ankle"):
            for axis in ("X", "Z"):
                value = cmds.getAttr(self.node + "." + level + "Radius" + axis)
                if value != 1.0:
                    raise AssertionError((level, axis, "radius default must be one", value))
        for level in ("thigh", "calf"):
            if cmds.getAttr(self.node + "." + level + "Position") != 0.5:
                raise AssertionError("Position profile default must be 0.5")
        entries = cmds.getAttr(self.node + ".bellScaleRamp", multiIndices=True) or []
        ramp = sorted(
            (
                cmds.getAttr(self.node + ".bellScaleRamp[%d].bellScaleRamp_Position" % index),
                cmds.getAttr(self.node + ".bellScaleRamp[%d].bellScaleRamp_FloatValue" % index),
            )
            for index in entries
        )
        if ramp != [(0.0, 1.0), (1.0, 1.0)]:
            raise AssertionError(("Bell ramp must be one throughout", ramp))
        for attribute in ("seams", "columnOffsetMatrix"):
            if cmds.getAttr(self.node + "." + attribute, multiIndices=True):
                raise AssertionError(attribute + " must be empty")
        self.reader = cmds.createNode("nurbsSurface")
        cmds.connectAttr(self.node + ".outputSurface", self.reader + ".create")

    def pose(self, theta, radius=2.0):
        cmds.setAttr(self.hips["left"] + ".rotateX", -theta)
        cmds.setAttr(self.node + ".ringScale", radius, 1.0, radius, type="double3")

    def read(self):
        surface = om.MFnNurbsSurface(om.MSelectionList().add(self.reader).getDagPath(0))
        dimensions = (surface.numCVsInU, surface.numCVsInV)
        # Copy the coordinates: elements of the returned array alias storage that later evaluations reuse.
        points = [om.MPoint(p.x, p.y, p.z, p.w) for p in surface.cvPositions(om.MSpace.kObject)]
        if len(points) != dimensions[0] * dimensions[1]:
            raise AssertionError(("CV count does not match dimensions", dimensions, len(points)))
        if not all(math.isfinite(value) for point in points for value in point_values(point)):
            raise AssertionError("Nonfinite outputSurface CV")
        if dimensions[0] - 3 != SUBDIVISION:
            raise AssertionError(("Independent periodic CV count changed", dimensions))
        return dimensions, points

    def base(self):
        cmds.setAttr(self.node + ".nodeState", 1)
        try:
            return self.read()
        finally:
            cmds.setAttr(self.node + ".nodeState", 0)

    @staticmethod
    def row(snapshot, v):
        dimensions, points = snapshot
        return [points[u * dimensions[1] + v] for u in range(1, SUBDIVISION + 1)]

    def select_rows(self, snapshot, heights):
        rows = []
        for height in heights:
            candidates = [
                v for v in range(snapshot[0][1]) if abs(sum(p.y for p in self.row(snapshot, v)) / SUBDIVISION - height) <= 1e-6
            ]
            if len(candidates) != 1:
                raise AssertionError(("Expected exactly one measured row at height", height, candidates))
            rows.append(candidates[0])
        return rows

    def left_ring(self):
        hip = om.MMatrix(cmds.getAttr(self.node + ".leftHipMatrix"))
        knee = matrix_origin(om.MMatrix(cmds.getAttr(self.node + ".leftKneeMatrix")))
        right_hip = matrix_origin(om.MMatrix(cmds.getAttr(self.node + ".rightHipMatrix")))
        right_knee = matrix_origin(om.MMatrix(cmds.getAttr(self.node + ".rightKneeMatrix")))
        translation = matrix_origin(hip)
        length = ((knee - translation).length() + (right_knee - right_hip).length()) * 0.5
        y = (-matrix_axis(hip, 1)).normal()
        ux = matrix_axis(hip, 2).normal()
        toward_knee = knee - translation
        if toward_knee.length() > 1e-6:
            rotation = y.rotateTo(toward_knee).asMatrix()
            y = y * rotation
            ux = ux * rotation
        x = (ux - y * (ux * y)).normal()
        z = (y ^ x).normal()
        scale = cmds.getAttr(self.node + ".ringScale")[0]
        return frame(x * scale[0], y * (length * scale[1]), z * scale[2], translation)

    def short_contact(self, base, height=10.0):
        v = self.select_rows(base, (height,))[0]
        waist = om.MMatrix(cmds.getAttr(self.node + ".bellMatrix"))
        origin = matrix_origin(waist)
        normal = (-matrix_axis(waist, 1)).normal()
        waist_x = matrix_axis(waist, 0).normal()
        x = (waist_x - normal * (waist_x * normal)).normal()
        z = (x ^ normal).normal()
        distance = origin.y - sum(p.y for p in self.row(base, v)) / SUBDIVISION
        scale = cmds.getAttr(self.node + ".bellScale")[0]
        bell = frame(x * scale[0], normal * (distance * scale[1]), z * scale[2], origin)
        ring = self.left_ring()
        return lift_geometry(bell, ring, self.row(base, v)), ring, v


KICK_BELL = (0.0, -0.03412, 0.99942, 0.0, 0.0, -0.99942, -0.03412, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 8.05221, -0.23446, 1.0)
KICK_RIGHT = dict(
    rightHip=(0.0, -0.88735, 0.46109, 0.0, 0.0, 0.46109, 0.88735, 0.0, -1.0, 0.0, 0.0, 0.0, -8.04797, -4.04666, -0.37385, 1.0),
    rightKnee=(0.0, -0.88919, 0.45754, 0.0, 0.0, 0.45754, 0.88919, 0.0, -1.0, 0.0, 0.0, 0.0, -8.04797, -39.38807, 17.99039, 1.0),
    rightHeel=(0.0, -0.88919, 0.45754, 0.0, 0.0, 0.45754, 0.88919, 0.0, -1.0, 0.0, 0.0, 0.0, -8.04797, -74.92938, 36.27826, 1.0),
)
# Left leg kicked forward until the thigh is close to horizontal.
KICK_LEFT = dict(
    leftHip=(
        0.30088,
        0.0271,
        -0.95328,
        0.0,
        0.0319,
        -0.99932,
        -0.01834,
        0.0,
        -0.95313,
        -0.02489,
        -0.30154,
        0.0,
        8.04797,
        -4.04666,
        -0.37385,
        1.0,
    ),
    leftKnee=(
        -0.02173,
        0.99967,
        -0.01384,
        0.0,
        0.43859,
        -0.00291,
        -0.89868,
        0.0,
        -0.89842,
        -0.0256,
        -0.43838,
        0.0,
        20.03142,
        -2.96723,
        -38.34085,
        1.0,
    ),
    leftHeel=(
        -0.02173,
        0.99967,
        -0.01384,
        0.0,
        0.43859,
        -0.00291,
        -0.89868,
        0.0,
        -0.89842,
        -0.0256,
        -0.43838,
        0.0,
        19.16295,
        36.98988,
        -38.89401,
        1.0,
    ),
)
KICK_RING_SCALE = (9.727, 1.0, 6.333)
KICK_THIGH_RADII = (0.996, 1.253)
KICK_THIGH_POSITION = 0.4125


class KickFixture:
    """Short skirt with a wide bell, hips off the axis and a thigh ring larger than the hip-height row gap."""

    def __init__(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.node = cmds.createNode("yddSkirtBellCollider")
        self.joints = {}
        bell = transform()
        cmds.xform(bell, matrix=KICK_BELL, worldSpace=True)
        cmds.connectAttr(bell + ".worldMatrix[0]", self.node + ".bellMatrix")
        for name, matrix in list(KICK_LEFT.items()) + list(KICK_RIGHT.items()):
            joint = transform()
            cmds.xform(joint, matrix=matrix, worldSpace=True)
            cmds.connectAttr(joint + ".worldMatrix[0]", self.node + "." + name + "Matrix")
            self.joints[name] = joint
        for attribute, value in dict(
            bellAxis=1,
            leftRingAxis=0,
            rightRingAxis=0,
            height=0.327,
            skirtType=0,
            bellSubdivision=SUBDIVISION,
            thighPosition=KICK_THIGH_POSITION,
            calfPosition=0.3125,
            thighRadiusX=KICK_THIGH_RADII[0],
            thighRadiusZ=KICK_THIGH_RADII[1],
            kneeRadiusX=0.512,
            kneeRadiusZ=0.627,
            calfRadiusX=0.642,
            calfRadiusZ=0.773,
            ankleRadiusX=0.327,
            ankleRadiusZ=0.412,
            falloff=-1.0,
            follow=0.2,
            smoothness=0.1,
            tightness=0.6,
        ).items():
            cmds.setAttr(self.node + "." + attribute, value)
        cmds.setAttr(self.node + ".bellScale", 15.65, 1.0, 15.65, type="double3")
        cmds.setAttr(self.node + ".ringScale", *KICK_RING_SCALE, type="double3")
        for index in cmds.getAttr(self.node + ".bellScaleRamp", multiIndices=True) or []:
            cmds.removeMultiInstance(self.node + ".bellScaleRamp[%d]" % index, b=True)
        for index, (position, value) in enumerate(((0.01365, 0.53984), (0.33858, 0.69433), (0.66352, 0.84899), (0.98059, 1.0))):
            plug = self.node + ".bellScaleRamp[%d]" % index
            cmds.setAttr(plug + ".bellScaleRamp_Position", position)
            cmds.setAttr(plug + ".bellScaleRamp_FloatValue", value)
            cmds.setAttr(plug + ".bellScaleRamp_Interp", 1)
        self.reader = cmds.createNode("nurbsSurface")
        cmds.connectAttr(self.node + ".outputSurface", self.reader + ".create")

    def read(self):
        surface = om.MFnNurbsSurface(om.MSelectionList().add(self.reader).getDagPath(0))
        dimensions = (surface.numCVsInU, surface.numCVsInV)
        points = [om.MPoint(p.x, p.y, p.z, p.w) for p in surface.cvPositions(om.MSpace.kObject)]
        if not all(math.isfinite(value) for point in points for value in point_values(point)):
            raise AssertionError("Nonfinite outputSurface CV")
        return dimensions, points

    def base(self):
        cmds.setAttr(self.node + ".nodeState", 1)
        try:
            return self.read()
        finally:
            cmds.setAttr(self.node + ".nodeState", 0)

    def matrix(self, attribute):
        return om.MMatrix(cmds.getAttr(self.node + "." + attribute))

    def left_thigh(self):
        """Frame and section radii of the left thigh ring (ring axis 0), matching the node's ring frame."""
        hip = self.matrix("leftHipMatrix")
        origin = matrix_origin(hip)
        knee = matrix_origin(self.matrix("leftKneeMatrix"))
        lengths = []
        for side in ("left", "right"):
            joints = [matrix_origin(self.matrix(side + part + "Matrix")) for part in ("Hip", "Knee")]
            lengths.append(((joints[1] - joints[0]).length(),))
        thigh = (lengths[0][0] + lengths[1][0]) * 0.5
        y = matrix_axis(hip, 0).normal()
        ux = matrix_axis(hip, 1).normal()
        toward_knee = knee - origin
        rotation = y.rotateTo(toward_knee).asMatrix()
        y = y * rotation
        ux = ux * rotation
        x = (ux - y * (ux * y)).normal()
        z = (y ^ x).normal()
        station = thigh * KICK_THIGH_POSITION

        def radii(axial):
            t = min(1.0, max(0.0, axial / station))
            return (
                KICK_RING_SCALE[0] * (1.0 + (KICK_THIGH_RADII[0] - 1.0) * t),
                KICK_RING_SCALE[2] * (1.0 + (KICK_THIGH_RADII[1] - 1.0) * t),
            )

        return origin, x, y, z, radii


class ContinuousLiftActivationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = LiftFixture()

    def assert_correspondence(self, snapshot, dimensions, rows):
        self.assertEqual(snapshot[0], dimensions)
        for row in rows:
            self.assertEqual(len(self.fixture.row(snapshot, row)), SUBDIVISION)

    def test_hip_height_forward_raise(self):
        fixture = self.fixture
        dimensions = fixture.read()[0]
        for theta in range(91):
            with self.subTest(theta=theta):
                fixture.pose(theta)
                output = fixture.read()
                self.assertEqual(output[0], dimensions)
                if theta != 89:
                    continue
                base = fixture.base()
                contact, ring, v = fixture.short_contact(base)
                self.assert_correspondence(base, dimensions, [v])
                self.assertTrue(contact["found"], contact)
                self.assertTrue(contact["extended"] or contact["R"] - contact["rOffset"] < 0.5 * contact["r"], contact)
                translation = matrix_origin(ring)
                normal = matrix_axis(ring, 1).normal()
                before, after = fixture.row(base, v), fixture.row(output, v)
                indices = [i for i, p in enumerate(before) if (p - translation) * normal > 0.0]
                self.assertTrue(indices, "No forward points at 89 degrees")
                upward = []
                for i in indices:
                    local = om.MVector(after[i] - translation) * ring.inverse()
                    self.assertGreaterEqual(math.hypot(local.x, local.z), 0.95, i)
                    upward.append((after[i] - before[i]) * om.MVector(0.0, 1.0, 0.0))
                self.assertGreaterEqual(min(upward), -1e-5 * contact["r"])
                self.assertGreater(max(upward), 1e-5 * contact["r"])

    def test_hip_row_forward_kick_lifts_over_thigh(self):
        fixture = KickFixture()
        base = fixture.base()
        dimensions, output = fixture.read()
        self.assertEqual(dimensions, base[0])
        count_v = dimensions[1]
        bell = fixture.matrix("bellMatrix")
        up = (-matrix_axis(bell, 1)).normal()
        waist = matrix_origin(bell)
        depths = []
        for v in range(count_v):
            row = [base[1][u * count_v + v] for u in range(1, SUBDIVISION + 1)]
            depths.append(sum((p - waist) * -up for p in row) / SUBDIVISION)
        hip = matrix_origin(fixture.matrix("leftHipMatrix"))
        hip_depth = (hip - waist) * -up
        v = min(range(1, count_v), key=lambda index: abs(depths[index] - hip_depth))
        self.assertLess(abs(depths[v] - hip_depth), 4.0, (depths, hip_depth))
        origin, x, y, z, radii = fixture.left_thigh()
        ring_radius = max(radii(0.0))
        checked = []
        for u in range(1, SUBDIVISION + 1):
            before = base[1][u * count_v + v]
            after = output[u * count_v + v]
            if before.x <= 0.5 or (before - origin) * y <= 0.0:
                continue
            offset = after - origin
            rx, rz = radii(offset * y)
            rho = math.hypot((offset * x) / rx, (offset * z) / rz)
            checked.append((u, rho, (after - before) * up))
        self.assertGreaterEqual(len(checked), 3, checked)
        for u, rho, rise in checked:
            self.assertGreaterEqual(rho, 0.95, (u, rho, rise))
            self.assertGreaterEqual(rise, -0.05 * ring_radius, (u, rho, rise))
        self.assertGreater(max(rise for _, _, rise in checked), 0.1 * ring_radius, checked)

    def assert_step_halving(self, pose, heights):
        fixture = self.fixture
        dimensions = None
        selected_rows = None

        def sample(steps):
            nonlocal dimensions, selected_rows
            samples = [[] for _ in heights]
            for k in range(steps + 1):
                theta, radius = pose(k / steps)
                fixture.pose(theta, radius)
                base = fixture.base()
                rows = fixture.select_rows(base, heights)
                output = fixture.read()
                if dimensions is None:
                    dimensions = output[0]
                    selected_rows = rows
                self.assertEqual(rows, selected_rows)
                self.assert_correspondence(base, dimensions, rows)
                self.assert_correspondence(output, dimensions, rows)
                for values, row in zip(samples, rows):
                    values.append(fixture.row(output, row))
            return samples

        coarse, fine = sample(64), sample(128)
        for height, coarse_row, fine_row in zip(heights, coarse, fine):

            def maximum_step(samples):
                return max((b - a).length() for first, second in zip(samples, samples[1:]) for a, b in zip(first, second))

            d_coarse, d_fine = maximum_step(coarse_row), maximum_step(fine_row)
            self.assertGreater(d_coarse, 1e-5, (height, d_coarse))
            self.assertLessEqual(d_fine, 0.75 * d_coarse, (height, d_coarse, d_fine))

    def test_generatrix_handover_step_halving(self):
        def pose(t):
            return 45.0, R0 + 0.4 * (t - 0.5) + 0.4 / 256.0

        margins = []
        for t in (0.0, 1.0):
            self.fixture.pose(*pose(t))
            contact, _, _ = self.fixture.short_contact(self.fixture.base())
            self.assertTrue(contact["found"], contact)
            margin = contact["R"] - contact["rOffset"]
            self.assertGreater(abs(margin), 1e-9, "Handover endpoint is too close to the root")
            margins.append(margin)
        self.assertLess(margins[0] * margins[1], 0.0, margins)
        self.assert_step_halving(pose, (10.0,))

    def test_return_to_vertical_step_halving(self):
        # The sweep stops at 1 degree: at exactly 0 degrees this row lies on the hip plane of the vertical leg, where the
        # relax cut at the ring origin switches the front half of the row on or off independently of the lift.
        self.assert_step_halving(lambda t: (44.0 * (1.0 - t) + 1.0, 2.0), (10.0,))

    def test_rest_matches_head_bits(self):
        if not REST_REFERENCE.is_file():
            self.skipTest("Rest reference is absent; run generate_continuous_lift_rest.py with the HEAD plugin first")
        reference = json.loads(REST_REFERENCE.read_text(encoding="utf-8"))
        self.assertEqual(reference["maya_version"], cmds.about(version=True), "Capture and test need the same Maya version")
        self.assertEqual(reference["maya_api_version"], cmds.about(apiVersion=True))
        self.assertTrue(reference["plugin_version"])
        self.assertTrue(reference["build_configuration"])
        self.fixture.pose(0.0)
        self.fixture.select_rows(self.fixture.base(), (10.0,))
        dimensions, points = self.fixture.read()
        self.assertEqual(list(dimensions), reference["dimensions"])
        self.assertEqual(packed_points(points), bytes.fromhex(reference["cv_xyzw_le_hex"]))

    def test_long_skirt_return_to_vertical_step_halving(self):
        cmds.setAttr(self.fixture.node + ".skirtType", 1)
        self.assert_step_halving(lambda t: (45.0 * (1.0 - t), 2.0), (7.5, 5.0))

    def test_margin_term_sample(self):
        outputs = []
        dimensions = None
        # The lower row of the short skirt is the one whose sphere radius exceeds the generatrix offset while the
        # section-top term is still small; on the hip-height row the height term always dominates the margin term.
        samples = ((1.75, False, 0.0742, 0.4153), (2.25, True, 0.0700, 0.0))
        for radius, extended, by_height, by_margin in samples:
            self.fixture.pose(10.0, radius)
            base = self.fixture.base()
            contact, _, v = self.fixture.short_contact(base, 8.75)
            self.assertTrue(contact["found"], contact)
            self.assertEqual(contact["extended"], extended, contact)
            self.assertLess(contact["collisionDelta"], 0.0)
            self.assertAlmostEqual(contact["byHeight"], by_height, delta=0.0005)
            self.assertEqual(contact["directionFactor"], 1.0)
            if extended:
                self.assertEqual(contact["byMargin"], 0.0)
                self.assertEqual(contact["a"], contact["byHeight"])
            else:
                self.assertAlmostEqual(contact["byMargin"], by_margin, delta=0.0005)
                self.assertGreater(contact["byMargin"], contact["byHeight"])
                self.assertEqual(contact["a"], contact["byMargin"])
            self.assertGreater(contact["a"], 0.0)
            self.assertLess(contact["a"], 1.0)
            output = self.fixture.read()
            if dimensions is None:
                dimensions = output[0]
            self.assert_correspondence(base, dimensions, [v])
            self.assert_correspondence(output, dimensions, [v])
            outputs.append(self.fixture.row(output, v))
        self.assertGreater(max((a - b).length() for a, b in zip(*outputs)), 1e-4)
