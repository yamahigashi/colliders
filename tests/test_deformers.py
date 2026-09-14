import math
import unittest

import maya.cmds as cmds
import maya.api.OpenMaya as om


def matrix(x=0, y=0, z=0):
    return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]


def distance(point, reference):
    return math.sqrt(sum((value - target) ** 2 for value, target in zip(point, reference)))


class DeformerTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def mesh(self, points):
        fn = om.MFnMesh()
        obj = fn.create(om.MPointArray([om.MPoint(*point) for point in points]), [len(points)], list(range(len(points))))
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
        node = cmds.deformer(mesh, type="yddSkirtCollideDeformer")[0]
        for side, x in [("left", 0), ("right", 10)]:
            for joint, y in [("Hip", 0), ("Knee", 1), ("Heel", 2)]:
                cmds.setAttr(node + "." + side + joint + "Matrix", *matrix(x, y), type="matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        cmds.setAttr(node + ".skirtType", int(long))
        cmds.setAttr(node + ".falloff", 0)
        return node

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

    def test_collision_endpoints(self):
        for long in [False, True]:
            for fade in [0, 0.1]:
                with self.subTest(long=long, fade=fade):
                    end = 2 if long else 1
                    outside_offset = end * fade * 0.5 if fade else 0.001
                    mesh = self.mesh(
                        [(0.25, end - 0.001, 0), (0.25, end, 0), (0.25, end + outside_offset, 0), (0.25, end + 0.3, 0)]
                    )
                    original = self.points(mesh)
                    node = self.collision(mesh, long)
                    cmds.setAttr(node + ".endFade", fade)
                    result = self.points(mesh)
                    self.assertAlmostEqual(result[0][0], 1, places=6)
                    self.assertAlmostEqual(result[1][0], 1, places=6)
                    if fade:
                        self.assertGreater(result[2][0], 0.25)
                        self.assertLess(result[2][0], 1)
                        # The straight Long leg has two overlapping full-length rings.
                        expected = 0.8125 if long else 0.625
                        self.assertAlmostEqual(result[2][0], expected, delta=1e-6)
                    else:
                        self.assertEqual(result[2], original[2])
                    self.assertEqual(result[3], original[3])
                    for attr in ["collision", "envelope"]:
                        cmds.setAttr(node + "." + attr, 0)
                        self.assertEqual(self.points(mesh), original)
                        cmds.setAttr(node + "." + attr, 1)
                    cmds.setAttr(node + ".weightList[0].weights[1]", 0)
                    self.assertEqual(self.points(mesh)[1], original[1])

    def set_profile(self, node, **values):
        for name, value in values.items():
            cmds.setAttr(node + "." + name, value)

    def bent_collision(self, mesh):
        node = self.collision(mesh, long=True)
        for side, x in (("left", 0), ("right", 10)):
            cmds.setAttr(node + "." + side + "KneeMatrix", *matrix(x, 1, 3), type="matrix")
        return node

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

    def test_leg_profile_knee_shrinks_push(self):
        mesh = self.mesh([(0.2, 0.1, 0), (0.2, 1, 0), (0.2, 0.5, 0.1)])
        node = self.collision(mesh)
        self.set_profile(node, kneeRadiusX=0.7, kneeRadiusZ=0.7)
        result = self.points(mesh)
        self.assertAlmostEqual(result[0][0], 1, delta=1e-6)
        self.assertAlmostEqual(result[1][0], 0.7, delta=1e-6)
        self.assertLess(result[1][0] - 0.2, result[0][0] - 0.2)
        self.set_profile(node, kneeRadiusX=1, kneeRadiusZ=1)
        self.assertAlmostEqual(self.points(mesh)[1][0], 1, delta=1e-6)

    def test_leg_profile_ankle_extension_with_end_fade(self):
        leg_length = 2 * math.sqrt(10)
        mesh = self.mesh([(0.25, 1.05 * leg_length, 0), (0.25, leg_length, 0), (0.25, 1.2 * leg_length, 0)])
        node = self.bent_collision(mesh)
        self.set_profile(node, ankleRadiusX=0.5, ankleRadiusZ=0.5, endFade=0.1)
        result = self.points(mesh)
        self.assertAlmostEqual(result[0][0], 0.375, delta=1e-6)
        self.assertAlmostEqual(result[1][0], 0.5, delta=1e-6)
        self.assertAlmostEqual(result[2][0], 0.25, delta=1e-6)

    def test_leg_profile_coincident_thigh_and_calf_stations(self):
        for thigh_position in (0, 1):
            for calf_position in (0, 1):
                with self.subTest(thigh_position=thigh_position, calf_position=calf_position):
                    mesh = self.mesh([(0.1, 1e-6, 0), (0.1, 1, 0), (0.1, 0.5, 0)])
                    node = self.collision(mesh)
                    self.set_profile(
                        node,
                        thighPosition=thigh_position,
                        calfPosition=calf_position,
                        thighRadiusX=1.2,
                        thighRadiusZ=1.2,
                        kneeRadiusX=0.7,
                        kneeRadiusZ=0.7,
                        calfRadiusX=0.8,
                        calfRadiusZ=0.8,
                        ankleRadiusX=0.6,
                        ankleRadiusZ=0.6,
                        endFade=0,
                    )
                    result = self.points(mesh)
                    expected_near_hip = (1.2 if thigh_position == 0 else 1) * (1 - 1e-6) + 0.7 * 1e-6
                    self.assertAlmostEqual(result[0][0], expected_near_hip, delta=1e-6)
                    self.assertAlmostEqual(result[1][0], 0.7, delta=1e-6)
                    self.assertTrue(all(math.isfinite(v) for point in result for v in point))
                    leg_length = 2 * math.sqrt(10)
                    station_y = leg_length * (0.5 + 0.5 * calf_position)
                    mesh = self.mesh([(0.1, station_y, 0), (0.2, station_y, 0), (0.1, station_y, 0.1)])
                    node = self.bent_collision(mesh)
                    self.set_profile(
                        node,
                        thighPosition=thigh_position,
                        calfPosition=calf_position,
                        thighRadiusX=1.2,
                        thighRadiusZ=1.2,
                        kneeRadiusX=0.7,
                        kneeRadiusZ=0.7,
                        calfRadiusX=0.8,
                        calfRadiusZ=0.8,
                        ankleRadiusX=0.6,
                        ankleRadiusZ=0.6,
                        endFade=0.1,
                    )
                    result = self.points(mesh)
                    self.assertAlmostEqual(result[0][0], 0.7 if calf_position == 0 else 0.6, delta=1e-6)
                    self.assertTrue(all(math.isfinite(v) for point in result for v in point))

    def test_leg_profile_anisotropic_chord_and_push_direction(self):
        original = [(0.2, 0.5, 0.15), (0.1, 0.5, 0.3), (0.3, 0.5, -0.2)]
        mesh = self.mesh(original)
        node = self.collision(mesh)
        self.set_profile(node, thighRadiusX=1.4, thighRadiusZ=0.6, collision=1, falloff=0)
        cmds.setAttr(node + ".bellMatrix", *matrix(-0.4, 0, 0.2), type="matrix")
        for point, source in zip(self.points(mesh), original):
            x, y, z = point
            self.assertAlmostEqual((z / 1.4) ** 2 + (x / 0.6) ** 2, 1, delta=1e-6)
            self.assertAlmostEqual(y, source[1], delta=1e-6)
            dx, dz = x - source[0], z - source[2]
            self.assertAlmostEqual(dx * (source[2] - 0.2) - dz * (source[0] + 0.4), 0, delta=1e-6)
            self.assertGreater(dx * (source[0] + 0.4) + dz * (source[2] - 0.2), 0)

    def test_leg_profile_calf_changes_long_hem(self):
        leg_length = 2 * math.sqrt(10)
        mesh = self.mesh([(0.2, s * leg_length, 0) for s in (0.7, 0.75, 0.8)])
        node = self.bent_collision(mesh)
        baseline = self.points(mesh)
        self.set_profile(node, calfRadiusX=0.8, calfRadiusZ=0.8)
        result = self.points(mesh)
        for point, previous, expected in zip(result, baseline, (0.84, 0.8, 0.84)):
            self.assertAlmostEqual(previous[0], 1, delta=1e-6)
            self.assertAlmostEqual(point[0], expected, delta=1e-6)
            self.assertLess(point[0], previous[0])

    def test_leg_profile_extended_keeps_knee_radius(self):
        mesh = self.mesh([(0.2, 1.5, 0), (0.2, 1.6, 0), (0.2, 1.7, 0)])
        node = self.collision(mesh, long=True)
        self.set_profile(
            node, kneeRadiusX=0.9, kneeRadiusZ=0.9, calfRadiusX=0.5, calfRadiusZ=0.5, ankleRadiusX=0.4, ankleRadiusZ=0.4
        )
        for point in self.points(mesh):
            self.assertAlmostEqual(point[0], 0.9, delta=1e-6)

    def test_leg_profile_degenerate_lengths_and_effective_length_guard(self):
        for calf_length in (0, 1):
            mesh = self.mesh([(0.1, 1, 0), (0.2, 1, 0), (0.1, 1, 0.1)])
            node = self.collision(mesh)
            for side, x in (("left", 0), ("right", 10)):
                cmds.setAttr(node + "." + side + "HeelMatrix", *matrix(x, 1 + calf_length), type="matrix")
            self.set_profile(node, kneeRadiusX=0.7, kneeRadiusZ=0.7, ankleRadiusX=1.3, ankleRadiusZ=1.3)
            self.assertAlmostEqual(self.points(mesh)[0][0], 1.3 if calf_length == 0 else 0.7, delta=1e-6)
        mesh = self.mesh([(0.1, 1e-5, 0), (0.1, 1, 0), (0.1, 2, 0)])
        node = self.collision(mesh, long=True)
        for side, x in (("left", 0), ("right", 10)):
            cmds.setAttr(node + "." + side + "KneeMatrix", *matrix(x, 0), type="matrix")
        self.set_profile(
            node,
            kneeRadiusX=0.7,
            kneeRadiusZ=0.7,
            calfRadiusX=1.2,
            calfRadiusZ=1.2,
            ankleRadiusX=0.9,
            ankleRadiusZ=0.9,
            endFade=0,
        )
        for point, expected in zip(self.points(mesh), (0.700005, 1.2, 0.9)):
            self.assertAlmostEqual(point[0], expected, delta=1e-6)
        mesh = self.mesh([(0.2, 0.5e-7, 0), (0.3, 0.5e-7, 0), (0.2, 0.5e-7, 0.1)])
        node = self.collision(mesh)
        self.set_profile(node, thighPosition=0, thighRadiusX=0.2, thighRadiusZ=0.2, kneeRadiusX=0.2, kneeRadiusZ=0.2, endFade=0)
        cmds.setAttr(node + ".ringScale1", 1e-7)
        self.assertAlmostEqual(self.points(mesh)[0][0], 1, delta=1e-6)

    def test_leg_profile_falloff_band_is_relative_per_axis(self):
        # With ring axis Y the ring frame X lies along world Z and the frame Z along
        # world X. Points at the same normalized radius 1.1 on both frame axes sit
        # inside the falloff band of their own axis and must receive the same
        # normalized shift; a point at absolute radius 1.1 along world X (frame Z,
        # multiplier 0.6) is far outside the band and must not move.
        mesh = self.mesh([(0, 0.5, 1.4 * 1.1), (0.6 * 1.1, 0.5, 0), (1.1, 0.5, 0)])
        node = self.collision(mesh)
        cmds.setAttr(node + ".falloff", 0.2)
        self.set_profile(node, thighRadiusX=1.4, thighRadiusZ=0.6, kneeRadiusX=1.4, kneeRadiusZ=0.6, endFade=0)
        points = self.points(mesh)
        shift_frame_x = (points[0][2] - 1.4 * 1.1) / 1.4
        shift_frame_z = (points[1][0] - 0.6 * 1.1) / 0.6
        self.assertGreater(shift_frame_x, 1e-4)
        self.assertAlmostEqual(shift_frame_x, shift_frame_z, delta=1e-6)
        self.assertAlmostEqual(points[2][0], 1.1, delta=1e-6)

    def test_leg_profile_zero_calf_extended_ring_uses_ankle(self):
        # With the heel on the knee the knee and ankle nodes coincide and the ankle
        # node wins, so every ring including the extended one reads the ankle value.
        mesh = self.mesh([(0.3, 1.0, 0), (0.3, 0.5, 0), (0.3, 0.9, 0.05)])
        node = self.collision(mesh, long=True)
        for side, x in (("left", 0), ("right", 10)):
            cmds.setAttr(node + "." + side + "HeelMatrix", *matrix(x, 1), type="matrix")
        self.set_profile(node, kneeRadiusX=0.7, kneeRadiusZ=0.7, ankleRadiusX=0.5, ankleRadiusZ=0.5, endFade=0)
        points = self.points(mesh)
        self.assertAlmostEqual(points[0][0], 0.5, delta=1e-6)
        self.assertTrue(all(math.isfinite(value) for point in points for value in point))

    def test_leg_profile_length_scale_moves_sampling(self):
        mesh = self.mesh([(0.2, 0.5, 0), (0.2, 0.75, 0), (0.2, 1, 0)])
        node = self.collision(mesh)
        self.set_profile(node, thighRadiusX=1.2, thighRadiusZ=1.2, kneeRadiusX=0.7, kneeRadiusZ=0.7)
        expected = self.points(mesh)
        cmds.setAttr(node + ".ringScale1", 2)
        self.assert_points_close(self.points(mesh), expected)

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
        for name, maximum in [("idleAmplitudeV", 2), ("idleAmplitudeU", 2), ("phaseSpread", 0.5), ("impulseAmount", 2)]:
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
                    bell = om.MMatrix([entry * scale if index < 12 else entry for index, entry in enumerate(list(rotation))])
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
            node, idleAmplitude=0, impulseAmount=1, impulseX=-0.5, impulseZ=0.25, impulsePosition=0.5, directionality=0.7
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
                self.assertAlmostEqual(distance(tilted[0], original[0]), reference_delta * abs(math.cos(angle)), delta=1e-6)
                self.assertEqual(tilted[2:], original[2:])
        self.assert_points_close(self.points(mesh), original)

    def test_wave_rotation_bell_local_ignores_rotation(self):
        mesh = self.mesh([(1, 1, 0), (-1, 1, 0), (0, 1, 1), (0, 1, -1)])
        node = self.directional_wave(mesh)
        self.set_wave_values(node, idleDirectionSpace=1, impulseSpace=1, impulseAmount=1, impulseX=-0.5, impulsePosition=0.5)
        local = self.points(mesh)
        for values in [self.rotation_matrix(y=1.0), [2 if i in (0, 5, 10) else 0 for i in range(16)], [float("nan")] * 16]:
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
                    self.assert_values_close(cmds.getAttr(node + ".evaluationToWorldRotation"), rotation, tolerance=1e-12)
                    self.assert_points_close(self.points(mesh), expected, tolerance=1e-12)
                    self.assertEqual(
                        cmds.listConnections(connected + ".evaluationToWorldRotation", source=True, destination=False), [source]
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
            combined, [tuple(v + u - p for v, u, p in zip(a, b, c)) for a, b, c in zip(vertical, around, original)]
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
                distance(spread[index], original[index]), distance(spread[index + 1], original[index + 1]), delta=1e-6
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
            doubled, [tuple(2 * a - b for a, b in zip(point, reference)) for point, reference in zip(single, original)]
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
                cmds.setKeyframe(node + "." + name, time=time, value=value, inTangentType="linear", outTangentType="linear")
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
