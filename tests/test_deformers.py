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
