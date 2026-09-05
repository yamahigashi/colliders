import math
import unittest

import maya.cmds as cmds
import maya.api.OpenMaya as om


def matrix(x=0, y=0, z=0):
    return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]


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
        node = cmds.deformer(mesh, type="skirtCollideDeformer")[0]
        for side, x in [("left", 0), ("right", 10)]:
            for joint, y in [("Hip", 0), ("Knee", 1), ("Heel", 2)]:
                cmds.setAttr(node + "." + side + joint + "Matrix", *matrix(x, y), type="matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        cmds.setAttr(node + ".skirtType", int(long))
        cmds.setAttr(node + ".falloff", 0)
        return node

    def wave(self, mesh):
        node = cmds.deformer(mesh, type="skirtWaveDeformer")[0]
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
