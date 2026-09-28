import math
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds


class WaveNoiseInputTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def mesh(self):
        fn = om.MFnMesh()
        points = om.MPointArray([
            om.MPoint(1.0, 0.2, 0.1),
            om.MPoint(0.7, 0.5, 0.5),
            om.MPoint(0.1, 1.0, 1.0),
        ])
        obj = fn.create(points, [3], [0, 1, 2])
        return om.MFnDagNode(obj).fullPathName()

    def evaluate(self, frequency_u, frequency_v, phase, noise_amplitude=0.2):
        mesh = self.mesh()
        node = cmds.deformer(mesh, type="yddSkirtWaveDeformer")[0]
        cmds.setAttr(node + ".noiseAmplitude", noise_amplitude)
        cmds.setAttr(node + ".idleAmplitude", 0.0)
        cmds.setAttr(node + ".impulseX", 0.0)
        cmds.setAttr(node + ".impulseZ", 0.0)
        cmds.setAttr(node + ".noiseFrequencyU", frequency_u)
        cmds.setAttr(node + ".noiseFrequencyV", frequency_v)
        cmds.setAttr(node + ".noisePhase", phase)

        shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
        result = om.MFnMesh(om.MSelectionList().add(shape).getDagPath(0)).getPoints()
        return [component for point in result for component in point]

    def test_extreme_finite_noise_inputs_remain_finite(self):
        baseline = self.evaluate(1.5, 1.5, 0.0, noise_amplitude=0.0)
        for phase in (1.0e300, -1.0e300):
            with self.subTest(phase=phase):
                result = self.evaluate(1.0e300, 1.0e300, phase)
                self.assertTrue(all(math.isfinite(value) for value in result))
                self.assertTrue(any(abs(a - b) > 1.0e-12 for a, b in zip(result, baseline)))

    def test_nonfinite_noise_input_is_a_noop_for_noise(self):
        baseline = self.evaluate(1.5, 1.5, 0.0, noise_amplitude=0.0)
        result = self.evaluate(float("nan"), float("inf"), float("-inf"))
        self.assertTrue(all(math.isfinite(value) for value in result))
        self.assertEqual(result, baseline)


if __name__ == "__main__":
    unittest.main()
