"""Cancellation continuity of the bell output curve."""

import math
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from test_solver import _bell, _points


def _ring_matrix(tx):
    matrix = om.MTransformationMatrix()
    matrix.setScale((2, 1, 2), om.MSpace.kTransform)
    matrix.setTranslation(om.MVector(tx, 0, 0), om.MSpace.kTransform)
    return list(matrix.asMatrix())


class MergeContinuityTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")

    def _sweep(self, node, divisions):
        frames = []
        for frame in range(divisions + 1):
            tx = 1.5 + frame / float(divisions)
            cmds.setAttr(node + ".ringMatrix[1]", *_ring_matrix(tx), type="matrix")
            points = _points(node, "outputCurve")
            self.assertEqual(len(points), 17)
            self.assertEqual(points[0], points[-1])
            self.assertTrue(all(math.isfinite(value) for point in points for value in point))
            for index, point in enumerate(points):
                angle = (index % 16) * 2.0 * math.pi / 16.0
                for center in (0.0, tx):
                    self.assertGreater(math.hypot(math.cos(angle) - center, math.sin(angle)), 1e-5)
                    self.assertGreater(math.hypot(point[0] - center, point[2]), 1e-5)
            frames.append(points)
        return frames

    @staticmethod
    def _steps(frames, index):
        return [
            tuple(b - a for a, b in zip(frames[i][index], frames[i + 1][index]))
            for i in range(len(frames) - 1)
        ]

    def test_bell_cancellation_sweep(self):
        node = _bell([0, 1], [_ring_matrix(0), _ring_matrix(1.5)], collision=1)
        cmds.setAttr(node + ".bellSubdivision", 16)
        cmds.setAttr(node + ".falloff", 0)
        coarse = self._sweep(node, 40)
        fine = self._sweep(node, 80)

        # A jump survives halving the sweep step; continuous motion halves its largest step.
        measured = 0
        for index in range(17):
            steps = self._steps(coarse, index)
            lengths = [math.sqrt(sum(value * value for value in step)) for step in steps]
            mean = sum(lengths) / len(lengths)
            if mean > 1e-3:
                measured += 1
                fine_lengths = [
                    math.sqrt(sum(value * value for value in step)) for step in self._steps(fine, index)
                ]
                self.assertLessEqual(max(fine_lengths), 0.75 * max(lengths),
                                     (index, max(fine_lengths), max(lengths)))
                for first, second in zip(steps, steps[1:]):
                    self.assertGreaterEqual(sum(a * b for a, b in zip(first, second)), 0.0,
                                            (index, first, second))
        self.assertGreater(measured, 0)
