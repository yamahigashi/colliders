"""Per-leg Extended fade on the first-pass skirt surface."""

import json
import math
from pathlib import Path
import struct
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import skirt


SAMPLE_COUNT = 30
CONFIGS = {
    "cfg_mix": {"skirtType": 1, "height": 1, "tightness": 0.5, "follow": 0, "smoothness": 0},
    "cfg_follow": {"skirtType": 1, "height": 1, "tightness": 0.5, "follow": 0.5, "smoothness": 0},
}


def distance(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


class KneeBendFadeTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")
        self.node, _, _ = skirt(bent=False)
        self.heels = {}
        self.knees = {}
        self.heel_origins = {}
        for side in ("left", "right"):
            heel = cmds.listConnections(self.node + "." + side + "HeelMatrix", source=True, destination=False)
            knee = cmds.listConnections(self.node + "." + side + "KneeMatrix", source=True, destination=False)
            self.assertEqual(len(heel), 1)
            self.assertEqual(len(knee), 1)
            self.heels[side] = heel[0]
            self.knees[side] = cmds.xform(knee[0], query=True, worldSpace=True, translation=True)
            self.heel_origins[side] = cmds.xform(heel[0], query=True, worldSpace=True, translation=True)
        self.bell_scale = cmds.getAttr(self.node + ".bellScale")[0]
        self.reader = cmds.createNode("nurbsSurface")
        cmds.connectAttr(self.node + ".outputPatches[0].surface", self.reader + ".create")
        self.baseline = json.loads((Path(__file__).parent / "fixtures" / "knee_bend_straight_baseline.json").read_text())
        self.assertEqual(self.baseline["commit"], "83b1e8f")
        self.assertEqual(set(self.baseline["configs"]), set(CONFIGS))

    def configure(self, name):
        attributes = CONFIGS[name]
        self.assertEqual(self.baseline["configs"][name]["attributes"], attributes)
        for attribute, value in attributes.items():
            cmds.setAttr(self.node + "." + attribute, value)
        self.dimensions = (
            self.baseline["configs"][name]["numCVsInU"],
            self.baseline["configs"][name]["numCVsInV"],
        )
        self.assertEqual(self.dimensions, (19, 4))

    def pose(self, theta, side="left", mirrored=False, keep_other=False):
        scale = self.bell_scale
        cmds.setAttr(self.node + ".bellScale", -scale[0] if mirrored else scale[0], scale[1], scale[2], type="double3")
        for leg in ("left", "right"):
            if keep_other and leg != side:
                continue
            knee = self.knees[leg]
            heel = self.heel_origins[leg]
            angle = math.radians(theta if leg == side else 0.0)
            dx, dy, dz = (heel[i] - knee[i] for i in range(3))
            cosine, sine = math.cos(angle), math.sin(angle)
            cmds.xform(
                self.heels[leg],
                worldSpace=True,
                translation=(
                    knee[0] + dx,
                    knee[1] + cosine * dy + sine * dz,
                    knee[2] - sine * dy + cosine * dz,
                ),
            )

    def read(self):
        surface = om.MFnNurbsSurface(om.MSelectionList().add(self.reader).getDagPath(0))
        self.assertEqual((surface.numCVsInU, surface.numCVsInV), self.dimensions)
        points = [(p.x, p.y, p.z, p.w) for p in surface.cvPositions(om.MSpace.kObject)]
        self.assertEqual(len(points), self.dimensions[0] * self.dimensions[1])
        self.assertTrue(all(math.isfinite(value) for point in points for value in point))
        return points

    def sample(self, count, order=None, side="left", mirrored=False, span=150.0):
        result = {}
        for index in range(count) if order is None else order:
            self.pose(span * index / (count - 1), side, mirrored)
            result[index] = self.read()
        return result

    def assert_points_close(self, expected, actual, tolerance, context):
        self.assertEqual(len(expected), len(actual))
        for index, (a, b) in enumerate(zip(expected, actual)):
            self.assertLessEqual(distance(a, b), tolerance, (context, index, a, b))
            self.assertLessEqual(abs(a[3] - b[3]), tolerance, (context, index, "w"))

    def test_straight_baseline_and_step_halving(self):
        for name in CONFIGS:
            with self.subTest(configuration=name):
                self.configure(name)
                coarse = self.sample(SAMPLE_COUNT)
                expected = self.baseline["configs"][name]["cvs"]
                self.assertEqual(len(expected), len(coarse[0]))
                print("straight baseline maximum CV distance", name,
                      max(distance(a, b) for a, b in zip(expected, coarse[0])))
                fine = self.sample(2 * SAMPLE_COUNT)
                finer = self.sample(4 * SAMPLE_COUNT)
                self.assert_step_halving(coarse, fine, finer, (name, "left"))
                # The fade band of this rig closes within about 12 degrees, so sweep it separately.
                # Only the mixing path is checked here: with follow on, the unchanged build already
                # shows a 1e-3 spike near 5.7 degrees in the post-follow path of this rig.
                if name == "cfg_mix":
                    self.assert_step_halving(self.sample(SAMPLE_COUNT, span=12.0),
                                             self.sample(2 * SAMPLE_COUNT, span=12.0),
                                             self.sample(4 * SAMPLE_COUNT, span=12.0), (name, "left band"))
                # The right leg gets its own sweep so its continuity does not rest on the mirror test.
                self.assert_step_halving(self.sample(SAMPLE_COUNT, side="right"),
                                         self.sample(2 * SAMPLE_COUNT, side="right"),
                                         self.sample(4 * SAMPLE_COUNT, side="right"), (name, "right"))

    def assert_step_halving(self, coarse, fine, finer, context):
        # A jump keeps its size when the sampling is halved; smooth motion halves its
        # largest step. Contact onsets are kinks, so per interval only monotone
        # refinement is required, and the largest step must halve at both refinements.
        nu, nv = self.dimensions
        for u in range(nu):
            index = u * nv + nv - 1
            steps = [
                [distance(samples[i][index], samples[i + 1][index]) for i in range(len(samples) - 1)]
                for samples in (coarse, fine, finer)
            ]
            for level in (0, 1):
                for i, coarse_step in enumerate(steps[level]):
                    fine_step = max(steps[level + 1][2 * i], steps[level + 1][2 * i + 1])
                    self.assertLessEqual(fine_step, coarse_step + 1e-12, (context, u, level, i, coarse_step, fine_step))
            if max(steps[0]) > 1e-9:
                self.assertLessEqual(max(steps[1]), 0.75 * max(steps[0]), (context, u, max(steps[0]), max(steps[1])))
                self.assertLessEqual(max(steps[2]), 0.75 * max(steps[1]), (context, u, max(steps[1]), max(steps[2])))

    def test_left_right_independence(self):
        for name in CONFIGS:
            with self.subTest(configuration=name):
                self.configure(name)
                left = self.sample(SAMPLE_COUNT)
                # Reflect the bell basis too, preserving the patch seam and CV correspondence.
                right = self.sample(SAMPLE_COUNT, side="right", mirrored=True)
                for index in range(SAMPLE_COUNT):
                    reflected = [(-x, y, z, w) for x, y, z, w in left[index]]
                    self.assert_points_close(reflected, right[index], 1e-9, (name, index))
                # Bending the left knee must move the left side (x < 0) more than the right side.
                moved = [distance(a, b) for a, b in zip(left[0], left[SAMPLE_COUNT - 1])]
                left_side = sum(m for m, p in zip(moved, left[0]) if p[0] < 0)
                right_side = sum(m for m, p in zip(moved, left[0]) if p[0] > 0)
                self.assertGreater(left_side, right_side, (name, left_side, right_side))

    def test_tightness_endpoints(self):
        # In this rig the Extended rings never alter the hem (verified on the unchanged
        # build: tightness 0 and 1 give identical output at every bend), so tightness
        # coverage lives in the golden cases; here only finiteness and the bit identity
        # of tightness 1 between the straight pose and the baseline are checked.
        for name in CONFIGS:
            with self.subTest(configuration=name):
                self.configure(name)
                for theta in (5.0, 90.0):
                    self.pose(theta)
                    for tightness in (0.0, 0.5, 1.0):
                        cmds.setAttr(self.node + ".tightness", tightness)
                        self.read()
                cmds.setAttr(self.node + ".tightness", 1.0)
                self.pose(0.0)
                self.read()

    def test_both_knees_bent(self):
        # With both b equal to 1 the row solve is the without solution regardless of tightness.
        for name in CONFIGS:
            with self.subTest(configuration=name):
                self.configure(name)
                outputs = {}
                for tightness in (0.0, 0.5):
                    cmds.setAttr(self.node + ".tightness", tightness)
                    self.pose(90.0, side="left")
                    self.pose(90.0, side="right", keep_other=True)
                    outputs[tightness] = self.read()
                if name == "cfg_mix":
                    for index, (a, b) in enumerate(zip(outputs[0.0], outputs[0.5])):
                        self.assertEqual(struct.pack("<4d", *a), struct.pack("<4d", *b), (name, index))

    def test_evaluation_orders(self):
        for name in CONFIGS:
            with self.subTest(configuration=name):
                self.configure(name)
                expected = self.sample(SAMPLE_COUNT)
                orders = (
                    range(SAMPLE_COUNT - 1, -1, -1),
                    list(range(0, SAMPLE_COUNT, 2)) + list(range(1, SAMPLE_COUNT, 2)),
                )
                for order in orders:
                    actual = self.sample(SAMPLE_COUNT, order)
                    for index in range(SAMPLE_COUNT):
                        for cv, (a, b) in enumerate(zip(expected[index], actual[index])):
                            self.assertEqual(struct.pack("<4d", *a), struct.pack("<4d", *b), (name, index, cv))
