"""Stable skirt output fixtures; golden double hashes target Maya 2026."""

import hashlib
import json
import math
from pathlib import Path
import struct
import unittest

import maya.api.OpenMaya as om
from maya import cmds

from helpers import output_object, skirt

FIXTURE = json.loads((Path(__file__).parent / "fixtures" / "compatible_outputs.json").read_text(encoding="utf-8"))


def build_surface(case):
    cmds.file(new=True, force=True)
    cmds.evaluationManager(mode="off")
    node, _, _ = skirt(case["bent"])
    for attribute, value in case["attributes"].items():
        cmds.setAttr(node + "." + attribute, value)
    data = output_object(node, "outputSurface")
    surface = om.MFnNurbsSurface(data)
    points = [tuple(point) for point in surface.cvPositions()]
    seam_distances = [
        (surface.cvPosition(column, row) - surface.cvPosition(surface.numCVsInU - 3 + column, row)).length()
        for column in range(3)
        for row in range(surface.numCVsInV)
    ]
    return points, surface.formInU, seam_distances


class CompatibilityTests(unittest.TestCase):
    def test_finite_periodic_surfaces_all_versions(self):
        for fixture in FIXTURE["cases"]:
            with self.subTest(case=fixture["case"]):
                points, form, seam_distances = build_surface(fixture["case"])
                self.assertEqual(len(points), fixture["point_count"])
                self.assertTrue(all(math.isfinite(value) for point in points for value in point))
                self.assertEqual(form, om.MFnNurbsSurface.kPeriodic)
                self.assertTrue(seam_distances)
                self.assertLessEqual(max(seam_distances), 1e-12)

    def test_maya2026_golden_double_hashes(self):
        if cmds.about(version=True) != FIXTURE["provenance"]["maya"]:
            self.skipTest("Golden double hashes cover Maya 2026 only; geometric checks run on every version")
        for fixture in FIXTURE["cases"]:
            with self.subTest(case=fixture["case"]):
                points, _, _ = build_surface(fixture["case"])
                packed = b"".join(struct.pack("<4d", *point) for point in points)
                self.assertEqual(hashlib.sha256(packed).hexdigest(), fixture["sha256"])
