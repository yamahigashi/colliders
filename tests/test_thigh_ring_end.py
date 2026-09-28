"""Distal ring fade regression on recorded matrix inputs."""

import json
import math
from pathlib import Path
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

import yddColliders
from helpers import transform


class ThighRingEndTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")

    def test_circular_vibration_first_pass(self):
        fixture = json.loads(
            (Path(__file__).parent / "fixtures" / "circular_vibration_rig.json").read_text()
        )
        joints = {}
        for side, x in (("left", -1), ("right", 1)):
            for joint, y in (("Hip", 10), ("Knee", 5), ("Heel", 0)):
                joints[side + joint + "Obj"] = transform((x, y, 0))
        node, _, _ = yddColliders.createSkirtBellCollider(**joints)
        matrices = ["bellMatrix"] + [
            side + joint + "Matrix"
            for side in ("left", "right")
            for joint in ("Hip", "Knee", "Heel")
        ]
        for attribute in matrices:
            plug = node + "." + attribute
            for source in cmds.listConnections(plug, source=True, destination=False, plugs=True) or []:
                cmds.disconnectAttr(source, plug)
        attributes = [
            "skirtType", "height", "bellSubdivision", "ringSubdivision",
            "leftRingAxis", "rightRingAxis", "bellAxis", "thighPosition", "calfPosition",
            "followRange", "referenceMaterialHeight",
        ] + [part + "Radius" + axis for part in ("thigh", "knee", "calf", "ankle") for axis in ("X", "Z")]
        for attribute in attributes:
            cmds.setAttr(node + "." + attribute, fixture["attrs"][attribute])
        for attribute in ("ringScale", "bellScale"):
            cmds.setAttr(node + "." + attribute, *fixture[attribute][0], type="double3")
        for attribute, value in (("falloff", -1), ("tightness", 0.6), ("smoothness", 0.1), ("follow", 0.2)):
            cmds.setAttr(node + "." + attribute, value)
        reader = cmds.createNode("nurbsSurface")
        cmds.connectAttr(node + ".outputPatches[0].surface", reader + ".create")
        frames = fixture["frames"]
        self.assertEqual(len(frames), 41)

        def evaluate(order):
            result = {}
            for index in order:
                frame = frames[index]
                cmds.currentTime(frame["frame"])
                for attribute in matrices:
                    cmds.setAttr(node + "." + attribute, *frame[attribute], type="matrix")
                surface = om.MFnNurbsSurface(om.MSelectionList().add(reader).getDagPath(0))
                self.assertEqual((surface.numCVsInU, surface.numCVsInV), (19, 4))
                points = surface.cvPositions(om.MSpace.kObject)
                result[index] = [(p.x, p.y, p.z) for p in points]
                self.assertTrue(all(math.isfinite(v) for point in result[index] for v in point))
            return result

        forward = evaluate(range(41))
        for order in (range(40, -1, -1), list(range(0, 41, 2)) + list(range(1, 41, 2))):
            actual = evaluate(order)
            for frame in range(41):
                for index, (expected, point) in enumerate(zip(forward[frame], actual[frame])):
                    for a, b in zip(expected, point):
                        self.assertLessEqual(abs(a - b), 1e-12, (frame, index))
        baseline = json.loads(
            (Path(__file__).parent / "fixtures" / "circular_vibration_baseline.json").read_text()
        )
        self.assertEqual(baseline["frames"], [i / 4.0 for i in range(41)])
        self.assertEqual([frame["frame"] for frame in frames], baseline["frames"])
        self.assertEqual(baseline["numCVsInU"], 19)
        self.assertEqual(baseline["numCVsInV"], 4)
        self.assertEqual(len(baseline["cvs"]), 41)
        for frame in baseline["cvs"]:
            self.assertEqual(len(frame), 19)
            for row in frame:
                self.assertEqual(len(row), 4)
                for point in row:
                    self.assertEqual(len(point), 3)
        for time in range(40):
            for u in range(19):
                for v in range(4):
                    index = u * 4 + v
                    a, b = forward[time][index], forward[time + 1][index]
                    old_a = baseline["cvs"][time][u][v]
                    old_b = baseline["cvs"][time + 1][u][v]
                    step_new = math.sqrt(sum((y - x) ** 2 for x, y in zip(a, b)))
                    step_baseline = math.sqrt(sum((y - x) ** 2 for x, y in zip(old_a, old_b)))
                    self.assertLessEqual(
                        step_new - step_baseline, 0.2,
                        (baseline["frames"][time], u, v, step_new, step_baseline),
                    )
