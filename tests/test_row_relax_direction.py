"""Direction-aware row relaxation on recorded matrix inputs."""

import json
import math
from pathlib import Path
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

import yddColliders
from helpers import transform


class RowRelaxDirectionTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")

    def test_circular_vibration_sweep(self):
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
        maximum_ratio = 0.0
        target_step = 0.0
        for u in range(16):
            for v in range(4):
                index = u * 4 + v
                steps = [tuple(b - a for a, b in zip(forward[i][index], forward[i + 1][index])) for i in range(40)]
                lengths = [math.sqrt(sum(x * x for x in step)) for step in steps]
                mean = sum(lengths) / len(lengths)
                if (u, v) == (12, 2):
                    target_step = max(lengths)
                    self.assertLessEqual(target_step, 1.2, (u, v, lengths))
                if mean > 1e-3:
                    ratio = max(lengths) / mean
                    maximum_ratio = max(maximum_ratio, ratio)
                    self.assertLessEqual(ratio, 3.0, (u, v, ratio))
                    for a, b in zip(steps, steps[1:]):
                        self.assertGreaterEqual(sum(x * y for x, y in zip(a, b)), 0.0, (u, v, a, b))
        print("row direction maximum ratio", maximum_ratio, "raw cv[12][2] maximum step", target_step)
