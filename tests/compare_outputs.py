"""Capture in fresh mayapy processes, then compare explicit double-coordinate fixtures."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

from helpers import coordinates, maya_session, skirt, transform, wave_settings


def sample(case, node, attribute, kind="mesh"):
    points = coordinates(node, attribute, kind)
    if not all(math.isfinite(value) for point in points for value in point):
        raise ValueError(f"Non-finite output in {case}")
    packed = b"".join(struct.pack("<4d", *p) for p in points)
    return dict(case=case, points=points, sha256=hashlib.sha256(packed).hexdigest())


def capture(plugin):
    with maya_session(plugin) as metadata:
        from maya import cmds
        import maya.api.OpenMaya as om

        rows = []
        defaults = [(0, 0.4, 0), (0, 1, 0.5), (1, 0.4, 0.5), (1, 1, 0), (1, 1, 0.5), (1, 1, 1)]
        for bent in (False, True):
            for skirt_type, height, tightness in defaults:
                cmds.file(new=True, force=True)
                node, _, _ = skirt(bent)
                attrs = dict(skirtType=skirt_type, height=height, tightness=tightness, smoothness=0, follow=0)
                for attr, value in attrs.items():
                    cmds.setAttr(node + "." + attr, value)
                rows.append(
                    sample(
                        dict(
                            group="review12",
                            bent=bent,
                            attributes=attrs,
                            known_fix="U1 translation correction can change legacy coordinates",
                        ),
                        node,
                        "outputSurface",
                        "surface",
                    )
                )
        for smoothness, follow in [(0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5), (1, 1)]:
            cmds.file(new=True, force=True)
            node, _, _ = skirt(True)
            attrs = dict(skirtType=1, height=1, tightness=1, smoothness=smoothness, follow=follow)
            for attr, value in attrs.items():
                cmds.setAttr(node + "." + attr, value)
            rows.append(sample(dict(group="smoothing", attributes=attrs), node, "outputSurface", "surface"))
        for mode in ("off", "serial", "parallel"):
            for condition in (
                "active",
                "amplitude0",
                "envelope0",
                "signals0",
                "impulse_only",
                "periodic_only",
                "noise_only",
            ):
                for complexity in (0, 0.35, 1):
                    cmds.file(new=True, force=True)
                    cmds.evaluationManager(mode=mode)
                    obj = cmds.polyCylinder(
                        radius=1, height=4, subdivisionsX=12, subdivisionsY=4, constructionHistory=False
                    )[0]
                    cmds.move(0, 2, 0, obj + ".vtx[*]", relative=True, objectSpace=True)
                    node = cmds.deformer(obj, type="yddSkirtWaveDeformer")[0]
                    attrs = wave_settings(node, condition)
                    attrs["idleComplexity"] = complexity
                    cmds.setAttr(node + ".idleComplexity", complexity)
                    for step, (frame, phase) in enumerate([(1, 0), (10, 0.7), (5, 0.2), (1, 0)]):
                        cmds.currentTime(frame)
                        cmds.setAttr(node + ".wavePhaseV", phase)
                        rows.append(
                            sample(
                                dict(
                                    group="wave",
                                    mode=mode,
                                    condition=condition,
                                    attributes=attrs,
                                    frame=frame,
                                    phase=phase,
                                    step=step,
                                ),
                                node,
                                "outputGeometry[0]",
                            )
                        )
        cmds.evaluationManager(mode="off")
        for zero_weight in (False, True):
            cmds.file(new=True, force=True)
            mesh = om.MFnMesh().create([om.MPoint(1, y, z) for y, z in [(1, 0), (2, 0), (4, 0.1)]], [3], [0, 1, 2])
            shape = om.MFnDagNode(mesh).fullPathName()
            node = cmds.deformer(shape, type="yddSkirtWaveDeformer")[0]
            attrs = wave_settings(node, "periodic_only")
            cmds.setAttr(node + ".waveCountU", 0)
            cmds.setAttr(node + ".waveCountV", 0)
            if zero_weight:
                cmds.percent(node, shape + ".vtx[2]", value=0)
            rows.append(
                sample(
                    dict(group="highest_vertex_weight", zero_weight=zero_weight, attributes=attrs),
                    node,
                    "outputGeometry[0]",
                )
            )
        for scale in ((0.4, 1.5, 0.4), (2, 1.5, 2)):
            for y in (0, 10):
                cmds.file(new=True, force=True)
                node = cmds.createNode("yddBellCollider")
                bell = transform((0, y, 0))
                ring = transform((0, y, 0), (0, 0, 45), scale)
                cmds.connectAttr(bell + ".worldMatrix[0]", node + ".bellMatrix")
                cmds.connectAttr(ring + ".worldMatrix[0]", node + ".ringMatrix[0]")
                rows.append(
                    sample(
                        dict(group="U1", translation=y, scale=scale, known_fix="Translation invariance"),
                        node,
                        "outputBellMesh",
                    )
                )
        for skirt_type in (0, 1):
            for falloff in (0, 0.2):
                cmds.file(new=True, force=True)
                points = [(0.25, y, z) for y, z in [(0.05, 0), (0.5, 0.2), (1, 0), (1.5, 0.1), (2.1, 0)]]
                points.extend([(1.1, 0.5, 0), (0.1, 0.75, 1.05)])
                mesh = om.MFnMesh().create([om.MPoint(*p) for p in points], [len(points)], list(range(len(points))))
                shape = om.MFnDagNode(mesh).fullPathName()
                node = cmds.deformer(shape, type="yddSkirtCollideDeformer")[0]
                for side, x in (("left", 0), ("right", 10)):
                    for joint, y in (("Hip", 0), ("Knee", 1), ("Heel", 2)):
                        guide = transform((x, y, 0))
                        cmds.connectAttr(guide + ".worldMatrix[0]", node + "." + side + joint + "Matrix")
                    cmds.setAttr(node + "." + side + "RingAxis", 1)
                cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
                attrs = dict(skirtType=skirt_type, falloff=falloff, endFade=0.1, collision=1)
                for name, value in attrs.items():
                    cmds.setAttr(node + "." + name, value)
                rows.append(
                    sample(
                        dict(group="leg_collision", attributes=attrs, input_points=points), node, "outputGeometry[0]"
                    )
                )
        return dict(metadata=metadata, cases=rows)


def compare(left, right, tolerance):
    if len(left["cases"]) != len(right["cases"]):
        raise ValueError("Fixture case counts differ")
    results = []
    for a, b in zip(left["cases"], right["cases"]):
        if a["case"] != b["case"]:
            raise ValueError("Fixture case definitions differ")
        if len(a["points"]) != len(b["points"]):
            results.append(dict(case=a["case"], exact=False, within_tolerance=False, error="Point counts differ"))
            continue
        if not all(math.isfinite(value) for fixture in (a, b) for point in fixture["points"] for value in point):
            raise ValueError(f"Non-finite coordinates in {a['case']}")
        delta = max((abs(x - y) for p, q in zip(a["points"], b["points"]) for x, y in zip(p, q)), default=0)
        results.append(
            dict(
                case=a["case"],
                exact=a["sha256"] == b["sha256"],
                max_absolute_delta=delta,
                within_tolerance=delta <= tolerance,
            )
        )
    return dict(
        tolerance=tolerance,
        all_exact=all(r["exact"] for r in results),
        all_within_tolerance=all(r["within_tolerance"] for r in results),
        results=results,
    )


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("capture")
    from helpers import add_plugin_argument

    add_plugin_argument(c)
    c.add_argument("--output", required=True)
    d = sub.add_parser("compare")
    d.add_argument("baseline")
    d.add_argument("candidate")
    d.add_argument("--tolerance", type=float, default=0)
    d.add_argument("--output", required=True)
    pair = sub.add_parser("pair")
    pair.add_argument("--mayapy", required=True)
    pair.add_argument("--baseline-plugin", required=True)
    pair.add_argument("--candidate-plugin", required=True)
    pair.add_argument("--output-dir", required=True)
    pair.add_argument("--tolerance", type=float, default=0)
    args = p.parse_args()
    if args.command == "capture":
        write(args.output, capture(args.plugin))
        return 0
    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        p.error("--tolerance must be nonnegative")
    if args.command == "pair":
        directory = Path(args.output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        for label, plugin in [("baseline", args.baseline_plugin), ("candidate", args.candidate_plugin)]:
            subprocess.run(
                [
                    args.mayapy,
                    "-B",
                    str(Path(__file__).resolve()),
                    "capture",
                    "--plugin",
                    plugin,
                    "--output",
                    str(directory / (label + ".json")),
                ],
                check=True,
            )
        args.baseline, args.candidate = directory / "baseline.json", directory / "candidate.json"
        args.output = directory / "comparison.json"
    result = compare(
        json.loads(Path(args.baseline).read_text(encoding="utf-8")),
        json.loads(Path(args.candidate).read_text(encoding="utf-8")),
        args.tolerance,
    )
    write(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != "results"}))
    return 0 if result["all_within_tolerance"] else 1


if __name__ == "__main__":
    sys.exit(main())
