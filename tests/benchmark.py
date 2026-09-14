"""Headless DG timing; run baseline and candidate separately without competing jobs."""

import json
from pathlib import Path
import platform
import statistics
import time

from helpers import cylinder, maya_session, output_object, parser, skirt, transform, wave_settings


def measure(node, attribute, dirty, values, iterations):
    import maya.api.OpenMaya as om

    driver = om.MSelectionList().add(node + "." + dirty).getPlug(0)
    output = om.MSelectionList().add(node + "." + attribute).getPlug(0)

    def evaluate(i):
        driver.setDouble(values[i % 2])
        output.asMObject()

    for i in range(20):
        evaluate(i)
    batches = []
    for batch in range(5):
        start = time.perf_counter_ns()
        for i in range(iterations):
            evaluate(batch * iterations + i)
        batches.append((time.perf_counter_ns() - start) / iterations / 1e6)
    return dict(batch_ms=batches, median_ms=statistics.median(batches))


def run(iterations, node_filter=None):
    from maya import cmds

    rows = []
    for subdivisions in (16, 64, 256):
        if node_filter not in (None, "yddSkirtBellCollider"):
            continue
        for tightness, smoothness, follow in ((0.5, 0, 0), (1, 0, 0), (0.5, 0.5, 0.5)):
            cmds.file(new=True, force=True)
            node, _, _ = skirt(bent=True)
            attrs = dict(bellSubdivision=subdivisions, tightness=tightness, smoothness=smoothness, follow=follow)
            for attr, value in attrs.items():
                cmds.setAttr(node + "." + attr, value)
            rows.append(
                dict(
                    node="yddSkirtBellCollider",
                    attributes=attrs,
                    **measure(node, "outputSurface", "height", (0.7, 0.7001), iterations),
                )
            )
    for size in (133, 1088, 16512):
        if node_filter not in (None, "yddSkirtWaveDeformer"):
            continue
        for condition in ("active", "amplitude0", "envelope0", "signals0", "impulse_only", "periodic_only"):
            cmds.file(new=True, force=True)
            if size == 133:
                _, locator, obj = skirt()
                cmds.delete(obj, constructionHistory=True)
            else:
                obj = cylinder(size)
                locator = transform()
            node = cmds.deformer(obj, type="yddSkirtWaveDeformer")[0]
            cmds.connectAttr(locator + ".worldMatrix[0]", node + ".bellMatrix")
            attrs = wave_settings(node, condition)
            import maya.api.OpenMaya as om

            data = output_object(node, "outputGeometry[0]")
            count = len(om.MFnNurbsSurface(data).cvPositions()) if size == 133 else len(om.MFnMesh(data).getPoints())
            if count != size:
                raise AssertionError(f"Expected {size} points, got {count}")
            rows.append(
                dict(
                    node="yddSkirtWaveDeformer",
                    points=size,
                    condition=condition,
                    attributes=attrs,
                    **measure(node, "outputGeometry[0]", "wavePhaseV", (0, 0.0001), iterations),
                )
            )
    for size in (80, 1088, 16512):
        if node_filter not in (None, "yddSkirtCollideDeformer"):
            continue
        for condition in ("active", "collision0", "envelope0"):
            cmds.file(new=True, force=True)
            obj = cylinder(size)
            node = cmds.deformer(obj, type="yddSkirtCollideDeformer")[0]
            cmds.setAttr(node + ".ringScale", 0.7, 1, 0.7)
            cmds.setAttr(node + ".falloff", 0.2)
            for side, x in [("left", -1), ("right", 1)]:
                for part, y in [("Hip", 10), ("Knee", 5), ("Heel", 0)]:
                    joint = transform((x, y, 0))
                    cmds.connectAttr(joint + ".worldMatrix[0]", node + "." + side + part + "Matrix")
                cmds.setAttr(node + "." + side + "RingAxis", 4)
            if condition != "active":
                cmds.setAttr(node + "." + condition[:-1], 0)
            rows.append(
                dict(
                    node="yddSkirtCollideDeformer",
                    points=size,
                    condition=condition,
                    **measure(node, "outputGeometry[0]", "endFade", (0.1, 0.1001), iterations),
                )
            )
    return rows


def main():
    p = parser(__doc__)
    p.add_argument("--output", required=True)
    p.add_argument("--iterations", type=int, default=40)
    p.add_argument("--node", choices=("yddSkirtBellCollider", "yddSkirtWaveDeformer", "yddSkirtCollideDeformer"))
    args = p.parse_args()
    if args.iterations < 2 or args.iterations % 2:
        p.error("--iterations must be a positive even number >= 2")
    with maya_session(args.plugin) as metadata:
        result = dict(
            metadata=metadata,
            cpu=platform.processor(),
            evaluation_mode="off",
            warmup=20,
            batches=5,
            iterations=args.iterations,
            results=run(args.iterations, args.node),
        )
        Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
