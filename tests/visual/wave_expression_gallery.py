"""Sample real Maya deformer surfaces and write a dependency-free review gallery.

Run with Maya's mayapy: wave_expression_gallery.py --plugin ABSOLUTE --output-dir DIR
No wave formula is implemented here; every displayed vertex is evaluated by Maya.
"""

import argparse
import json
import math
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helpers import maya_session, output_object, skirt  # noqa: E402

NU, NV = 32, 17
PHASES = [i / 16 for i in range(17)]
BASE = dict(
    amplitude=1.5,
    envelope=1,
    idleAmplitude=1,
    idleAmplitudeV=2,
    idleAmplitudeU=0,
    idleDirectionality=1,
    idleDirectionX=-1,
    idleDirectionZ=0,
    idleDirectionSpace=1,
    idleComplexity=0,
    waveCountV=0.75,
    waveCountU=2,
    wavePhaseU=0,
    wavePhaseV=0,
    phaseSpread=0,
    noiseAmplitude=0,
    noiseFrequencyU=3,
    noiseFrequencyV=1,
    noisePhase=0.37,
    impulseAmount=0,
    impulseX=-1,
    impulseZ=0,
    impulsePosition=0.5,
    impulseWidth=0.35,
    impulseSpace=1,
    directionality=0.8,
    skew=0.3,
    sharpness=0.3,
)
PRESETS = [
    ("Axisymmetric V", dict(idleDirectionality=0)),
    ("Directional V / spread 0", {}),
    ("U only", dict(idleAmplitudeV=0, idleAmplitudeU=2)),
    ("Directional V / spread 0.15", dict(phaseSpread=0.15)),
    ("Send amount 0 / moving position", dict(idleAmplitude=0, impulseAmount=0)),
    ("Send amount 1 / moving position", dict(idleAmplitude=0, impulseAmount=1)),
]


def distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def capture():
    from maya import cmds
    import maya.api.OpenMaya as om

    _, bell, surface = skirt()
    cmds.delete(surface, constructionHistory=True)
    # More vertical CV rows make the travelling send visible between samples.
    cmds.rebuildSurface(
        surface,
        constructionHistory=False,
        replaceOriginal=True,
        rebuildType=0,
        spansU=32,
        spansV=16,
        degreeU=3,
        degreeV=3,
        keepRange=0,
        keepCorners=True,
    )
    node = cmds.deformer(surface, type="yddSkirtWaveDeformer")[0]
    cmds.connectAttr(bell + ".worldMatrix[0]", node + ".bellMatrix")
    inverse = om.MMatrix(cmds.getAttr(bell + ".worldMatrix[0]")).inverse()
    shape = cmds.listRelatives(surface, shapes=True, noIntermediate=True, fullPath=True)[0]
    world = om.MMatrix(cmds.xform(surface, query=True, worldSpace=True, matrix=True))
    for index, value in enumerate((0.0, 1.0)):
        prefix = node + ".amplitudeRamp[%d].amplitudeRamp_" % index
        cmds.setAttr(prefix + "Position", value)
        cmds.setAttr(prefix + "FloatValue", value)
        cmds.setAttr(prefix + "Interp", 1)

    def sample():
        obj = output_object(shape, "local")
        fn = om.MFnNurbsSurface(obj)
        u0, u1 = fn.knotDomainInU
        v0, v1 = fn.knotDomainInV
        points = []
        for u in range(NU + 1):
            for v in range(NV):
                p = fn.getPointAtParam(u0 + (u1 - u0) * u / NU, v0 + (v1 - v0) * v / (NV - 1))
                p = p * world * inverse
                points.append([p.x, p.y, p.z])
        cvs = fn.cvPositions()
        cv_points = []
        for i in range(len(cvs)):
            p = om.MPoint(cvs[i]) * world * inverse
            cv_points.append([p.x, p.y, p.z])
        return points, cv_points, [fn.numCVsInU, fn.numCVsInV]

    cmds.setAttr(node + ".envelope", 0)
    rest, rest_cvs, cv_counts = sample()
    cases = []
    for title, overrides in PRESETS:
        settings = dict(BASE, **overrides)
        for attr, value in settings.items():
            cmds.setAttr(node + "." + attr, value)
        frames, cvs, inputs = [], [], []
        send = title.startswith("Send")
        for phase in PHASES:
            varying = (
                dict(impulsePosition=phase * (1 + settings["impulseWidth"]))
                if send
                else dict(wavePhaseV=phase, wavePhaseU=phase)
            )
            for attr, value in varying.items():
                cmds.setAttr(node + "." + attr, value)
            points, cv_points, _ = sample()
            frames.append(points)
            cvs.append(cv_points)
            inputs.append(varying)
        endpoint = max(distance(a, b) for a, b in zip(frames[0], frames[-1]))
        cases.append(
            dict(
                title=title,
                settings=settings,
                frameInputs=inputs,
                frames=frames,
                cvFrames=cvs,
                endpointMaxDistance=endpoint,
                loopExpected=not send,
                finite=all(math.isfinite(x) for frame in frames + cvs for p in frame for x in p),
            )
        )
    return dict(
        coordinateSpace="bell local",
        bellMatrix=list(inverse.inverse()),
        objectMatrix=list(world),
        nu=NU,
        nv=NV,
        phases=PHASES,
        cvCounts=cv_counts,
        cvCount=cv_counts[0] * cv_counts[1],
        sampleCount=(NU + 1) * NV,
        ramp=[[0, 0, "linear"], [1, 1, "linear"]],
        rest=rest,
        restCVs=rest_cvs,
        cases=cases,
    )


def write_gallery(data, output):
    payload = json.dumps(data, allow_nan=False, separators=(",", ":"))
    (output / "samples.json").write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")
    template = Path(__file__).with_name("gallery.html").read_text(encoding="utf-8")
    (output / "index.html").write_text(template.replace("__DATA__", payload), encoding="utf-8")
    # Orthographic X/Y view makes opposite-side signs easy to compare.
    all_points = data["rest"] + [p for case in data["cases"] for frame in case["frames"] for p in frame]
    xmin, xmax = min(p[0] for p in all_points), max(p[0] for p in all_points)
    ymin, ymax = min(p[1] for p in all_points), max(p[1] for p in all_points)
    scale = min(210 / (xmax - xmin), 270 / (ymax - ymin))
    rows = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="2220" viewBox="0 0 1200 2220">',
        '<rect width="1200" height="2220" fill="#101721"/>',
        '<text x="20" y="25" fill="white" font-family="sans-serif">Actual Maya NURBS samples · bell X/Y view · gray = rest · waist above, hem below (+Y)</text>',
    ]
    for row, case in enumerate(data["cases"]):
        for col, frame_index in enumerate((0, 4, 8, 12, 16)):
            cx, top = col * 240 + 120, row * 360 + 55
            rows.append(
                f'<text x="{col * 240 + 8}" y="{top}" fill="white" font-size="11" font-family="sans-serif">{case["title"]}</text>'
            )
            frame_inputs = case["frameInputs"][frame_index]
            key = "wavePhaseV" if case["loopExpected"] else "impulsePosition"
            label = "phase" if case["loopExpected"] else "sendPosition"
            rows.append(
                f'<text x="{col * 240 + 8}" y="{top + 16}" fill="#abbccc" font-size="11">{label} = {frame_inputs[key]:.3f}</text>'
            )
            for points, color in ((data["rest"], "#48515c"), (case["frames"][frame_index], "#69d7e8")):
                for u in range(0, NU + 1, 2):
                    coords = " ".join(
                        f"{cx + (p[0] - (xmin + xmax) / 2) * scale:.2f},{top + 35 + (p[1] - ymin) * scale:.2f}"
                        for p in points[u * NV : (u + 1) * NV]
                    )
                    rows.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="1"/>')
    rows.append("</svg>")
    (output / "contact-sheet.svg").write_text("\n".join(rows), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plugin", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not Path(args.plugin).is_absolute():
        parser.error("--plugin must be absolute")
    # Resolve the native Maya plug-in folder without changing user preferences.
    os.environ["MAYA_PLUG_IN_PATH"] = (
        str(Path(sys.executable).parent / "plug-ins") + os.pathsep + os.environ.get("MAYA_PLUG_IN_PATH", "")
    )
    with maya_session(args.plugin) as metadata:
        data = capture()
        data["metadata"] = metadata
    if not all(case["finite"] for case in data["cases"]):
        raise ValueError("Maya produced non-finite samples")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_gallery(data, args.output_dir)
    print(
        json.dumps(
            dict(
                output=str(args.output_dir),
                metadata=metadata,
                cases=[{k: c[k] for k in ("title", "finite", "endpointMaxDistance", "loopExpected")} for c in data["cases"]],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
