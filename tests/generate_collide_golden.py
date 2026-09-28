"""Capture the golden second-pass surface outputs with Maya 2026."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

from helpers import ROOT, add_plugin_argument, maya_session

EXPECTED_PLUGIN_SHA256 = "3da88a4b3a848e92943616e74d35347aaca8cde0fdd474ae82aac906aa6d4faf"
EXPECTED_MANIFEST_SHA256 = "852619ac24baaf873924aeb38cd9f085232714bd0c935d967aeb2010bb022fb0"
OWNED_SOURCE_PATHS = (
    "sources/skirtLegModel.h",
    "sources/skirtLegBuild.h",
    "sources/skirtCollideDeformer.cpp",
    "sources/skirtCollideDeformer.h",
    "tests/skirt_leg_model_test.cpp",
    "sources/skirtSurfaceCache.h",
)
INPUT_DIR = ROOT / "tests/fixtures/long_skirt_poses"
OUTPUT = ROOT / "tests/fixtures/collide_golden_maya2026.json"
SOURCE_MANIFEST = ROOT / "tests/fixtures/collide_golden_source_manifest.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinate_hash(points):
    return hashlib.sha256(b"".join(struct.pack("<3d", *point) for point in points)).hexdigest()


def validate_source_manifest(manifest_path):
    manifest_path = Path(manifest_path).resolve(strict=True)
    manifest_hash = sha256(manifest_path)
    if manifest_hash != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("golden source manifest SHA-256 mismatch: {}".format(manifest_hash))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    variant = manifest.get("golden", {})
    sources = variant.get("source_sha256", {})
    if variant.get("dll_sha256") != EXPECTED_PLUGIN_SHA256:
        raise RuntimeError("golden manifest does not identify the gauss3_boundary_reuse plugin")
    if set(sources) != set(OWNED_SOURCE_PATHS):
        raise RuntimeError("golden manifest source list does not match the six owned files")
    for relative_path in OWNED_SOURCE_PATHS:
        source_path = (ROOT / relative_path).resolve(strict=True)
        if ROOT not in source_path.parents:
            raise RuntimeError("source manifest path escapes repository: {}".format(relative_path))
        actual_hash = sha256(source_path)
        if actual_hash != sources[relative_path]:
            raise RuntimeError("golden source hash mismatch for {}: {}".format(relative_path, actual_hash))
    plugin_hash = sha256(ROOT / "plugins/2026/yddColliders.mll")
    if plugin_hash != variant["dll_sha256"]:
        raise RuntimeError("current Maya 2026 plugin does not match the golden manifest")
    return {
        "manifestPath": str(manifest_path),
        "manifestSha256": manifest_hash,
        "variant": "gauss3_boundary_reuse",
        "sourceSha256": sources,
        "dllSha256": variant["dll_sha256"],
    }


def capture(plugin, source_manifest, output):
    import maya.cmds as cmds

    sys.path.insert(0, str(ROOT / "tests"))
    from test_collide_projection import surface_snapshot
    from test_deformers import DeformerTests

    plugin_path = Path(plugin).resolve(strict=True)
    canonical_plugin = (ROOT / "plugins/2026/yddColliders.mll").resolve(strict=True)
    if plugin_path != canonical_plugin:
        raise RuntimeError("generator only accepts the canonical repository Maya 2026 plugin")
    plugin_hash = sha256(plugin_path)
    if plugin_hash != EXPECTED_PLUGIN_SHA256:
        raise RuntimeError("golden plugin SHA-256 mismatch: {}".format(plugin_hash))
    source_metadata = validate_source_manifest(source_manifest)
    output_path = Path(output).resolve()
    if output_path.exists():
        raise FileExistsError("refusing to overwrite existing golden: {}".format(output_path))
    fixtures = {
        frame: json.loads((INPUT_DIR / "frame_{}.json".format(frame)).read_text(encoding="utf-8"))
        for frame in (0, 1, 2, 10)
    }
    builder = DeformerTests()
    cases = []
    with maya_session(str(plugin_path)) as session:
        loaded_path = Path(cmds.pluginInfo("yddColliders", query=True, path=True)).resolve(strict=True)
        if loaded_path != plugin_path:
            raise RuntimeError("Maya loaded {}, expected {}".format(loaded_path, plugin_path))
        cmds.evaluationManager(mode="off")
        for frame in (0, 1, 2, 10):
            data = fixtures[frame]
            input_hash = sha256(INPUT_DIR / "frame_{}.json".format(frame))
            for closed in (False, True):
                cmds.file(new=True, force=True)
                surface = builder.fixture_surface(data["input"])
                rest = builder.fixture_surface(data["rest"])
                node = cmds.deformer(surface, type="yddSkirtCollideDeformer")[0]
                cmds.setAttr(node + ".closedU", closed)
                rest_shape = cmds.listRelatives(rest, shapes=True, noIntermediate=True, fullPath=True)[0]
                cmds.connectAttr(rest_shape + ".local", node + ".restGeometry", force=True)
                builder.configure_fixture(node, data)
                cmds.currentTime(frame)
                output = surface_snapshot(node + ".outputGeometry[0]")["cvs"]
                if len(output) != len(data["input"]["cvs"]):
                    raise RuntimeError("output CV count mismatch for frame {}, closedU={}".format(frame, closed))
                if not all(math.isfinite(component) for point in output for component in point):
                    raise RuntimeError("nonfinite output for frame {}, closedU={}".format(frame, closed))
                cases.append(
                    {
                        "frame": frame,
                        "closedU": closed,
                        "inputSha256": input_hash,
                        "cvCount": len(output),
                        "sha256LittleEndian3d": coordinate_hash(output),
                        "cvs": output,
                    }
                )
        payload = {
            "schema": "yddColliders.collideGolden.gauss3.v1",
            "mayaVersion": cmds.about(version=True),
            "pluginPath": str(loaded_path),
            "pluginSha256": plugin_hash,
            "source": source_metadata,
            "generatorPath": str(Path(__file__).resolve()),
            "generatorSha256": sha256(Path(__file__).resolve()),
            "evaluationMode": "DG (Evaluation Manager off)",
            "captureMode": "first explicit output pull after fixture configuration in a fresh scene; one contact generation with three Gauss points per nonzero V span and three Gauss-Seidel sweeps",
            "nodeType": "yddSkirtCollideDeformer",
            "nodeParameters": "serialized from each input JSON; closedU is the case value",
            "coordinateEncoding": "finite little-endian IEEE-754 doubles, packed as x,y,z per CV in Maya surface CV order",
            "inputs": {
                str(frame): {
                    "path": str(INPUT_DIR / "frame_{}.json".format(frame)),
                    "sha256": sha256(INPUT_DIR / "frame_{}.json".format(frame)),
                }
                for frame in (0, 1, 2, 10)
            },
            "cases": cases,
            "mayaSession": {"maya": session["maya"], "pluginSha256": session["sha256"]},
        }
    if not payload["mayaVersion"].startswith("2026"):
        raise RuntimeError("golden capture requires Maya 2026, got {}".format(payload["mayaVersion"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_plugin_argument(parser)
    parser.add_argument(
        "--source-manifest",
        default=str(SOURCE_MANIFEST),
        help="frozen golden gauss3_boundary_reuse manifest",
    )
    parser.add_argument("--output", default=str(OUTPUT), help="new output path; existing files are never overwritten")
    args = parser.parse_args()
    capture(args.plugin, args.source_manifest, args.output)


if __name__ == "__main__":
    main()
