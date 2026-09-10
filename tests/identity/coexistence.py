"""Run one isolated Maya load-order, serialization, and unload acceptance check."""

import argparse
from contextlib import contextmanager
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import traceback

FORK = {
    "yddBellCollider": 0x7DD00,
    "yddPlaneCollider": 0x7DD01,
    "yddSkirtBellCollider": 0x7DD02,
    "yddSkirtCollideDeformer": 0x7DD03,
    "yddSkirtWaveDeformer": 0x7DD04,
}
UPSTREAM = {"bellCollider": 1274434, "planeCollider": 1274435, "skirtBellCollider": 1274436}
ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def absolute_file(value):
    path = Path(value)
    if not path.is_absolute() or not path.is_file():
        raise argparse.ArgumentTypeError("Expected an existing absolute file: " + value)
    return path.resolve()


@contextmanager
def session():
    os.environ["MAYA_SKIP_USERSETUP_PY"] = "1"
    os.environ["MAYA_PLUG_IN_PATH"] = (
        str(Path(sys.executable).parent / "plug-ins") + os.pathsep + os.environ.get("MAYA_PLUG_IN_PATH", "")
    )
    sys.dont_write_bytecode = True
    import maya.standalone

    maya.standalone.initialize(name="python")
    try:
        yield
    finally:
        maya.standalone.uninitialize()


def identity(node):
    from maya import cmds
    import maya.api.OpenMaya as om

    obj = om.MSelectionList().add(node).getDependNode(0)
    return {"type": cmds.nodeType(node), "id": om.MFnDependencyNode(obj).typeId.id()}


def assert_identity(node, node_type):
    expected = dict(FORK, **UPSTREAM)[node_type]
    require(identity(node) == {"type": node_type, "id": expected}, "Wrong node identity: " + node)


def registration(path, expected):
    from maya import cmds

    name = path.stem
    require(cmds.pluginInfo(name, query=True, loaded=True), "Plugin is not loaded: " + name)
    actual = Path(cmds.pluginInfo(name, query=True, path=True)).resolve()
    require(actual == path, "Unexpected loaded path: " + str(actual))
    types = cmds.pluginInfo(name, query=True, dependNode=True) or []
    require(set(types) == set(expected), "Unexpected registrations: " + repr(types))
    return dict(path=str(actual), nodes=sorted(types), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def defaults_and_creation(types):
    from maya import cmds

    defaults = {}
    for node_type in types:
        if node_type.endswith("Deformer"):
            geometry = cmds.polyCylinder(constructionHistory=False)[0]
            node = cmds.deformer(geometry, type=node_type)[0]
            attr, expected = ("idleAmplitude", 0.0) if "Wave" in node_type else ("collision", 1.0)
        else:
            node = cmds.createNode(node_type)
            attr, expected = (
                ("normalAxis", 1) if "Plane" in node_type or node_type == "planeCollider" else ("bellSubdivision", 16)
            )
        assert_identity(node, node_type)
        value = cmds.getAttr(node + "." + attr)
        require(value == expected, "Default changed: " + node_type + "." + attr)
        defaults[node_type] = {attr: value}
    return defaults


def import_backend(name, directory):
    require((directory / (name + ".py")).is_file(), "Missing backend module in " + str(directory))
    sys.path.insert(0, str(directory))
    module = importlib.import_module(name)
    require(Path(module.__file__).resolve() == (directory / (name + ".py")).resolve(), "Wrong Python module: " + name)
    return module


def scene_fixture(fork_module, upstream_module):
    from maya import cmds

    records, connections, outputs = {}, [], {}
    for module, prefix, expected in (
        (fork_module, "fork_", ("yddBellCollider", "yddPlaneCollider", "yddSkirtBellCollider")),
        (upstream_module, "upstream_", ("bellCollider", "planeCollider", "skirtBellCollider")),
    ):
        bell, locator, _, _ = module.createBellCollider(prefix=prefix)
        assert_identity(bell, expected[0])
        cmds.setAttr(bell + ".bellBottomRadius", 0.625 if prefix == "fork_" else 0.875)
        records[bell] = ["bellBottomRadius"]
        outputs[bell + ".outputCurve"] = "curve"
        connections.append((locator + ".worldMatrix[0]", bell + ".bellMatrix"))
        joints = []
        for i, position in enumerate(((-1, 10, 0), (-1, 5, 0), (-1, 0, 0), (1, 10, 0), (1, 5, 0), (1, 0, 0))):
            joint = cmds.createNode("transform", name=prefix + "fixtureJoint" + str(i))
            cmds.setAttr(joint + ".translate", *position)
            joints.append(joint)
        skirt, skirt_locator, _ = module.createSkirtBellCollider(
            prefix=prefix,
            **dict(zip(("leftHipObj", "leftKneeObj", "leftHeelObj", "rightHipObj", "rightKneeObj", "rightHeelObj"), joints)),
        )
        assert_identity(skirt, expected[2])
        cmds.setAttr(skirt + ".tightness", 0.375 if prefix == "fork_" else 0.625)
        records[skirt] = ["tightness"]
        outputs[skirt + ".outputSurface"] = "surface"
        connections.append((skirt_locator + ".worldMatrix[0]", skirt + ".bellMatrix"))
        plane = cmds.createNode(expected[1], name=prefix + "plane")
        cmds.setAttr(plane + ".inputPosition", 0.25, -0.5, 0.75)
        records[plane] = ["normalAxis", "inputPosition"]
        cmds.connectAttr(locator + ".worldMatrix[0]", plane + ".planeMatrix")
        connections.append((locator + ".worldMatrix[0]", plane + ".planeMatrix"))
        outputs[plane + ".outputPosition"] = "numeric"
    for node_type in ("yddSkirtCollideDeformer", "yddSkirtWaveDeformer"):
        geometry = cmds.polyCylinder(radius=1, height=4, subdivisionsX=16, subdivisionsY=4, constructionHistory=False)[0]
        node = cmds.deformer(geometry, type=node_type)[0]
        locator = cmds.createNode("transform", name=node + "Driver")
        cmds.connectAttr(locator + ".worldMatrix[0]", node + ".bellMatrix")
        connections.append((locator + ".worldMatrix[0]", node + ".bellMatrix"))
        if "Wave" in node_type:
            attrs = dict(idleAmplitude=0.375, wavePhaseV=0.25, idleDirectionality=1, phaseSpread=0.125)
        else:
            attrs = dict(collision=0.625, falloff=0.25)
        for attr, value in attrs.items():
            cmds.setAttr(node + "." + attr, value)
        records[node] = list(attrs)
        outputs[node + ".outputGeometry[0]"] = "mesh"
    return records, connections, outputs


def snapshot(records, connections, outputs):
    from maya import cmds
    import maya.api.OpenMaya as om

    result = {"nodes": {}, "connections": [], "outputs": {}}
    for node, attrs in records.items():
        result["nodes"][node] = dict(identity(node), attributes={attr: cmds.getAttr(node + "." + attr) for attr in attrs})
    for source, destination in connections:
        require(cmds.isConnected(source, destination), "Connection lost: " + source + " -> " + destination)
        result["connections"].append([source, destination])
    for plug, kind in outputs.items():
        if kind == "numeric":
            values = cmds.getAttr(plug)[0]
        else:
            obj = om.MSelectionList().add(plug).getPlug(0).asMObject()
            if kind == "mesh":
                points = om.MFnMesh(obj).getPoints()
            elif kind == "curve":
                points = om.MFnNurbsCurve(obj).cvPositions()
            else:
                points = om.MFnNurbsSurface(obj).cvPositions()
            # Keep the owning MPointArray alive while reading borrowed points.
            values = [coordinate for p in points for coordinate in (p.x, p.y, p.z, p.w)]
        require(values and all(math.isfinite(x) for x in values), "Invalid geometry output: " + plug)
        result["outputs"][plug] = list(values)
    return result


def serialization(records, connections, outputs, directory):
    from maya import cmds

    before = snapshot(records, connections, outputs)
    results = {}
    for extension, file_type in (("ma", "mayaAscii"), ("mb", "mayaBinary")):
        path = directory / ("coexistence." + extension)
        cmds.file(rename=str(path))
        cmds.file(save=True, type=file_type, force=True)
        cmds.file(new=True, force=True)
        cmds.file(str(path), open=True, force=True, prompt=False)
        after = snapshot(records, connections, outputs)
        require(before["nodes"] == after["nodes"], extension + ": node identities or attributes changed")
        require(before["connections"] == after["connections"], extension + ": connections changed")
        maximum = 0.0
        for plug, values in before["outputs"].items():
            restored = after["outputs"][plug]
            require(len(values) == len(restored), extension + ": output size changed")
            maximum = max(maximum, max(abs(a - b) for a, b in zip(values, restored)))
        require(maximum <= 1e-10, extension + ": output changed by " + str(maximum))
        results[extension] = dict(
            nodeCount=len(records),
            connectionCount=len(connections),
            maxOutputDifference=maximum,
            exactOutputs=before["outputs"] == after["outputs"],
            tolerance=1e-10,
        )
    return results


def run(args, result):
    from maya import cmds, mel

    paths = [args.upstream_plugin, args.plugin] if args.load_order == "upstream-first" else [args.plugin, args.upstream_plugin]
    for path in paths:
        cmds.loadPlugin(str(path))
    result["plugins"] = dict(fork=registration(args.plugin, FORK), upstream=registration(args.upstream_plugin, UPSTREAM))
    result["maya"] = cmds.about(version=True)
    result["python"] = sys.version
    cmds.evaluationManager(mode="off")
    result["defaults"] = defaults_and_creation(dict(FORK, **UPSTREAM))
    result["classifications"] = {}
    for node_type in list(FORK)[:3]:
        values = cmds.getClassification(node_type)
        expected = "drawdb/geometry/" + node_type
        require(expected in values, "Wrong draw classification for " + node_type + ": " + repr(values))
        result["classifications"][node_type] = values
    cmds.file(new=True, force=True)
    fork_module = import_backend("yddColliders", ROOT / "scripts")
    upstream_module = import_backend("colliders", args.upstream_scripts.resolve())
    result["modules"] = dict(fork=fork_module.__file__, upstream=upstream_module.__file__)
    result["aeTemplates"] = {}
    for name in ("AEyddBellColliderTemplate", "AEyddSkirtBellColliderTemplate"):
        path = ROOT / "scripts" / (name + ".mel")
        require(path.is_file(), "Missing AE template: " + str(path))
        mel.eval('source "' + path.as_posix() + '";')
        location = mel.eval('whatIs "' + name + '"')
        require(location.startswith("Mel procedure found in:"), "Missing AE procedure: " + location)
        result["aeTemplates"][name] = location
    records, connections, outputs = scene_fixture(fork_module, upstream_module)
    with tempfile.TemporaryDirectory(prefix="ydd-identity-") as temporary:
        result["serialization"] = serialization(records, connections, outputs, Path(temporary))
        cmds.file(new=True, force=True)
    cmds.flushUndo()
    cmds.unloadPlugin(args.plugin.stem)
    registration(args.upstream_plugin, UPSTREAM)
    defaults_and_creation(UPSTREAM)
    require(not cmds.pluginInfo(args.plugin.stem, query=True, loaded=True), "Fork remained loaded")
    result["unloadForkLeavesUpstream"] = True
    cmds.file(new=True, force=True)
    cmds.flushUndo()
    cmds.loadPlugin(str(args.plugin))
    cmds.unloadPlugin(args.upstream_plugin.stem)
    registration(args.plugin, FORK)
    defaults_and_creation(FORK)
    require(not cmds.pluginInfo(args.upstream_plugin.stem, query=True, loaded=True), "Upstream remained loaded")
    result["unloadUpstreamLeavesFork"] = True
    cmds.file(new=True, force=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plugin", type=absolute_file, required=True)
    parser.add_argument("--upstream-plugin", type=absolute_file, required=True)
    parser.add_argument("--upstream-scripts", type=Path, required=True)
    parser.add_argument("--load-order", choices=("upstream-first", "fork-first"), required=True)
    parser.add_argument("--output", type=Path, help="Optional JSON result file (also printed to stdout)")
    args = parser.parse_args()
    require(args.plugin.stem == "yddColliders", "Expected yddColliders plugin")
    require(args.upstream_plugin.stem == "colliders", "Expected upstream colliders plugin")
    result = dict(loadOrder=args.load_order, passed=False)
    try:
        with session():
            run(args, result)
        result["passed"] = True
    except Exception:
        result["error"] = traceback.format_exc()
    payload = json.dumps(result, allow_nan=False, separators=(",", ":"))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
