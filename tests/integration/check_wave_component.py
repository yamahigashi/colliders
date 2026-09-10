"""Build real mGear skirt rigs and verify the directional-wave integration contract."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


# member, host name, node input, default, hard range (None means unbounded)
HOST_ATTRIBUTES = (
    ("sway_vertical_att", "swayVertical", "idleAmplitudeV", 1.0, (0.0, 2.0)),
    ("sway_around_att", "swayAround", "idleAmplitudeU", 0.0, (0.0, 2.0)),
    ("sway_directional_att", "swayDirectional", "idleDirectionality", 1.0, (0.0, 1.0)),
    ("sway_dir_x_att", "swayDirX", "idleDirectionX", -1.0, None),
    ("sway_dir_z_att", "swayDirZ", "idleDirectionZ", 0.0, None),
    ("sway_spread_att", "swaySpread", "phaseSpread", 0.0, (0.0, 0.5)),
    ("send_amount_att", "sendAmount", "impulseAmount", 0.0, (0.0, 2.0)),
)
EXISTING_HOST_NAMES = (
    "waveAmplitude",
    "swayAmount",
    "swayPhase",
    "swaySpin",
    "swayNoise",
    "swayNoisePhase",
    "sendDirX",
    "sendDirZ",
    "sendPos",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def absolute_path(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise argparse.ArgumentTypeError("Expected an absolute path: %s" % value)
    try:
        return path.resolve(strict=True)
    except OSError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry_hash(points) -> str:
    return hashlib.sha256(json.dumps(points, separators=(",", ":"), allow_nan=False).encode("ascii")).hexdigest()


def host_name(component, name: str) -> str:
    if component.options["classicChannelNames"]:
        return component.getName(name)
    return component.getAttrName(name)


def build_fixture(wave: bool, post_collision: bool):
    import maya.cmds as cmds
    from mgear import shifter
    from mgear.core import transform
    from mgear.shifter import guide as rigguide
    from ymt_components.ymt_skirt_01.guide import Guide

    cmds.file(new=True, force=True)
    parent = rigguide.Rig()
    parent.initialHierarchy()
    guide = Guide()
    for row in range(5):
        v = row / 4.0
        radius = 1.8 + 0.8 * v
        for column in range(8):
            theta = 2 * math.pi * column / 8.0
            guide.tra["skirt_%d_%d_loc" % (row, column)] = transform.getTransformFromPos(
                [radius * math.cos(theta), -9 * v, radius * math.sin(theta)]
            )
    guide.draw(parent.model)
    cmds.setAttr(str(guide.root) + ".wave", wave)
    cmds.setAttr(str(guide.root) + ".postCollision", post_collision)
    for name in ("data_collector", "data_collector_embedded"):
        cmds.setAttr(str(parent.model) + "." + name, False)
    cmds.select(str(parent.model), replace=True)
    rig = shifter.Rig()
    rig.buildFromSelection()
    require(len(rig.components) == 1, "Expected exactly one built skirt component")
    component = next(iter(rig.components.values()))
    require(len(component.fk_ctls) == 40, "Expected 40 skirt FK controls")
    require(len(component.ring_ctls) == 2, "Expected 2 ring controls")
    return component


def check_build(wave: bool, post_collision: bool):
    import maya.cmds as cmds
    from helpers import coordinates

    component = build_fixture(wave, post_collision)
    host = str(component.uihost)
    expected_order = [component.ring_skin_cluster]
    result = {"wave": wave, "post_collision": post_collision, "fk_controls": 40, "ring_controls": 2}
    if wave:
        node = component.wave_deformer
        expected_order.insert(0, node)
        checked = {}
        for member, name, node_input, default, limits in HOST_ATTRIBUTES:
            plug = str(getattr(component, member))
            attribute = plug.rsplit(".", 1)[1]
            require(attribute == host_name(component, name), "Unexpected host name: " + plug)
            require(cmds.getAttr(plug) == default, "Incorrect initial value: " + plug)
            require(cmds.attributeQuery(attribute, node=host, listDefault=True) == [default], "Incorrect default: " + plug)
            require(cmds.getAttr(plug, keyable=True), "Non-keyable host attribute: " + plug)
            require(cmds.getAttr(plug, type=True) == "double", "Incorrect host type: " + plug)
            require(cmds.isConnected(plug, node + "." + node_input), "Missing direct connection: " + plug)
            for exists_flag, value_flag, expected in (
                ("minExists", "minimum", limits[0] if limits else None),
                ("maxExists", "maximum", limits[1] if limits else None),
            ):
                exists = cmds.attributeQuery(attribute, node=host, **{exists_flag: True})
                require(exists == (expected is not None), "Incorrect range constraint: " + plug)
                if exists:
                    require(
                        cmds.attributeQuery(attribute, node=host, **{value_flag: True}) == [expected], "Incorrect range: " + plug
                    )
            if limits is None:
                require(
                    cmds.attributeQuery(attribute, node=host, softRange=True) == [-1, 1],
                    "Incorrect direction soft range: " + plug,
                )
            checked[name] = {"plug": plug, "input": node_input, "default": default, "range": limits}
        require(cmds.getAttr(node + ".idleDirectionSpace") == 0, "Expected World sway direction")
        require(not cmds.getAttr(node + ".idleDirectionSpace", keyable=True), "Direction space must be non-keyable")
        require(cmds.getAttr(node + ".idleDirectionSpace", channelBox=True), "Direction space must be visible")
        initial = coordinates(node, "outputGeometry[0]", kind="surface")
        upstream = coordinates(node, "input[0].inputGeometry", kind="surface")
        require(initial == upstream, "New host defaults must leave wave geometry exactly unchanged")
        surface_before = coordinates(component.collider_surface_shape, "local", kind="surface")
        cmds.setAttr(str(component.sway_amount_att), 1)
        cmds.setAttr(str(component.sway_phase_att), 0.23)
        active = coordinates(node, "outputGeometry[0]", kind="surface")
        surface_after = coordinates(component.collider_surface_shape, "local", kind="surface")
        require(active != initial, "Raising sway must deform the wave output")
        require(surface_after != surface_before, "Raising sway must deform the final rig surface")
        result.update(
            host_attributes=checked,
            initial_wave_sha256=geometry_hash(initial),
            active_wave_sha256=geometry_hash(active),
            initial_surface_sha256=geometry_hash(surface_before),
            active_surface_sha256=geometry_hash(surface_after),
        )
    else:
        require(not cmds.ls(type="yddSkirtWaveDeformer"), "wave=False must not create a wave node")
        for name in EXISTING_HOST_NAMES + tuple(item[1] for item in HOST_ATTRIBUTES):
            require(
                not cmds.attributeQuery(host_name(component, name), node=host, exists=True),
                "wave=False created host attribute: " + name,
            )
    if post_collision:
        expected_order.insert(0, component.post_collide_deformer)
    else:
        require(not cmds.ls(type="yddSkirtCollideDeformer"), "postCollision=False created a corrective deformer")
    history = cmds.listHistory(component.collider_surface_shape, pruneDagObjects=True) or []
    actual_order = [
        node for node in history if cmds.nodeType(node) in ("skinCluster", "yddSkirtWaveDeformer", "yddSkirtCollideDeformer")
    ]
    require(actual_order == expected_order, "Incorrect downstream-first deformer order: %s" % actual_order)
    result["deformer_order"] = list(reversed(actual_order))
    return result


def check_loaded_plugins(candidate: Path, upstream: Path):
    import maya.cmds as cmds

    result = {}
    for name, expected_path, expected_types in (
        (
            "yddColliders",
            candidate,
            {"yddBellCollider", "yddPlaneCollider", "yddSkirtBellCollider", "yddSkirtCollideDeformer", "yddSkirtWaveDeformer"},
        ),
        ("colliders", upstream, {"bellCollider", "planeCollider", "skirtBellCollider"}),
    ):
        require(cmds.pluginInfo(name, query=True, loaded=True), "Expected loaded plugin: " + name)
        loaded_path = Path(cmds.pluginInfo(name, query=True, path=True)).resolve(strict=True)
        require(loaded_path.samefile(expected_path), "Unexpected loaded plugin path: " + str(loaded_path))
        node_types = set(cmds.pluginInfo(name, query=True, dependNode=True) or [])
        require(node_types == expected_types, "Unexpected registered node types for " + name)
        result[name] = {"path": str(loaded_path), "sha256": digest(loaded_path), "node_types": sorted(node_types)}
    return result


def check_component_autoload(candidate: Path, upstream: Path):
    import maya.cmds as cmds
    from ymt_components.ymt_skirt_01 import Component

    cmds.file(new=True, force=True)
    cmds.unloadPlugin("yddColliders")
    cmds.loadPlugin(str(upstream), quiet=True)
    require(not cmds.pluginInfo("yddColliders", query=True, loaded=True), "Candidate must be unloaded before the guard")
    require(cmds.pluginInfo("colliders", query=True, loaded=True), "Upstream must be loaded before the guard")
    Component.__new__(Component)._ensure_ydd_colliders_plugin()
    return {"upstream_only_before_guard": True, "loaded_plugins": check_loaded_plugins(candidate, upstream)}


def main() -> int:
    from helpers import maya_session

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plugin", required=True, type=absolute_path)
    parser.add_argument(
        "--upstream-plugin", type=absolute_path, help="Original colliders plugin for consumer autoload/coexistence checks"
    )
    parser.add_argument(
        "--mgear-release", required=True, type=absolute_path, help="mGear release directory containing scripts and platforms"
    )
    parser.add_argument(
        "--components-root",
        required=True,
        type=absolute_path,
        help="mgear_shifter_components repository directory containing python",
    )
    args = parser.parse_args()
    scripts = args.mgear_release / "scripts"
    component_python = args.components_root / "python"
    component_package = component_python / "ymt_components" / "ymt_skirt_01"
    for path in (
        args.plugin,
        scripts / "mgear" / "__init__.py",
        component_package / "__init__.py",
        component_package / "guide.py",
    ):
        require(path.is_file(), "Required input file missing: " + str(path))
    if args.upstream_plugin is not None:
        require(args.upstream_plugin.is_file(), "Required upstream plugin missing: " + str(args.upstream_plugin))
    sys.path[:0] = [str(component_python), str(scripts)]
    os.environ["MGEAR_SHIFTER_COMPONENT_PATH"] = str(component_python / "ymt_components")
    native_plugins = Path(sys.executable).resolve().parent / "plug-ins"
    os.environ["MAYA_PLUG_IN_PATH"] = os.pathsep.join(
        (str(args.plugin.parent), str(native_plugins), os.environ.get("MAYA_PLUG_IN_PATH", ""))
    )
    captured = io.StringIO()
    try:
        with redirect_stdout(captured), maya_session(str(args.plugin)) as metadata:
            import maya.cmds as cmds

            platform, extension = {
                "win32": ("windows", ".mll"),
                "linux": ("linux", ".so"),
                "darwin": ("osx", ".bundle"),
            }[sys.platform]
            for name in ("matrixNodes", "quatNodes", "cacheEvaluator"):
                if not cmds.pluginInfo(name, query=True, loaded=True):
                    cmds.loadPlugin(str(native_plugins / (name + extension)), quiet=True)
            solver_dir = args.mgear_release / "platforms" / metadata["maya"] / platform / "x64" / "plug-ins"
            os.environ["MAYA_PLUG_IN_PATH"] = str(solver_dir) + os.pathsep + os.environ["MAYA_PLUG_IN_PATH"]
            cmds.loadPlugin(str(solver_dir / ("mgear_solvers" + extension)), quiet=True)
            import mgear
            from ymt_components import ymt_skirt_01
            from ymt_components.ymt_skirt_01 import guide

            require(ymt_skirt_01.VERSION == guide.VERSION == [3, 0, 0], "Expected Component and Guide version 3.0.0")
            result = dict(metadata)
            result.update(
                mgear_release=str(args.mgear_release),
                components_root=str(args.components_root),
                component_version=ymt_skirt_01.VERSION,
                mgear_source_sha256=digest(Path(mgear.__file__)),
                component_source_sha256=digest(component_package / "__init__.py"),
                guide_source_sha256=digest(component_package / "guide.py"),
                test_source_sha256=digest(Path(__file__)),
            )
            if args.upstream_plugin is not None:
                result["consumer_autoload"] = check_component_autoload(args.plugin, args.upstream_plugin)
            result["checks"] = []
            for wave, post in ((False, True), (True, False), (True, True)):
                check = check_build(wave, post)
                if args.upstream_plugin is not None:
                    check_loaded_plugins(args.plugin, args.upstream_plugin)
                    check["consumer_coexistence"] = True
                result["checks"].append(check)
            result["passed"] = True
    except Exception:
        print(captured.getvalue()[-12000:], file=sys.stderr)
        raise
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
