"""Shared standalone lifecycle and deterministic Maya scene builders."""

import argparse
from contextlib import contextmanager
import hashlib
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def add_plugin_argument(parser):
    parser.add_argument(
        "--plugin", default=os.environ.get("COLLIDERS_PLUGIN_PATH"), help="Absolute plugin path (or COLLIDERS_PLUGIN_PATH)"
    )


@contextmanager
def maya_session(plugin):
    if not plugin:
        raise ValueError("Specify --plugin or COLLIDERS_PLUGIN_PATH")
    plugin = Path(plugin).resolve(strict=True)
    os.environ["MAYA_SKIP_USERSETUP_PY"] = "1"
    sys.path.insert(0, str(ROOT / "scripts"))
    import maya.standalone

    maya.standalone.initialize(name="python")
    try:
        import maya.cmds as cmds

        cmds.loadPlugin(str(plugin))
        cmds.evaluationManager(mode="off")
        yield {
            "plugin": str(plugin),
            "sha256": hashlib.sha256(plugin.read_bytes()).hexdigest(),
            "maya": cmds.about(version=True),
            "python": sys.version,
        }
    finally:
        maya.standalone.uninitialize()


def transform(position=(0, 0, 0), rotation=(0, 0, 0), scale=(1, 1, 1)):
    from maya import cmds

    node = cmds.createNode("transform")
    for attr, value in [("translate", position), ("rotate", rotation), ("scale", scale)]:
        cmds.setAttr(node + "." + attr, *value)
    return node


def skirt(bent=False):
    import colliders

    positions = [
        (-1, 10, 0),
        (-1, 6, 3) if bent else (-1, 5, 0),
        (-1, 2, 0) if bent else (-1, 0, 0),
        (1, 10, 0),
        (1, 5, 0),
        (1, 0, 0),
    ]
    joints = [transform(p) for p in positions]
    return colliders.createSkirtBellCollider(
        leftHipObj=joints[0],
        leftKneeObj=joints[1],
        leftHeelObj=joints[2],
        rightHipObj=joints[3],
        rightKneeObj=joints[4],
        rightHeelObj=joints[5],
    )


def output_object(node, attribute):
    import maya.api.OpenMaya as om

    return om.MSelectionList().add(node + "." + attribute).getPlug(0).asMObject()


def coordinates(node, attribute, kind="mesh"):
    import maya.api.OpenMaya as om

    obj = output_object(node, attribute)
    points = om.MFnNurbsSurface(obj).cvPositions() if kind == "surface" else om.MFnMesh(obj).getPoints()
    return [[p.x, p.y, p.z, p.w] for p in points]


def cylinder(size):
    from maya import cmds

    sx, sy = {80: (16, 4), 1088: (64, 16), 16512: (128, 128)}[size]
    obj = cmds.polyCylinder(
        radius=1, height=4, subdivisionsX=sx, subdivisionsY=sy, subdivisionsCaps=0, constructionHistory=False
    )[0]
    cmds.move(0, 2, 0, obj + ".vtx[*]", relative=True, objectSpace=True)
    actual = cmds.polyEvaluate(obj, vertex=True)
    if actual != size:
        raise AssertionError(f"Expected {size} vertices, got {actual}")
    return obj


def wave_settings(node, condition="active"):
    from maya import cmds

    attrs = dict(amplitude=1, envelope=1, idleAmplitude=0.3, impulseX=-0.5, impulseZ=0, noiseAmplitude=0.2, impulsePosition=0.5)
    if condition == "amplitude0":
        attrs["amplitude"] = 0
    elif condition == "envelope0":
        attrs["envelope"] = 0
    elif condition == "signals0":
        attrs.update(idleAmplitude=0, impulseX=0, impulseZ=0, noiseAmplitude=0)
    elif condition == "impulse_only":
        attrs.update(idleAmplitude=0, noiseAmplitude=0)
    elif condition == "periodic_only":
        attrs.update(impulseX=0, noiseAmplitude=0)
    elif condition == "noise_only":
        attrs.update(idleAmplitude=0, impulseX=0)
    elif condition != "active":
        raise ValueError(condition)
    for attr, value in attrs.items():
        cmds.setAttr(node + "." + attr, value)
    return attrs


def parser(description):
    result = argparse.ArgumentParser(description=description)
    add_plugin_argument(result)
    return result
