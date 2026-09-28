"""Destructive smoke test for a NEW, DEDICATED Maya GUI process only.

Call main(plugin_path, output_dir) after the GUI has initialized (executeDeferred).
It replaces the scene and viewport preferences. The caller owns process shutdown.
Pixel checks establish visible drawing and refresh behavior, not visual equivalence
with a previous plugin build or Cached Playback coverage.
"""

import ctypes
import hashlib
import json
from pathlib import Path
import traceback

import maya.api.OpenMaya as om
import maya.cmds as cmds


def _matrix(position, scale=(1, 1, 1)):
    transform = om.MTransformationMatrix()
    transform.setTranslation(om.MVector(*position), om.MSpace.kTransform)
    transform.setScale(scale, om.MSpace.kTransform)
    return list(transform.asMatrix())


def main(plugin_path, output_dir):
    """Save PNG evidence and result.json; return its JSON-serializable dictionary."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    result = {"passed": False, "captures": {}, "checks": {}, "errors": []}
    try:
        if cmds.about(batch=True):
            raise RuntimeError("This test requires a dedicated Maya GUI process")
        cmds.loadPlugin(str(plugin_path), quiet=True)
        cmds.file(new=True, force=True)
        cmds.currentTime(1)
        cmds.evaluationManager(mode="parallel")
        cmds.displayRGBColor("background", 0.08, 0.08, 0.08)
        cmds.displayRGBColor("backgroundTop", 0.08, 0.08, 0.08)
        cmds.displayRGBColor("backgroundBottom", 0.08, 0.08, 0.08)
        cmds.displayPref(displayGradient=False)

        bell = cmds.createNode("yddBellCollider", name="smokeBellShape")
        skirt = cmds.createNode("yddSkirtBellCollider", name="smokeSkirtShape")
        cmds.setAttr(bell + ".bellMatrix", *_matrix((-3, -1, 0), (1.5, 3, 1.5)), type="matrix")
        cmds.setAttr(bell + ".ringMatrix[0]", *_matrix((-2.5, -1, 0), (0.7, 3, 0.7)), type="matrix")
        cmds.setAttr(bell + ".drawColor", 0.1, 0.8, 1, type="double3")
        cmds.setAttr(bell + ".drawOpacity", 0.6)
        cmds.setAttr(skirt + ".bellMatrix", *_matrix((3, 2, 0)), type="matrix")
        cmds.setAttr(skirt + ".bellAxis", 4)
        cmds.setAttr(skirt + ".skirtType", 1)
        cmds.setAttr(skirt + ".height", 0.9)
        cmds.setAttr(skirt + ".ringScale", 0.55, 1, 0.55, type="double3")
        for attribute, position in {
            "leftHipMatrix": (2.3, 1, 0),
            "leftKneeMatrix": (2.0, -1, 0.3),
            "leftHeelMatrix": (1.7, -3, 0.7),
            "rightHipMatrix": (3.7, 1, 0),
            "rightKneeMatrix": (4.1, -1.8, -0.4),
            "rightHeelMatrix": (4.5, -4.5, 0.3),
        }.items():
            cmds.setAttr(skirt + "." + attribute, *_matrix(position), type="matrix")

        camera, camera_shape = cmds.camera(name="smokeCamera", orthographic=True, orthographicWidth=15)
        cmds.setAttr(camera + ".translate", 12, 8, 20, type="double3")
        target = cmds.spaceLocator(name="smokeAim")[0]
        cmds.setAttr(target + ".translateY", -1)
        constraint = cmds.aimConstraint(target, camera, aimVector=(0, 0, -1), upVector=(0, 1, 0))[0]
        rotation = cmds.getAttr(camera + ".rotate")[0]
        cmds.delete(constraint, target)
        cmds.setAttr(camera + ".rotate", *rotation, type="double3")
        cmds.setAttr(camera_shape + ".nearClipPlane", 0.1)
        cmds.setAttr(camera_shape + ".farClipPlane", 1000)
        window = cmds.window(title="Collider viewport smoke", widthHeight=(800, 600))
        cmds.paneLayout()
        panel = cmds.modelPanel(camera=camera)
        cmds.modelEditor(
            panel,
            edit=True,
            rendererName="vp2Renderer",
            allObjects=False,
            locators=True,
            grid=False,
            hud=False,
            displayAppearance="smoothShaded",
            displayLights="default",
            selectionHiliteDisplay=False,
        )
        cmds.showWindow(window)
        cmds.setFocus(panel)
        cmds.select(clear=True)

        def capture(name):
            cmds.refresh(force=True)
            path = directory / (name + ".png")
            cmds.playblast(
                completeFilename=str(path),
                format="image",
                compression="png",
                frame=cmds.currentTime(query=True),
                viewer=False,
                showOrnaments=False,
                forceOverwrite=True,
                offScreen=False,
                percent=100,
                widthHeight=(800, 600),
                editorPanelName=panel,
            )
            image = om.MImage()
            image.readFromFile(str(path))
            width, height = image.getSize()
            pixels = ctypes.string_at(image.pixels(), width * height * 4)
            if len(pixels) != width * height * 4:
                raise RuntimeError("Unexpected MImage pixel buffer size")
            colored = sum(max(pixels[i : i + 3]) - min(pixels[i : i + 3]) > 30 for i in range(0, len(pixels), 4))
            record = {
                "path": str(path),
                "sha256": hashlib.sha256(pixels).hexdigest(),
                "colored_pixels": colored,
                "size": [width, height],
            }
            result["captures"][name] = record
            return record["sha256"]

        capture("warmup")
        baseline = capture("before")
        result["checks"]["visible_colored_drawing"] = result["captures"]["before"]["colored_pixels"] > 100
        edits = [
            (skirt + ".ringScale", (0.9, 0.8, 0.7), "double3"),
            (skirt + ".ringSubdivision", 5, None),
            (skirt + ".height", 0.2, None),
            (bell + ".ringSubdivision", 5, None),
            (bell + ".drawColor", (1, 0.1, 0.2), "double3"),
            (bell + ".drawOpacity", 0.1, None),
            (skirt + ".leftKneeMatrix", _matrix((1, -2, 1.5)), "matrix"),
        ]
        for index, (plug, value, kind) in enumerate(edits):
            original = cmds.getAttr(plug)
            if kind == "double3":
                original = original[0]
            name = "edit_{}_{}".format(index, plug.rsplit(".", 1)[1])
            if kind:
                cmds.setAttr(plug, *value, type=kind)
            else:
                cmds.setAttr(plug, value)
            result["checks"][name + "_changed"] = capture(name) != baseline
            if kind:
                cmds.setAttr(plug, *original, type=kind)
            else:
                cmds.setAttr(plug, original)
            result["checks"][name + "_restored"] = capture(name + "_restored") == baseline

        for node in (bell, skirt):
            cmds.setAttr(node + ".visibility", False)
        result["checks"]["hidden_differs"] = capture("hidden") != baseline
        result["checks"]["hidden_no_colored_drawing"] = result["captures"]["hidden"]["colored_pixels"] < 10
        for node in (bell, skirt):
            cmds.setAttr(node + ".visibility", True)
        result["checks"]["visibility_restored"] = capture("visible_restored") == baseline
        cmds.setAttr(camera_shape + ".orthographicWidth", 19)
        result["checks"]["camera_changed"] = capture("camera_changed") != baseline
        cmds.setAttr(camera_shape + ".orthographicWidth", 15)
        result["checks"]["camera_restored"] = capture("camera_restored") == baseline
        cmds.currentTime(12)
        capture("time_forward")
        cmds.currentTime(1)
        result["checks"]["time_restored"] = capture("time_restored") == baseline
        result["passed"] = all(result["checks"].values())
    except Exception:
        result["errors"].append(traceback.format_exc())
    (directory / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
