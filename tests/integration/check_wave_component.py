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
GRID_ROWS = 5
GRID_COLS = 8
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


def grid_loc_name(column: int, row: int) -> str:
    return "skirt_%d_%d_loc" % (column, row)


def draw_guide(row_offset: float = 0.0, radius_scale: float = 1.0):
    import maya.cmds as cmds
    from mgear.core import transform
    from mgear.shifter import guide as rigguide
    from ymt_components.ymt_skirt_01.guide import Guide

    cmds.file(new=True, force=True)
    parent = rigguide.Rig()
    parent.initialHierarchy()
    guide = Guide()
    for row in range(GRID_ROWS):
        v = row / float(GRID_ROWS - 1)
        radius = (1.8 + 0.8 * v) * radius_scale
        for column in range(GRID_COLS):
            theta = 2 * math.pi * column / float(GRID_COLS)
            guide.tra[grid_loc_name(column, row)] = transform.getTransformFromPos(
                [radius * math.cos(theta), row_offset - 9 * v, radius * math.sin(theta)]
            )
    guide.draw(parent.model)
    return parent, guide


def build_fixture(
    wave: bool,
    post_collision: bool,
    row_offset: float = 0.0,
    ring_controls: int | None = 2,
    collision: float | None = None,
    radius_scale: float = 1.0,
    add_joints: bool = True,
):
    parent, guide = draw_guide(row_offset, radius_scale)
    return build_guide_fixture(parent, guide, wave, post_collision, ring_controls, collision, add_joints)


def build_guide_fixture(parent, guide, wave, post_collision, ring_controls=2, collision=None, add_joints=True):
    """Build an authored guide, sharing fixture validation across guide variants."""
    import maya.cmds as cmds
    from mgear import shifter

    cmds.setAttr(str(guide.root) + ".addJoints", add_joints)
    cmds.setAttr(str(guide.root) + ".wave", wave)
    cmds.setAttr(str(guide.root) + ".postCollision", post_collision)
    if collision is not None:
        cmds.setAttr(str(guide.root) + ".collision", collision)
    for name in ("data_collector", "data_collector_embedded"):
        cmds.setAttr(str(parent.model) + "." + name, False)
    cmds.select(str(parent.model), replace=True)
    rig = shifter.Rig()
    rig.buildFromSelection()
    require(len(rig.components) == 1, "Expected exactly one built skirt component")
    component = next(iter(rig.components.values()))
    require(len(component.fk_ctls) == GRID_ROWS * GRID_COLS, "Expected %d skirt FK controls" % (GRID_ROWS * GRID_COLS))
    if ring_controls is not None:
        require(len(component.ring_ctls) == ring_controls, "Expected %d ring controls" % ring_controls)
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


def check_rows_relative_to_waist(row_offset: float):
    # Regression: rows above the waist reference mapped to V < 0 and failed with
    # "singular driver matrix"; a hem above the hips failed the height check. The
    # component now seats the bell origin on the top row and clamps to the
    # collider surface length, so any offset must build with drivers on the guide.
    import maya.cmds as cmds

    # Auto ring stations legitimately collapse to one ring when the hem sits near
    # the knee, so only require that at least one ring built. Wider rows keep the
    # rest pose clear of the leg rings (extent 2.0), so no rest-pose collision can
    # move the drivers off the guide positions.
    component = build_fixture(False, True, row_offset=row_offset, ring_controls=None, radius_scale=1.25)
    require(len(component.ring_ctls) >= 1, "Expected at least one ring control")
    expected_waist_y = max(0.0, row_offset)
    require(
        abs(component.reference_positions["waist"][1] - expected_waist_y) < 1.0e-6,
        "Bell origin must sit on the top row when it is above the waist reference",
    )
    require(component.row_axial_projections[0] >= 0.0, "Row projections must start at or below the bell origin")
    for cell, expected in component.grid_positions.items():
        actual = cmds.xform(str(component.npos[cell]), query=True, worldSpace=True, translation=True)
        require(
            max(abs(a - e) for a, e in zip(actual, expected)) < 1.0e-3,
            "Driver %s must sit on its guide position: %s vs %s" % (cell, actual, expected),
        )
    return {
        "row_offset": row_offset,
        "fk_controls": len(component.fk_ctls),
        "ring_controls": len(component.ring_ctls),
        "waist_y": component.reference_positions["waist"][1],
        "height": component.height,
        "surface_length": component.surface_length,
    }


def short_name(node) -> str:
    return str(node).split("|")[-1]


def parent_short_name(node) -> str:
    import maya.cmds as cmds

    parents = cmds.listRelatives(str(node), parent=True, fullPath=True) or []
    require(len(parents) == 1, "Expected exactly one parent for %s" % node)
    return short_name(parents[0])


def check_grid_hierarchy():
    # Column-major naming: skirt_<col>_<row>. Guide locators chain down each
    # column from the guide root; joints chain the same way from the
    # parent-relative joint; jointRelatives index the chain col * rows + row.
    import maya.cmds as cmds

    component = build_fixture(False, True)
    guide_prefix = component.guide.fullName + "_"
    for column in range(GRID_COLS):
        for row in range(GRID_ROWS):
            locator = guide_prefix + grid_loc_name(column, row)
            require(cmds.objExists(locator), "Missing guide locator " + locator)
            expected_parent = guide_prefix + ("root" if row == 0 else grid_loc_name(column, row - 1))
            actual_parent = parent_short_name(locator)
            require(
                actual_parent == expected_parent,
                "Guide locator %s parent %s, expected %s" % (locator, actual_parent, expected_parent),
            )

    ctls_group = component.fullName + "_skirtCtls_grp"
    for column in range(GRID_COLS):
        group = component.fullName + "_skirtCol%d_grp" % column
        require(cmds.objExists(group), "Missing column group " + group)
        require(parent_short_name(group) == ctls_group, "Column group %s must sit under %s" % (group, ctls_group))
        children = [short_name(node) for node in cmds.listRelatives(group, children=True, type="transform", fullPath=True) or []]
        expected_children = [component.fullName + "_skirt_%d_%d_npo" % (column, row) for row in range(GRID_ROWS)]
        require(children == expected_children, "Column group %s children %s, expected %s" % (group, children, expected_children))
    stale_groups = cmds.ls(component.fullName + "_skirtRow*_grp")
    require(not stale_groups, "Per-row control groups must no longer exist: %s" % stale_groups)
    expected_ctls = {"skirt_%d_%d_ctl" % (column, row) for column in range(GRID_COLS) for row in range(GRID_ROWS)}
    actual_ctls = set()
    for ctl in component.fk_ctls:
        name = short_name(ctl)
        require(name.startswith(component.fullName + "_"), "Control outside component namespace: " + name)
        actual_ctls.add(name[len(component.fullName) + 1 :])
    require(actual_ctls == expected_ctls, "Unexpected control names: %s" % sorted(actual_ctls ^ expected_ctls))
    for key in list(component.relatives):
        if key.startswith("skirt_"):
            require(
                key in {grid_loc_name(c, r) for c in range(GRID_COLS) for r in range(GRID_ROWS)}, "Unexpected relative " + key
            )
            column, row = (int(part) for part in key[len("skirt_") : -len("_loc")].split("_"))
            require(component.aliasRelatives[key] == "%d_%d" % (column, row), "Unexpected alias for " + key)
            require(
                short_name(component.relatives[key]).endswith("_skirt_%d_%d_ctl" % (column, row)),
                "Relative %s points at the wrong control" % key,
            )

    joints = list(component.jointList)
    require(len(joints) == GRID_ROWS * GRID_COLS, "Expected %d grid joints, got %d" % (GRID_ROWS * GRID_COLS, len(joints)))
    row0_parents = set()
    for column in range(GRID_COLS):
        for row in range(GRID_ROWS):
            index = column * GRID_ROWS + row
            joint = joints[index]
            name = short_name(joint)
            require(name.endswith("_%d_%d_jnt" % (column, row)), "Joint %d has unexpected name %s" % (index, name))
            require(
                component.jointRelatives.get(grid_loc_name(column, row)) == index,
                "jointRelatives missing or wrong for %s" % grid_loc_name(column, row),
            )
            actual_parent = parent_short_name(joint)
            if row == 0:
                row0_parents.add(actual_parent)
            else:
                expected_parent = short_name(joints[index - 1])
                require(
                    actual_parent == expected_parent, "Joint %s parent %s, expected %s" % (name, actual_parent, expected_parent)
                )
    require(len(row0_parents) == 1, "Row 0 joints must share one parent, found %s" % sorted(row0_parents))
    require(row0_parents == {short_name(component.parent_relative_jnt)}, "Row 0 joints must hang from the parent-relative joint")
    require(component.jointRelatives.get("root") == 0, "jointRelatives['root'] must index the first grid joint")
    return {
        "controls": len(actual_ctls),
        "joints": len(joints),
        "row0_parent": sorted(row0_parents)[0],
        "joint_names": [short_name(joint) for joint in joints[: GRID_ROWS + 1]],
    }


def check_no_joints():
    # addJoints off: the grid still builds, registers no joints and no joint relatives.
    component = build_fixture(False, True, add_joints=False)
    require(len(component.jointList) == 0, "addJoints off must not create joints, got %d" % len(component.jointList))
    grid_relatives = [key for key in component.jointRelatives if key.startswith("skirt_") or key == "root"]
    require(not grid_relatives, "addJoints off must not register joint relatives: %s" % grid_relatives)
    return {"controls": len(component.fk_ctls), "joints": 0}


def build_expecting_failure(parent, needle: str, label: str) -> str:
    import maya.cmds as cmds
    from mgear import shifter

    cmds.select(str(parent.model), replace=True)
    rig = shifter.Rig()
    log = io.StringIO()
    with redirect_stdout(log):
        rig.buildFromSelection()
    text = log.getvalue()
    require(needle in text, "%s was not rejected with %r; build log:\n%s" % (label, needle, text[-3000:]))
    component = next(iter(rig.components.values()))
    require(not getattr(component, "fk_ctls", []), "%s must not build grid controls" % label)
    return next(line for line in text.splitlines() if needle in line).strip()


def check_flat_guide_rejected():
    # A guide whose grid locators all sit directly under the guide root (the
    # layout of versions before 4.0.0) must stop the build with a message that
    # names the first offending locator and its parent. A second case breaks
    # column 0 row 3 and column 7 row 1 only: a column-major scan reports
    # skirt_0_3_loc, a row-major scan would report skirt_7_1_loc.
    import maya.cmds as cmds

    parent, guide = draw_guide()
    guide_prefix = guide.fullName + "_"
    root = guide_prefix + "root"
    for column in range(GRID_COLS):
        for row in range(1, GRID_ROWS):
            cmds.parent(guide_prefix + grid_loc_name(column, row), root)
    flat = build_expecting_failure(parent, "skirt_0_1_loc must be parented under skirt_0_0_loc; found root", "Flat guide")

    parent, guide = draw_guide()
    guide_prefix = guide.fullName + "_"
    root = guide_prefix + "root"
    cmds.parent(guide_prefix + grid_loc_name(GRID_COLS - 1, 1), root)
    cmds.parent(guide_prefix + grid_loc_name(0, 3), root)
    partial = build_expecting_failure(
        parent, "skirt_0_3_loc must be parented under skirt_0_2_loc; found root", "Partially flat guide"
    )

    parent, guide = draw_guide()
    cmds.setAttr(str(guide.root) + ".rows", GRID_ROWS + 1)
    mismatch = build_expecting_failure(parent, "grid locator names do not match rows/cols", "rows mismatch")
    require(
        "skirt_0_%d_loc" % GRID_ROWS in mismatch and "must be parented" not in mismatch,
        "rows mismatch must be reported before the hierarchy check: " + mismatch,
    )
    return {"rejected": True, "flat": flat, "partial": partial, "rows_mismatch": mismatch}


PROFILE_RADIUS_NAMES = (
    "thighRadiusX",
    "thighRadiusZ",
    "kneeRadiusX",
    "kneeRadiusZ",
    "calfRadiusX",
    "calfRadiusZ",
    "ankleRadiusX",
    "ankleRadiusZ",
)
PROFILE_POSITION_NAMES = ("thighPosition", "calfPosition")


def check_leg_profile_plumbing():
    # The ten guide parameters reach both collider nodes unchanged; defaults are 1.0 / 0.5.
    import maya.cmds as cmds

    component = build_fixture(False, True)
    nodes = (component.collider_node, component.post_collide_deformer)
    for node in nodes:
        for name in PROFILE_RADIUS_NAMES:
            require(abs(cmds.getAttr(node + "." + name) - 1.0) < 1e-12, "Default %s on %s must be 1.0" % (name, node))
        for name in PROFILE_POSITION_NAMES:
            require(abs(cmds.getAttr(node + "." + name) - 0.5) < 1e-12, "Default %s on %s must be 0.5" % (name, node))

    values = dict(zip(PROFILE_RADIUS_NAMES, (1.2, 0.8, 0.7, 0.75, 1.1, 1.05, 0.6, 0.55)))
    values.update(thighPosition=0.4, calfPosition=0.6)
    parent, guide = draw_guide()
    for name, value in values.items():
        cmds.setAttr(str(guide.root) + "." + name, value)
    cmds.setAttr(str(guide.root) + ".postCollision", True)
    for name in ("data_collector", "data_collector_embedded"):
        cmds.setAttr(str(parent.model) + "." + name, False)
    cmds.select(str(parent.model), replace=True)
    from mgear import shifter

    rig = shifter.Rig()
    rig.buildFromSelection()
    component = next(iter(rig.components.values()))
    require(len(component.fk_ctls) == GRID_ROWS * GRID_COLS, "Profile build must still create the grid controls")
    for node in (component.collider_node, component.post_collide_deformer):
        for name, value in values.items():
            actual = cmds.getAttr(node + "." + name)
            require(abs(actual - value) < 1e-9, "%s.%s is %s, expected %s" % (node, name, actual, value))
    return {"nodes": 2, "values": values}


def check_leg_profile_upgrade():
    # A guide root from 4.0.0 has none of the profile parameters; opening the settings
    # adds them with defaults and the rig then builds with the default profile.
    import maya.cmds as cmds
    from ymt_components.ymt_skirt_01.guide import componentSettings

    parent, guide = draw_guide()
    root = str(guide.root)
    names = PROFILE_RADIUS_NAMES + PROFILE_POSITION_NAMES + ("profileMesh",)
    for name in names:
        cmds.deleteAttr(root + "." + name)
    require(not cmds.attributeQuery("thighRadiusX", node=root, exists=True), "Attribute removal failed")

    class Settings:
        _ensure_leg_profile_parameters = componentSettings._ensure_leg_profile_parameters

    settings = Settings()
    settings.root = guide.root
    added = settings._ensure_leg_profile_parameters()
    require(sorted(added) == sorted(names), "Unexpected added parameters: %s" % added)
    for name in PROFILE_RADIUS_NAMES:
        require(abs(cmds.getAttr(root + "." + name) - 1.0) < 1e-12, "Restored %s must default to 1.0" % name)
    for name in PROFILE_POSITION_NAMES:
        require(abs(cmds.getAttr(root + "." + name) - 0.5) < 1e-12, "Restored %s must default to 0.5" % name)
    require(cmds.getAttr(root + ".profileMesh") in ("", None), "Restored profileMesh must be empty")
    require(settings._ensure_leg_profile_parameters() == [], "Second call must add nothing")
    for name in ("data_collector", "data_collector_embedded"):
        cmds.setAttr(str(parent.model) + "." + name, False)
    cmds.select(str(parent.model), replace=True)
    from mgear import shifter

    rig = shifter.Rig()
    rig.buildFromSelection()
    component = next(iter(rig.components.values()))
    require(len(component.fk_ctls) == GRID_ROWS * GRID_COLS, "Upgraded guide must build the grid")
    require(abs(cmds.getAttr(component.collider_node + ".kneeRadiusX") - 1.0) < 1e-12, "Upgraded guide must send defaults")
    return {"added": sorted(added), "controls": len(component.fk_ctls)}


def check_leg_profile_validation():
    # The guide parameters carry the node limits, so out-of-range values are refused
    # before a build can see them; the rig-side validation is the second line.
    import maya.cmds as cmds

    results = {}
    parent, guide = draw_guide()
    for name, value in (("kneeRadiusX", 0.0), ("calfPosition", 1.5), ("thighPosition", -0.1)):
        try:
            cmds.setAttr(str(guide.root) + "." + name, value)
        except RuntimeError as exc:
            results[name] = str(exc).strip()
        else:
            require(False, "Guide must refuse %s=%s" % (name, value))
    return results


def _leg_loft(points, radii, sides=64):
    # Sections perpendicular to each segment with a circular profile; returns (vertices, counts, connects).
    import math as _math

    vertices = []
    for (start, end), station_radii in zip(zip(points[:-1], points[1:]), radii):
        axis = [e - s for s, e in zip(start, end)]
        length = _math.sqrt(sum(a * a for a in axis))
        axis = [a / length for a in axis]
        # basis perpendicular to the segment axis
        up = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 0.0, 1.0)
        u = [axis[1] * up[2] - axis[2] * up[1], axis[2] * up[0] - axis[0] * up[2], axis[0] * up[1] - axis[1] * up[0]]
        norm = _math.sqrt(sum(c * c for c in u))
        u = [c / norm for c in u]
        v = [axis[1] * u[2] - axis[2] * u[1], axis[2] * u[0] - axis[0] * u[2], axis[0] * u[1] - axis[1] * u[0]]
        for t, radius in station_radii:
            center = [s + (e - s) * t for s, e in zip(start, end)]
            for k in range(sides):
                angle = 2.0 * _math.pi * k / sides
                vertices.append(
                    tuple(center[i] + u[i] * radius * _math.cos(angle) + v[i] * radius * _math.sin(angle) for i in range(3))
                )
    rings = len(vertices) // sides
    counts, connects = [], []
    for ring in range(rings - 1):
        for k in range(sides):
            a, b = ring * sides + k, ring * sides + (k + 1) % sides
            counts.append(4)
            connects.extend([a, b, b + sides, a + sides])
    return vertices, counts, connects


def _profile_radius(t, stations):
    for (t0, r0), (t1, r1) in zip(stations[:-1], stations[1:]):
        if t0 <= t <= t1:
            return r0 + (r1 - r0) * (t - t0) / (t1 - t0) if t1 > t0 else r1
    return stations[-1][1]


def check_leg_profile_fitting():
    # Synthetic legs: hip 0.6, thigh peak 0.75 at 40 percent, knee 0.55, calf peak 0.7 at 50
    # percent, ankle 0.4. Fitting must recover the ratios within 5 percent and the peak
    # positions within 0.1, and reject an empty or non-mesh profileMesh.
    import maya.cmds as cmds
    import maya.api.OpenMaya as om
    from ymt_components.ymt_skirt_01.guide import componentSettings

    parent, guide = draw_guide()
    prefix = guide.fullName + "_"
    thigh_stations = ((0.0, 0.6), (0.4, 0.75), (1.0, 0.55))
    calf_stations = ((0.0, 0.55), (0.5, 0.7), (1.0, 0.4))
    vertices, counts, connects = [], [], []
    for side in ("L", "R"):
        joints = [
            tuple(cmds.xform(prefix + name + "_" + side, query=True, worldSpace=True, translation=True))
            for name in ("hip", "knee", "heel")
        ]
        section_radii = []
        for stations in (thigh_stations, calf_stations):
            section_radii.append([(i / 32.0, _profile_radius(i / 32.0, stations)) for i in range(33)])
        leg_vertices, leg_counts, leg_connects = _leg_loft(joints, section_radii)
        offset = len(vertices)
        vertices.extend(leg_vertices)
        counts.extend(leg_counts)
        connects.extend(index + offset for index in leg_connects)
    fn = om.MFnMesh()
    mesh_obj = fn.create(om.MPointArray([om.MPoint(*v) for v in vertices]), counts, connects)
    mesh = om.MFnDagNode(mesh_obj).fullPathName()
    mesh = cmds.rename(mesh, "profileBody")
    # Second body: the same legs split into two separate meshes under a group, as
    # production characters are; the group name must fit identically.
    half = len(counts) // 2
    parts = []
    for index, (count_slice, connect_slice) in enumerate(
        ((counts[:half], connects[: sum(counts[:half])]), (counts[half:], connects[sum(counts[:half]) :]))
    ):
        part_fn = om.MFnMesh()
        part_obj = part_fn.create(om.MPointArray([om.MPoint(*v) for v in vertices]), count_slice, connect_slice)
        parts.append(cmds.rename(om.MFnDagNode(part_obj).fullPathName(), "profilePart%d" % index))
    parts_group = cmds.group(parts, name="profileParts")

    class Fitter:
        fit_profile_from_mesh = componentSettings.fit_profile_from_mesh
        update_ring_preview = componentSettings.update_ring_preview
        _create_ring_preview_circle = componentSettings._create_ring_preview_circle
        _delete_ring_preview = componentSettings._delete_ring_preview
        _guide_position = componentSettings._guide_position
        _root_front = componentSettings._root_front
        _normalized = componentSettings._normalized

    fitter = Fitter()
    fitter.root = guide.root
    result = {}
    cmds.setAttr(str(guide.root) + ".thighRadiusX", 2.5)
    for bad in ("", "nonexistent_mesh_xyz", prefix + "root"):
        cmds.setAttr(str(guide.root) + ".profileMesh", bad, type="string")
        try:
            fitter.fit_profile_from_mesh()
        except RuntimeError as exc:
            require("profileMesh" in str(exc), "Unexpected fitting error for %r: %s" % (bad, exc))
        else:
            require(False, "Fitting must reject profileMesh=%r" % bad)
        require(
            abs(cmds.getAttr(str(guide.root) + ".thighRadiusX") - 2.5) < 1e-12, "Rejected fit must leave parameters unchanged"
        )
    cmds.setAttr(str(guide.root) + ".thighRadiusX", 1.0)
    result["rejected"] = ["", "nonexistent", "non-mesh"]

    cmds.setAttr(str(guide.root) + ".profileMesh", mesh, type="string")
    fitter.fit_profile_from_mesh()
    fitted = {
        name: cmds.getAttr(str(guide.root) + "." + name)
        for name in PROFILE_RADIUS_NAMES + PROFILE_POSITION_NAMES + ("ringScaleX", "ringScaleZ")
    }
    expected = {
        "thighRadiusX": 0.75 / 0.6,
        "thighRadiusZ": 0.75 / 0.6,
        "kneeRadiusX": 0.55 / 0.6,
        "kneeRadiusZ": 0.55 / 0.6,
        "calfRadiusX": 0.7 / 0.6,
        "calfRadiusZ": 0.7 / 0.6,
        "ankleRadiusX": 0.4 / 0.6,
        "ankleRadiusZ": 0.4 / 0.6,
        "ringScaleX": 0.6 / 1.0,
        "ringScaleZ": 0.6 / 1.0,
    }
    for name, value in expected.items():
        require(abs(fitted[name] - value) <= 0.05 * value, "Fitted %s = %s, expected %s" % (name, fitted[name], value))
    require(abs(fitted["thighPosition"] - 0.4) <= 0.1, "Fitted thighPosition %s, expected about 0.4" % fitted["thighPosition"])
    require(abs(fitted["calfPosition"] - 0.5) <= 0.1, "Fitted calfPosition %s, expected about 0.5" % fitted["calfPosition"])
    previews = cmds.ls(prefix + "ringPreview*_Crv")
    require(len(previews) == 10, "Expected 10 preview ellipses after fitting, found %d" % len(previews))
    result["fitted"] = fitted

    cmds.setAttr(str(guide.root) + ".profileMesh", parts_group, type="string")
    fitter.fit_profile_from_mesh()
    for name in expected:
        again = cmds.getAttr(str(guide.root) + "." + name)
        require(
            abs(again - fitted[name]) < 1e-6, "Group fit %s = %s differs from single-mesh fit %s" % (name, again, fitted[name])
        )
    empty_group = cmds.group(empty=True, name="profileEmptyGroup")
    cmds.setAttr(str(guide.root) + ".profileMesh", empty_group, type="string")
    try:
        fitter.fit_profile_from_mesh()
    except RuntimeError as exc:
        require("profileMesh" in str(exc), "Unexpected error for an empty group: %s" % exc)
    else:
        require(False, "Fitting must reject a group without meshes")
    result["group_fit"] = parts_group
    return result


def check_rebuild_grid():
    # The settings UI rebuild path (componentSettings.rebuild_grid_locators)
    # must create the same col-major chain, and rebuilding over an existing
    # chain must delete parents and children without raising. The Qt dialog
    # is not constructed; the rebuild methods run on a stand-in object.
    import maya.cmds as cmds
    from ymt_components.ymt_skirt_01 import guide as guide_module
    from ymt_components.ymt_skirt_01.guide import componentSettings

    pm = guide_module.pm

    class SpinBox:
        def __init__(self, value: int) -> None:
            self._value = value

        def value(self) -> int:
            return self._value

    class Tab:
        def __init__(self, rows: int, cols: int) -> None:
            self.rows_spinBox = SpinBox(rows)
            self.cols_spinBox = SpinBox(cols)

    class Rebuilder:
        rebuild_grid_locators = componentSettings.rebuild_grid_locators
        _guide_position = componentSettings._guide_position
        _root_front = componentSettings._root_front
        _normalized = componentSettings._normalized
        _delete_existing_grid = componentSettings._delete_existing_grid
        _create_grid_locator = componentSettings._create_grid_locator
        _create_row_display_curves = componentSettings._create_row_display_curves

    def grid_state(guide_prefix: str, rows: int, cols: int):
        top = cmds.listRelatives(guide_prefix + "root", parent=True, fullPath=True)[0]
        names = [
            short_name(node) for node in cmds.listRelatives(top, allDescendents=True, type="transform", fullPath=True) or []
        ]
        grid = sorted(name for name in names if name.startswith(guide_prefix + "skirt_") and name.endswith("_loc"))
        curves = sorted(name for name in names if name.startswith(guide_prefix + "skirtRow") and name.endswith("Crv"))
        require(len(grid) == rows * cols, "Expected %d grid locators after rebuild, found %d" % (rows * cols, len(grid)))
        require(len(curves) == rows, "Expected %d row curves after rebuild, found %d" % (rows, len(curves)))
        positions = {}
        for column in range(cols):
            for row in range(rows):
                locator = guide_prefix + grid_loc_name(column, row)
                require(locator in grid, "Rebuild did not create " + locator)
                expected_parent = guide_prefix + ("root" if row == 0 else grid_loc_name(column, row - 1))
                actual_parent = parent_short_name(locator)
                require(
                    actual_parent == expected_parent,
                    "Rebuilt %s parent %s, expected %s" % (locator, actual_parent, expected_parent),
                )
                positions[(column, row)] = cmds.xform(locator, query=True, worldSpace=True, translation=True)
        return positions

    parent, guide = draw_guide()
    guide_prefix = guide.fullName + "_"
    rebuilder = Rebuilder()
    rebuilder.root = pm.PyNode(guide_prefix + "root")
    rebuilder.settingsTab = Tab(GRID_ROWS, GRID_COLS)
    rebuilder.rebuild_grid_locators()
    first = grid_state(guide_prefix, GRID_ROWS, GRID_COLS)
    # Rows share a center and radius, so world positions must be preserved through
    # parenting: every locator in a row keeps the same distance to the row centroid.
    for row in range(GRID_ROWS):
        ring = [first[(column, row)] for column in range(GRID_COLS)]
        center = [sum(p[axis] for p in ring) / GRID_COLS for axis in range(3)]
        radii = [math.sqrt(sum((p[axis] - center[axis]) ** 2 for axis in range(3))) for p in ring]
        require(max(radii) - min(radii) < 1.0e-6 and min(radii) > 0.1, "Rebuilt row %d is not a ring: %s" % (row, radii))
    rebuilder.settingsTab = Tab(GRID_ROWS - 1, GRID_COLS + 2)
    rebuilder.rebuild_grid_locators()
    second = grid_state(guide_prefix, GRID_ROWS - 1, GRID_COLS + 2)
    rebuilder.settingsTab = Tab(GRID_ROWS, GRID_COLS)
    rebuilder.rebuild_grid_locators()
    third = grid_state(guide_prefix, GRID_ROWS, GRID_COLS)
    require(
        all(max(abs(a - b) for a, b in zip(first[cell], third[cell])) < 1.0e-6 for cell in first),
        "Rebuilding twice must reproduce the same world positions",
    )
    return {"first": len(first), "second": len(second), "third": len(third)}


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

            require(ymt_skirt_01.VERSION == guide.VERSION == [4, 1, 0], "Expected Component and Guide version 4.1.0")
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
            result["rows_relative_to_waist"] = [check_rows_relative_to_waist(offset) for offset in (1.0, 20.0, -3.0)]
            result["grid_hierarchy"] = check_grid_hierarchy()
            result["no_joints"] = check_no_joints()
            result["rebuild_grid"] = check_rebuild_grid()
            result["leg_profile_plumbing"] = check_leg_profile_plumbing()
            result["leg_profile_validation"] = check_leg_profile_validation()
            result["leg_profile_upgrade"] = check_leg_profile_upgrade()
            result["leg_profile_fitting"] = check_leg_profile_fitting()
            result["flat_guide_rejected"] = check_flat_guide_rejected()
            from check_local_evaluation import check_local_evaluation

            result["local_evaluation"] = check_local_evaluation()
            result["passed"] = True
    except Exception:
        print(captured.getvalue()[-12000:], file=sys.stderr)
        raise
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
