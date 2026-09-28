# yddColliders

Maya collider and skirt deformer nodes, maintained by yamahigashi as a
fork of [azagoruyko/colliders](https://github.com/azagoruyko/colliders)
by Alexander Zagoruyko. The upstream project provides the Bell Collider,
Plane Collider, Skirt Bell Collider, and the UI; this fork keeps them and
adds skirt post deformers, a surface fit node, and Attribute Editor
templates. Both projects use the Apache 2.0 license (see
[LICENSE](LICENSE)).

## Relation to upstream

Since version 3.0.0 this fork builds as the `yddColliders` plugin and
registers its nodes with `ydd` names and its own IDs, so it loads alongside
the upstream `colliders` plugin without conflicts. The upstream repository,
its `colliders.mll`, and its node names are not changed by this fork.

| Node type | Local internal-use ID |
| --- | --- |
| `yddBellCollider` | `0x0007DD00` |
| `yddPlaneCollider` | `0x0007DD01` |
| `yddSkirtBellCollider` | `0x0007DD02` |
| `yddSkirtCollideDeformer` | `0x0007DD03` |
| `yddSkirtWaveDeformer` | `0x0007DD04` |
| `yddSkirtSurfaceFit` | `0x0007DD05` |

These IDs are local assignments for internal use, not Autodesk global
allocations. Before distributing outside the organization, obtain an
allocated range and plan scene migration.

Screenshot and video from the upstream project:

<img width="977" height="550" alt="image" src="https://github.com/user-attachments/assets/922bf5db-6b6c-40f4-8ed8-0e4748e0cf30" />

Youtube: https://www.youtube.com/watch?v=89Dbd-n8EzY

## Features
From upstream:
- **Bell Collider**: Multi-ring bell collider.
- **Skirt Bell Collider**: Mixture of Bell Colliders specifically for long/shirt skirts in a single node.
- **Simple UI**

Added in this fork:
- **Skirt Collide Deformer**: post pass that pushes the skirt surface out of the collider surface and smooths the result.
- **Skirt Wave Deformer**: directional waves on the skirt surface.
- **Skirt Surface Fit**: rebuilds a NURBS surface with a chosen V span count while keeping selected V parameters fixed.
- **Attribute Editor templates** for the collider nodes.

## Installation
#### Compile C++ plugin for Maya version you want.
Use Visual Studio 2022 and CMake. Run `buildAll.bat 2026` for one Maya
version, or `buildAll.bat` for Maya 2022 through 2027. The build target is
`yddColliders`; the script places `yddColliders.mll` under
`plugins/<Maya version>/`.

Add that version's plugin directory to Maya's plugin path and load
`yddColliders`.

#### Using Python
Use the script `scripts/yddColliders.py` to run the UI.

```python
import sys
sys.path.insert(0, r"<path_to_colliders>/scripts")
import maya.cmds as cmds
cmds.loadPlugin("yddColliders")
import yddColliders
yddColliders.show()
```

Keep the `AEyddBellColliderTemplate.mel` and
`AEyddSkirtBellColliderTemplate.mel` files on Maya's script path. See
[tests/README.md](tests/README.md) for registration, geometry, and rig checks.
