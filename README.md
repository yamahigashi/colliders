# yddColliders

Maya collider and wave nodes for skirts, maintained as the yamahigashi
fork. Version 3.0.0 uses the `yddColliders` plugin and `ydd` node names
so you can load it alongside the upstream `colliders` plugin.
See [ADR-0003](docs/adr/0003-ydd-plugin-identity.md) for the identity contract.
See the [validation record](docs/validation/ydd-identity.md) for Maya 2022
through 2027 coexistence checks, scene round trips, and output comparisons.

| Node type | Local internal-use ID |
| --- | --- |
| `yddBellCollider` | `0x0007DD00` |
| `yddPlaneCollider` | `0x0007DD01` |
| `yddSkirtBellCollider` | `0x0007DD02` |
| `yddSkirtCollideDeformer` | `0x0007DD03` |
| `yddSkirtWaveDeformer` | `0x0007DD04` |

These IDs are local assignments for internal use, not Autodesk global
allocations. Before distributing outside the organization, obtain an
allocated range and plan scene migration.

Open existing fork scenes with the archived old plugin, then rebuild rigs
or migrate their data. The new names and IDs do not preserve the old node
identity in `.ma` or `.mb` scenes.

<img width="977" height="550" alt="image" src="https://github.com/user-attachments/assets/922bf5db-6b6c-40f4-8ed8-0e4748e0cf30" />

Youtube: https://www.youtube.com/watch?v=89Dbd-n8EzY

## Features
- **Bell Collider**: Multi-ring bell collider.
- **Skirt Bell Collider**: Mixture of Bell Colliders specifically for long/shirt skirts in a single node.
- **Simple UI**

## Installation
#### Compile C++ plugin for Maya version you want.
Use Visual Studio 2022 and CMake. Run `buildAll.bat 2026` for one Maya
version, or `buildAll.bat` for Maya 2022 through 2027. The build target is
`yddColliders`; the script places `yddColliders.mll` under
`plugins/<Maya version>/`.

Add that version's plugin directory to Maya's plugin path and load
`yddColliders`. Keep the archived fork `colliders.mll` outside the active
search path. The upstream repository and its `colliders` plugin remain
separate.

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
