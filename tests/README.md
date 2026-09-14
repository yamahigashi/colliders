Run these commands from the repository root with Maya's `mayapy` executable. Use an absolute plugin path to select the binary under test. You can set `YDD_COLLIDERS_PLUGIN_PATH` instead of passing `--plugin`.

```powershell
$mayapy = 'C:\Program Files\Autodesk\Maya2026\bin\mayapy.exe'
& $mayapy -B tests/run_maya_tests.py --plugin C:\build\yddColliders.mll
& $mayapy -B tests/run_maya_tests.py --plugin C:\build\yddColliders.mll --pattern test_solver.py
& $mayapy -B tests/benchmark.py --plugin C:\build\yddColliders.mll --output candidate-timing.json
```

The runner uses standard-library `unittest`, skips Python user setup, adds this checkout's `scripts` directory to `sys.path`, and closes Maya in a `finally` block. A failed test, load error, or empty discovery returns a nonzero exit status. Run each binary in a separate process.

`test_input_validation.py` drives invalid values through connected plugs and
checks that computation reports an error and recovers after restoring a valid
value. Bell and ring subdivision accept 3 through 4096, skirt type accepts 0
or 1, and axis enums accept 0 through 5. Invalid ring drawing inputs clear the
draw cache. `test_noise_inputs.py` covers extreme finite noise coordinates and
checks that nonfinite noise coordinates contribute zero noise. Run hostile
input tests only against a binary containing the input validation changes, in
a separate process with an external timeout; older binaries can hang or crash.

Current capture tools create only `ydd` node types. For comparisons between
two current binaries, capture both with this checkout and the same Maya
version:

```powershell
& $mayapy -B tests/compare_outputs.py capture --plugin C:\baseline\yddColliders.mll --output baseline.json
& $mayapy -B tests/compare_outputs.py capture --plugin C:\build\yddColliders.mll --output candidate.json
python tests/compare_outputs.py compare baseline.json candidate.json --output comparison.json
```

For a pre-identity baseline, run `capture` from the matching archived old
checkout with its old `colliders.mll`, then capture the candidate with this
checkout and `yddColliders.mll`. Compare the two JSON files with the current
`compare` command. The current capture code has no old-name fallback. Keep
historical result files and their recorded plugin names unchanged.

You can run two current-identity captures through one command:

```powershell
python tests/compare_outputs.py pair --mayapy $mayapy --baseline-plugin C:\baseline\yddColliders.mll --candidate-plugin C:\build\yddColliders.mll --output-dir comparison
```

Use native Windows Python for `pair` with Windows mayapy so both processes understand the same paths. The `compare` command needs no Maya installation. Its default tolerance is zero; use `--tolerance` to set an absolute coordinate tolerance. The report includes bitwise equality and the largest absolute component difference for each case. A comparison outside tolerance returns exit status 1, including known fixes; inspect those differences instead of treating legacy results as expected behavior.

Capture covers the 12 `REVIEW.md` skirt conditions, five smoothing/follow combinations, wave amplitude and individual signal bypasses, isolated signal contributions, complexity endpoints, DG/Serial/Parallel evaluation, and backward time changes. Additional cases cover a zero-weight highest vertex and the U1 translation setup. Each fixture stores case attributes, double-precision `(x,y,z,w)` coordinates, and a SHA-256 over little-endian doubles. Review-case metadata identifies the known U1 correction. Coordinates support numerical comparisons when a hash differs.

Benchmark each binary without another Maya job running. The benchmark uses DG evaluation, 20 warmups, five batches, and the median batch mean in milliseconds. Each evaluation changes an input and requests the output data object. Set `--iterations` to an even value of at least 2 (default 40). Timing includes Python calls and DG propagation. It excludes GUI drawing and downstream solver surface rebuilding; the 133-CV wave fixture uses a baked rebuilt surface. Mesh sizes and surface CV counts are checked at runtime. JSON includes the plugin hash, Maya/Python versions, processor description, settings, and every batch result. Compare timings only after inspecting output differences.

`test_compatibility.py` checks the 12 default and five smoothing cases against `fixtures/compatible_outputs.json`. That file records the baseline source commit, binary SHA-256, Maya version, and double encoding. Maya 2026 runs the golden-hash checks; other Maya versions skip that check and still run finite-coordinate, point-count, and periodic-seam assertions. The golden set excludes U1 probes. The benchmark resolves its input/output plugs before warmup and uses `MPlug.setDouble()` and `MPlug.asMObject()` during timing.

Build and run the C++ draw-buffer tests from the repository root. Use a devkit that matches your installed Maya version, and add Maya's `bin` directory to `PATH` so Windows can resolve the Maya DLLs:

```powershell
cmake -S tests -B build/draw-tests -A x64 -DMAYA_DEVKIT_DIR=C:/devkits/Maya2026
cmake --build build/draw-tests --config Release
$env:PATH = 'C:\Program Files\Autodesk\Maya2026\bin;' + $env:PATH
ctest --test-dir build/draw-tests -C Release --output-on-failure
```

With these tests you can check triangle connectivity, wire closure, transformed coordinates, cache invalidation, and the shared Skirt ring frames for unequal left/right leg lengths. The `solver_input` test also checks invalid subdivision, mismatched point arrays, matrix row bounds, and recovery at the C++ solver boundary. It also compares `BellCircleTable` bell points against an inline copy of the historical per-point trigonometry for subdivisions 3 through 4096, three axes, and a skewed matrix, requiring bitwise equality.

Run the viewport smoke test in a new, dedicated Maya GUI process. You will replace its scene and change its display preferences. Paste this example into the Python tab of the Script Editor; use absolute paths for the checkout, plugin, and output directory:

```python
import sys
import maya.utils

sys.path.insert(0, r"C:\src\colliders\tests")
import viewport_smoke

maya.utils.executeDeferred(
    lambda: viewport_smoke.main(
        r"C:\build\yddColliders.mll",
        r"C:\validation\viewport",
    )
)
```

Inspect `result.json` and the PNG files in the output directory. Require `passed: true`, an empty `errors` list, and all pixel checks to pass. Check the images for Bell and Skirt drawing, including the asymmetric legs. The test compares raw pixel hashes after same-frame edits and restoration, visibility changes, camera changes, and a time round trip. It also checks for colored drawing before hiding the colliders and its absence afterward. You must assess Cached Playback and playback FPS in separate tests. Close the dedicated Maya process after reviewing the results; `main()` leaves it open.


`test_deformers.py` covers the wave deformer's hidden `evaluationToWorldRotation`
(`etwr`) matrix. The tests check the identity default and attribute flags, that
an explicit identity leaves the output bitwise unchanged, that a World
direction under a rotation equals the same direction pre-rotated by the inverse
with the input strength kept for impulse, that a tilt keeps the transformed Y
component and attenuates by the projection without renormalizing, that Bell
Local directions ignore the matrix, that scale, shear, reflection,
translation, perspective, and nonfinite matrices pass the geometry through and
recover on the next valid value, and that the value and its connection survive
`.ma` and `.mb` round trips.

`test_registration.py` checks the plugin vendor/version and an independent
list of five node names, IDs, API kinds, and locator draw classifications.
It also checks the public Python module and prefixed custom node names.
The identity contract is [ADR-0003](../docs/adr/0003-ydd-plugin-identity.md).

For a real mGear component build, supply the mGear release directory and
the component repository root:

```powershell
& $mayapy -B tests/integration/check_wave_component.py --plugin C:\build\yddColliders.mll --mgear-release C:\mgear\release --components-root C:\mgear_shifter_components
```

The integration check builds three wave/post-collision configurations and
checks the host attributes, connections, deformer order, initial zero wave,
and active surface deformation. It then rebuilds with the grid rows offset
above the waist reference, entirely above the hips, and past the heels, and
requires every driver to sit on its guide position. Its JSON includes source and plugin hashes.
The component and guide must both report version 4.1.0. The check also verifies the column-major grid naming, the guide locator chain, the joint chain, that a flat pre-4.0.0 guide is rejected, that the ten leg profile parameters reach both collider nodes, and that fitting the profile from a synthetic body mesh recovers known station ratios.

The local component requires plugin major 4. Wave and post collision operate in
geometry object coordinates. Their references must use that same space; the
surface transform inherits the component root. World signal directions use
`evaluationToWorldRotation`. There is no World/Object evaluation switch.

`check_local_evaluation.py`, called by the integration runner, verifies the
surface's inherited transform, local geometry and matrix connections, skin
influence and final-surface weight queries through both Maya APIs, and zero
collider/rebuild computes on root edits. Maya's standard skin geometry connection
is retained; skin may compute once on a root edit. Six guide configurations cover translation, rotation,
scale, short skirts, and different rebuild/ring settings. Keyed forward and
backward time changes exercise a nonempty Serial/Parallel evaluation graph and
reject fallback. Maya 2026 also compares positions and basis directions against
the recorded world-space rig in `reference/skirt_world_2026.json.gz` and checks
74 contact poses at strengths 0, 1, and 2 against
`reference/skirt_secondary_world_2026.json.gz`. Five authored guide variants
also compare against `reference/skirt_guides_world_2026.json.gz`, including
ring edits after moving and rotating the guide. All three references include
source hashes and tolerances. Contact toggles must return to their original output.

To time the complete component, including final surface, controls, and joints:

```powershell
& $mayapy -B tests/integration/benchmark_skirt_component.py --plugin C:\build\yddColliders.mll --mgear-release C:\mgear\release --components-root C:\mgear_shifter_components --output rig-timing.json --mode off --mode serial --mode parallel --warmup 20 --batches 5 --iterations 60 --profile
```

Run baseline and candidate in separate processes with no concurrent Maya jobs.
Numerical capture uses the surface DAG path in world space and runs outside the
timed interval. `--profile` adds a separate DG-only `dgtimer` pass; it does not
time Serial/Parallel execution. The JSON retains every batch and per-node
compute/dirty/fetch/callback counters. Inclusive node times overlap and must not
be added together. Inputs change through `MPlug`, followed by explicit output
pulls; selecting Serial/Parallel here does not establish EM playback throughput.
Use the keyed integration check for EM correctness. These measurements exclude
GUI drawing and Cached Playback.


Add `--upstream-plugin C:\upstream\colliders.mll` to test consumer autoload
and coexistence. The CLI adds the candidate directory to the plugin search
path before initializing Maya, unloads the initial candidate in an empty
scene without forcing, and loads upstream alone. It invokes the real
component plugin guard and checks both loaded paths and registered node
sets. It then builds all three rig configurations with both plugins loaded.
The JSON records the autoload check, both binary hashes, and coexistence
checks after each rig build.

SIMD 緩和処理の実測と再現手順は [simd-relaxation.md](../docs/validation/simd-relaxation.md) を参照してください。
