Run these commands from the repository root with Maya's `mayapy` executable. Use an absolute plugin path to select the binary under test. You can set `YDD_COLLIDERS_PLUGIN_PATH` instead of passing `--plugin`.

```powershell
$mayapy = 'C:\Program Files\Autodesk\Maya2026\bin\mayapy.exe'
& $mayapy -B tests/run_maya_tests.py --plugin C:\build\yddColliders.mll
& $mayapy -B tests/run_maya_tests.py --plugin C:\build\yddColliders.mll --pattern test_solver.py
& $mayapy -B tests/benchmark.py --plugin C:\build\yddColliders.mll --output candidate-timing.json
```

The runner uses standard-library `unittest`, skips Python user setup, adds this checkout's `scripts` directory to `sys.path`, and closes Maya in a `finally` block. A failed test, load error, or empty discovery returns a nonzero exit status. Run each binary in a separate process.

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

With these tests you can check triangle connectivity, wire closure, transformed coordinates, cache invalidation, and the shared Skirt ring frames for unequal left/right leg lengths.

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
and active surface deformation. Its JSON includes source and plugin hashes.
The component and guide must both report version 3.0.0.


Add `--upstream-plugin C:\upstream\colliders.mll` to test consumer autoload
and coexistence. The CLI adds the candidate directory to the plugin search
path before initializing Maya, unloads the initial candidate in an empty
scene without forcing, and loads upstream alone. It invokes the real
component plugin guard and checks both loaded paths and registered node
sets. It then builds all three rig configurations with both plugins loaded.
The JSON records the autoload check, both binary hashes, and coexistence
checks after each rig build.
