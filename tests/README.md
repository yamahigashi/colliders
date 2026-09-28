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

Capture covers the 12 default skirt conditions, five smoothing/follow combinations, wave amplitude and individual signal bypasses, isolated signal contributions, complexity endpoints, DG/Serial/Parallel evaluation, and backward time changes. Additional cases cover a zero-weight highest vertex and the U1 translation setup. Each fixture stores case attributes, double-precision `(x,y,z,w)` coordinates, and a SHA-256 over little-endian doubles. Review-case metadata identifies the known U1 correction. Coordinates support numerical comparisons when a hash differs.

Benchmark each binary without another Maya job running. The benchmark uses DG evaluation, 20 warmups, five batches, and the median batch mean in milliseconds. Each evaluation changes an input and requests the output data object. Set `--iterations` to an even value of at least 2 (default 40). Timing includes Python calls and DG propagation. It excludes GUI drawing and downstream solver surface rebuilding; the 133-CV wave fixture uses a baked rebuilt surface. Mesh sizes and surface CV counts are checked at runtime. JSON includes the plugin hash, Maya/Python versions, processor description, settings, and every batch result. Compare timings only after inspecting output differences.

`test_compatibility.py` checks the 12 default and five smoothing cases against `fixtures/compatible_outputs.json`. That file records the baseline source commit, binary SHA-256, Maya version, and double encoding. Maya 2026 runs the golden-hash checks; other Maya versions skip that check and still run finite-coordinate, point-count, and periodic-seam assertions. The golden set excludes U1 probes. The benchmark resolves its input/output plugs before warmup and uses `MPlug.setDouble()` and `MPlug.asMObject()` during timing.

The second-pass golden capture is `fixtures/collide_golden_maya2026.json`.
It stores yddSkirtCollideDeformer outputs for the four long-skirt poses in
`fixtures/long_skirt_poses/` (frames 0, 1, 2, and 10) with both `closedU`
values, captured from a fresh scene with one contact generation followed by
three Gauss-Seidel sweeps and three Gauss points per span. Its metadata
records the plugin, the source manifest, the inputs, the evaluation mode, and
the CV packing. Regenerate it only from a plugin built from the sources listed
in `fixtures/collide_golden_source_manifest.json`:

```powershell
& $mayapy -B tests/generate_collide_golden.py --plugin (Resolve-Path plugins/2026/yddColliders.mll) --source-manifest tests/fixtures/collide_golden_source_manifest.json --output C:\validation\collide_golden.json
```

The generator verifies the manifest and its source hashes before loading the
plugin, refuses to overwrite an existing capture, and checks CV counts and
finite coordinates before writing. `fixtures/compatible_outputs.json` is the
separate first-pass skirt-surface compatibility golden. The Python collision
oracle in `collide_reference.py` uses the same single-contact, three-sweep
schedule when no schedule is supplied.

The first-pass fixtures are `fixtures/circular_vibration_rig.json` (a 19 by 4
bell rig with per-frame leg matrices) with `circular_vibration_baseline.json`,
`knee_bend_straight_baseline.json`, and `row_suspension_rest_baseline.json`.
They are read by test_thigh_ring_end.py, test_row_relax_direction.py,
test_knee_bend_fade.py, and test_row_suspension.py.

Build and run the C++ draw-buffer tests from the repository root. Use a devkit that matches your installed Maya version, and add Maya's `bin` directory to `PATH` so Windows can resolve the Maya DLLs:

```powershell
cmake -S tests -B build/draw-tests -A x64 -DMAYA_DEVKIT_DIR=C:/devkits/Maya2026
cmake --build build/draw-tests --config Release
$env:PATH = 'C:\Program Files\Autodesk\Maya2026\bin;' + $env:PATH
ctest --test-dir build/draw-tests -C Release --output-on-failure
```

With these tests you can check triangle connectivity, wire closure, transformed coordinates, cache invalidation, and the shared Skirt ring frames for unequal left/right leg lengths. The `solver_input` test also checks invalid subdivision, mismatched point arrays, matrix row bounds, and recovery at the C++ solver boundary. It also compares `BellCircleTable` bell points against an inline copy of the historical per-point trigonometry for subdivisions 3 through 4096, three axes, and a skewed matrix, requiring bitwise equality.

Run the nurbs_refit_kernel CTest target to check basis functions, uniform
knots and Greville abscissae, linear precision through de Boor evaluation,
convex bounds, repeated knots, and the clamped right endpoint. The executable
links no Maya libraries. You can also compile tests/nurbs_refit_kernel_test.cpp
with a C++11 compiler and sources as the include directory; a failed assertion
returns a nonzero exit status.

With test_surface_fit.py, you can check yddSkirtSurfaceFit on long and short
skirts at spansV values 1, 4, and 9: topology, knots, convex bounds, endpoints,
periodicity, invalid inputs and recovery, and degenerate height. The test
prints the maximum shape difference and checks it against 10 percent of the
long skirt's mean leg length. You can also check CV hashes across DG, Serial,
and Parallel time changes and restore spansV and connections from a .ma file.

With test_collide_projection.py, you can check the second pass,
yddSkirtCollideDeformer, on a radius-2.6 mesh with 32 circumference segments
and 16 height segments. Supply a separate rest mesh through restGeometry and
capture restBellMatrix plus the six rest hip, knee, and heel matrices. You can
also connect NURBS surface data to restGeometry. Keep rest points, current
points, and every matrix in the geometry's object coordinates.

For each point, use its outward direction in the rest shape to choose the
contact side, then rotate that direction with each leg. Build each constraint
from a leg's exterior support plane and its activation weight, and include an
upward displacement constraint when contact requires it.
Short skirts use both thighs; Long skirts add both lower legs. Resolve both
legs together to find the minimum displacement. Use falloff (0 through 1,
default 0.2) to set an activation band outside each cylinder, measured as a
fraction of its section radius. Outside that band, you get no correction from
that cylinder. Through the band, you get a smooth increase in constraint
weight, reaching full weight at the surface and throughout the interior.
Remove endFade from second-pass wiring; endpoint handling now uses the
support normal. Set envelope and vertex weights to 1 to apply the
full solution; smaller strengths scale the correction. With full activation
and feasible constraints, you get the exterior support condition for that
cylinder. Check vertices and surface CVs; their interpolated faces or surface
patches can still intersect a leg.

The projection tests cover unchanged rest points with positive clearance,
interior points projected out even at rest, horizontal legs, forward/reverse
sweeps, bilateral symmetry, exchanged leg inputs, and DG/Serial/Parallel time
changes. You can check every support inequality and cross-section interior
with the independent Python geometry oracle. Read the printed infeasible
vertex count, maximum displacement per degree over 0 through 90 degrees,
and maximum difference in displacement across mesh edges with the test results.
For the local sweep from 60 through 66 degrees in 0.1-degree steps, require a
maximum step displacement of 0.2. Check reverse-sweep hashes in a separate
subtest; the per-degree and mesh-edge values are reports without thresholds.
Missing or mismatched rest geometry must preserve the input and produce one
warning per node; reconnecting valid rest data must restore the correction.
The suite also covers NURBS rest data and half-strength displacement.

Two expected-failure regressions document limits of partial activation.
In the horizontal-leg fixture, you can satisfy the weighted constraints
without reaching every full support plane. With a current point across the
leg from its rest contact side, you can also end inside the cylinder after a
partial correction: at radius 1, falloff 0.2, rest point (2, 0.5, 0), and
current point (-1.1, 0.5, 0), the result is about (-0.05, 0.5, 0).
The sweep test has no expected-failure marker.

With the collision cases in test_deformers.py, you can check endpoint planes,
station radii and positions through their effect on support, coincident
stations, degenerate segments, singular matrices, strength bypasses, and
object-coordinate invariance. For falloff, check unchanged points outside the
band, increasing displacement as you widen the band, continuity at both band
boundaries, and full correction inside the cylinder even at falloff zero.
The capture tool includes static and horizontal leg poses with rest inputs.
For the correction benchmark, toggle nodeState between 0 and 1 in active
and zero-envelope cases.

With tests/test_waist_ring.py, you can check upper-ring contact on
yddSkirtBellCollider through nodeState and smoothness. Coverage includes
contact against both legs, unchanged CVs behind the hip plane when smoothness
is zero, a higher noncontact waist, forward/reverse leg sweeps,
DG/Serial/Parallel time changes, bilateral symmetry, root tilt, static poses
at subdivisions 3, 16, and 64, nodeState bypass, and upper-row follow isolation.
The nodeState test checks every raw ring against analytic circle coordinates,
retains the knots and periodic CVs, and restores the normal output after a toggle.
The geometric assertions use nodeState = 1 as the baseline and an independent
Python reconstruction of the leg cylinders. Only points on the knee side of
the hip plane need to clear the cylinder radius. These tests check CVs;
cubic interpolation and the surface between rows can pass inside the cylinders.

Run the native relax_solver_reference CTest target to compare the shared
relaxation kernel's SSE2 and scalar results with the Maya reference, bit for
bit, with hip-plane clipping enabled and disabled. The boundary cases include
negative, zero, and positive plane distances: clipping preserves points at or
behind the plane, while the default unclipped path pushes interior points at
all three distances. The comparison also covers odd counts, inactive SIMD
lanes, sequential rings, homogeneous coordinates, and nonfinite inputs.

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

`visual/wave_expression_gallery.py` samples the wave deformer with mayapy and
writes a self-contained review gallery (`index.html`, `samples.json`, and a
contact sheet) to `--output-dir`. `visual/rasterize_contact_sheet.py` turns the
SVG contact sheet into a PNG with Pillow. The outputs are generated files and
are not committed; `tests/visual/output/` is ignored.


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
list of six node names, IDs, API kinds, and locator draw classifications.
It also checks the public Python module and prefixed custom node names.
The identity contract is [ADR-0003](../docs/adr/0003-ydd-plugin-identity.md).
