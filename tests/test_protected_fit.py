"""Protected-band surface fit tests for yddSkirtSurfaceFit."""

import hashlib
import math
import struct
import unittest

import maya.api.OpenMaya as om
import maya.cmds as cmds

from helpers import coordinates, output_object, skirt, transform


def distance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def point_hash(points):
    digest = hashlib.sha256()
    for point in points:
        digest.update(struct.pack("<4d", *point))
    return digest.hexdigest()


def set_double_array(node, values):
    plug = om.MSelectionList().add(node + ".protectedVParameters").getPlug(0)
    plug.setMObject(om.MFnDoubleArrayData().create(om.MDoubleArray(values)))


def tangent_angle(first, second):
    a, b = om.MVector(*first), om.MVector(*second)
    if a.length() == 0.0 or b.length() == 0.0:
        raise AssertionError("A tangent angle requires two nonzero vectors")
    return math.degrees(math.atan2((a ^ b).length(), a * b))


def local_full_knots(spans):
    return [0.0] * 4 + [float(j) / spans for j in range(1, spans)] + [1.0] * 4


def old_sample_rows(input_rows, bounds, spans):
    rows = list(input_rows)
    for boundary in bounds[1:-1]:
        matched = [j for j, value in enumerate(rows) if abs(value - boundary) <= 1e-9]
        if len(matched) != 1:
            raise AssertionError("Expected one input row per protected boundary")
        rows[matched[0]] = boundary
    matrix = []
    for band, count in enumerate(spans):
        full = local_full_knots(count)
        for local in range(0 if band == 0 else 1, count + 3):
            greville = 0.0
            for value in full[local + 1:local + 4]:
                greville += value
            greville /= 3.0
            sample = bounds[band + 1] if local == count + 2 else (
                bounds[band] + (bounds[band + 1] - bounds[band]) * greville
            )
            coefficients = [0.0] * len(rows)
            in_band = [j for j, value in enumerate(rows) if bounds[band] <= value <= bounds[band + 1]]
            exact = [j for j in in_band if abs(rows[j] - sample) <= 1e-12]
            if exact:
                coefficients[exact[0]] = 1.0
            else:
                for left, right in zip(in_band, in_band[1:]):
                    if rows[left] <= sample <= rows[right]:
                        weight = (sample - rows[left]) / (rows[right] - rows[left])
                        coefficients[left], coefficients[right] = 1.0 - weight, weight
                        break
                else:
                    raise AssertionError("Unbracketed sample")
            matrix.append(coefficients)
    return matrix


def sampled_columns(columns, matrix):
    output = []
    for column in columns:
        for row in matrix:
            point = [0.0, 0.0, 0.0, 1.0]
            for j, coefficient in enumerate(row):
                for c in range(3):
                    point[c] += coefficient * column[j][c]
            output.append(point)
    return output


def de_boor(points, knots, parameter, degree=3):
    count = len(points)
    span = count - 1
    for i in range(degree, count):
        if knots[i] <= parameter < knots[i + 1]:
            span = i
            break
    values = [list(points[span - degree + j][:3]) for j in range(degree + 1)]
    for level in range(1, degree + 1):
        for j in range(degree, level - 1, -1):
            index = span - degree + j
            alpha = (parameter - knots[index]) / (knots[index + degree + 1 - level] - knots[index])
            values[j] = [(1.0 - alpha) * a + alpha * b for a, b in zip(values[j - 1], values[j])]
    return values[degree]


class ProtectedFitTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)
        cmds.evaluationManager(mode="off")
        self.addCleanup(cmds.evaluationManager, mode="off")
        self.errors = []
        self.surface_data = []

        def capture(message, message_type, _client_data):
            if message_type == om.MCommandMessage.kError:
                self.errors.append(message)

        callback = om.MCommandMessage.addCommandOutputCallback(capture)
        self.addCleanup(om.MMessage.removeCallback, callback)

    def points(self, node):
        return coordinates(node, "outputSurface", "surface")

    def shape_points(self, fit):
        shapes = cmds.listConnections(
            fit + ".outputSurface", source=False, destination=True, shapes=True, type="nurbsSurface"
        ) or []
        self.assertEqual(len(shapes), 1)
        dag_path = om.MSelectionList().add(shapes[0]).getDagPath(0)
        points = om.MFnNurbsSurface(dag_path).cvPositions(om.MSpace.kObject)
        return [[p.x, p.y, p.z, p.w] for p in points]

    def surface(self, node):
        data = output_object(node, "outputSurface")
        self.surface_data.append(data)
        return om.MFnNurbsSurface(data)

    def sampler(self, node):
        sampler = cmds.createNode("pointOnSurfaceInfo")
        cmds.connectAttr(node + ".outputSurface", sampler + ".inputSurface")
        cmds.setAttr(sampler + ".turnOnPercentage", False)
        return sampler

    def sample(self, sampler, u, v, attribute="position"):
        cmds.setAttr(sampler + ".parameterU", u)
        cmds.setAttr(sampler + ".parameterV", v)
        return cmds.getAttr(sampler + "." + attribute)[0]

    def assert_finite(self, points):
        self.assertTrue(points)
        self.assertTrue(all(math.isfinite(value) for point in points for value in point))

    def assert_invalid_preserves(self, fit, expected, reason, repeats=3):
        nodes = set(cmds.ls(long=True))
        connections = cmds.listConnections(fit, connections=True, plugs=True) or []
        spans = cmds.getAttr(fit + ".spansV")
        for attempt in range(repeats):
            with self.subTest(attempt=attempt):
                self.errors.clear()
                try:
                    output_object(fit, "outputSurface")
                except RuntimeError:
                    pass
                points = self.shape_points(fit)
                self.assertTrue(
                    any("yddSkirtSurfaceFit:" in message and reason in message for message in self.errors),
                    self.errors,
                )
                self.assertEqual(points, expected)
                self.assertEqual(cmds.getAttr(fit + ".spansV"), spans)
                self.assertEqual(set(cmds.ls(long=True)), nodes)
                self.assertEqual(cmds.listConnections(fit, connections=True, plugs=True) or [], connections)

    def make_linear_v(self, columns, v_knots):
        nu = len(columns)
        nv = len(columns[0])
        points = om.MPointArray()
        for u in range(nu):
            for v in range(nv):
                x, y, z = columns[u][v]
                points.append(om.MPoint(x, y, z))
        u_knots = om.MDoubleArray([float(i) / float(nu - 1) for i in range(nu)])
        parent = om.MFnTransform().create()
        fn = om.MFnNurbsSurface()
        obj = fn.create(
            points,
            u_knots,
            om.MDoubleArray(v_knots),
            1,
            1,
            om.MFnNurbsSurface.kOpen,
            om.MFnNurbsSurface.kOpen,
            False,
            parent,
        )
        return om.MFnDagNode(obj).fullPathName()

    def fit_surface(self, source_shape, spans=2):
        fit = cmds.createNode("yddSkirtSurfaceFit")
        cmds.setAttr(fit + ".spansV", spans)
        cmds.connectAttr(source_shape + ".local", fit + ".inputSurface")
        plane = cmds.nurbsPlane()[0]
        shape = cmds.listRelatives(plane, shapes=True, fullPath=True)[0]
        cmds.connectAttr(fit + ".outputSurface", shape + ".create", force=True)
        return fit

    def slit_columns(self):
        return [
            [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0), (0.0, 1.0, 0.0)],
            [(1.0, 0.0, 0.0), (1.0, 0.5, 0.0), (1.0, 1.0, 0.5)],
        ]

    def band_layout(self, surface, protected):
        bounds = [0.0] + list(protected) + [1.0]
        knots = list(surface.knotsInV())
        spans = [1 + sum(lo < value < hi for value in knots) for lo, hi in zip(bounds, bounds[1:])]
        starts = [0]
        for count in spans:
            starts.append(starts[-1] + count + 2)
        self.assertEqual(starts[-1] + 1, surface.numCVsInV)
        return bounds, spans, starts

    def local_band(self, surface, start, spans, lo, hi):
        cvs = surface.cvPositions()
        points = om.MPointArray()
        for u in range(surface.numCVsInU):
            for j in range(spans + 3):
                points.append(cvs[u * surface.numCVsInV + start + j])
        knots = [lo] * 3 + [lo + (hi - lo) * j / spans for j in range(1, spans)] + [hi] * 3
        owned = om.MFnNurbsSurfaceData().create()
        self.surface_data.append(owned)
        fn = om.MFnNurbsSurface()
        fn.create(points, surface.knotsInU(), om.MDoubleArray(knots), surface.degreeInU, 3,
                  surface.formInU, om.MFnNurbsSurface.kOpen, False, owned)
        return fn

    def assert_junctions(self, fit, protected):
        surface = self.surface(fit)
        bounds, spans, starts = self.band_layout(surface, protected)
        cvs = self.points(fit)
        sampler = self.sampler(fit)
        for band, boundary in enumerate(protected, 1):
            upper_full = local_full_knots(spans[band - 1])
            lower_full = local_full_knots(spans[band])
            last = spans[band - 1] + 2
            up_width = (bounds[band] - bounds[band - 1]) * (upper_full[last + 3] - upper_full[last])
            down_width = (bounds[band + 1] - bounds[band]) * (lower_full[4] - lower_full[1])
            self.assertLessEqual(down_width / up_width, 4.0)
            derivatives = []
            for u in range(surface.numCVsInU):
                row = u * surface.numCVsInV + starts[band]
                up = [3.0 * (cvs[row][c] - cvs[row - 1][c]) / up_width for c in range(3)]
                down = [3.0 * (cvs[row + 1][c] - cvs[row][c]) / down_width for c in range(3)]
                self.assertLess(tangent_angle(up, down), 1e-6)
                self.assertLessEqual(distance(up, down), 1e-9)
                derivatives.append(up)
            upper = self.local_band(surface, starts[band - 1], spans[band - 1], bounds[band - 1], boundary)
            lower = self.local_band(surface, starts[band], spans[band], boundary, bounds[band + 1])
            umin, umax = surface.knotDomainInU
            u_knots = list(surface.knotsInU())
            u_full = [u_knots[0]] + u_knots + [u_knots[-1]]
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
                u = umin + (umax - umin) * fraction
                expected = de_boor(derivatives, u_full, u, surface.degreeInU)
                up = upper.tangents(u, boundary)[1]
                down = lower.tangents(u, boundary)[1]
                common = surface.tangents(u, boundary)[1]
                posi = self.sample(sampler, u, boundary, "tangentV")
                self.assertLess(tangent_angle(up, down), 1e-6)
                self.assertLess(tangent_angle(up, expected), 1e-6)
                self.assertLess(tangent_angle(up, common), 1e-6)
                self.assertLess(tangent_angle(up, posi), 1e-6)

    def assert_preserved_rows(self, fit, source, protected):
        result = self.surface(fit)
        bounds, spans, starts = self.band_layout(result, protected)
        source_cvs = source.cvPositions()
        columns = [
            [tuple(source_cvs[u * source.numCVsInV + j])[:3] for j in range(source.numCVsInV)]
            for u in range(source.numCVsInU)
        ]
        input_knots = list(source.knotsInV())
        lo, hi = source.knotDomainInV
        rows = [(value - lo) / (hi - lo) for value in input_knots]
        old = sampled_columns(columns, old_sample_rows(rows, bounds, spans))
        actual = self.points(fit)
        m = result.numCVsInV
        for u in range(result.numCVsInU):
            first = slice(u * m, u * m + starts[1] + 1)
            self.assertEqual(point_hash(actual[first]), point_hash(old[first]))
            for band, boundary in enumerate(protected, 1):
                j = next(j for j, value in enumerate(rows) if abs(value - boundary) <= 1e-9)
                self.assertEqual(actual[u * m + starts[band]][:3], list(columns[u][j]))
        umin, umax = result.knotDomainInU
        for boundary in protected:
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
                u = umin + fraction * (umax - umin)
                self.assertLessEqual(distance(tuple(result.getPointAtParam(u, boundary))[:3],
                                              tuple(source.getPointAtParam(u, lo + boundary * (hi - lo)))[:3]), 1e-9)

    def test_bent_four_row_junction(self):
        rows = [0.0, 0.25, 0.5, 1.0]
        columns = [[(float(u), v, 0.5 if v == 1.0 else 0.0) for v in rows] for u in range(2)]
        self.assertGreaterEqual(tangent_angle((0, 0.25, 0), (0, 0.5, 0.5)), 45.0)
        shape = self.make_linear_v(columns, rows)
        fit = self.fit_surface(shape, spans=2)
        set_double_array(fit, [0.5])
        self.assert_junctions(fit, [0.5])
        source = om.MFnNurbsSurface(om.MSelectionList().add(shape).getDagPath(0))
        self.assert_preserved_rows(fit, source, [0.5])
        old_matrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.33333333333333337, 0.66666666666666663, 0.0, 0.0],
            [0.0, 0.66666666666666674, 0.33333333333333326, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.66666666666666674, 0.33333333333333326],
            [0.0, 0.0, 0.33333333333333348, 0.66666666666666652],
            [0.0, 0.0, 0.0, 1.0],
        ]
        old_cvs = sampled_columns(columns, old_matrix)
        actual = self.points(fit)
        for u in range(2):
            self.assertEqual(point_hash(actual[u * 7:u * 7 + 4]), point_hash(old_cvs[u * 7:u * 7 + 4]))
        matrix = [list(row) for row in old_matrix]
        matrix[4] = [2.0 * a - b for a, b in zip(matrix[3], matrix[2])]
        expected = sampled_columns(columns, matrix)
        result = self.surface(fit)
        for u in range(2):
            value = de_boor(expected[u * 7 + 3:u * 7 + 7], local_full_knots(1), 2.0 / 3.0)
            self.assertLessEqual(distance(tuple(result.getPointAtParam(float(u), 5.0 / 6.0))[:3], value), 1e-9)

    def test_empty_array_matches_existing_output(self):
        collider, _, surface = skirt()
        shape = cmds.listRelatives(surface, shapes=True, fullPath=True)[0]
        fit = cmds.createNode("yddSkirtSurfaceFit")
        cmds.setAttr(fit + ".spansV", 4)
        cmds.connectAttr(collider + ".outputSurface", fit + ".inputSurface")
        dest = cmds.listRelatives(cmds.nurbsPlane()[0], shapes=True, fullPath=True)[0]
        cmds.connectAttr(fit + ".outputSurface", dest + ".create", force=True)
        before = self.points(fit)
        self.assert_finite(before)
        set_double_array(fit, [])
        after = self.points(fit)
        self.assertEqual(point_hash(after), point_hash(before))
        self.assertEqual(after, before)

    def test_protected_differences_only_below_position(self):
        source = self.make_linear_v(self.slit_columns(), [0.0, 0.5, 1.0])
        fit = self.fit_surface(source, spans=2)
        set_double_array(fit, [0.5])
        result = self.surface(fit)
        self.assertEqual(result.numCVsInV, 7)
        self.assertEqual(result.degreeInV, 3)
        self.assertEqual(result.formInV, om.MFnNurbsSurface.kOpen)
        self.assertEqual(list(result.knotsInV()), [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0])
        cvs = self.points(fit)
        m = 7
        for k in range(4):
            left = cvs[k]
            right = cvs[m + k]
            self.assertAlmostEqual(left[2], right[2], delta=1e-9)
            self.assertAlmostEqual(left[2], 0.0, delta=1e-9)
        self.assertLessEqual(distance(cvs[3][:3], (0.0, 0.5, 0.0)), 1e-9)
        self.assertLessEqual(distance(cvs[m + 3][:3], (1.0, 0.5, 0.0)), 1e-9)
        self.assertGreater(abs(cvs[m + 6][2] - cvs[6][2]), 1e-9)
        self.assertEqual(cvs[m + 4][2] - cvs[4][2], 0.0)
        left = result.getPointAtParam(0.0, 0.75)
        right = result.getPointAtParam(1.0, 0.75)
        self.assertAlmostEqual(right.z - left.z, 3.0 / 16.0, delta=1e-12)
        self.assert_junctions(fit, [0.5])

    def test_plain_fit_quarter_bleed(self):
        source = self.make_linear_v(self.slit_columns(), [0.0, 0.5, 1.0])
        fit = self.fit_surface(source, spans=1)
        sampler = self.sampler(fit)
        left = self.sample(sampler, 0.0, 0.25)
        right = self.sample(sampler, 1.0, 0.25)
        self.assertAlmostEqual(right[2] - left[2], 1.0 / 32.0, delta=1e-12)

    def connection_collider(self, skirt_type, seams):
        collider = cmds.createNode("yddSkirtBellCollider")
        waist = transform((0, 11.5, 0), (180, 0, 0))
        cmds.connectAttr(waist + ".worldMatrix[0]", collider + ".bellMatrix")
        for side, x in (("left", -1), ("right", 1)):
            for part, y in (("Hip", 10), ("Knee", 5), ("Heel", 0)):
                joint = transform((x, y, 0))
                cmds.connectAttr(joint + ".worldMatrix[0]", collider + "." + side + part + "Matrix")
            cmds.setAttr(collider + "." + side + "RingAxis", 4)
        for name, value in (("skirtType", skirt_type), ("bellAxis", 1), ("bellSubdivision", 16)):
            cmds.setAttr(collider + "." + name, value)
        cmds.setAttr(collider + ".bellScale", 2, 1, 2, type="double3")
        cmds.setAttr(collider + ".ringScale", 2, 1, 2, type="double3")
        for index, u, height in seams:
            base = collider + ".seams[{}].".format(index)
            cmds.setAttr(base + "enabled", True)
            cmds.setAttr(base + "materialU", u)
            cmds.setAttr(base + "startHeight", height)
        return collider

    def test_fit_connections(self):
        cases = (
            (),
            ((0, 0.25, 0.4), (1, 0.75, 1.0)),
            ((0, 0.2, 0.25), (1, 0.5, 0.55), (2, 0.8, 0.75)),
        )
        for skirt_type in (1, 0):
            for seams in cases:
                with self.subTest(skirt_type=skirt_type, seams=seams):
                    collider = self.connection_collider(skirt_type, seams)
                    patches = om.MSelectionList().add(collider + ".outputPatches").getPlug(0)
                    patches.evaluateNumElements()
                    indices = sorted(patches.getExistingArrayAttributeIndices())
                    self.assertEqual(indices, list(range(1, len(seams) + 1)) if seams else [0])
                    for index in indices:
                        base = collider + ".outputPatches[{}].".format(index)
                        fit = cmds.createNode("yddSkirtSurfaceFit")
                        cmds.setAttr(fit + ".spansV", 4)
                        cmds.connectAttr(base + "surface", fit + ".inputSurface")
                        cmds.connectAttr(base + "vBreaks", fit + ".protectedVParameters")
                        self.assertTrue(cmds.isConnected(base + "surface", fit + ".inputSurface"))
                        self.assertTrue(cmds.isConnected(base + "vBreaks", fit + ".protectedVParameters"))
                        data = output_object(collider, "outputPatches[{}].vBreaks".format(index))
                        protected = list(om.MFnDoubleArrayData(data).array())
                        self.assertTrue(protected)
                        if not seams:
                            self.assertEqual(len(protected), 1)
                        else:
                            self.assertEqual(protected, sorted({height for _, _, height in seams if 0 < height < 1}))
                        owned = om.MFnNurbsSurfaceData().create()
                        om.MFnNurbsSurface().copy(output_object(collider, "outputPatches[{}].surface".format(index)), owned)
                        self.surface_data.append(owned)
                        source = om.MFnNurbsSurface(owned)
                        self.assert_junctions(fit, protected)
                        self.assert_preserved_rows(fit, source, protected)
                        baseline = self.points(fit)
                        result = self.surface(fit)
                        _, _, starts = self.band_layout(result, protected)
                        rows = list(source.knotsInV())
                        lo, hi = source.knotDomainInV
                        for band, boundary in enumerate(protected, 1):
                            injected = om.MFnNurbsSurfaceData().create()
                            om.MFnNurbsSurface().copy(owned, injected)
                            self.surface_data.append(injected)
                            changed = om.MFnNurbsSurface(injected)
                            cvs = changed.cvPositions()
                            for u in range(source.numCVsInU):
                                for j, value in enumerate(rows):
                                    if (value - lo) / (hi - lo) > boundary + 1e-9:
                                        k = u * source.numCVsInV + j
                                        point = cvs[k]
                                        point.z += (u + 1) * 0.125
                                        cvs[k] = point
                            changed.setCVPositions(cvs)
                            cmds.disconnectAttr(base + "surface", fit + ".inputSurface")
                            om.MSelectionList().add(fit + ".inputSurface").getPlug(0).setMObject(injected)
                            after = self.points(fit)
                            for u in range(result.numCVsInU):
                                first = u * result.numCVsInV
                                last = first + starts[band] + 1
                                self.assertEqual(point_hash(after[first:last]), point_hash(baseline[first:last]))
                            self.assertNotEqual(point_hash(after), point_hash(baseline))
                            self.assert_junctions(fit, protected)
                            cmds.connectAttr(base + "surface", fit + ".inputSurface", force=True)
                            self.assertEqual(point_hash(self.points(fit)), point_hash(baseline))

    def test_two_patch_input_isolation(self):
        collider = self.connection_collider(1, ((0, 0.25, 0.4), (1, 0.75, 1.0)))
        patches = om.MSelectionList().add(collider + ".outputPatches").getPlug(0)
        patches.evaluateNumElements()
        self.assertEqual(sorted(patches.getExistingArrayAttributeIndices()), [1, 2])
        fits, inputs, sources, breaks = [], [], [], []
        for index in (1, 2):
            base = collider + ".outputPatches[{}].".format(index)
            fit = cmds.createNode("yddSkirtSurfaceFit")
            cmds.setAttr(fit + ".spansV", 4)
            cmds.connectAttr(base + "surface", fit + ".inputSurface")
            cmds.connectAttr(base + "vBreaks", fit + ".protectedVParameters")
            shape = cmds.createNode("nurbsSurface")
            cmds.connectAttr(fit + ".outputSurface", shape + ".create")
            self.assertTrue(cmds.isConnected(base + "surface", fit + ".inputSurface"))
            self.assertTrue(cmds.isConnected(base + "vBreaks", fit + ".protectedVParameters"))
            owned = om.MFnNurbsSurfaceData().create()
            om.MFnNurbsSurface().copy(output_object(collider, "outputPatches[{}].surface".format(index)), owned)
            self.surface_data.append(owned)
            inputs.append(om.MFnNurbsSurface(owned))
            sources.append(base + "surface")
            fits.append(fit)
            data = output_object(collider, "outputPatches[{}].vBreaks".format(index))
            breaks.append(tuple(om.MFnDoubleArrayData(data).array()))
        self.assertEqual(breaks[0], breaks[1])
        self.assertTrue(any(abs(value - 0.4) <= 1e-9 for value in breaks[0]))
        baseline = [self.points(fit) for fit in fits]
        surfaces = [self.surface(fit) for fit in fits]
        nv = surfaces[0].numCVsInV
        self.assertEqual(nv, surfaces[1].numCVsInV)
        knots = list(surfaces[0].knotsInV())
        full = [knots[0]] + knots + [knots[-1]]
        rows = [sum(full[j + 1:j + 4]) / 3.0 for j in range(nv)]
        upper = [j for j, v in enumerate(rows) if v <= 0.4 + 1e-12]
        tip = min(range(nv), key=lambda j: abs(rows[j] - 0.4))
        self.assertAlmostEqual(rows[tip], 0.4, delta=1e-12)
        far_start = (surfaces[1].numCVsInU - 1) * nv
        for j in upper:
            self.assertLessEqual(distance(baseline[0][j][:3], baseline[1][far_start + j][:3]), 1e-9)
        samplers = [self.sampler(fit) for fit in fits]
        for i, u in ((0, 0.0), (1, 1.0)):
            expected = tuple(inputs[i].getPointAtParam(u, 0.4))[:3]
            self.assertLessEqual(distance(self.sample(samplers[i], u, 0.4), expected), 1e-9)
            evaluated = tuple(surfaces[i].getPointAtParam(u, 0.2))[:3]
            self.assertLessEqual(distance(self.sample(samplers[i], u, 0.2), evaluated), 1e-6)

        plain = cmds.createNode("yddSkirtSurfaceFit")
        cmds.setAttr(plain + ".spansV", 1)
        cmds.connectAttr(sources[0], plain + ".inputSurface")
        plain_sampler = self.sampler(plain)
        plain_before = self.sample(plain_sampler, 0.0, 0.2)
        injected = om.MFnNurbsSurfaceData().create()
        fresh_source = om.MSelectionList().add(sources[0]).getPlug(0).asMObject()
        om.MFnNurbsSurface().copy(fresh_source, injected)
        self.surface_data.append(injected)
        changed = om.MFnNurbsSurface(injected)
        cvs = changed.cvPositions()
        for j, v in enumerate(inputs[0].knotsInV()):
            if v > 0.4 + 1e-9:
                point = cvs[j]
                point.z += 1.0
                cvs[j] = point
        changed.setCVPositions(cvs)
        for fit in (fits[0], plain):
            cmds.disconnectAttr(sources[0], fit + ".inputSurface")
            om.MSelectionList().add(fit + ".inputSurface").getPlug(0).setMObject(injected)
        after = [self.points(fit) for fit in fits]
        for u in range(surfaces[0].numCVsInU):
            indices = [u * nv + j for j in upper]
            self.assertEqual(
                point_hash([baseline[0][j] for j in indices]),
                point_hash([after[0][j] for j in indices]),
            )
        self.assertEqual(point_hash(after[1]), point_hash(baseline[1]))
        self.assertNotEqual(point_hash(after[0]), point_hash(baseline[0]))
        self.assertGreater(distance(after[0][nv - 1][:3], after[1][far_start + nv - 1][:3]), 1e-9)
        self.assertGreater(distance(self.sample(plain_sampler, 0.0, 0.2), plain_before), 1e-9)
        self.assertLessEqual(
            distance(self.sample(samplers[0], 0.0, 0.4), tuple(inputs[0].getPointAtParam(0.0, 0.4))[:3]),
            1e-9,
        )
        for fit in (fits[0], plain):
            cmds.connectAttr(sources[0], fit + ".inputSurface", force=True)
        self.assertEqual(point_hash(self.points(fits[0])), point_hash(baseline[0]))

    def test_two_fits_plain_quarter_bleed(self):
        fits = []
        for hem_offset in (0.0, 0.5):
            columns = [
                [(float(u), 0.0, 0.0), (float(u), 0.5, 0.0), (float(u), 1.0, hem_offset)]
                for u in range(2)
            ]
            source = self.make_linear_v(columns, [0.0, 0.5, 1.0])
            fits.append(self.fit_surface(source, spans=1))
        samplers = [self.sampler(fit) for fit in fits]
        positions = [self.sample(sampler, 0.0, 0.25) for sampler in samplers]
        self.assertAlmostEqual(positions[1][2] - positions[0][2], 1.0 / 32.0, delta=1e-12)
        for fit in fits:
            cmds.setAttr(fit + ".spansV", 2)
            set_double_array(fit, [0.5])
        points = [self.points(fit) for fit in fits]
        for u in range(2):
            self.assertEqual(
                point_hash(points[0][u * 7:u * 7 + 4]),
                point_hash(points[1][u * 7:u * 7 + 4]),
            )
        self.assertNotEqual(point_hash(points[0]), point_hash(points[1]))
        positions = [self.sample(sampler, 0.0, 0.25) for sampler in samplers]
        self.assertEqual(positions[0], positions[1])
        for u in range(2):
            self.assertEqual(points[1][u * 7 + 4][2] - points[0][u * 7 + 4][2], 0.0)
        positions = [self.sample(sampler, 0.0, 0.75) for sampler in samplers]
        self.assertAlmostEqual(positions[1][2] - positions[0][2], 3.0 / 16.0, delta=1e-12)
        for fit in fits:
            self.assert_junctions(fit, [0.5])

    def test_rounded_protected_rows_keep_exact_coefficients(self):
        rows = [0.0, 0.1, 0.4, 0.7, 1.0]
        columns = [
            [(float(u), float(j * j), float((u + 2) * (j + 1))) for j in range(len(rows))]
            for u in range(2)
        ]
        source = self.make_linear_v(columns, rows)
        fit = self.fit_surface(source, spans=4)
        set_double_array(fit, [0.1, 0.4, 0.7])
        result = self.points(fit)
        for u in range(2):
            for boundary, output_row in enumerate((3, 6, 9), 1):
                self.assertEqual(result[u * 13 + output_row][:3], list(columns[u][boundary]))

    def test_validation_and_recovery_conditions_10_to_14(self):
        source = self.make_linear_v(
            [
                [(0.0, v, 0.0) for v in (0.0, 0.25, 0.5, 0.75, 1.0)],
                [(1.0, v, 0.0) for v in (0.0, 0.25, 0.5, 0.75, 1.0)],
            ],
            [0.0, 0.25, 0.5, 0.75, 1.0],
        )
        fit = self.fit_surface(source, spans=4)
        set_double_array(fit, [0.5])
        valid = self.shape_points(fit)
        self.assert_finite(valid)

        cases = [
            ([0.7, 0.3], "strictly increasing in (0, 1)"),
            ([0.0], "strictly increasing in (0, 1)"),
            ([1.0], "strictly increasing in (0, 1)"),
            ([0.5, 0.5], "strictly increasing in (0, 1)"),
            ([float("nan")], "strictly increasing in (0, 1)"),
            ([float("inf")], "strictly increasing in (0, 1)"),
        ]
        for values, reason in cases:
            with self.subTest(values=values):
                set_double_array(fit, [0.5])
                before = self.shape_points(fit)
                set_double_array(fit, values)
                self.assert_invalid_preserves(fit, before, reason)

        set_double_array(fit, [1e-12])
        self.assert_invalid_preserves(fit, valid, "band length times V domain length")

        set_double_array(fit, [0.25, 0.5, 0.75])
        cmds.setAttr(fit + ".spansV", 2)
        self.assert_invalid_preserves(fit, valid, "at least the number of protected bands")
        cmds.setAttr(fit + ".spansV", 4)
        set_double_array(fit, [0.5])
        self.assertEqual(self.shape_points(fit), valid)

        plane = cmds.nurbsPlane(degree=3)[0]
        plane_shape = cmds.listRelatives(plane, shapes=True, fullPath=True)[0]
        set_double_array(fit, [0.5])
        cmds.connectAttr(plane_shape + ".local", fit + ".inputSurface", force=True)
        self.assert_invalid_preserves(fit, valid, "V degree must be 1")
        cmds.connectAttr(source + ".local", fit + ".inputSurface", force=True)

        set_double_array(fit, [0.5 - 0.25e-9, 0.5 + 0.25e-9])
        self.assert_invalid_preserves(fit, valid, "each input V row must match at most one")

        set_double_array(fit, [i / 258.0 for i in range(1, 258)])
        self.assert_invalid_preserves(fit, valid, "at least the number of protected bands")

        set_double_array(fit, [0.3])
        self.assert_invalid_preserves(fit, valid, "match exactly one input V row")

        set_double_array(fit, [0.7, 0.3])
        cmds.setAttr(fit + ".spansV", 1)
        self.assert_invalid_preserves(fit, valid, "strictly increasing in (0, 1)")
        cmds.setAttr(fit + ".spansV", 4)
        set_double_array(fit, [0.5])
        recovered = self.shape_points(fit)
        self.assert_finite(recovered)
        self.assertEqual(recovered, valid)
