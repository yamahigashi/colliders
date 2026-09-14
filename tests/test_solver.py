"""Regression coverage for the shared Bell point solver and sparse ring inputs."""

import math
import unittest

import maya.cmds as cmds
import maya.api.OpenMaya as om


def _distance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _plug(node, attribute):
    selection = om.MSelectionList()
    selection.add(node)
    return om.MFnDependencyNode(selection.getDependNode(0)).findPlug(attribute, False)


def _points(node, attribute="outputBellMesh"):
    data = _plug(node, attribute).asMObject()
    points = om.MFnMesh(data).getPoints() if attribute == "outputBellMesh" else om.MFnNurbsCurve(data).cvPositions()
    return [tuple(p)[:3] for p in points]


def _matrix(translation=(0, 0, 0), rotation=(0, 0, 45), radius=0.4):
    transform = om.MTransformationMatrix()
    transform.setScale((radius, 1.5, radius), om.MSpace.kTransform)
    transform.setRotation(om.MEulerRotation(*(math.radians(x) for x in rotation)))
    transform.setTranslation(om.MVector(*translation), om.MSpace.kTransform)
    return list(transform.asMatrix())


def _bell(indices, matrices, translation=(0, 0, 0), collision=0):
    node = cmds.createNode("yddBellCollider")
    matrix = om.MTransformationMatrix()
    matrix.setTranslation(om.MVector(*translation), om.MSpace.kTransform)
    cmds.setAttr(node + ".bellMatrix", *list(matrix.asMatrix()), type="matrix")
    cmds.setAttr(node + ".collision", collision)
    for index, ring in zip(indices, matrices):
        cmds.setAttr("{}.ringMatrix[{}]".format(node, index), *ring, type="matrix")
    return node


class BellSolverTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def assertPointsClose(self, first, second, tolerance=2e-6):
        self.assertEqual(len(first), len(second))
        self.assertGreater(len(first), 0)
        self.assertLessEqual(max(_distance(a, b) for a, b in zip(first, second)), tolerance)

    def test_valid_and_missing_intersections_are_translation_invariant(self):
        for radius in (0.4, 2.0):
            for collision in (0, 1):
                with self.subTest(radius=radius, collision=collision):
                    origin = _bell([0], [_matrix(radius=radius)], collision=collision)
                    moved = _bell([0], [_matrix((0, 10, 0), radius=radius)], (0, 10, 0), collision)
                    for attribute in ("outputBellMesh", "outputCurve"):
                        expected = _points(origin, attribute)
                        actual = [(x, y - 10, z) for x, y, z in _points(moved, attribute)]
                        self.assertPointsClose(expected, actual)

    def test_missing_ring_intersection_preserves_relaxation(self):
        ring = _matrix(radius=2)
        baseline = _bell([], [])
        unchanged = _bell([0], [ring])
        relaxed = _bell([0], [ring], collision=1)
        self.assertPointsClose(_points(baseline), _points(unchanged), 0)
        first, second = _points(unchanged), _points(relaxed)
        self.assertEqual(len(first), len(second))
        self.assertGreater(max(_distance(a, b) for a, b in zip(first, second)), 0.1)
        self.assertPointsClose(first[:17], second[:17], 0)

    def test_tangent_and_missed_bell_intersections_leave_shape_unchanged(self):
        baseline = _bell([], [])
        for z in (1.001, 2.0):
            with self.subTest(z=z):
                node = _bell([4], [_matrix((0, 0, z))])
                self.assertPointsClose(_points(baseline), _points(node), 0)
                moved = _bell([4], [_matrix((0, 10, z))], (0, 10, 0))
                actual = [(x, y - 10, z) for x, y, z in _points(moved)]
                self.assertPointsClose(_points(node), actual)

    def test_curve_retains_double_solution_and_mesh_rounds_to_float(self):
        import struct

        node = _bell([0], [_matrix(rotation=(23, 11, 47))])
        mesh = _points(node)
        curve = _points(node, "outputCurve")
        self.assertEqual(len(mesh), 33)
        self.assertEqual(len(curve), 17)
        self.assertEqual(curve[0], curve[-1])
        rounded = [tuple(struct.unpack("f", struct.pack("f", value))[0] for value in point) for point in curve[:-1]]
        self.assertEqual(mesh[17:], rounded)
        self.assertNotEqual(mesh[17:], curve[:-1])

    def test_sparse_indices_match_dense_outputs(self):
        rings = [_matrix(), _matrix(rotation=(40, 0, -30))]
        for mode in ("off", "serial", "parallel"):
            cmds.evaluationManager(mode=mode)
            dense = _bell([0, 1], rings)
            for indices in ([0, 5], [3, 9]):
                with self.subTest(mode=mode, indices=indices):
                    sparse = _bell(indices, rings)
                    for attribute in ("outputCurve", "outputBellMesh"):
                        self.assertPointsClose(_points(dense, attribute), _points(sparse, attribute), 0)

    def test_sparse_remove_and_reinsert_and_empty(self):
        rings = [_matrix(), _matrix(rotation=(40, 0, -30))]
        sparse = _bell([3, 9], rings)
        _points(sparse)
        cmds.removeMultiInstance(sparse + ".ringMatrix[3]", b=True)
        single = _bell([0], [rings[1]])
        self.assertPointsClose(_points(single), _points(sparse), 0)
        cmds.setAttr(sparse + ".ringMatrix[3]", *rings[0], type="matrix")
        dense = _bell([0, 1], rings)
        self.assertPointsClose(_points(dense), _points(sparse), 0)
        for index in (3, 9):
            cmds.removeMultiInstance(sparse + ".ringMatrix[{}]".format(index), b=True)
        empty = _bell([], [])
        for attribute in ("outputBellMesh", "outputCurve"):
            self.assertPointsClose(_points(empty, attribute), _points(sparse, attribute), 0)


class SkirtLegProfileSolverTests(unittest.TestCase):
    def setUp(self):
        cmds.file(new=True, force=True)

    def skirt(self, bent=False):
        from test_deformers import matrix

        node = cmds.createNode("yddSkirtBellCollider")
        for side in ("left", "right"):
            for joint, y in (("Hip", 0), ("Knee", 2), ("Heel", 4)):
                cmds.setAttr(node + "." + side + joint + "Matrix", *matrix(0, y, y if bent else 0), type="matrix")
            cmds.setAttr(node + "." + side + "RingAxis", 1)
        cmds.setAttr(node + ".bellMatrix", *matrix(0, 0 if bent else -2, 0), type="matrix")
        cmds.setAttr(node + ".ringScale", 1, 1, 1, type="double3")
        cmds.setAttr(node + ".bellScale", 2 if bent else 0.1, 1, 2 if bent else 0.1, type="double3")
        for name, value in dict(skirtType=0, height=1, tightness=1, follow=0, collision=1, bellAxis=1).items():
            cmds.setAttr(node + "." + name, value)
        return node

    def surface_rows(self, node):
        from helpers import output_object

        # The data object must stay referenced while the function set reads it.
        surface = output_object(node, "outputSurface")
        fn = om.MFnNurbsSurface(surface)
        points = fn.cvPositions()
        return [[tuple(points[u * fn.numCVsInV + v])[:3] for u in range(fn.numCVsInU - 3)] for v in range(fn.numCVsInV)]

    def test_knee_profile_shrinks_row_with_and_without_follow(self):
        for follow in (0, 0.5):
            with self.subTest(follow=follow):
                node = self.skirt(bent=True)
                cmds.setAttr(node + ".follow", follow)
                baseline = self.surface_rows(node)
                for axis in ("X", "Z"):
                    cmds.setAttr(node + ".kneeRadius" + axis, 0.7)
                rows = self.surface_rows(node)
                front = max(range(len(rows[0])), key=lambda i: rows[0][i][2])

                def leg_radius(point):
                    return math.hypot(point[0], (point[2] - point[1]) / math.sqrt(2))

                self.assertEqual(rows[0], baseline[0])
                self.assertLess(leg_radius(rows[2][front]), leg_radius(rows[0][front]))
                self.assertLess(leg_radius(rows[2][front]), leg_radius(baseline[2][front]))
                self.assertTrue(all(math.isfinite(value) for row in rows for point in row for value in point))

    def test_collision_rotation_moves_line_point_toward_hip(self):
        radii = []
        for radius in (1, 0.7):
            node = _bell([0], [_matrix(rotation=(0, 0, 45), radius=radius)], collision=0)
            point = min(_points(node, "outputCurve"), key=lambda p: p[0])
            projected_radius = abs((point[0] + point[1]) / math.sqrt(2))
            self.assertAlmostEqual(projected_radius, radius, delta=0.002)
            radii.append(projected_radius)
        self.assertLess(radii[1], radii[0])

    def test_both_nodes_agree_at_profile_stations(self):
        from test_deformers import DeformerTests, matrix

        fixture = DeformerTests()
        values = dict(
            thighRadiusX=1.2,
            thighRadiusZ=0.8,
            kneeRadiusX=0.5,
            kneeRadiusZ=0.4,
            calfRadiusX=0.9,
            calfRadiusZ=0.7,
            ankleRadiusX=0.8,
            ankleRadiusZ=0.6,
        )
        for station, height, row, x_radius, z_radius in (
            ("thigh", 0.5, 1, 1.2, 0.8),
            ("knee", 0.5, 2, 0.5, 0.4),
            ("calf", 0.5, 3, 0.9, 0.7),
            ("ankle", 1, 3, 0.8, 0.6),
        ):
            with self.subTest(station=station):
                node = self.skirt()
                cmds.setAttr(node + ".skirtType", 1)
                cmds.setAttr(node + ".height", height)
                cmds.setAttr(node + ".ringScale", 1.3, 1, 0.9, type="double3")
                fixture.set_profile(node, **values)
                surface_row = self.surface_rows(node)[row]
                y = dict(thigh=1, knee=2, calf=3, ankle=4)[station]
                mesh = fixture.mesh([(0.05, y, 0), (0, y, 0.05), (0.05, y, 0.05)])
                deformer = fixture.collision(mesh, long=station in ("calf", "ankle"))
                for side in ("left", "right"):
                    for joint, joint_y in (("Hip", 0), ("Knee", 2), ("Heel", 4)):
                        cmds.setAttr(deformer + "." + side + joint + "Matrix", *matrix(0, joint_y), type="matrix")
                cmds.setAttr(deformer + ".ringScale", 1.3, 1, 0.9, type="double3")
                fixture.set_profile(deformer, **values)
                points = fixture.points(mesh)
                self.assertAlmostEqual(points[0][0], 0.9 * z_radius, delta=2e-6)
                self.assertAlmostEqual(points[1][2], 1.3 * x_radius, delta=2e-6)
                self.assertAlmostEqual(max(p[0] for p in surface_row), points[0][0], delta=2e-6)
                self.assertAlmostEqual(max(p[2] for p in surface_row), points[1][2], delta=2e-6)

    def test_rows_above_hip_ignore_profile_and_follow_uses_row_profile(self):
        node = self.skirt()
        # At the minimum height the hem level sits just below the hip, so the waist
        # and mid rows are above the hip and must ignore the profile; only the hem
        # row may move.
        cmds.setAttr(node + ".height", 0.01)
        baseline = self.surface_rows(node)
        for station in ("thigh", "knee", "calf", "ankle"):
            for axis in ("X", "Z"):
                cmds.setAttr(node + "." + station + "Radius" + axis, 1.7)
        scaled = self.surface_rows(node)
        self.assertEqual(scaled[:2], baseline[:2])
        self.assertTrue(all(math.isfinite(value) for row in scaled for point in row for value in point))
        cmds.setAttr(node + ".height", 1)
        cmds.setAttr(node + ".follow", 0.5)
        for axis in ("X", "Z"):
            cmds.setAttr(node + ".kneeRadius" + axis, 0.7)
        rows = self.surface_rows(node)
        self.assertAlmostEqual(max(p[0] for p in rows[1]), 1.7, delta=2e-6)
        self.assertAlmostEqual(max(p[0] for p in rows[2]), 0.7, delta=2e-6)
