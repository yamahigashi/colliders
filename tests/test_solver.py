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
    node = cmds.createNode("bellCollider")
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
