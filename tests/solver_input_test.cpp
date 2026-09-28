#include "bellColliderSolver.h"
#include "utils.hpp"

#include <climits>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <maya/MLibrary.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

MPointArray pointsFor(int sides) {
  return BellColliderSolver::makeBellPoints(MMatrix(), 1, sides);
}

void requireEmptyPoints(const MPointArray &points) {
  require(points.length() == 0, "invalid bell point input was accepted");
}

void requireNullObject(const MObject &object) {
  require(object.isNull(), "invalid mesh/curve input was accepted");
}

void testPointInputs() {
  for (int invalid : {-1, 0, 2, 4097, INT_MAX})
    requireEmptyPoints(pointsFor(invalid));

  requireEmptyPoints(BellColliderSolver::makeBellPoints(MMatrix(), 3, 3));

  const MPointArray minimum = pointsFor(3);
  require(minimum.length() == 7, "minimum bell point count is wrong");
  const MPointArray maximum = pointsFor(4096);
  require(maximum.length() == 8193, "maximum bell point count is wrong");
}

void testMeshAndCurveInputs() {
  for (int invalid : {-1, 0, 2, 4097, INT_MAX}) {
    requireNullObject(BellColliderSolver::makeBellMesh(MPointArray(), invalid));
    requireNullObject(BellColliderSolver::makeBellCurve(MPointArray(), invalid));
  }

  const MPointArray valid = pointsFor(3);
  MPointArray tooShort = valid;
  tooShort.setLength(6);
  MPointArray tooLong = valid;
  tooLong.append(MPoint());
  requireNullObject(BellColliderSolver::makeBellMesh(tooShort, 3));
  requireNullObject(BellColliderSolver::makeBellMesh(tooLong, 3));
  requireNullObject(BellColliderSolver::makeBellCurve(tooShort, 3));
  requireNullObject(BellColliderSolver::makeBellCurve(tooLong, 3));

  require(!BellColliderSolver::makeBellMesh(valid, 3).isNull(),
          "valid bell mesh was rejected");
  require(!BellColliderSolver::makeBellCurve(valid, 3).isNull(),
          "valid bell curve was rejected");
}

void testSolverInputs() {
  BellColliderInputs inputs;
  inputs.bellMatrix = MMatrix();
  inputs.rings.clear();
  BellColliderOutputs outputs;
  outputs.points.append(MPoint(99, 99, 99));
  outputs.meanDisplacement = MVector(99, 99, 99);

  const MPointArray valid = pointsFor(3);
  for (int invalid : {-1, 0, 2, 4097, INT_MAX}) {
    inputs.bellSubdivision = invalid;
    require(BellColliderSolver::solve(inputs, valid, outputs) ==
                MS::kInvalidParameter,
            "invalid subdivision was accepted by solve");
    require(outputs.points.length() == 0, "invalid solve retained output points");
    require(outputs.meanDisplacement == MVector(0, 0, 0),
            "invalid solve retained mean displacement");
  }

  inputs.bellSubdivision = 3;
  MPointArray wrongLength = valid;
  wrongLength.setLength(6);
  require(BellColliderSolver::solve(inputs, wrongLength, outputs) ==
              MS::kInvalidParameter,
          "mismatched base points were accepted by solve");
  require(outputs.points.length() == 0, "mismatched solve retained output points");

  require(BellColliderSolver::solve(inputs, valid, outputs) == MS::kSuccess,
          "valid empty-ring solve failed after invalid inputs");
  require(outputs.points.length() == valid.length(),
          "valid empty-ring solve returned wrong point count");
}

// Historical per-point generation, kept verbatim as the reference for the
// shared circle table.
MPointArray referenceBellPoints(const MMatrix &matrix, unsigned int axis, int numSides,
                                double height, double bottomRadius, double topRadius) {
  MPointArray points;
  points.append(MPoint(0, 0, 0) * matrix);
  for (int row = 0; row < 2; ++row) {
    const double radius = row == 0 ? bottomRadius : topRadius;
    const double offset = row == 0 ? 0.0 : height;
    for (int i = 0; i < numSides; ++i) {
      const double rad = (double)i / numSides * 2 * M_PI;
      const double x = radius * cos(rad);
      const double z = radius * sin(rad);
      MPoint p;
      switch (axis) {
      case 0:
        p = MPoint(offset, x, z);
        break;
      case 1:
        p = MPoint(x, offset, z);
        break;
      case 2:
        p = MPoint(x, z, offset);
        break;
      }
      points.append(p * matrix);
    }
  }
  return points;
}

bool bitwiseEqual(const MPointArray &a, const MPointArray &b) {
  if (a.length() != b.length())
    return false;
  for (unsigned int i = 0; i < a.length(); ++i) {
    if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z || a[i].w != b[i].w)
      return false;
  }
  return true;
}

void testCircleTable() {
  for (int invalid : {-1, 0, 2, 4097, INT_MAX}) {
    const BellCircleTable table(invalid);
    require(!table.valid() && table.numSides() == 0, "invalid circle table was accepted");
    requireEmptyPoints(BellColliderSolver::makeBellPoints(MMatrix(), 1, table));
  }
  requireEmptyPoints(BellColliderSolver::makeBellPoints(MMatrix(), 3, BellCircleTable(3)));

  const double values[4][4] = {{0.7, 0.2, -0.1, 0.0},
                               {0.3, 1.9, 0.4, 0.0},
                               {-0.5, 0.1, 1.2, 0.0},
                               {3.5, -2.25, 10.125, 1.0}};
  const MMatrix skewed(values);
  const MMatrix matrices[] = {MMatrix(), skewed};
  for (int sides : {3, 4, 16, 33, 4096}) {
    const BellCircleTable table(sides);
    require(table.valid() && table.numSides() == sides, "circle table size is wrong");
    for (const MMatrix &matrix : matrices) {
      for (unsigned int axis = 0; axis < 3; ++axis) {
        const MPointArray reference = referenceBellPoints(matrix, axis, sides, 1.5, 0.8, 1.1);
        const MPointArray shared = BellColliderSolver::makeBellPoints(matrix, axis, table, 1.5, 0.8, 1.1);
        const MPointArray direct = BellColliderSolver::makeBellPoints(matrix, axis, sides, 1.5, 0.8, 1.1);
        require(bitwiseEqual(reference, shared), "circle table changed bell points");
        require(bitwiseEqual(reference, direct), "subdivision overload changed bell points");
      }
    }
    // The same table serves several bells in one evaluation.
    const MPointArray first = BellColliderSolver::makeBellPoints(skewed, 1, table, 1.0, 0.5, 1.0);
    const MPointArray second = BellColliderSolver::makeBellPoints(skewed, 1, table, 1.0, 0.5, 1.0);
    require(bitwiseEqual(first, second), "circle table reuse is not deterministic");
  }
}

void testMatrixInputs() {
  const MMatrix matrix;
  require(maxis(matrix, 0).isEquivalent(MVector(1, 0, 0), 1e-12),
          "valid matrix axis lookup failed");
  for (unsigned int invalid : {4u, 5u, UINT_MAX})
    require(maxis(matrix, invalid).length() == 0,
            "invalid matrix row was read");
}

} // namespace

int main(int argc, char **argv) {
  MStatus status = MLibrary::initialize(argv[0], true);
  if (status != MS::kSuccess)
    return 2;
  try {
    testPointInputs();
    testMeshAndCurveInputs();
    testSolverInputs();
    testCircleTable();
    testMatrixInputs();
    std::cout << "solver input tests passed\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    MLibrary::cleanup(1, false);
    return 1;
  }
  MLibrary::cleanup(0, false);
  return 0;
}
