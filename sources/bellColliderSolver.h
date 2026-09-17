#pragma once

#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <vector>

#include "utils.hpp"

struct PreparedBellRing
{
    explicit PreparedBellRing(const MMatrix &matrix)
        : inverse(matrix.inverse()), direction(maxis(matrix, 1)), translation(taxis(matrix)),
          normal(direction.normal()), plane(translation, normal)
    {
    }

    MMatrix inverse;
    MVector direction;
    MPoint translation;
    MVector normal;
    Plane plane;
};

// Unit circle samples shared by every bell built with the same subdivision in
// one evaluation. The angle expression matches the historical per-point code so
// the resulting coordinates stay bitwise identical.
struct BellCircleTable
{
    explicit BellCircleTable(int numSides);

    bool valid() const { return !cosines.empty(); }
    int numSides() const { return static_cast<int>(cosines.size()); }

    std::vector<double> cosines;
    std::vector<double> sines;
};

struct BellColliderInputs
{
    MMatrix bellMatrix;
    std::vector<PreparedBellRing> rings;
    int bellSubdivision = 16;
    float falloff = 0.0f;
    float collision = 0.0f;
    bool capAtRingOrigin = false;
    double smoothness = 0.0;
    double followGain = 0.0;
};

struct BellColliderOutputs
{
    MPointArray points;
    MVector meanDisplacement = MVector(0, 0, 0);
};

class BellColliderSolver
{
public:
  static MPointArray makeBellPoints(const MMatrix &matrix, unsigned int axis, int numSides, double height = 1,
                                    double bottomRadius = 1, double topRadius = 1);
  static MPointArray makeBellPoints(const MMatrix &matrix, unsigned int axis, const BellCircleTable &circle,
                                    double height = 1, double bottomRadius = 1, double topRadius = 1);
  static MObject makeBellMesh(const MPointArray &points, int numSides);
  static MObject makeBellCurve(const MPointArray &points, int bellSubdivision);
  static void roundMeshPoints(MPointArray &points);
  static bool collisionPoints(const MMatrix &bellMatrix, const MMatrix &bellInverse, const Plane &bellPlane,
                              const PreparedBellRing &ring, MPoint &bellPoint, MPoint &ringPoint, MPoint &linePoint);
  static void relaxTowardRingBoundary(MPointArray &points, const PreparedBellRing &ring, double collision,
                                      int startIndex, int count, bool capAtRingOrigin = false);
  static void smoothDisplacements(std::vector<MVector> &displacements, double smoothness);
  static MStatus solve(const BellColliderInputs &inputs, const MPointArray &baseBellPoints,
                       BellColliderOutputs &outputs);

private:
    static void deformPoints(const BellColliderInputs& inputs, const std::vector<PreparedBellRing>& rings, const MPointArray& baseBellPoints, const Plane& bellPlane, std::vector<MPointArray>& bellPointsList);
    static void averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const std::vector<MPointArray>& bellPointsList, MPointArray& outBellPoints, bool useUnnormalized);
};
