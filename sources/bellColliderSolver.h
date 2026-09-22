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

struct BellRowSide
{
    int seamIndex = -1;
    int bank = 0;
};

struct BellRowComponent
{
    double startU = 0.0;
    double endU = 1.0;
    bool closed = false;
};

struct BellRowVertex
{
    double materialU = 0.0;
    BellRowSide side;
    int previous = -1;
    int next = -1;
    int componentId = 0;
    int panelId = 0;
    std::vector<unsigned int> outputDuplicates;
};

struct BellRowTopology
{
    std::vector<BellRowVertex> vertices;
    std::vector<BellRowComponent> components;
    unsigned int outputCount = 0;
};

struct BellRowInputs
{
    MMatrix bellMatrix;
    std::vector<PreparedBellRing> rings;
    float falloff = 0.0f;
    float collision = 0.0f;
    bool capAtRingOrigin = false;
    double smoothness = 0.0;
    double followGain = 0.0;
    double followRange = 0.1875;
    // Material width over which a row component fades out of a ring rotation
    // once the ring points past its end (0 keeps only the components that
    // contain the contact position).
    double contactBlendWidth = 0.0;
    int physicalLevel = 0;
};

struct BellDirectField
{
    BellRowTopology topology;
    std::vector<MVector> values;
};

struct BellRowTransfer
{
    std::vector<MVector> values;
};

struct BellRowOutputs
{
    // Points are ordered by independent CV; the caller expands output duplicates.
    MPointArray points;
    std::vector<MVector> directDisplacements;
    BellDirectField directField;
    // One weight per component for each input ring, in ring order.
    std::vector<std::vector<double>> componentWeights;
};
class BellColliderSolver
{
public:
  static MStatus solveRow(const BellRowInputs &inputs, const MPointArray &baseRow,
                         const BellRowTopology &topology, BellRowOutputs &outputs);
  static MStatus smoothDisplacements(std::vector<MVector> &displacements, double smoothness,
                                    const BellRowTopology &topology);
  static MStatus computeLocalFollow(const std::vector<MVector> &directDisplacements,
                                   const BellRowTopology &topology, const MMatrix &bellMatrix,
                                   double followRange, std::vector<MVector> &follow);
  static MStatus transferDirectField(const BellDirectField &source, const BellRowTopology &destination,
                                    BellRowTransfer &transfer);
  // componentWeights, when given, receives one weight per component for each
  // ring (1 for rings without a rotation).
  // Relax toward one ring, faded for the points of components that ring lifts
  // with a partial weight (see the row contract); complete for weights 0 and 1.
  static void relaxRowFaded(MPointArray &points, const PreparedBellRing &ring, double collision,
                            bool capAtRingOrigin, const BellRowTopology &topology,
                            const std::vector<double> &componentWeights);
  static MStatus deformPoints(const BellRowInputs &inputs, const MPointArray &baseRow,
                              const BellRowTopology &topology,
                              const std::vector<std::vector<unsigned int>> &ringSeeds,
                              std::vector<MPointArray> &ringPoints,
                              std::vector<std::vector<double>> *componentWeights = nullptr);
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

    static void deformPoints(const BellColliderInputs& inputs, const std::vector<PreparedBellRing>& rings, const MPointArray& baseBellPoints, const Plane& bellPlane, std::vector<MPointArray>& bellPointsList);
private:
    static void averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const std::vector<MPointArray>& bellPointsList, MPointArray& outBellPoints, bool useUnnormalized);
};
