#pragma once

#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <cstddef>
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
    bool distalEnd = false;
    double distalLength = 0.0;
    double distalWidth = 0.0;
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

class BellRowCorrespondence
{
    friend class BellColliderSolver;
    struct Entry
    {
        std::size_t left = 0;
        std::size_t right = 0;
        double lambda = 0.0;
        bool matched = false;
        bool interpolate = false;
    };
    std::size_t sourceCount = 0;
    bool valid = false;
    std::vector<Entry> entries;
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
    // Material width over which an open row component fades out of a ring
    // rotation as the contact position approaches one of its free ends (0
    // keeps the full rotation up to the end).
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
  static void prepareDistalEnd(PreparedBellRing &ring, const MMatrix &matrix, bool distal);
  static MStatus solveRow(const BellRowInputs &inputs, const MPointArray &baseRow,
                         const BellRowTopology &topology, BellRowOutputs &outputs);
  static MStatus smoothDisplacements(std::vector<MVector> &displacements, double smoothness,
                                    const BellRowTopology &topology);
  static MStatus computeLocalFollow(const std::vector<MVector> &directDisplacements,
                                   const BellRowTopology &topology, const MMatrix &bellMatrix,
                                   double followRange, std::vector<MVector> &follow);
  static MStatus transferDirectField(const BellDirectField &source, const BellRowTopology &destination,
                                    BellRowTransfer &transfer);
  static MStatus buildRowCorrespondence(const BellRowTopology &source, const BellRowTopology &destination,
                                        BellRowCorrespondence &correspondence);
  static MStatus transferRowValues(const std::vector<MVector> &values, const BellRowTopology &source,
                                   const BellRowTopology &destination, std::vector<MVector> &out,
                                   std::vector<bool> &matched);
  static MStatus transferRowValues(const std::vector<MVector> &values, const BellRowCorrespondence &correspondence,
                                   std::vector<MVector> &out, std::vector<bool> &matched);
  static MStatus projectSuspendedRow(const MPointArray &anchorFinal, const MPointArray &anchorBase,
                                     const BellRowTopology &anchorTopology, const MPointArray *aboveFinal,
                                     const BellRowTopology *aboveTopology, const MPointArray &base,
                                     const MPointArray &current, const BellRowTopology &topology, MPointArray &out);
  static MStatus projectSuspendedRow(const MPointArray &anchorFinal, const MPointArray &anchorBase,
                                     const BellRowCorrespondence &anchor, const MPointArray *aboveFinal,
                                     const BellRowCorrespondence *above, const MPointArray &base,
                                     const MPointArray &current, MPointArray &out);
  static MStatus smoothRowDisplacements(const MPointArray &upBase, const MPointArray &upCurrent,
                                        const BellRowTopology &upTopology, const MPointArray &base,
                                        const MPointArray &current, const BellRowTopology &topology,
                                        const MPointArray &downBase, const MPointArray &downCurrent,
                                        const BellRowTopology &downTopology, MPointArray &out);
  static MStatus smoothRowDisplacements(const MPointArray &upBase, const MPointArray &upCurrent,
                                        const BellRowCorrespondence &up, const MPointArray &base,
                                        const MPointArray &current, const MPointArray &downBase,
                                        const MPointArray &downCurrent, const BellRowCorrespondence &down,
                                        MPointArray &out);
  // componentWeights, when given, receives one weight per component for each
  // ring (1 for rings without a rotation).
  // Relax toward one ring, faded for the points of components that ring lifts
  // with a partial weight; complete for weights 0 and 1.
  static void relaxRowFaded(MPointArray &points, const PreparedBellRing &ring, double collision,
                            bool capAtRingOrigin, const BellRowTopology &topology,
                            const std::vector<double> &componentWeights);
  static void relaxRowFaded(MPointArray &points, const PreparedBellRing &ring, double collision,
                            bool capAtRingOrigin, const BellRowTopology &topology,
                            const std::vector<double> &componentWeights, const std::vector<MVector> &directions);
  static std::vector<MVector> rowDirections(const MPointArray &points, const MPointArray &base);
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
  static void relaxTowardRingBoundary(MPointArray &points, const PreparedBellRing &ring, double collision,
                                      int startIndex, int count, bool capAtRingOrigin,
                                      const std::vector<MVector> &directions, const std::vector<double> &betaScales);
  static void smoothDisplacements(std::vector<MVector> &displacements, double smoothness);
  static MStatus solve(const BellColliderInputs &inputs, const MPointArray &baseBellPoints,
                       BellColliderOutputs &outputs);

    static void deformPoints(const BellColliderInputs& inputs, const std::vector<PreparedBellRing>& rings, const MPointArray& baseBellPoints, const Plane& bellPlane, std::vector<MPointArray>& bellPointsList);
    static MVector mergeDisplacement(const MPointArray& basePoints,
                                     const std::vector<MPointArray>& ringPoints, unsigned int index);
private:
  static void transferRowValuesImpl(const BellRowTopology &source, const BellRowTopology &destination,
                                    BellRowCorrespondence &correspondence);
  static void applyRowCorrespondence(const std::vector<MVector> &values, const BellRowCorrespondence &correspondence,
                                     std::vector<MVector> &out, std::vector<bool> &matched);
  static void averageDisplacements(int bellSubdivision, const MPointArray &baseBellPoints,
                                   const std::vector<MPointArray> &bellPointsList, MPointArray &outBellPoints);
};
