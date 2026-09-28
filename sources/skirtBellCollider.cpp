#include <maya/MArrayDataBuilder.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MPlug.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MRampAttribute.h>
#include <maya/MUIDrawManager.h>
#include <maya/MVector.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "bellColliderSolver.h"
#include "colliderInputValidation.h"
#include "pluginIdentity.h"
#include "skirtBellCollider.h"
#include "skirtLegBuild.h"
#include "skirtLegProfile.h"
#include "skirtRingFrames.h"
#include "utils.hpp"

using namespace std;

MTypeId SkirtBellCollider::typeId(PluginIdentity::kSkirtBellTypeId);
MString SkirtBellCollider::typeName(PluginIdentity::kSkirtBellNodeName);

MString SkirtBellCollider::drawDbClassification = MString("drawdb/geometry/") + PluginIdentity::kSkirtBellNodeName;
MString SkirtBellCollider::drawRegistrantId = PluginIdentity::kDrawRegistrant;

MObject SkirtBellCollider::attr_bellMatrix;
MObject SkirtBellCollider::attr_leftHipMatrix;
MObject SkirtBellCollider::attr_leftKneeMatrix;
MObject SkirtBellCollider::attr_leftHeelMatrix;
MObject SkirtBellCollider::attr_rightHipMatrix;
MObject SkirtBellCollider::attr_rightKneeMatrix;
MObject SkirtBellCollider::attr_rightHeelMatrix;

MObject SkirtBellCollider::attr_skirtType;
MObject SkirtBellCollider::attr_height;
MObject SkirtBellCollider::attr_ringScale;
MObject SkirtBellCollider::attr_thighRadiusX;
MObject SkirtBellCollider::attr_thighRadiusZ;
MObject SkirtBellCollider::attr_kneeRadiusX;
MObject SkirtBellCollider::attr_kneeRadiusZ;
MObject SkirtBellCollider::attr_calfRadiusX;
MObject SkirtBellCollider::attr_calfRadiusZ;
MObject SkirtBellCollider::attr_ankleRadiusX;
MObject SkirtBellCollider::attr_ankleRadiusZ;
MObject SkirtBellCollider::attr_thighPosition;
MObject SkirtBellCollider::attr_calfPosition;

MObject SkirtBellCollider::attr_bellScale;
MObject SkirtBellCollider::attr_bellSubdivision;
MObject SkirtBellCollider::attr_ringSubdivision;
MObject SkirtBellCollider::attr_falloff;
MObject SkirtBellCollider::attr_tightness;
MObject SkirtBellCollider::attr_smoothness;
MObject SkirtBellCollider::attr_follow;
MObject SkirtBellCollider::attr_bellScaleRamp;
MObject SkirtBellCollider::attr_leftRingAxis;
MObject SkirtBellCollider::attr_rightRingAxis;
MObject SkirtBellCollider::attr_bellAxis;

MObject SkirtBellCollider::attr_outputSurface;
MObject SkirtBellCollider::attr_seams;
MObject SkirtBellCollider::attr_seamEnabled;
MObject SkirtBellCollider::attr_seamMaterialU;
MObject SkirtBellCollider::attr_seamStartHeight;
MObject SkirtBellCollider::attr_panelHems;
MObject SkirtBellCollider::attr_hemUSamples;
MObject SkirtBellCollider::attr_hemHeightSamples;
MObject SkirtBellCollider::attr_followRange;
MObject SkirtBellCollider::attr_referenceMaterialHeight;
MObject SkirtBellCollider::attr_columnMaterialU;
MObject SkirtBellCollider::attr_columnOffsetMatrix;
MObject SkirtBellCollider::attr_outputReferenceHeight;
MObject SkirtBellCollider::attr_outputPatches;
MObject SkirtBellCollider::attr_patchSurface;
MObject SkirtBellCollider::attr_patchMaterialUStart;
MObject SkirtBellCollider::attr_patchMaterialUEnd;
MObject SkirtBellCollider::attr_patchVBreaks;
MObject SkirtBellCollider::attr_patchHemUSamples;
MObject SkirtBellCollider::attr_patchHemHeightSamples;

namespace
{
constexpr double materialTolerance = 1e-9;
constexpr unsigned int maxPanelIndex = 2147483646u;

struct SkirtSeam
{
    unsigned int index;
    double u;
    double height;
};

struct SkirtPanel
{
    unsigned long long id;
    double start;
    double end;
    MDoubleArray hemU;
    MDoubleArray hemHeight;
};

struct SkirtCutSettings
{
    std::vector<SkirtSeam> seams;
    std::vector<SkirtPanel> panels;
    std::vector<unsigned int> hemIndices;
    bool panelMode = false;
    bool ignoredHem = false;
};

struct ColumnOffset
{
    unsigned int index;
    double u;
    MMatrix matrix;
    bool unit;
};

MStatus skirtError(const MObject &node, int condition, const std::string &reason);
MStatus readDoubleArray(const MObject &object, MDoubleArray &values);

bool exactUnit(const MMatrix &matrix)
{
    const MMatrix identity;
    for (unsigned int r = 0; r < 4; ++r)
        for (unsigned int c = 0; c < 4; ++c)
            if (matrix[r][c] != identity[r][c])
                return false;
    return true;
}

MStatus readColumnOffsets(MDataBlock &block, const MObject &node, std::vector<ColumnOffset> &columns)
{
    MStatus status;
    MArrayDataHandle array = block.inputArrayValue(SkirtBellCollider::attr_columnOffsetMatrix, &status);
    if (!status)
        return status;
    const unsigned int count = array.elementCount(&status);
    if (!status)
        return status;
    std::vector<std::pair<unsigned int, unsigned int>> elements;
    for (unsigned int i = 0; i < count; ++i)
    {
        status = array.jumpToArrayElement(i);
        if (!status)
            return status;
        const unsigned int index = array.elementIndex(&status);
        if (!status)
            return status;
        elements.emplace_back(index, i);
    }
    std::sort(elements.begin(), elements.end());
    MDataHandle materialData = block.inputValue(SkirtBellCollider::attr_columnMaterialU, &status);
    if (!status)
        return status;
    MDoubleArray materialU;
    status = readDoubleArray(materialData.data(), materialU);
    if (!status)
        return status;
    for (const auto &entry : elements)
    {
        status = array.jumpToArrayElement(entry.second);
        if (!status)
            return status;
        MDataHandle element = array.inputValue(&status);
        if (!status)
            return status;
        if (entry.first >= materialU.length())
            return skirtError(node, 17, "columnMaterialU[" + std::to_string(entry.first) +
                                             "] is missing for columnOffsetMatrix[" + std::to_string(entry.first) + "].");
        columns.push_back({entry.first, materialU[entry.first], element.asMatrix(), false});
    }
    for (const auto &column : columns)
        if (!std::isfinite(column.u) || column.u < 0.0 || column.u >= 1.0)
            return skirtError(node, 17, "columnMaterialU[" + std::to_string(column.index) +
                                             "] for columnOffsetMatrix[" + std::to_string(column.index) +
                                             "] must be finite in [0,1).");
    for (auto &column : columns)
        column.unit = exactUnit(column.matrix);
    std::vector<size_t> order(columns.size());
    for (size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return columns[a].u == columns[b].u ? columns[a].index < columns[b].index : columns[a].u < columns[b].u;
    });
    for (size_t i = 0; i < order.size(); ++i)
    {
        const size_t next = (i + 1) % order.size();
        const double gap = columns[order[next]].u + (next == 0 ? 1.0 : 0.0) - columns[order[i]].u;
        if (order.size() > 1 && !(gap > materialTolerance))
            return skirtError(node, 18, "columnOffsetMatrix[" + std::to_string(columns[order[i]].index) +
                                             "] and columnOffsetMatrix[" +
                                             std::to_string(columns[order[next]].index) + "] spacing " +
                                             std::to_string(gap) + " must exceed 1e-9.");
    }
    std::sort(columns.begin(), columns.end(), [](const ColumnOffset &a, const ColumnOffset &b) {
        return a.u == b.u ? a.index < b.index : a.u < b.u;
    });
    return MS::kSuccess;
}

MStatus skirtError(const MObject &node, int condition, const std::string &reason)
{
    MGlobal::displayError(MFnDependencyNode(node).name() + ": condition " + MString(std::to_string(condition).c_str()) +
                          ": " + MString(reason.c_str()));
    return MS::kInvalidParameter;
}

MStatus readDoubleArray(const MObject &object, MDoubleArray &values)
{
    if (object.isNull())
    {
        values.clear();
        return MS::kSuccess;
    }
    MStatus status;
    MFnDoubleArrayData data(object, &status);
    if (!status)
        return status;
    values = data.array(&status);
    return status;
}

bool validHem(const MDoubleArray &u, const MDoubleArray &height)
{
    if (u.length() == 0 && height.length() == 0)
        return true;
    if (u.length() != height.length() || u.length() < 2)
        return false;
    for (unsigned int i = 0; i < u.length(); ++i)
        if (!std::isfinite(u[i]) || !std::isfinite(height[i]) || height[i] <= 0.0 || height[i] > 1.0 ||
            (i == 0 && u[i] != 0.0) || (i + 1 == u.length() && u[i] != 1.0) || (i != 0 && u[i] <= u[i - 1]))
            return false;
    return true;
}

bool defaultHem(const MDoubleArray &heights)
{
    for (unsigned int i = 0; i < heights.length(); ++i)
        if (heights[i] != 1.0)
            return false;
    return true;
}

MStatus readCutSettings(MDataBlock &block, const MObject &node, SkirtCutSettings &settings)
{
    MStatus status;
    MArrayDataHandle seams = block.inputArrayValue(SkirtBellCollider::attr_seams, &status);
    if (!status)
        return status;
    const unsigned int seamCount = seams.elementCount(&status);
    if (!status)
        return status;
    for (unsigned int i = 0; i < seamCount; ++i)
    {
        status = seams.jumpToArrayElement(i);
        if (!status)
            return status;
        const unsigned int index = seams.elementIndex(&status);
        if (!status)
            return status;
        MDataHandle element = seams.inputValue(&status);
        if (!status)
            return status;
        if (!element.child(SkirtBellCollider::attr_seamEnabled).asBool())
            continue;
        settings.seams.push_back({index, element.child(SkirtBellCollider::attr_seamMaterialU).asDouble(),
                                  element.child(SkirtBellCollider::attr_seamStartHeight).asDouble()});
    }
    std::sort(settings.seams.begin(), settings.seams.end(),
              [](const SkirtSeam &a, const SkirtSeam &b) { return a.index < b.index; });
    for (const auto &seam : settings.seams)
        if (!std::isfinite(seam.u) || seam.u < 0.0 || seam.u >= 1.0 || !std::isfinite(seam.height) ||
            seam.height < 0.0 || seam.height > 1.0)
            return skirtError(node, 4,
                              "seams[" + std::to_string(seam.index) +
                                  "]: materialU must be finite in [0,1), startHeight "
                                  "finite in [0,1].");
    std::sort(settings.seams.begin(), settings.seams.end(),
              [](const SkirtSeam &a, const SkirtSeam &b) { return a.u < b.u; });
    for (size_t i = 0; i < settings.seams.size(); ++i)
    {
        const double end = i + 1 < settings.seams.size() ? settings.seams[i + 1].u : settings.seams[0].u + 1.0;
        if (end - settings.seams[i].u <= materialTolerance)
            return skirtError(node, 5,
                              "seams[" + std::to_string(settings.seams[i].index) +
                                  "]: periodic materialU spacing must exceed 1e-9.");
        settings.panels.push_back(
            {static_cast<unsigned long long>(settings.seams[i].index) + 1u, settings.seams[i].u, end, {}, {}});
    }
    if (settings.seams.empty())
        settings.panels.push_back({0, 0.0, 1.0, {}, {}});
    settings.panelMode = !settings.seams.empty();
    MArrayDataHandle hems = block.inputArrayValue(SkirtBellCollider::attr_panelHems, &status);
    if (!status)
        return status;
    const unsigned int hemCount = hems.elementCount(&status);
    if (!status)
        return status;
    std::map<unsigned int, unsigned int> physicalIndices;
    for (unsigned int i = 0; i < hemCount; ++i)
    {
        status = hems.jumpToArrayElement(i);
        if (!status)
            return status;
        const unsigned int index = hems.elementIndex(&status);
        if (!status)
            return status;
        physicalIndices.emplace(index, i);
        settings.hemIndices.push_back(index);
    }
    for (const auto &entry : physicalIndices)
    {
        auto panel = std::find_if(settings.panels.begin(), settings.panels.end(), [&](const SkirtPanel &p) {
            return p.id == entry.first && (settings.seams.empty() || p.id != 0);
        });
        if (panel == settings.panels.end())
        {
            settings.ignoredHem = true;
            continue;
        }
        status = hems.jumpToArrayElement(entry.second);
        if (!status)
            return status;
        MDataHandle element = hems.inputValue(&status);
        if (!status)
            return status;
        status = readDoubleArray(element.child(SkirtBellCollider::attr_hemUSamples).data(), panel->hemU);
        if (!status)
            return status;
        status = readDoubleArray(element.child(SkirtBellCollider::attr_hemHeightSamples).data(), panel->hemHeight);
        if (!status)
            return status;
        if (!validHem(panel->hemU, panel->hemHeight))
            return skirtError(node, 6,
                              "panelHems[" + std::to_string(entry.first) +
                                  "]: expected empty arrays or matching finite samples, increasing "
                                  "U from 0 to 1 and 0 < height <= 1.");
        settings.panelMode = settings.panelMode || !defaultHem(panel->hemHeight);
    }
    return MS::kSuccess;
}

struct SkirtHeight
{
    double t;
    int physical = -1;
};

std::vector<SkirtHeight> clusterHeights(const std::vector<double> &physical, SkirtCutSettings &settings,
                                        bool &collapsed)
{
    struct Candidate
    {
        double t;
        int physical;
        int seam;
    };
    std::vector<Candidate> candidates;
    for (size_t i = 0; i < physical.size(); ++i)
        candidates.push_back({physical[i], static_cast<int>(i), -1});
    for (size_t i = 0; i < settings.seams.size(); ++i)
        candidates.push_back({settings.seams[i].height, -1, static_cast<int>(i)});
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &a, const Candidate &b) { return a.t != b.t ? a.t < b.t : a.physical < b.physical; });
    std::vector<SkirtHeight> rows;
    for (size_t first = 0; first < candidates.size();)
    {
        size_t end = first + 1;
        while (end < candidates.size() && candidates[end].t - candidates[first].t <= materialTolerance)
            ++end;
        SkirtHeight row{candidates[first].t, -1};
        for (size_t i = first; i < end; ++i)
            if (candidates[i].physical >= 0)
            {
                if (row.physical >= 0)
                    collapsed = true;
                else
                    row = {candidates[i].t, candidates[i].physical};
            }
        for (size_t i = first; i < end; ++i)
            if (candidates[i].seam >= 0)
                settings.seams[candidates[i].seam].height = row.t;
        rows.push_back(row);
        first = end;
    }
    return rows;
}

double hemEndpoint(const SkirtPanel &panel, bool end)
{
    return panel.hemHeight.length() == 0 ? 1.0 : panel.hemHeight[end ? panel.hemHeight.length() - 1 : 0];
}

MStatus validateHemBoundaries(const SkirtCutSettings &settings, double beta, const MObject &node)
{
    for (const auto &panel : settings.panels)
    {
        if (panel.hemHeight.length() == 0 && !(1.0 > beta))
            return skirtError(node, 8, "panelHems[" + std::to_string(panel.id) + "]: default height must exceed beta.");
        for (unsigned int i = 0; i < panel.hemHeight.length(); ++i)
            if (!(panel.hemHeight[i] > beta))
                return skirtError(node, 8,
                                  "panelHems[" + std::to_string(panel.id) + "]: every height must exceed beta.");
    }
    for (size_t i = 0; i < settings.panels.size(); ++i)
    {
        if (!settings.seams.empty() && settings.seams[i].height != 1.0)
            continue;
        const auto &start = settings.panels[i];
        const auto &end = settings.panels[(i + settings.panels.size() - 1) % settings.panels.size()];
        if (std::abs(hemEndpoint(start, false) - hemEndpoint(end, true)) > materialTolerance)
            return skirtError(node, 8,
                              "panelHems[" + std::to_string(start.id) +
                                  "]: uncut boundary input heights differ by more than 1e-9.");
    }
    return MS::kSuccess;
}

struct SkirtCurve
{
    std::vector<MPoint> points;
    std::vector<double> knots;
};

bool finitePoint(const MPoint &point)
{
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && std::isfinite(point.w);
}

bool validCurve(const SkirtCurve &curve)
{
    if (curve.points.size() < 4 || curve.knots.size() != curve.points.size() + 4)
        return false;
    for (const auto &point : curve.points)
        if (!finitePoint(point))
            return false;
    size_t multiplicity = 0;
    for (size_t i = 0; i < curve.knots.size(); ++i)
    {
        const double knot = curve.knots[i];
        if (!std::isfinite(knot) || (i != 0 && knot < curve.knots[i - 1]))
            return false;
        multiplicity = i != 0 && knot == curve.knots[i - 1] ? multiplicity + 1 : 1;
        if (knot > curve.knots[3] && knot < curve.knots[curve.points.size()] && multiplicity > 3)
            return false;
    }
    return curve.knots[curve.points.size()] > curve.knots[3];
}

bool insertKnot(SkirtCurve &curve, double u)
{
    const int n = static_cast<int>(curve.points.size()) - 1;
    const int k =
        static_cast<int>(std::upper_bound(curve.knots.begin(), curve.knots.end(), u) - curve.knots.begin()) - 1;
    const int multiplicity = static_cast<int>(std::count(curve.knots.begin(), curve.knots.end(), u));
    if (multiplicity >= 3)
        return true;
    if (k < 3 || k > n)
        return false;
    std::vector<MPoint> points(curve.points.size() + 1);
    for (int i = 0; i <= k - 3; ++i)
        points[i] = curve.points[i];
    for (int i = k - multiplicity; i <= n; ++i)
        points[i + 1] = curve.points[i];
    for (int i = k - 2; i <= k - multiplicity; ++i)
    {
        const double denominator = curve.knots[i + 3] - curve.knots[i];
        if (!(denominator > 0.0))
            return false;
        const double alpha = (u - curve.knots[i]) / denominator;
        points[i] = curve.points[i - 1] * (1.0 - alpha) + curve.points[i] * alpha;
    }
    curve.points.swap(points);
    curve.knots.insert(curve.knots.begin() + k + 1, u);
    return true;
}

bool clampAt(SkirtCurve &curve, double u)
{
    while (std::count(curve.knots.begin(), curve.knots.end(), u) < 3)
        if (!insertKnot(curve, u))
            return false;
    return true;
}

bool extractCurve(const SkirtCurve &curve, double start, double end, SkirtCurve &result)
{
    if (!(end > start))
        return false;
    const size_t firstKnot = std::lower_bound(curve.knots.begin(), curve.knots.end(), start) - curve.knots.begin();
    const size_t lastKnot = std::lower_bound(curve.knots.begin(), curve.knots.end(), end) - curve.knots.begin();
    const size_t firstPoint = start == curve.knots[3] && firstKnot == 0 ? 0 : firstKnot - 1;
    const size_t lastPoint = end == curve.knots[curve.points.size()] ? curve.points.size() - 1 : lastKnot - 1;
    if (firstPoint > lastPoint || lastPoint >= curve.points.size())
        return false;
    result.points.assign(curve.points.begin() + firstPoint, curve.points.begin() + lastPoint + 1);
    result.knots.assign(4, start);
    for (double knot : curve.knots)
        if (knot > start && knot < end)
            result.knots.push_back(knot);
    result.knots.insert(result.knots.end(), 4, end);
    return validCurve(result);
}

SkirtCurve extendPeriodic(const SkirtCurve &curve, int count)
{
    SkirtCurve result;
    result.points.reserve(3 * count + 3);
    for (int i = 0; i < 3 * count + 3; ++i)
        result.points.push_back(curve.points[i % count]);
    for (int i = 0; i < 3 * count + 7; ++i)
    {
        const int cycle = i / count;
        const int index = i % count;
        result.knots.push_back(i < static_cast<int>(curve.knots.size()) ? curve.knots[i] : curve.knots[index] + cycle);
    }
    return result;
}

bool restrictPeriodic(const SkirtCurve &periodic, int count, double start, double end, SkirtCurve &result)
{
    SkirtCurve extended = extendPeriodic(periodic, count);
    return clampAt(extended, start) && clampAt(extended, end) && extractCurve(extended, start, end, result);
}

struct SkirtRowComponent
{
    SkirtCurve base;
    std::vector<unsigned int> vertices;
    bool periodic = false;
};

struct SkirtRow
{
    SkirtHeight height;
    int parent = 0;
    BellRowTopology topology;
    std::vector<SkirtRowComponent> components;
    MPointArray base;
    MPointArray points;
    BellDirectField direct;
    MMatrix matrix;
    std::vector<PreparedBellRing> normalRings;
    std::vector<PreparedBellRing> without;
    std::vector<PreparedBellRing> with;
};

struct SkirtSegment
{
    SkirtCurve bottom;
    SkirtCurve top;
    MMatrix matrix;
    MVector x;
    MVector z;
};

struct PreparedColumn
{
    double u;
    bool unit;
    MPoint origin;
    MMatrix delta;
};

bool prepareColumns(const std::vector<ColumnOffset> &columns, const SkirtCurve &waist, int count,
                    const MMatrix &bellFrame, MMatrix &inverseBell, std::vector<PreparedColumn> &prepared)
{
    if (columns.empty() || std::all_of(columns.begin(), columns.end(), [](const ColumnOffset &c) { return c.unit; }))
        return true;
    inverseBell = bellFrame.inverse();
    for (const auto &column : columns)
    {
        PreparedColumn item;
        item.u = column.u;
        item.unit = column.unit;
        if (!item.unit)
        {
            SkirtCurve extended = extendPeriodic(waist, count);
            if (!clampAt(extended, column.u))
                return false;
            const auto first = std::lower_bound(extended.knots.begin(), extended.knots.end(), column.u);
            if (first == extended.knots.begin() || first == extended.knots.end())
                return false;
            const size_t f = static_cast<size_t>(first - extended.knots.begin());
            item.origin = extended.points[f - 1] * inverseBell;
            const MVector derivative = (extended.points[f] - extended.points[f - 1]) *
                                       (3.0 / (extended.knots[f + 3] - column.u));
            const MVector tangentRaw = derivative * inverseBell;
            const MVector axis(0.0, 1.0, 0.0);
            const MVector tangent = (tangentRaw - axis * (tangentRaw * axis)).normal();
            const MVector radialRaw(item.origin.x, 0.0, item.origin.z);
            const MVector radial = (radialRaw - tangent * (radialRaw * tangent)).normal();
            MMatrix frame;
            for (unsigned int c = 0; c < 3; ++c)
            {
                frame[0][c] = tangent[c];
                frame[1][c] = radial[c];
                frame[2][c] = axis[c];
            }
            frame[0][3] = frame[1][3] = frame[2][3] = 0.0;
            frame[3][0] = frame[3][1] = frame[3][2] = 0.0;
            frame[3][3] = 1.0;
            item.delta = frame.transpose() * column.matrix * frame;
        }
        prepared.push_back(item);
    }
    return true;
}

MVector columnDisplacement(const MPoint &point, const PreparedColumn &column)
{
    if (column.unit)
        return MVector(0.0, 0.0, 0.0);
    const MPoint relative(point.x - column.origin.x, point.y - column.origin.y, point.z - column.origin.z, 1.0);
    MPoint transformed = relative * column.delta;
    transformed.x += column.origin.x;
    transformed.y += column.origin.y;
    transformed.z += column.origin.z;
    return transformed - point;
}

bool applyColumnOffsets(SkirtRow &row, const std::vector<PreparedColumn> &columns, const MMatrix &bellFrame,
                        const MMatrix &inverseBell)
{
    if (columns.empty() || std::all_of(columns.begin(), columns.end(), [](const PreparedColumn &c) { return c.unit; }))
        return true;
    for (const auto &component : row.topology.components)
    {
        std::vector<std::pair<double, size_t>> samples;
        for (size_t j = 0; j < columns.size(); ++j)
        {
            double s = columns[j].u;
            if (!component.closed)
            {
                s += std::ceil(component.startU - s);
                if (!(s > component.startU && s < component.endU))
                    continue;
            }
            samples.emplace_back(s, j);
        }
        if (!component.closed && samples.empty())
        {
            for (size_t j = 0; j < columns.size(); ++j)
            {
                const auto seamDistance = [](double a, double b) {
                    const double d = std::abs(a - b);
                    return std::min(d, 1.0 - d);
                };
                if (seamDistance(columns[j].u, component.startU - std::floor(component.startU)) <= materialTolerance)
                    samples.emplace_back(component.startU, j);
                if (seamDistance(columns[j].u, component.endU - std::floor(component.endU)) <= materialTolerance)
                    samples.emplace_back(component.endU, j);
            }
        }
        std::sort(samples.begin(), samples.end());
        for (unsigned int i = 0; i < row.topology.vertices.size(); ++i)
        {
            const auto &vertex = row.topology.vertices[i];
            if (&component != &row.topology.components[vertex.componentId])
                continue;
            MVector delta(0.0, 0.0, 0.0);
            if (samples.size() == 1)
                delta = columnDisplacement(row.base[i] * inverseBell, columns[samples[0].second]);
            else if (samples.size() > 1)
            {
                double s = vertex.materialU;
                if (component.closed)
                    s += std::ceil(samples.front().first - s);
                else if (vertex.side.seamIndex >= 0)
                    s = vertex.side.bank < 0 ? component.startU : component.endU;
                auto upper = std::lower_bound(samples.begin(), samples.end(), std::make_pair(s, size_t(0)));
                if (upper != samples.end() && upper->first == s)
                    delta = columnDisplacement(row.base[i] * inverseBell, columns[upper->second]);
                else if (!component.closed && upper == samples.begin())
                    delta = columnDisplacement(row.base[i] * inverseBell, columns[samples.front().second]);
                else if (!component.closed && upper == samples.end())
                    delta = columnDisplacement(row.base[i] * inverseBell, columns[samples.back().second]);
                else
                {
                    size_t left, right;
                    double lo, hi;
                    if (upper == samples.end())
                    {
                        left = samples.back().second;
                        right = samples.front().second;
                        lo = samples.back().first;
                        hi = samples.front().first + 1.0;
                        if (s < lo)
                            s += 1.0;
                    }
                    else
                    {
                        right = upper->second;
                        left = (upper - 1)->second;
                        lo = (upper - 1)->first;
                        hi = upper->first;
                    }
                    const double lambda = (s - lo) / (hi - lo);
                    const MPoint pointB = row.base[i] * inverseBell;
                    delta = columnDisplacement(pointB, columns[left]) * (1.0 - lambda) +
                            columnDisplacement(pointB, columns[right]) * lambda;
                }
            }
            const MVector deltaG = delta * bellFrame;
            row.base[i] += deltaG;
            const MPoint &p = row.base[i];
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || !std::isfinite(p.w))
                return false;
        }
    }
    return true;
}

size_t panelForMaterial(const SkirtCutSettings &settings, double u, int bank)
{
    for (size_t i = 0; i < settings.panels.size(); ++i)
    {
        const auto &panel = settings.panels[i];
        const double endpoint = bank == 1 ? panel.end : panel.start;
        if (u == endpoint + std::round(u - endpoint))
            return i;
    }
    for (size_t i = 0; i < settings.panels.size(); ++i)
    {
        const auto &panel = settings.panels[i];
        double s = u;
        if (s < panel.start || s > panel.end)
            s += std::ceil(panel.start - s);
        if (s > panel.start && s < panel.end)
            return i;
        if (s == panel.start && bank != 1)
            return i;
        if (s == panel.end && bank == 1)
            return i;
    }
    return 0;
}

double sharedHemEndpoint(const SkirtCutSettings &settings, size_t panelIndex, bool atEnd)
{
    const size_t boundary = atEnd ? (panelIndex + 1) % settings.panels.size() : panelIndex;
    if (!settings.seams.empty() && settings.seams[boundary].height != 1.0)
        return hemEndpoint(settings.panels[panelIndex], atEnd);
    const auto &start = settings.panels[boundary];
    const auto &end = settings.panels[(boundary + settings.panels.size() - 1) % settings.panels.size()];
    return start.id <= end.id ? hemEndpoint(start, false) : hemEndpoint(end, true);
}

double materialHem(const SkirtCutSettings &settings, double u, int bank)
{
    const size_t selected = panelForMaterial(settings, u, bank);
    const auto &panel = settings.panels[selected];
    if (u == panel.start + std::round(u - panel.start) && bank != 1)
        return sharedHemEndpoint(settings, selected, false);
    if (u == panel.end + std::round(u - panel.end))
        return sharedHemEndpoint(settings, selected, true);
    double s = u;
    if (s < panel.start || s > panel.end)
        s += std::ceil(panel.start - s);
    const double x = (s - panel.start) / (panel.end - panel.start);
    const unsigned int sampleCount = panel.hemU.length() == 0 ? 2 : panel.hemU.length();
    const auto sampleU = [&](unsigned int i) {
        return panel.hemU.length() == 0 ? static_cast<double>(i) : panel.hemU[i];
    };
    const auto height = [&](unsigned int i) {
        if (i == 0)
            return sharedHemEndpoint(settings, selected, false);
        if (i + 1 == sampleCount)
            return sharedHemEndpoint(settings, selected, true);
        return panel.hemHeight[i];
    };
    for (unsigned int i = 0; i < sampleCount; ++i)
    {
        if (x == sampleU(i))
            return height(i);
        if (x < sampleU(i) && i != 0)
        {
            const double alpha = (x - sampleU(i - 1)) / (sampleU(i) - sampleU(i - 1));
            const double lower = height(i - 1);
            const double upper = height(i);
            return lower == upper ? lower : lower * (1.0 - alpha) + upper * alpha;
        }
    }
    return height(sampleCount - 1);
}

bool buildRow(SkirtRow &row, const SkirtSegment &segment, const SkirtCutSettings &settings,
              const std::vector<double> &distances, double height, double beta, int count)
{
    std::vector<size_t> cuts;
    for (size_t i = 0; i < settings.seams.size(); ++i)
        if (settings.seams[i].height == 0.0 || row.height.t > settings.seams[i].height)
            cuts.push_back(i);
    struct Entry
    {
        BellRowVertex vertex;
        MPoint point;
        size_t component;
        size_t cv;
    };
    std::vector<Entry> entries;
    unsigned int output = 0;
    const size_t componentCount = cuts.empty() ? 1 : cuts.size();
    for (size_t c = 0; c < componentCount; ++c)
    {
        const bool closed = cuts.empty();
        const size_t begin = closed ? 0 : cuts[c];
        const size_t finish = closed ? 0 : cuts[(c + 1) % cuts.size()];
        const bool mappedHem = row.height.t > beta &&
                               std::any_of(settings.panels.begin(), settings.panels.end(),
                                           [](const SkirtPanel &panel) { return !defaultHem(panel.hemHeight); });
        const bool periodic = closed && !mappedHem;
        const double a = closed ? (mappedHem ? settings.panels.front().start : 0.0) : settings.seams[begin].u;
        const double b = closed ? a + 1.0 : settings.seams[finish].u + (c + 1 == cuts.size() ? 1.0 : 0.0);
        SkirtCurve bottom, top;
        if (periodic)
        {
            bottom = segment.bottom;
            top = segment.top;
        }
        else if (!restrictPeriodic(segment.bottom, count, a, b, bottom) ||
                 !restrictPeriodic(segment.top, count, a, b, top))
            return false;
        if (mappedHem)
            for (const auto &panel : settings.panels)
            {
                const double boundary = panel.start + std::ceil(a - panel.start);
                if (boundary > a && boundary < b &&
                    (!clampAt(bottom, boundary) || !clampAt(top, boundary)))
                    return false;
            }
        if (!validCurve(bottom) || !validCurve(top) || bottom.knots != top.knots)
            return false;
        SkirtRowComponent component;
        component.base = top;
        component.periodic = periodic;
        component.vertices.resize(top.points.size());
        row.topology.components.push_back({a, b, closed});
        const size_t first = periodic ? 1 : 0;
        const size_t end = periodic ? static_cast<size_t>(count + 1) : top.points.size() - (closed ? 1 : 0);
        double previous = -std::numeric_limits<double>::infinity();
        for (size_t i = first; i < end; ++i)
        {
            BellRowVertex vertex;
            if (periodic)
                vertex.materialU = (top.knots[i + 1] + top.knots[i + 2] + top.knots[i + 3]) / 3.0;
            else
            {
                const double width = b - a;
                const double greville =
                    ((top.knots[i + 1] - a) / width + (top.knots[i + 2] - a) / width + (top.knots[i + 3] - a) / width) /
                    3.0;
                vertex.materialU = a + greville * width;
            }
            if (!periodic && top.knots[i + 1] == top.knots[i + 3])
                vertex.materialU = top.knots[i + 1];
            if (!closed && i + 1 == end)
                vertex.materialU = b;
            if (!std::isfinite(vertex.materialU) || vertex.materialU <= previous)
                return false;
            previous = vertex.materialU;
            vertex.componentId = static_cast<int>(c);
            if (!closed && i == 0)
                vertex.side = {static_cast<int>(settings.seams[begin].index & 0x7fffffffu), -1};
            if (!closed && i + 1 == end)
                vertex.side = {static_cast<int>(settings.seams[finish].index & 0x7fffffffu), 1};
            size_t panelIndex = 0;
            if (vertex.side.bank != 0)
                panelIndex = panelForMaterial(settings, vertex.materialU, vertex.side.bank);
            else
                for (size_t p = 0; p < settings.panels.size(); ++p)
                {
                    const auto &panel = settings.panels[p];
                    double u = vertex.materialU;
                    if (u < panel.start || u >= panel.end)
                        u += std::ceil(panel.start - u);
                    if (panel.start <= u && u < panel.end)
                    {
                        panelIndex = p;
                        break;
                    }
                }
            vertex.panelId = static_cast<int>(settings.panels[panelIndex].id & 0x7fffffffu);
            double t = row.height.t;
            const double hem = t > beta ? materialHem(settings, vertex.materialU, vertex.side.bank) : 1.0;
            if (t > beta && hem != 1.0)
                t = beta + (t - beta) * (hem - beta) / (1.0 - beta);
            const double length = distances[row.parent + 1] - distances[row.parent];
            const double alpha = row.height.physical == 0 ? 0.0
                                 : row.height.physical > 0 && t == row.height.t
                                     ? 1.0
                                     : (length > 0.0 ? (height * t - distances[row.parent]) / length : 1.0);
            MPoint point = alpha == 0.0   ? bottom.points[i]
                           : alpha == 1.0 ? top.points[i]
                                          : bottom.points[i] * (1.0 - alpha) + top.points[i] * alpha;
            if (!finitePoint(point))
                return false;
            if (periodic)
            {
                for (size_t j = 0; j < top.points.size(); ++j)
                    if ((j + count - 1) % count == i - 1)
                        vertex.outputDuplicates.push_back(output + static_cast<unsigned int>(j));
            }
            else
            {
                vertex.outputDuplicates.push_back(output + static_cast<unsigned int>(i));
                if (closed && i == 0)
                    vertex.outputDuplicates.push_back(output + static_cast<unsigned int>(top.points.size() - 1));
            }
            entries.push_back({vertex, point, c, i});
        }
        output += static_cast<unsigned int>(top.points.size());
        row.components.push_back(std::move(component));
    }
    std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
        if (a.vertex.materialU != b.vertex.materialU)
            return a.vertex.materialU < b.vertex.materialU;
        if (a.vertex.panelId != b.vertex.panelId)
            return a.vertex.panelId < b.vertex.panelId;
        return a.vertex.side.bank < b.vertex.side.bank;
    });
    row.topology.outputCount = output;
    row.base.setLength(static_cast<unsigned int>(entries.size()));
    for (size_t i = 0; i < entries.size(); ++i)
    {
        row.topology.vertices.push_back(entries[i].vertex);
        row.base[static_cast<unsigned int>(i)] = entries[i].point;
        auto &component = row.components[entries[i].component];
        if (component.periodic)
        {
            for (size_t j = 0; j < component.vertices.size(); ++j)
                if ((j + count - 1) % count == entries[i].cv - 1)
                    component.vertices[j] = static_cast<unsigned int>(i);
        }
        else
        {
            component.vertices[entries[i].cv] = static_cast<unsigned int>(i);
            if (row.topology.components[entries[i].component].closed && entries[i].cv == 0)
                component.vertices.back() = static_cast<unsigned int>(i);
        }
    }
    for (size_t c = 0; c < row.components.size(); ++c)
    {
        std::vector<int> members;
        for (size_t i = 0; i < entries.size(); ++i)
            if (entries[i].component == c)
                members.push_back(static_cast<int>(i));
        for (size_t j = 0; j < members.size(); ++j)
        {
            auto &vertex = row.topology.vertices[members[j]];
            vertex.previous = j ? members[j - 1] : row.topology.components[c].closed ? members.back() : -1;
            vertex.next = j + 1 < members.size()              ? members[j + 1]
                          : row.topology.components[c].closed ? members.front()
                                                              : -1;
        }
    }
    return true;
}

bool patchCurves(const SkirtRow &row, const SkirtCutSettings &settings, int count, std::vector<SkirtCurve> &patches)
{
    patches.resize(settings.panels.size());
    for (size_t c = 0; c < row.components.size(); ++c)
    {
        SkirtCurve solved = row.components[c].base;
        for (size_t i = 0; i < solved.points.size(); ++i)
            solved.points[i] = row.points[row.components[c].vertices[i]];
        const auto &component = row.topology.components[c];
        if (row.components[c].periodic)
            solved = extendPeriodic(solved, count);
        struct Range
        {
            size_t panel;
            double a;
            double b;
        };
        std::vector<Range> ranges;
        for (size_t p = 0; p < settings.panels.size(); ++p)
        {
            double a = settings.panels[p].start;
            double b = settings.panels[p].end;
            double shift = 0.0;
            if (!component.closed)
            {
                shift = std::ceil(component.startU - a);
                a += shift;
                b += shift;
                if (a < component.startU || b > component.endU)
                    continue;
            }
            if (!clampAt(solved, a) || !clampAt(solved, b))
                return false;
            ranges.push_back({p, a, b});
        }
        for (const auto &range : ranges)
        {
            auto &patch = patches[range.panel];
            if (!extractCurve(solved, range.a, range.b, patch))
                return false;
            const double width = range.b - range.a;
            for (double &knot : patch.knots)
                knot = (knot - range.a) / width;
        }
    }
    for (const auto &patch : patches)
        if (!validCurve(patch))
            return false;
    return true;
}

void sharePatchBoundaries(std::vector<std::vector<SkirtCurve>> &surfaces, const std::vector<SkirtHeight> &heights,
                          const SkirtCutSettings &settings)
{
    for (size_t i = 0; i < surfaces.size(); ++i)
    {
        const size_t previous = (i + surfaces.size() - 1) % surfaces.size();
        for (size_t v = 0; v < heights.size(); ++v)
        {
            if (!settings.seams.empty() &&
                (settings.seams[i].height == 0.0 || heights[v].t > settings.seams[i].height))
                continue;
            const MPoint shared = settings.panels[i].id <= settings.panels[previous].id
                                      ? surfaces[i][v].points.front()
                                      : surfaces[previous][v].points.back();
            surfaces[i][v].points.front() = shared;
            surfaces[previous][v].points.back() = shared;
        }
    }
}

bool unifyPatchBasis(std::vector<SkirtCurve> &curves)
{
    std::vector<double> representatives;
    std::map<double, size_t> knots;
    for (auto &curve : curves)
    {
        std::map<double, size_t> local;
        for (double &knot : curve.knots)
        {
            const auto representative = std::find_if(
                representatives.begin(), representatives.end(), [knot](double value) {
                    return knot == value || std::abs(knot - value) <= 1e-12 * std::max(1.0, std::abs(value));
                });
            if (representative == representatives.end())
                representatives.push_back(knot);
            else
                knot = *representative;
            ++local[knot];
        }
        for (const auto &item : local)
            knots[item.first] = std::max(knots[item.first], item.second);
    }
    for (auto &curve : curves)
        for (const auto &item : knots)
            while (static_cast<size_t>(std::count(curve.knots.begin(), curve.knots.end(), item.first)) < item.second)
                if (!insertKnot(curve, item.first))
                    return false;
    return true;
}

MStatus makeSkirtSurface(const std::vector<SkirtCurve> &rows, const std::vector<SkirtHeight> &heights, bool periodic,
                         MObject &data)
{
    MStatus status;
    MPointArray points;
    const unsigned int numU = static_cast<unsigned int>(rows.front().points.size());
    const unsigned int numV = static_cast<unsigned int>(rows.size());
    status = points.setLength(numU * numV);
    if (!status)
        return status;
    for (unsigned int u = 0; u < numU; ++u)
        for (unsigned int v = 0; v < numV; ++v)
        {
            status = points.set(rows[v].points[u], u * numV + v);
            if (!status)
                return status;
        }
    MDoubleArray uKnots, vKnots;
    for (size_t i = 1; i + 1 < rows.front().knots.size(); ++i)
        uKnots.append(rows.front().knots[i]);
    for (const auto &height : heights)
        vKnots.append(height.t);
    MFnNurbsSurfaceData dataFn;
    data = dataFn.create(&status);
    if (!status)
        return status;
    MFnNurbsSurface surface;
    surface.create(points, uKnots, vKnots, 3, 1, periodic ? MFnNurbsSurface::kPeriodic : MFnNurbsSurface::kOpen,
                   MFnNurbsSurface::kOpen, false, data, &status);
    return status;
}
}

void SkirtBellCollider::postConstructor()
{
    MObject thisObj = thisMObject();
    MRampAttribute rampAttr(thisObj, attr_bellScaleRamp);

    MFloatArray positions;
    MFloatArray values;
    MIntArray interpolations;

    positions.append(0.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    positions.append(1.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    rampAttr.addEntries(positions, values, interpolations);
}

MStatus SkirtBellCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;
    MStatus stat;

    // Joint matrices
    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_bellMatrix);

    attr_leftHipMatrix = mAttr.create("leftHipMatrix", "leftHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHipMatrix);

    attr_leftKneeMatrix = mAttr.create("leftKneeMatrix", "leftKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftKneeMatrix);

    attr_leftHeelMatrix = mAttr.create("leftHeelMatrix", "leftHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHeelMatrix);

    attr_rightHipMatrix = mAttr.create("rightHipMatrix", "rightHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHipMatrix);

    attr_rightKneeMatrix = mAttr.create("rightKneeMatrix", "rightKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightKneeMatrix);

    attr_rightHeelMatrix = mAttr.create("rightHeelMatrix", "rightHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHeelMatrix);

    // Skirt Type
    attr_skirtType = eAttr.create("skirtType", "skirtType", 1);
    eAttr.addField("Short", 0);
    eAttr.addField("Long", 1);
    eAttr.setKeyable(true);
    addAttribute(attr_skirtType);

    // Height
    attr_height = nAttr.create("height", "height", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.01f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_height);

    // Ring Scale
    attr_ringScale = nAttr.create("ringScale", "ringScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.5, 1.0, 0.5);
    nAttr.setKeyable(true);
    addAttribute(attr_ringScale);

    // Bell Scale
    attr_bellScale = nAttr.create("bellScale", "bellScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.8, 1.0, 0.8);
    nAttr.setKeyable(true);
    addAttribute(attr_bellScale);

    // Bell Subdivision
    attr_bellSubdivision = nAttr.create("bellSubdivision", "bellSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_bellSubdivision);

    // Ring Subdivision
    attr_ringSubdivision = nAttr.create("ringSubdivision", "ringSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_ringSubdivision);

    // Falloff
    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(-1.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    // Tightness
    attr_tightness = nAttr.create("tightness", "tightness", MFnNumericData::kFloat, 0.5f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_tightness);

    attr_smoothness = nAttr.create("smoothness", "smoothness", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_smoothness);

    attr_follow = nAttr.create("follow", "follow", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_follow);

    // Bell Scale Ramp (Curve)
    attr_bellScaleRamp = MRampAttribute::createCurveRamp("bellScaleRamp", "bellScaleRamp");
    addAttribute(attr_bellScaleRamp);

    // Left Leg Axis
    attr_leftRingAxis = eAttr.create("leftRingAxis", "leftRingAxis", 0);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_leftRingAxis);

    // Right Leg Axis
    attr_rightRingAxis = eAttr.create("rightRingAxis", "rightRingAxis", 3);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_rightRingAxis);

    // Skirt Axis
    attr_bellAxis = eAttr.create("bellAxis", "bellAxis", 1);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_bellAxis);

    // Output NURBS Surface (Visible to users)
    attr_outputSurface = tAttr.create("outputSurface", "outputSurface", MFnData::kNurbsSurface);
    tAttr.setWritable(false);
    tAttr.setStorable(false);
    addAttribute(attr_outputSurface);

    attr_thighRadiusX = nAttr.create("thighRadiusX", "thighRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_thighRadiusX);

    attr_thighRadiusZ = nAttr.create("thighRadiusZ", "thighRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_thighRadiusZ);

    attr_kneeRadiusX = nAttr.create("kneeRadiusX", "kneeRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_kneeRadiusX);

    attr_kneeRadiusZ = nAttr.create("kneeRadiusZ", "kneeRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_kneeRadiusZ);

    attr_calfRadiusX = nAttr.create("calfRadiusX", "calfRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_calfRadiusX);

    attr_calfRadiusZ = nAttr.create("calfRadiusZ", "calfRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_calfRadiusZ);

    attr_ankleRadiusX = nAttr.create("ankleRadiusX", "ankleRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_ankleRadiusX);

    attr_ankleRadiusZ = nAttr.create("ankleRadiusZ", "ankleRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_ankleRadiusZ);

    attr_thighPosition = nAttr.create("thighPosition", "thighPosition", MFnNumericData::kDouble, 0.5);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    addAttribute(attr_thighPosition);

    attr_calfPosition = nAttr.create("calfPosition", "calfPosition", MFnNumericData::kDouble, 0.5);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    addAttribute(attr_calfPosition);

    auto flags = [](MFnAttribute &attribute, bool output, bool multi) {
        attribute.setWritable(!output);
        attribute.setStorable(!output);
        attribute.setKeyable(false);
        attribute.setConnectable(true);
        attribute.setReadable(true);
        attribute.setArray(multi);
        if (multi)
            attribute.setIndexMatters(true);
    };
    auto number = [&](const char *name, const char *shortName, double value, bool output) {
        MObject attribute = nAttr.create(name, shortName, MFnNumericData::kDouble, value);
        flags(nAttr, output, false);
        return attribute;
    };
    auto array = [&](const char *name, const char *shortName, bool output) {
        MFnDoubleArrayData arrayData;
        MObject empty = arrayData.create(MDoubleArray());
        MObject attribute = tAttr.create(name, shortName, MFnData::kDoubleArray, empty);
        flags(tAttr, output, false);
        return attribute;
    };
    MFnCompoundAttribute cAttr;
    attr_seamEnabled = nAttr.create("enabled", "sen", MFnNumericData::kBoolean, true);
    flags(nAttr, false, false);
    attr_seamMaterialU = number("materialU", "smu", 0.0, false);
    attr_seamStartHeight = number("startHeight", "ssh", 0.0, false);
    attr_seams = cAttr.create("seams", "sms");
    flags(cAttr, false, true);
    cAttr.addChild(attr_seamEnabled);
    cAttr.addChild(attr_seamMaterialU);
    cAttr.addChild(attr_seamStartHeight);
    addAttribute(attr_seams);

    attr_hemUSamples = array("uSamples", "phu", false);
    attr_hemHeightSamples = array("heightSamples", "phh", false);
    attr_panelHems = cAttr.create("panelHems", "phs");
    flags(cAttr, false, true);
    cAttr.addChild(attr_hemUSamples);
    cAttr.addChild(attr_hemHeightSamples);
    addAttribute(attr_panelHems);
    attr_followRange = number("followRange", "flr", 0.1875, false);
    addAttribute(attr_followRange);
    attr_referenceMaterialHeight = number("referenceMaterialHeight", "rmh", 1.0, false);
    addAttribute(attr_referenceMaterialHeight);
    attr_columnMaterialU = array("columnMaterialU", "comu", false);
    attr_columnOffsetMatrix = mAttr.create("columnOffsetMatrix", "comx", MFnMatrixAttribute::kDouble);
    flags(mAttr, false, true);
    mAttr.setArray(true);
    mAttr.setIndexMatters(true);
    addAttribute(attr_columnOffsetMatrix);
    addAttribute(attr_columnMaterialU);
    attr_outputReferenceHeight = number("outputReferenceHeight", "orh", 0.0, true);
    addAttribute(attr_outputReferenceHeight);

    attr_patchSurface = tAttr.create("surface", "opsf", MFnData::kNurbsSurface);
    flags(tAttr, true, false);
    attr_patchMaterialUStart = number("materialUStart", "opua", 0.0, true);
    attr_patchMaterialUEnd = number("materialUEnd", "opub", 1.0, true);
    attr_patchVBreaks = array("vBreaks", "opvb", true);
    attr_patchHemUSamples = array("hemUSamples", "ophu", true);
    attr_patchHemHeightSamples = array("hemHeightSamples", "ophh", true);
    attr_outputPatches = cAttr.create("outputPatches", "ops");
    flags(cAttr, true, true);
    cAttr.setUsesArrayDataBuilder(true);
    for (const MObject &child : {attr_patchSurface, attr_patchMaterialUStart, attr_patchMaterialUEnd, attr_patchVBreaks,
                                 attr_patchHemUSamples, attr_patchHemHeightSamples})
        cAttr.addChild(child);
    addAttribute(attr_outputPatches);

    // Set up attribute affects relationships
    const MObject affects[] = {attr_thighRadiusX,
                               attr_thighRadiusZ,
                               attr_kneeRadiusX,
                               attr_kneeRadiusZ,
                               attr_calfRadiusX,
                               attr_calfRadiusZ,
                               attr_ankleRadiusX,
                               attr_ankleRadiusZ,
                               attr_thighPosition,
                               attr_calfPosition,
                               attr_bellMatrix,
                               attr_leftHipMatrix,
                               attr_leftKneeMatrix,
                               attr_leftHeelMatrix,
                               attr_rightHipMatrix,
                               attr_rightKneeMatrix,
                               attr_rightHeelMatrix,
                               attr_skirtType,
                               attr_height,
                               attr_ringScale,
                               attr_bellScale,
                               attr_bellSubdivision,
                               attr_falloff,
                               attr_tightness,
                               attr_smoothness,
                               attr_follow,
                               attr_bellScaleRamp,
                               attr_leftRingAxis,
                               attr_rightRingAxis,
                               attr_bellAxis,
                               attr_seams,
                               attr_seamEnabled,
                               attr_seamMaterialU,
                               attr_seamStartHeight,
                               attr_panelHems,
                               attr_hemUSamples,
                               attr_hemHeightSamples,
                               attr_followRange,
                               attr_referenceMaterialHeight,
                               attr_columnOffsetMatrix,
                               attr_columnMaterialU,
                               MPxNode::state};
    for (const MObject &attr : affects)
    {
        for (const MObject &output : {attr_outputSurface, attr_outputReferenceHeight, attr_outputPatches,
                                      attr_patchSurface, attr_patchMaterialUStart, attr_patchMaterialUEnd,
                                      attr_patchVBreaks, attr_patchHemUSamples, attr_patchHemHeightSamples})
            attributeAffects(attr, output);
    }

    return MS::kSuccess;
}

// Maya propagates dirtiness to the outputPatches parent plug but not to the existing
// elements and their children; consumers connected to a child would otherwise keep
// reading the previous data object after it has been released by the next evaluation.
MStatus SkirtBellCollider::setDependentsDirty(const MPlug &plugBeingDirtied, MPlugArray &affectedPlugs)
{
    const unsigned int before = affectedPlugs.length();
    MStatus stat = MPxLocatorNode::setDependentsDirty(plugBeingDirtied, affectedPlugs);
    const MObject attribute = plugBeingDirtied.attribute();
    const bool output = attribute == attr_outputSurface || attribute == attr_outputPatches ||
                        attribute == attr_outputReferenceHeight || attribute == attr_patchSurface ||
                        attribute == attr_patchMaterialUStart || attribute == attr_patchMaterialUEnd ||
                        attribute == attr_patchVBreaks || attribute == attr_patchHemUSamples ||
                        attribute == attr_patchHemHeightSamples;
    if (output || (affectedPlugs.length() == before && attribute != MPxNode::state))
        return stat;
    MPlug patches(thisMObject(), attr_outputPatches);
    affectedPlugs.append(patches);
    const unsigned int count = patches.numElements();
    for (unsigned int i = 0; i < count; ++i)
    {
        MPlug element = patches.elementByPhysicalIndex(i);
        affectedPlugs.append(element);
        for (unsigned int c = 0; c < element.numChildren(); ++c)
            affectedPlugs.append(element.child(c));
    }
    return stat;
}

MStatus SkirtBellCollider::compute(const MPlug &plug, MDataBlock &dataBlock)
{
    const MObject attribute = plug.attribute();
    const bool referenceRequest = attribute == attr_outputReferenceHeight;
    const bool surfaceRequest = attribute == attr_outputSurface;
    const bool patchRequest = attribute == attr_outputPatches || attribute == attr_patchSurface ||
                              attribute == attr_patchMaterialUStart || attribute == attr_patchMaterialUEnd ||
                              attribute == attr_patchVBreaks || attribute == attr_patchHemUSamples ||
                              attribute == attr_patchHemHeightSamples;
    if (!referenceRequest && !surfaceRequest && !patchRequest)
        return MS::kUnknownParameter;
    std::lock_guard<std::mutex> lock(computeMutex_);
    MStatus stat;
    if (referenceRequest)
    {
        const double value = dataBlock.inputValue(attr_referenceMaterialHeight, &stat).asDouble();
        if (!stat)
            return stat;
        if (!std::isfinite(value) || value <= 0.0)
            return skirtError(thisMObject(), 10, "referenceMaterialHeight must be finite and positive.");
        MDataHandle output = dataBlock.outputValue(attr_outputReferenceHeight, &stat);
        if (!stat)
            return stat;
        output.setDouble(value);
        return dataBlock.setClean(plug);
    }
    // Get input values
    const MMatrix inputBellMatrix = dataBlock.inputValue(attr_bellMatrix).asMatrix();
    const MMatrix leftHipMatrix = dataBlock.inputValue(attr_leftHipMatrix).asMatrix();
    const MMatrix leftKneeMatrix = dataBlock.inputValue(attr_leftKneeMatrix).asMatrix();
    const MMatrix leftHeelMatrix = dataBlock.inputValue(attr_leftHeelMatrix).asMatrix();
    const MMatrix rightHipMatrix = dataBlock.inputValue(attr_rightHipMatrix).asMatrix();
    const MMatrix rightKneeMatrix = dataBlock.inputValue(attr_rightKneeMatrix).asMatrix();
    const MMatrix rightHeelMatrix = dataBlock.inputValue(attr_rightHeelMatrix).asMatrix();

    const short skirtType = dataBlock.inputValue(attr_skirtType).asShort();
    const float height = dataBlock.inputValue(attr_height).asFloat();
    const MVector ringScale = dataBlock.inputValue(attr_ringScale).asVector();
    const MVector bellScale = dataBlock.inputValue(attr_bellScale).asVector();
    const int bellSubdivision = dataBlock.inputValue(attr_bellSubdivision).asInt();
    const float falloff = dataBlock.inputValue(attr_falloff).asFloat();
    const bool hasNoEffect = dataBlock.inputValue(MPxNode::state).asShort() == 1;
    const float tightness = dataBlock.inputValue(attr_tightness).asFloat();
    const float smoothness = dataBlock.inputValue(attr_smoothness).asFloat();
    const float follow = dataBlock.inputValue(attr_follow).asFloat();
    const short leftRingAxis = dataBlock.inputValue(attr_leftRingAxis).asShort();
    const short rightRingAxis = dataBlock.inputValue(attr_rightRingAxis).asShort();
    const short bellAxis = dataBlock.inputValue(attr_bellAxis).asShort();

    const double followRange = dataBlock.inputValue(attr_followRange).asDouble();
    if (!ColliderInput::validSubdivision(bellSubdivision))
        return skirtError(thisMObject(), 1, "bellSubdivision must be between 3 and 4096.");
    if (!ColliderInput::validSkirtType(skirtType))
        return skirtError(thisMObject(), 2, "skirtType must be 0 or 1.");
    if (!ColliderInput::validAxis(leftRingAxis))
        return skirtError(thisMObject(), 3, "leftRingAxis must be between 0 and 5.");
    if (!ColliderInput::validAxis(rightRingAxis))
        return skirtError(thisMObject(), 3, "rightRingAxis must be between 0 and 5.");
    if (!ColliderInput::validAxis(bellAxis))
        return skirtError(thisMObject(), 3, "bellAxis must be between 0 and 5.");
    SkirtCutSettings settings;
    stat = readCutSettings(dataBlock, thisMObject(), settings);
    if (!stat)
        return stat;
    if (settings.ignoredHem && !warnedIgnoredHem_.exchange(true))
        MGlobal::displayWarning(MFnDependencyNode(thisMObject()).name() +
                                ": panelHems contains an ignored logical index "
                                "without an active panel.");

    // Ramp attribute
    MRampAttribute rampAttr(thisMObject(), attr_bellScaleRamp);

    // Midpoints
    const MPoint LH = taxis(leftHipMatrix);
    const MPoint RH = taxis(rightHipMatrix);

    const MPoint H = (LH + RH) * 0.5;

    const MPoint W = taxis(inputBellMatrix);

    const SkirtRingFrames ringFrames(leftHipMatrix, leftKneeMatrix, leftHeelMatrix, rightHipMatrix, rightKneeMatrix,
                                     rightHeelMatrix, ringScale, leftRingAxis, rightRingAxis, skirtType == 1);
    const double L_thigh = ringFrames.thighLength;
    const double L_calf = ringFrames.calfLength;

    const SkirtLegProfile profile(
        dataBlock.inputValue(attr_thighRadiusX).asDouble(), dataBlock.inputValue(attr_thighRadiusZ).asDouble(),
        dataBlock.inputValue(attr_kneeRadiusX).asDouble(), dataBlock.inputValue(attr_kneeRadiusZ).asDouble(),
        dataBlock.inputValue(attr_calfRadiusX).asDouble(), dataBlock.inputValue(attr_calfRadiusZ).asDouble(),
        dataBlock.inputValue(attr_ankleRadiusX).asDouble(), dataBlock.inputValue(attr_ankleRadiusZ).asDouble(),
        dataBlock.inputValue(attr_thighPosition).asDouble(), dataBlock.inputValue(attr_calfPosition).asDouble(),
        L_thigh, L_calf);

    const double d_hip = (H - W).length();
    const double d_knee = d_hip + L_thigh;
    const double d_heel = d_knee + L_calf;
    const double d_mid = (d_hip + d_knee) * 0.5;

    double h_param = (double)height;
    if (h_param < 0.0)
        h_param = 0.0;
    if (h_param > 1.0)
        h_param = 1.0;

    double h_val = 0.0;
    if (skirtType == 1) // Long skirt: bottom interpolates from Knee to Heel
    {
        h_val = d_knee + (d_heel - d_knee) * h_param;
    }
    else // Short skirt: bottom interpolates from Hip to Knee
    {
        h_val = d_hip + (d_knee - d_hip) * h_param;
    }

    const double raw_defaultHeight = (skirtType == 1) ? d_heel : d_knee;
    const double defaultHeight = raw_defaultHeight < 1e-5 ? 1.0 : raw_defaultHeight;
    const double s = h_val / defaultHeight;

    // Get direction vector based on the selected bell matrix axis
    const MVector raw_dir_vector = getAxis(inputBellMatrix, bellAxis);
    const MVector dir_vector = raw_dir_vector.length() < 1e-4 ? MVector(0, -1, 0) : raw_dir_vector.normal();

    // Start point of the skirt aligned rigidly along the chosen waist axis
    const MPoint P_start = W;

    const int N = (skirtType == 0) ? 2 : 3;

    // Setup level distances along the skirt axis
    vector<double> levelDistances;
    levelDistances.push_back(0.0);
    if (skirtType == 1) // Long skirt: height controls the heel level only
    {
        levelDistances.push_back(d_mid);
        levelDistances.push_back(d_knee);
        levelDistances.push_back(h_val);
    }
    else // Short skirt: height scales all levels
    {
        levelDistances.push_back(d_mid * s);
        levelDistances.push_back(d_knee * s);
    }

    bool invalidMaterial = !std::isfinite(h_val) || h_val <= 0.0;
    bool invalidSegments = false;
    std::vector<double> physicalHeights;
    for (size_t i = 0; i < levelDistances.size(); ++i)
    {
        const double t = i == 0 ? 0.0 : i + 1 == levelDistances.size() ? 1.0 : levelDistances[i] / h_val;
        if (!std::isfinite(t))
            invalidMaterial = true;
        physicalHeights.push_back(std::isfinite(t) ? t : static_cast<double>(i) / N);
        if (!std::isfinite(levelDistances[i]))
            invalidMaterial = true;
        if (i != 0 && !(levelDistances[i] > levelDistances[i - 1]))
            invalidSegments = true;
    }
    bool collapsed = false;
    std::vector<SkirtHeight> heights = clusterHeights(physicalHeights, settings, collapsed);
    double beta = physicalHeights[physicalHeights.size() - 2];
    MDoubleArray vBreaks;
    bool internalStart = false;
    for (const auto &seam : settings.seams)
        if (seam.height > 0.0 && seam.height < 1.0)
        {
            beta = internalStart ? std::max(beta, seam.height) : seam.height;
            internalStart = true;
            vBreaks.append(seam.height);
        }
    if (!internalStart && beta > 0.0 && beta < 1.0)
        vBreaks.append(beta);
    std::vector<double> breaks;
    for (unsigned int i = 0; i < vBreaks.length(); ++i)
        breaks.push_back(vBreaks[i]);
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
    vBreaks.clear();
    for (double value : breaks)
        vBreaks.append(value);
    stat = validateHemBoundaries(settings, beta, thisMObject());
    if (!stat)
        return stat;
    std::vector<ColumnOffset> columnOffsets;
    stat = readColumnOffsets(dataBlock, thisMObject(), columnOffsets);
    if (!stat)
        return stat;
    if (invalidMaterial || !std::isfinite(follow) || !std::isfinite(smoothness) || !std::isfinite(falloff) ||
        !std::isfinite(tightness))
        return skirtError(thisMObject(), 9,
                          "height and row construction must produce finite "
                          "material coordinates and solver inputs.");
    const double h_safe = std::max(h_val, 1e-5);
    const MVector dir_y = dir_vector;

    const MVector raw_waist_X = xaxis(inputBellMatrix);
    const MVector waist_X = raw_waist_X.length() < 1e-4 ? MVector(1, 0, 0) : raw_waist_X.normal();

    const MVector raw_X = (waist_X - (waist_X * dir_y) * dir_y);
    const MVector X = [&]() {
        if (raw_X.length() < 1e-4)
        {
            const MVector raw_waist_Z = zaxis(inputBellMatrix);
            const MVector waist_Z = raw_waist_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_waist_Z.normal();
            const MVector crossed = dir_y ^ waist_Z;
            if (crossed.length() < 1e-4)
            {
                return MVector(1, 0, 0);
            }
            return crossed.normal();
        }
        return raw_X.normal();
    }();

    const MVector raw_Z = X ^ dir_y;
    const MVector Z = raw_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_Z.normal();
    MMatrix bellFrame;
    const bool needsColumnFrame = std::any_of(columnOffsets.begin(), columnOffsets.end(),
                                               [](const ColumnOffset &column) { return !column.unit; });
    if (needsColumnFrame)
    {
        const double bellFrameValues[4][4] = {{X.x, X.y, X.z, 0.0},
                                              {dir_y.x, dir_y.y, dir_y.z, 0.0},
                                              {Z.x, Z.y, Z.z, 0.0},
                                              {P_start.x, P_start.y, P_start.z, 1.0}};
        bellFrame = MMatrix(bellFrameValues);
    }

    const BellCircleTable circle(bellSubdivision);

    std::vector<SkirtSegment> segments;
    for (int i = 0; i < N; ++i)
    {
        const double distance = levelDistances[i + 1] - levelDistances[i];
        const double safeDist = std::max(distance, 1e-4);
        const MPoint origin = P_start + dir_vector * levelDistances[i];
        float bottomScale = 1.0f, topScale = 1.0f;
        rampAttr.getValueAtPosition(static_cast<float>(levelDistances[i] / h_safe), bottomScale, &stat);
        if (!stat)
            return stat;
        rampAttr.getValueAtPosition(static_cast<float>(levelDistances[i + 1] / h_safe), topScale, &stat);
        if (!stat)
            return stat;
        topScale = std::max(topScale, 1e-5f);
        SkirtSegment segment;
        segment.x = X * (bellScale.x * topScale);
        segment.z = Z * (bellScale.z * topScale);
        const MVector y = dir_y * (safeDist * bellScale.y);
        const double matrix[4][4] = {{segment.x.x, segment.x.y, segment.x.z, 0.0},
                                     {y.x, y.y, y.z, 0.0},
                                     {segment.z.x, segment.z.y, segment.z.z, 0.0},
                                     {origin.x, origin.y, origin.z, 1.0}};
        segment.matrix = MMatrix(matrix);
        MPointArray points =
            BellColliderSolver::makeBellPoints(segment.matrix, 1, circle, 1, bottomScale / topScale, 1);
        BellColliderSolver::roundMeshPoints(points);
        for (int u = 0; u < bellSubdivision + 3; ++u)
        {
            const int original = (bellSubdivision + 2 - u) % bellSubdivision;
            segment.bottom.points.push_back(points[1 + original]);
            segment.top.points.push_back(points[bellSubdivision + 1 + original]);
        }
        for (int k = 0; k < bellSubdivision + 7; ++k)
        {
            segment.bottom.knots.push_back(static_cast<double>(k - 3));
            segment.top.knots.push_back(static_cast<double>(k - 3));
        }
        if (!validCurve(segment.bottom) || !validCurve(segment.top))
            return skirtError(thisMObject(), 9,
                              "bellMatrix, bellScale or bellScaleRamp produced a "
                              "nonfinite reference curve.");
        segments.push_back(std::move(segment));
    }
    MObject referenceSurface;
    stat = makeSkirtSurface({segments[0].bottom, segments[0].top}, {{0.0, 0}, {1.0, 1}}, true, referenceSurface);
    if (!stat)
        return stat;
    MFnNurbsSurface reference(referenceSurface, &stat);
    if (!stat)
        return stat;
    double umin, umax, vmin, vmax;
    stat = reference.getKnotDomain(umin, umax, vmin, vmax);
    if (!stat)
        return stat;
    if (!std::isfinite(umin) || !std::isfinite(umax) || !(umax > umin))
        return skirtError(thisMObject(), 9, "reference surface has an invalid material U domain.");
    for (auto &segment : segments)
    {
        for (double &knot : segment.bottom.knots)
            knot = (knot - umin) / (umax - umin);
        segment.top.knots = segment.bottom.knots;
    }
    std::vector<PreparedColumn> preparedColumns;
    MMatrix inverseBell;
    if (!prepareColumns(columnOffsets, segments[0].bottom, bellSubdivision, bellFrame, inverseBell, preparedColumns))
        return skirtError(thisMObject(), 9, "column material could not be clamped on the waist curve.");
    std::vector<SkirtRow> rows;
    for (const auto &heightRow : heights)
    {
        SkirtRow row;
        row.height = heightRow;
        if (heightRow.physical >= 0)
            row.parent = std::max(0, heightRow.physical - 1);
        else
        {
            const auto upper = std::upper_bound(physicalHeights.begin(), physicalHeights.end(), heightRow.t);
            row.parent = std::max(0, std::min(N - 1, static_cast<int>(upper - physicalHeights.begin()) - 1));
        }
        const auto &segment = segments[row.parent];
        row.matrix = segment.matrix;
        if (heightRow.physical < 0)
        {
            const MPoint origin = P_start + dir_vector * (h_val * heightRow.t);
            const double matrix[4][4] = {{segment.x.x, segment.x.y, segment.x.z, 0.0},
                                         {dir_vector.x, dir_vector.y, dir_vector.z, 0.0},
                                         {segment.z.x, segment.z.y, segment.z.z, 0.0},
                                         {origin.x, origin.y, origin.z, 1.0}};
            row.matrix = MMatrix(matrix);
        }
        if (!buildRow(row, segment, settings, levelDistances, h_val, beta, bellSubdivision))
            return skirtError(thisMObject(), 9,
                              "seams or material coordinates cannot construct a "
                              "finite cubic row with distinct Greville samples.");
        if (!applyColumnOffsets(row, preparedColumns, bellFrame, inverseBell))
            return skirtError(thisMObject(), 9, "column offsets produced a nonfinite row base.");
        const double distance = h_val * row.height.t;
        const double level = profile.levelParameter(distance, d_hip);
        auto prepare = [&](const MMatrix &base, SkirtLegProfile::Ring kind) {
            const auto radius = profile.forRing(level, kind, ringScale.y);
            MMatrix matrix = base;
            for (unsigned int c = 0; c < 4; ++c)
            {
                matrix[0][c] *= radius.x;
                matrix[2][c] *= radius.z;
            }
            PreparedBellRing ring(matrix);
            BellColliderSolver::prepareDistalEnd(
                ring, matrix, kind == SkirtLegProfile::Ring::Knee || kind == SkirtLegProfile::Ring::Extended);
            return ring;
        };
        row.normalRings.push_back(prepare(ringFrames.leftKnee, SkirtLegProfile::Ring::Knee));
        row.normalRings.push_back(prepare(ringFrames.rightKnee, SkirtLegProfile::Ring::Knee));
        if (skirtType == 1)
        {
            row.normalRings.push_back(prepare(ringFrames.leftHeel, SkirtLegProfile::Ring::Heel));
            row.normalRings.push_back(prepare(ringFrames.rightHeel, SkirtLegProfile::Ring::Heel));
            row.without.assign(row.normalRings.begin() + 2, row.normalRings.end());
            row.with = row.without;
            row.with.push_back(prepare(ringFrames.leftExtended, SkirtLegProfile::Ring::Extended));
            row.with.push_back(prepare(ringFrames.rightExtended, SkirtLegProfile::Ring::Extended));
        }
        for (const auto *rings : {&row.normalRings, &row.with})
            for (const auto &ring : *rings)
                for (unsigned int r = 0; r < 4; ++r)
                    for (unsigned int c = 0; c < 4; ++c)
                        if (!std::isfinite(ring.inverse[r][c]))
                            return skirtError(thisMObject(), 9,
                                              "leg matrices or radius inputs produced a "
                                              "nonfinite row ring.");
        rows.push_back(std::move(row));
    }
    const bool validFollowRange = std::isfinite(followRange) && followRange >= 0.0;
    for (auto &row : rows)
    {
        const double distance = h_val * row.height.t;
        const auto &normalRings = row.normalRings;
        const auto &without = row.without;
        const auto &with = row.with;
        row.points = row.base;
        if (hasNoEffect)
            continue;
        if (row.height.physical == 0)
        {
            const int count = static_cast<int>(row.points.length());
            for (const auto &ring : normalRings)
                BellColliderSolver::relaxTowardRingBoundary(
                    row.points, ring, 1.0, 0, count, true, std::vector<MVector>(count, MVector(0, 0, 0)), {});
            if (smoothness > 0.0f)
            {
                std::vector<MVector> displacements(count);
                for (int i = 0; i < count; ++i)
                    displacements[i] = row.points[i] - row.base[i];
                stat = BellColliderSolver::smoothDisplacements(displacements, smoothness, row.topology);
                if (!stat)
                    return skirtError(thisMObject(), 9, "waist row topology failed smoothing.");
                for (int i = 0; i < count; ++i)
                    row.points[i] = row.base[i] + displacements[i];
                for (const auto &ring : normalRings)
                    BellColliderSolver::relaxTowardRingBoundary(
                        row.points, ring, 1.0, 0, count, true,
                        BellColliderSolver::rowDirections(row.points, row.base), {});
            }
            continue;
        }
        const bool belowKnee = skirtType == 1 && row.height.t > physicalHeights[2];
        auto solve = [&](const std::vector<PreparedBellRing> &rings, BellRowOutputs &output) {
            BellRowInputs input;
            input.bellMatrix = row.matrix;
            input.rings = rings;
            input.falloff = falloff;
            input.collision = 1.0f;
            input.capAtRingOrigin = true;
            input.smoothness = smoothness;
            input.followGain = follow * (distance / h_safe) * 0.5;
            input.followRange = validFollowRange ? followRange : 0.0;
            input.contactBlendWidth = 1.0 / bellSubdivision;
            input.physicalLevel = row.parent;
            return BellColliderSolver::solveRow(input, row.base, row.topology, output);
        };
        BellRowOutputs solution;
        if (belowKnee && tightness > 0.0f && tightness < 1.0f)
        {
            BellRowOutputs other;
            stat = solve(with, solution);
            if (stat)
                stat = solve(without, other);
            if (!stat)
                return skirtError(thisMObject(), 9, "row solve failed for a knee-to-hem component.");
            for (unsigned int i = 0; i < solution.points.length(); ++i)
            {
                solution.points[i] = solution.points[i] * (1.0 - tightness) + other.points[i] * tightness;
                solution.directField.values[i] =
                    solution.directField.values[i] * (1.0 - tightness) + other.directField.values[i] * tightness;
            }
        }
        else
        {
            stat = solve(belowKnee ? (tightness >= 1.0f ? without : with) : normalRings, solution);
            if (!stat)
                return skirtError(thisMObject(), 9, "row solve failed for a material component.");
        }
        row.points = solution.points;
        row.direct = solution.directField;
        BellColliderSolver::roundMeshPoints(row.points);
        bool receives = false;
        if (follow > 0.0f)
        {
            const double rowTightness = belowKnee ? tightness : 0.0;
            std::vector<MVector> propagation(row.points.length(), MVector(0, 0, 0));
            for (const auto &source : rows)
            {
                const int k = source.height.physical;
                if (k <= 0 || k >= N || source.height.t >= row.height.t)
                    continue;
                const double denominator = physicalHeights[k + 1] - physicalHeights[k];
                if (!(denominator > 0.0))
                    continue;
                const double coefficient =
                    std::max(0.0, std::min(1.0, (row.height.t - physicalHeights[k]) / denominator));
                if (coefficient == 0.0)
                    continue;
                BellRowTransfer transfer;
                stat = BellColliderSolver::transferDirectField(source.direct, row.topology, transfer);
                if (!stat)
                    return skirtError(thisMObject(), 9, "material correspondence failed during row transfer.");
                for (unsigned int i = 0; i < row.points.length(); ++i)
                    propagation[i] += transfer.values[i] * coefficient;
                receives = true;
            }
            if (receives)
                for (unsigned int i = 0; i < row.points.length(); ++i)
                    row.points[i] += propagation[i] * (follow * (1.0 - rowTightness));
        }
        if (receives)
        {
            // The weights of a ring depend only on the ring and the base row, so
            // the weights of the primary solve serve every ring list: with[i] and
            // without[i] are the same ring for i below without.size().
            const auto &weights = solution.componentWeights;
            const auto weightsFor = [&](size_t index) -> const std::vector<double> & {
                static const std::vector<double> none;
                return index < weights.size() ? weights[index] : none;
            };
            if (belowKnee)
            {
                for (size_t i = 0; i < without.size(); ++i)
                    BellColliderSolver::relaxRowFaded(row.points, without[i], 1.0, true, row.topology, weightsFor(i),
                                                      BellColliderSolver::rowDirections(row.points, row.base));
                if (tightness < 1.0f)
                    for (size_t i = 2; i < with.size(); ++i)
                        BellColliderSolver::relaxRowFaded(row.points, with[i], 1.0 - tightness, true, row.topology,
                                                          weightsFor(i),
                                                          BellColliderSolver::rowDirections(row.points, row.base));
            }
            else
                for (size_t i = 0; i < normalRings.size(); ++i)
                    BellColliderSolver::relaxRowFaded(row.points, normalRings[i], 1.0, true, row.topology,
                                                      weightsFor(i),
                                                      BellColliderSolver::rowDirections(row.points, row.base));
        }
    }

    if (settings.panelMode && surfaceRequest)
        return skirtError(thisMObject(), 11,
                          "outputSurface is unavailable with active seams or a "
                          "nondefault panel hem; use outputPatches.");

    std::vector<std::vector<SkirtCurve>> surfaces(settings.panels.size());
    std::vector<SkirtCurve> periodicRows;
    for (const auto &row : rows)
    {
        std::vector<SkirtCurve> curves;
        if (!patchCurves(row, settings, bellSubdivision, curves))
            return skirtError(thisMObject(), 12, "outputPatches cannot represent finite cubic restrictions.");
        for (size_t i = 0; i < curves.size(); ++i)
            surfaces[i].push_back(std::move(curves[i]));
        if (!settings.panelMode)
        {
            SkirtCurve curve = row.components[0].base;
            for (size_t i = 0; i < curve.points.size(); ++i)
                curve.points[i] = row.points[row.components[0].vertices[i]];
            periodicRows.push_back(std::move(curve));
        }
    }
    for (auto &surface : surfaces)
        if (!unifyPatchBasis(surface))
            return skirtError(thisMObject(), 12, "outputPatches cannot form a common U basis.");
    sharePatchBoundaries(surfaces, heights, settings);
    if (!settings.panelMode)
        surfaces.push_back(std::move(periodicRows));
    for (const auto &surface : surfaces)
        for (const auto &curve : surface)
            if (!validCurve(curve))
                return skirtError(thisMObject(), 12, "output CVs and knots must be finite cubic data.");
    for (size_t i = 0; i < heights.size(); ++i)
        if (!std::isfinite(heights[i].t) || (i && heights[i].t <= heights[i - 1].t))
            return skirtError(thisMObject(), 12, "output V knots must be finite and strictly increasing.");
    if (!validFollowRange)
        return skirtError(thisMObject(), 13, "followRange must be finite and nonnegative.");
    if (invalidSegments || collapsed)
        return skirtError(thisMObject(), 14,
                          "physical levels contain a zero-length segment or share "
                          "a height cluster.");
    for (const auto &seam : settings.seams)
        if (seam.index > maxPanelIndex)
            return skirtError(thisMObject(), 15,
                              "seams[" + std::to_string(seam.index) + "]: logical index exceeds 2147483646.");
    for (unsigned int index : settings.hemIndices)
        if (index > maxPanelIndex)
            return skirtError(thisMObject(), 15,
                              "panelHems[" + std::to_string(index) + "]: logical index exceeds 2147483646.");
    for (const auto &surface : surfaces)
        if (surface.size() > static_cast<size_t>(std::numeric_limits<int>::max()) / surface.front().points.size())
            return skirtError(thisMObject(), 16, "output CV count exceeds the supported integer index range.");

    std::vector<MObject> surfaceData(surfaces.size());
    for (size_t i = 0; i < surfaces.size(); ++i)
    {
        stat = makeSkirtSurface(surfaces[i], heights, i == settings.panels.size(), surfaceData[i]);
        if (!stat)
            return stat;
    }
    MArrayDataBuilder builder(&dataBlock, attr_outputPatches,
                              static_cast<unsigned int>(settings.panels.size()), &stat);
    if (!stat)
        return stat;
    MFnDoubleArrayData arrayData;
    const MDoubleArray emptyHem;
    for (size_t i = 0; i < settings.panels.size(); ++i)
    {
        const auto &panel = settings.panels[i];
        // Each element owns its data object; sharing one across elements frees it twice.
        MObject breakData = arrayData.create(vBreaks, &stat);
        if (!stat)
            return stat;
        MObject hemUData = arrayData.create(settings.panelMode ? panel.hemU : emptyHem, &stat);
        if (!stat)
            return stat;
        MObject hemHeightData = arrayData.create(settings.panelMode ? panel.hemHeight : emptyHem, &stat);
        if (!stat)
            return stat;
        MDataHandle element = builder.addElement(static_cast<unsigned int>(panel.id), &stat);
        if (!stat)
            return stat;
        element.child(attr_patchSurface).setMObject(surfaceData[i]);
        element.child(attr_patchMaterialUStart).setDouble(panel.start);
        element.child(attr_patchMaterialUEnd).setDouble(panel.end);
        element.child(attr_patchVBreaks).setMObject(breakData);
        element.child(attr_patchHemUSamples).setMObject(hemUData);
        element.child(attr_patchHemHeightSamples).setMObject(hemHeightData);
    }
    MDataHandle surfaceOutput = dataBlock.outputValue(attr_outputSurface, &stat);
    if (!stat)
        return stat;
    MArrayDataHandle patchesOutput = dataBlock.outputArrayValue(attr_outputPatches, &stat);
    if (!stat)
        return stat;
    stat = patchesOutput.set(builder);
    if (!stat)
        return stat;
    if (!settings.panelMode)
        surfaceOutput.setMObject(surfaceData.back());
    patchesOutput.setAllClean();
    const MPlug patchesPlug(thisMObject(), attr_outputPatches);
    dataBlock.setClean(patchesPlug);
    for (const auto &panel : settings.panels)
    {
        MPlug element = patchesPlug.elementByLogicalIndex(static_cast<unsigned int>(panel.id));
        for (unsigned int i = 0; i < element.numChildren(); ++i)
            dataBlock.setClean(element.child(i));
        dataBlock.setClean(element);
    }
    if (!settings.panelMode)
        dataBlock.setClean(MPlug(thisMObject(), attr_outputSurface));
    dataBlock.setClean(plug);
    return MS::kSuccess;
}

MUserData *SkirtBellColliderDrawOverride::prepareForDraw(const MDagPath &objPath, const MDagPath &cameraPath,
                                                         const MHWRender::MFrameContext &frameContext,
                                                         MUserData *oldData)
{
    MStatus stat;
    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto *data = dynamic_cast<SkirtBellColliderDrawData *>(oldData);
    if (!data)
        data = new SkirtBellColliderDrawData();

    // Extract attributes
    auto getMatrix = [&obj](const MObject &attr, MMatrix &outMat) {
        MPlug plug(obj, attr);
        MObject tempObj;
        if (plug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
            outMat = MFnMatrixData(tempObj).matrix();
    };

    MMatrix leftHipMatrix, leftKneeMatrix, leftHeelMatrix;
    MMatrix rightHipMatrix, rightKneeMatrix, rightHeelMatrix;

    getMatrix(SkirtBellCollider::attr_leftHipMatrix, leftHipMatrix);
    getMatrix(SkirtBellCollider::attr_leftKneeMatrix, leftKneeMatrix);
    getMatrix(SkirtBellCollider::attr_leftHeelMatrix, leftHeelMatrix);
    getMatrix(SkirtBellCollider::attr_rightHipMatrix, rightHipMatrix);
    getMatrix(SkirtBellCollider::attr_rightKneeMatrix, rightKneeMatrix);
    getMatrix(SkirtBellCollider::attr_rightHeelMatrix, rightHeelMatrix);

    short skirtType = 1;
    MPlug(obj, SkirtBellCollider::attr_skirtType).getValue(skirtType);

    MVector ringScale(0.5, 1.0, 0.5);
    MPlug ringScalePlug(obj, SkirtBellCollider::attr_ringScale);
    if (ringScalePlug.numChildren() == 3)
    {
        double sx = 0.5, sy = 1.0, sz = 0.5;
        ringScalePlug.child(0).getValue(sx);
        ringScalePlug.child(1).getValue(sy);
        ringScalePlug.child(2).getValue(sz);
        ringScale = MVector(sx, sy, sz);
    }

    int ringSubdivision = 16;
    MPlug(obj, SkirtBellCollider::attr_ringSubdivision).getValue(ringSubdivision);

    short leftRingAxis = 0, rightRingAxis = 0;
    MPlug(obj, SkirtBellCollider::attr_leftRingAxis).getValue(leftRingAxis);
    MPlug(obj, SkirtBellCollider::attr_rightRingAxis).getValue(rightRingAxis);

    if (!ColliderInput::validSkirtType(skirtType) || !ColliderInput::validAxis(leftRingAxis) ||
        !ColliderInput::validAxis(rightRingAxis))
    {
        data->drawData.rings.update({}, 0);
        data->drawData.curves.update(MPointArray(), 0, 0);
        return data;
    }

    auto getDouble = [&obj](const MObject &attr, double value) {
        MPlug(obj, attr).getValue(value);
        return value;
    };
    const double radii[8] = {
        getDouble(SkirtBellCollider::attr_thighRadiusX, 1.0), getDouble(SkirtBellCollider::attr_thighRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_kneeRadiusX, 1.0),  getDouble(SkirtBellCollider::attr_kneeRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_calfRadiusX, 1.0),  getDouble(SkirtBellCollider::attr_calfRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_ankleRadiusX, 1.0), getDouble(SkirtBellCollider::attr_ankleRadiusZ, 1.0)};
    const MMatrix current[2][3] = {{leftHipMatrix, leftKneeMatrix, leftHeelMatrix},
                                   {rightHipMatrix, rightKneeMatrix, rightHeelMatrix}};
    const short axes[2] = {leftRingAxis, rightRingAxis};
    Leg legs[2];
    std::vector<LegSegment> restSegments;
    buildLegs(current, current, axes, ringScale, radii, getDouble(SkirtBellCollider::attr_thighPosition, 0.5),
              getDouble(SkirtBellCollider::attr_calfPosition, 0.5), skirtType == 1, legs, restSegments);
    std::vector<std::array<LegVec3, 4>> rows;
    std::vector<std::array<double, 2>> farMultipliers;
    for (const auto &leg : legs)
        for (const auto &segment : leg)
            legRingMatrices(segment.current, rows, farMultipliers);
    std::vector<MMatrix> matrices;
    for (const auto &ring : rows)
        matrices.push_back(legRingMatrix(ring));
    data->drawData.rings.update(matrices, ringSubdivision, farMultipliers);

    bool panelMode = false;
    MPlug seams(obj, SkirtBellCollider::attr_seams);
    const unsigned int seamCount = seams.numElements(&stat);
    if (!stat)
        panelMode = true;
    for (unsigned int i = 0; i < seamCount && !panelMode; ++i)
    {
        MPlug seam = seams.elementByPhysicalIndex(i, &stat);
        if (!stat)
        {
            panelMode = true;
            break;
        }
        bool enabled = true;
        stat = seam.child(SkirtBellCollider::attr_seamEnabled).getValue(enabled);
        panelMode = !stat || enabled;
    }
    if (!panelMode)
    {
        MPlug hems(obj, SkirtBellCollider::attr_panelHems);
        const unsigned int hemCount = hems.numElements(&stat);
        if (!stat)
            panelMode = true;
        for (unsigned int i = 0; i < hemCount && !panelMode; ++i)
        {
            MPlug hem = hems.elementByPhysicalIndex(i, &stat);
            if (!stat)
            {
                panelMode = true;
                break;
            }
            if (hem.logicalIndex() != 0)
                continue;
            MObject uData, heightData;
            MDoubleArray u, height;
            stat = hem.child(SkirtBellCollider::attr_hemUSamples).getValue(uData);
            if (stat)
                stat = hem.child(SkirtBellCollider::attr_hemHeightSamples).getValue(heightData);
            if (stat)
                stat = readDoubleArray(uData, u);
            if (stat)
                stat = readDoubleArray(heightData, height);
            panelMode = !stat || !validHem(u, height) || !defaultHem(height);
        }
    }
    std::vector<ColliderDraw::Curves::Element> elements;
    auto readSurface = [&](MPlug plug, unsigned int panelId)
    {
        MObject surfaceData;
        MStatus status = plug.getValue(surfaceData);
        if (!status || surfaceData.isNull())
            return false;
        MFnNurbsSurface surfaceFn(surfaceData, &status);
        if (!status)
            return false;
        ColliderDraw::Curves::Element element;
        element.panelId = panelId;
        status = surfaceFn.getCVs(element.points, MSpace::kObject);
        if (!status)
            return false;
        element.numU = surfaceFn.numCVsInU(&status);
        if (!status)
            return false;
        element.numV = surfaceFn.numCVsInV(&status);
        if (!status)
            return false;
        elements.push_back(element);
        return true;
    };
    if (panelMode)
    {
        MPlug patches(obj, SkirtBellCollider::attr_outputPatches);
        patches.evaluateNumElements(&stat);
        MIntArray indices;
        if (stat)
            patches.getExistingArrayAttributeIndices(indices, &stat);
        if (stat)
            for (unsigned int i = 0; i < indices.length(); ++i)
            {
                MPlug patch = patches.elementByLogicalIndex(indices[i], &stat);
                if (!stat)
                    break;
                MPlug surface = patch.child(SkirtBellCollider::attr_patchSurface, &stat);
                if (!stat)
                    break;
                if (!readSurface(surface, static_cast<unsigned int>(indices[i])))
                {
                    stat = MS::kFailure;
                    break;
                }
            }
        if (!stat)
            elements.clear();
    }
    else
        readSurface(MPlug(obj, SkirtBellCollider::attr_outputSurface), 0);
    data->drawData.curves.update(elements);
    // Default transparent cyan color for drawing collider rings
    data->drawData.color = MColor(0.0f, 0.6f, 1.0f, 0.25f);

    return data;
}
void SkirtBellColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath,
    MHWRender::MUIDrawManager& drawManager,
    const MHWRender::MFrameContext& frameContext,
    const MUserData* data)
{
    auto* skirtDrawData = dynamic_cast<const SkirtBellColliderDrawData*>(data);
    if (!skirtDrawData)
        return;

    const auto& drawData = skirtDrawData->drawData;

    drawManager.beginDrawable();

    drawData.rings.geometry.draw(drawManager, drawData.color, MColor(0.0f, 0.1f, 0.2f, 1.0f));
    drawManager.setColor(MColor(1.0f, 0.75f, 0.0f, 1.0f));
    drawManager.setLineWidth(2.0f);
    for (const auto &element : drawData.curves.elements)
        if (element.lines.length())
            drawManager.lineList(element.lines, false);

    drawManager.endDrawable();
}
