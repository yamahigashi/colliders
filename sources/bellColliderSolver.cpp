#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

#include "bellColliderSolver.h"
#include "colliderInputValidation.h"
#include "utils.hpp"
#include "bellColliderRelaxKernel.h"

using namespace std;

void BellColliderSolver::roundMeshPoints(MPointArray &points)
{
    // MFnMesh create/setPoints stores float coordinates; retain that boundary
    // when a caller consumes mesh-equivalent points without creating a mesh.
    for (unsigned int i = 0; i < points.length(); ++i)
        points.set(
            MPoint(static_cast<float>(points[i].x), static_cast<float>(points[i].y), static_cast<float>(points[i].z)),
            i);
}

BellCircleTable::BellCircleTable(int numSides)
{
    if (!ColliderInput::validSubdivision(numSides))
        return;
    cosines.reserve(numSides);
    sines.reserve(numSides);
    for (int i = 0; i < numSides; ++i)
    {
        const double rad = (double)i / numSides * 2 * M_PI;
        cosines.push_back(cos(rad));
        sines.push_back(sin(rad));
    }
}

MPointArray BellColliderSolver::makeBellPoints(const MMatrix &matrix, unsigned int axis, int numSides,
                                               double height, double bottomRadius, double topRadius)
{
    return makeBellPoints(matrix, axis, BellCircleTable(numSides), height, bottomRadius, topRadius);
}

MPointArray BellColliderSolver::makeBellPoints(const MMatrix &matrix, unsigned int axis, const BellCircleTable &circle,
                                               double height, double bottomRadius, double topRadius)
{
    MPointArray points;
    if (!circle.valid() || axis > 2)
        return points;
    const int numSides = circle.numSides();
    points.append(MPoint(0, 0, 0) * matrix);
    for (int row = 0; row < 2; ++row)
    {
        const double radius = row == 0 ? bottomRadius : topRadius;
        const double offset = row == 0 ? 0.0 : height;
        for (int i = 0; i < numSides; ++i)
        {
            const double x = radius * circle.cosines[i];
            const double z = radius * circle.sines[i];
            MPoint p;
            switch (axis)
            {
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

MObject BellColliderSolver::makeBellMesh(const MPointArray &points, int numSides)
{
    if (!ColliderInput::validSubdivision(numSides) ||
        points.length() != static_cast<unsigned int>(2 * numSides + 1))
        return MObject::kNullObj;
    MIntArray polygonCounts, polygonConnects;
    for (int i = 0; i < numSides; ++i)
    {
        polygonCounts.append(3);
        polygonConnects.append(0);
        polygonConnects.append(i + 1);
        polygonConnects.append(i == numSides - 1 ? 1 : i + 2);
    }
    for (int i = 0; i < numSides; ++i)
    {
        polygonCounts.append(4);
        polygonConnects.append(i + 1);
        polygonConnects.append(numSides + i + 1);
        polygonConnects.append(i == numSides - 1 ? numSides + 1 : numSides + i + 2);
        polygonConnects.append(i == numSides - 1 ? 1 : i + 2);
    }
    MFnMeshData meshData;
    MObject meshObject = meshData.create();
    MFnMesh meshFn;
    meshFn.create(points.length(), numSides * 2, points, polygonCounts, polygonConnects, meshObject);
    return meshObject;
}

MObject BellColliderSolver::makeBellCurve(const MPointArray &points, int bellSubdivision)
{
    if (!ColliderInput::validSubdivision(bellSubdivision) ||
        points.length() != static_cast<unsigned int>(2 * bellSubdivision + 1))
        return MObject::kNullObj;
    const int START = bellSubdivision + 1;
    const int END = points.length();

    MPointArray cvs;
    MDoubleArray knots;
    for (int i = START; i < END; i++)
    {
        knots.append(cvs.length());
        cvs.append(points[i]);
    }

    knots.append(knots[knots.length() - 1] + 1);
    cvs.append(cvs[0]);

    MFnNurbsCurveData curveDataFn;
    MObject curveData = curveDataFn.create();

    MFnNurbsCurve curveFn;
    curveFn.create(cvs, knots, 1, MFnNurbsCurve::kPeriodic, false, false, curveData);
    return curveData;
}

void BellColliderSolver::prepareDistalEnd(PreparedBellRing &ring, const MMatrix &matrix, bool distal)
{
    const double L = ring.direction.length();
    const double a = MVector(matrix[0][0], matrix[0][1], matrix[0][2]).length();
    const double b = MVector(matrix[2][0], matrix[2][1], matrix[2][2]).length();
    const double w = (std::min)(L, (std::max)(a, b));
    ring.distalLength = L;
    ring.distalWidth = w;
    ring.distalEnd = distal && std::isfinite(L) && std::isfinite(w);
}

namespace
{
template <class Ops>
BellColliderRelax::Vector<Ops> broadcast(double x, double y, double z)
{
    return {Ops::splat(x), Ops::splat(y), Ops::splat(z)};
}

template <class Ops>
BellColliderRelax::Ring<Ops> prepareRelaxRing(const PreparedBellRing &ring, double collision, bool capAtRingOrigin)
{
    MPoint origin = ring.plane.orig;
    MPoint translation = ring.translation;
    if (origin.w != 1.0)
        origin.cartesianize();
    if (translation.w != 1.0)
        translation.cartesianize();
    BellColliderRelax::Ring<Ops> result;
    result.origin = broadcast<Ops>(origin.x, origin.y, origin.z);
    result.normal = broadcast<Ops>(ring.plane.normal.x, ring.plane.normal.y, ring.plane.normal.z);
    result.translation = broadcast<Ops>(translation.x, translation.y, translation.z);
    for (int column = 0; column != 3; ++column)
        result.inverseColumns[column] = broadcast<Ops>(ring.inverse[0][column], ring.inverse[1][column],
                                                       ring.inverse[2][column]);
    result.collision = Ops::splat(collision);
    result.capAtRingOrigin = capAtRingOrigin;
    result.distalEnd = ring.distalEnd;
    result.distalLength = Ops::splat(ring.distalLength);
    result.distalWidth = Ops::splat(ring.distalWidth);
    return result;
}

MPoint cartesianPoint(const MPoint &point)
{
    MPoint result = point;
    if (result.w != 1.0)
        result.cartesianize();
    return result;
}
} // namespace

std::vector<MVector> BellColliderSolver::rowDirections(const MPointArray &points, const MPointArray &base)
{
    std::vector<MVector> directions(points.length());
    for (unsigned int i = 0; i < points.length(); ++i)
        directions[i] = cartesianPoint(points[i]) - cartesianPoint(base[i]);
    return directions;
}

void BellColliderSolver::relaxTowardRingBoundary(MPointArray &points, const PreparedBellRing &ring, double collision,
                                                 int startIndex, int count, bool capAtRingOrigin)
{
    if (!(collision > 1e-5))
        return;

#ifdef YDD_RELAX_SSE2
    const auto prepared = prepareRelaxRing<BellColliderRelax::Pair>(ring, collision, capAtRingOrigin);
    for (int j = startIndex; j < startIndex + count; j += 2)
    {
        const int second = j + 1 < startIndex + count ? j + 1 : j;
        MPoint &a = points[j];
        MPoint &b = points[second];
        const MPoint ca = cartesianPoint(a);
        const MPoint cb = cartesianPoint(b);
        const BellColliderRelax::Vector<BellColliderRelax::Pair> raw = {
            _mm_set_pd(b.x, a.x), _mm_set_pd(b.y, a.y), _mm_set_pd(b.z, a.z)};
        const BellColliderRelax::Vector<BellColliderRelax::Pair> cartesian = {
            _mm_set_pd(cb.x, ca.x), _mm_set_pd(cb.y, ca.y), _mm_set_pd(cb.z, ca.z)};
        const auto output = BellColliderRelax::relax(raw, cartesian, prepared);
        a.x = _mm_cvtsd_f64(output.x);
        a.y = _mm_cvtsd_f64(output.y);
        a.z = _mm_cvtsd_f64(output.z);
        if (second != j)
        {
            b.x = _mm_cvtsd_f64(_mm_unpackhi_pd(output.x, output.x));
            b.y = _mm_cvtsd_f64(_mm_unpackhi_pd(output.y, output.y));
            b.z = _mm_cvtsd_f64(_mm_unpackhi_pd(output.z, output.z));
        }
    }
#else
    const auto prepared = prepareRelaxRing<BellColliderRelax::Scalar>(ring, collision, capAtRingOrigin);
    for (int j = startIndex; j < startIndex + count; ++j)
    {
        MPoint &point = points[j];
        const MPoint cartesian = cartesianPoint(point);
        const BellColliderRelax::Vector<BellColliderRelax::Scalar> raw = {point.x, point.y, point.z};
        const BellColliderRelax::Vector<BellColliderRelax::Scalar> input = {cartesian.x, cartesian.y, cartesian.z};
        const auto output = BellColliderRelax::relax(raw, input, prepared);
        point.x = output.x;
        point.y = output.y;
        point.z = output.z;
    }
#endif
}

void BellColliderSolver::relaxTowardRingBoundary(MPointArray &points, const PreparedBellRing &ring, double collision,
                                                 int startIndex, int count, bool capAtRingOrigin,
                                                 const std::vector<MVector> &directions,
                                                 const std::vector<double> &betaScales)
{
    if (!(collision > 1e-5))
        return;
    // A missing direction means no lift (radial exit); a missing scale means no fade of the centre shift.
    const auto directionAt = [&](int index) -> MVector {
        return index >= 0 && static_cast<size_t>(index) < directions.size() ? directions[index] : MVector(0, 0, 0);
    };
    const auto scaleAt = [&](int index) -> double {
        return index >= 0 && static_cast<size_t>(index) < betaScales.size() ? betaScales[index] : 1.0;
    };

#ifdef YDD_RELAX_SSE2
    const auto prepared = prepareRelaxRing<BellColliderRelax::Pair>(ring, collision, capAtRingOrigin);
    for (int j = startIndex; j < startIndex + count; j += 2)
    {
        const int second = j + 1 < startIndex + count ? j + 1 : j;
        MPoint &a = points[j];
        MPoint &b = points[second];
        const MPoint ca = cartesianPoint(a);
        const MPoint cb = cartesianPoint(b);
        const BellColliderRelax::Vector<BellColliderRelax::Pair> raw = {_mm_set_pd(b.x, a.x), _mm_set_pd(b.y, a.y),
                                                                        _mm_set_pd(b.z, a.z)};
        const BellColliderRelax::Vector<BellColliderRelax::Pair> cartesian = {
            _mm_set_pd(cb.x, ca.x), _mm_set_pd(cb.y, ca.y), _mm_set_pd(cb.z, ca.z)};
        const MVector da = directionAt(j);
        const MVector db = directionAt(second);
        const BellColliderRelax::Vector<BellColliderRelax::Pair> direction = {
            _mm_set_pd(db.x, da.x), _mm_set_pd(db.y, da.y), _mm_set_pd(db.z, da.z)};
        const auto scales = _mm_set_pd(scaleAt(second), scaleAt(j));
        const auto output = BellColliderRelax::relax(raw, cartesian, prepared, direction, scales);
        a.x = _mm_cvtsd_f64(output.x);
        a.y = _mm_cvtsd_f64(output.y);
        a.z = _mm_cvtsd_f64(output.z);
        if (second != j)
        {
            b.x = _mm_cvtsd_f64(_mm_unpackhi_pd(output.x, output.x));
            b.y = _mm_cvtsd_f64(_mm_unpackhi_pd(output.y, output.y));
            b.z = _mm_cvtsd_f64(_mm_unpackhi_pd(output.z, output.z));
        }
    }
#else
    const auto prepared = prepareRelaxRing<BellColliderRelax::Scalar>(ring, collision, capAtRingOrigin);
    for (int j = startIndex; j < startIndex + count; ++j)
    {
        MPoint &point = points[j];
        const MPoint cartesian = cartesianPoint(point);
        const BellColliderRelax::Vector<BellColliderRelax::Scalar> raw = {point.x, point.y, point.z};
        const BellColliderRelax::Vector<BellColliderRelax::Scalar> input = {cartesian.x, cartesian.y, cartesian.z};
        const MVector d = directionAt(j);
        const BellColliderRelax::Vector<BellColliderRelax::Scalar> direction = {d.x, d.y, d.z};
        const auto output = BellColliderRelax::relax(raw, input, prepared, direction, scaleAt(j));
        point.x = output.x;
        point.y = output.y;
        point.z = output.z;
    }
#endif
}

bool BellColliderSolver::collisionPoints(const MMatrix &bellMatrix, const MMatrix &bellMatrixInverse,
                                         const Plane &bellPlane, const PreparedBellRing &ring,
                                         MPoint &collisionPointBell, MPoint &collisionPointRing, MPoint &linePoint,
                                         bool *extended)
{
    if (extended)
        *extended = false;
    const MMatrix &ringMatrixInverse = ring.inverse;
    const MVector &ringDirection = ring.direction;
    const MPoint &ring_translate = ring.translation;
    const MVector &ringNormal = ring.normal;
    const Plane &ringPlane = ring.plane;
    const MVector bellAxis = maxis(bellMatrix, 1);
    const MVector bellNormal = bellAxis.normal();
    const MPoint ring_translate_proj = bellPlane.projectPoint(ring_translate);
    const MVector ringDirection_proj = bellPlane.projectVector(ringDirection);
    if (ringDirection_proj.length() <= 1e-3)
        return false;
    bool found = false;
    const MPointArray hitPoints = findSphereLineIntersection(
        ring_translate_proj * bellMatrixInverse, ringDirection_proj * bellMatrixInverse, MPoint(0, 0, 0), 1.001);

    if (hitPoints.length() > 0)
    {
        collisionPointBell = hitPoints[0] * bellMatrix + bellAxis;

        const double linePointCoeff = ringNormal * bellNormal > 0 ? 1 : -1;

        const MVector ring_proj = ringPlane.projectVector(ringDirection_proj * linePointCoeff);
        double delta = 1.0;
        MVector ring_proj_scaled(0, 0, 0);
        double ring_proj_len = ring_proj.length();
        if (ring_proj_len > 1e-5)
        {
            double local_len = (ring_proj * ringMatrixInverse).length();
            if (local_len > 1e-5)
                delta = ring_proj_len / local_len;
            ring_proj_scaled = ring_proj.normal() * delta;
        }
        linePoint = ring_translate + ring_proj_scaled;

        const MPointArray sphereLinePoints = findSphereLineIntersection(linePoint, ringDirection, ring_translate,
                                                                        (collisionPointBell - ring_translate).length());

        for (int k = 0; k < sphereLinePoints.length(); k++)
        {
            if ((sphereLinePoints[k] - ring_translate) * ringDirection > 0)
            {
                collisionPointRing = sphereLinePoints[k];
                found = true;
            }
        }
        if (!found && extended && ring_proj_len > 1e-5)
        {
            collisionPointRing = linePoint;
            found = true;
            *extended = true;
        }
    }

    return found;
}

void BellColliderSolver::deformPoints(const BellColliderInputs& inputs, const std::vector<PreparedBellRing>& rings, const MPointArray& baseBellPoints, const Plane& bellPlane, vector<MPointArray>& bellPointsList)
{
    const MMatrix bellMatrix = inputs.bellMatrix;
    const MMatrix bellMatrixInverse = bellMatrix.inverse();
    const int bellSubdivision = inputs.bellSubdivision;
    const float falloff = inputs.falloff;

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();

    bellPointsList.clear();

    for (const auto &ring : rings)
    {
        const MPoint &ring_translate = ring.translation;
        const MPoint ring_translate_proj = bellPlane.projectPoint(ring_translate);
        const MVector ringDirection_proj = bellPlane.projectVector(ring.direction);

        MPointArray bellPoints = baseBellPoints;

        MPoint collisionPointBell, collisionPointRing, linePoint;
        if (collisionPoints(bellMatrix, bellMatrixInverse, bellPlane, ring, collisionPointBell, collisionPointRing,
                            linePoint))
        {
            double bellAxisLen = bellAxis.length();
            const double collisionDelta = bellAxisLen > 1e-5 ? (bellPlane.distance(collisionPointRing) - bellPlane.distance(collisionPointBell)) / bellAxisLen : 0.0;
            if (collisionDelta < 0)
            {
                MTransformationMatrix rotationMatrixFn;
                rotationMatrixFn.setTranslation(ring_translate, MSpace::kWorld);
                const MMatrix rotateMatrixInverse = rotationMatrixFn.asMatrixInverse();

                const MQuaternion quat(collisionPointBell - ring_translate, collisionPointRing - ring_translate); // rotate X to final point
                rotationMatrixFn.rotateBy(quat, MSpace::kTransform);
                const MMatrix rotateMatrix = rotationMatrixFn.asMatrix();

                const Plane upperBellPlane(bell_translate + bellAxis, bellNormal);

                // Ring terms in bell space do not change across the top vertices of this ring.
                const MPoint ring_translate_proj_bell = ring_translate_proj * bellMatrixInverse;
                const MVector ringDirection_proj_bell = (ringDirection_proj * bellMatrixInverse).normal();

                // bell top deformation
                for (int j = bellSubdivision + 1; j < (int)bellPoints.length(); j++)
                {
                    const MPoint bellPoint_proj = bellPlane.projectPoint(bellPoints[j]);
                    const MVector offset_proj = bellPoint_proj * bellMatrixInverse - ring_translate_proj_bell;

                    double weight = offset_proj.normal() * ringDirection_proj_bell; // -1..1

                    if (weight > falloff)
                    {
                        double divisor = 1.0 - falloff;
                        weight = divisor > 1e-5 ? (weight - falloff) / divisor : 1.0;
                        if (inputs.smoothness > 0.0)
                            weight = weight * weight * (3.0 - 2.0 * weight);

                        const MPoint rp = bellPoints[j] * rotateMatrixInverse * rotateMatrix;

                        // constrain by upper plane
                        MPoint p = rp;
                        if (upperBellPlane.distance(rp) > 0)
                            p = upperBellPlane.projectPoint(rp);

                        bellPoints[j] = p * weight + bellPoints[j] * (1.0 - weight);                               
                    }
                }
            }
        }

        relaxTowardRingBoundary(bellPoints, ring, inputs.collision, bellSubdivision + 1,
                                (int)bellPoints.length() - bellSubdivision - 1, inputs.capAtRingOrigin);

        bellPointsList.push_back(bellPoints);
    }
}

MVector BellColliderSolver::mergeDisplacement(const MPointArray& basePoints,
                                              const std::vector<MPointArray>& ringPoints, unsigned int index)
{
    // Squared-length weighted mean, rescaled so that parallel pushes keep the largest single length.
    // The rescale factor max / (sum |d|^3 / sum |d|^2) is bounded, so the merge stays continuous
    // through cancellation (the mean reaches zero before the factor can matter).
    double sum = 0;
    double cubic = 0;
    double longest = 0;
    for (const auto& points : ringPoints)
    {
        const double length = (points[index] - basePoints[index]).length();
        sum += pow(length, 2);
        cubic += length * length * length;
        if (length > longest)
            longest = length;
    }
    if (!(sum > 0))
        return MVector(0, 0, 0);
    MVector merged(0, 0, 0);
    for (const auto& points : ringPoints)
    {
        const MVector displacement = points[index] - basePoints[index];
        const double squaredLength = pow(displacement.length(), 2);
        const double weight = squaredLength / sum;
        merged += displacement * weight;
    }
    const double average = cubic / sum;
    if (average > 0.0)
        merged = merged * (longest / average);
    return merged;
}

void BellColliderSolver::averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const vector<MPointArray>& bellPointsList, MPointArray& outBellPoints)
{
    outBellPoints = baseBellPoints;
    for (unsigned int i = bellSubdivision + 1; i < baseBellPoints.length(); i++)
        outBellPoints[i] += mergeDisplacement(baseBellPoints, bellPointsList, i);
}

void BellColliderSolver::smoothDisplacements(vector<MVector> &displacements, double smoothness)
{
    const int count = static_cast<int>(displacements.size());
    const double alpha = smoothness * 0.5;
    if (alpha != 0.0)
    {
        vector<MVector> smoothedDisplacements(count);
        for (int iteration = 0; iteration < 3; iteration++)
        {
            for (int i = 0; i < count; i++)
            {
                const int previous = (i + count - 1) % count;
                const int next = (i + 1) % count;
                smoothedDisplacements[i] = displacements[i] * (1.0 - alpha)
                    + (displacements[previous] + displacements[next]) * (alpha * 0.5);
            }
            displacements.swap(smoothedDisplacements);
        }
    }
}

MStatus BellColliderSolver::solve(const BellColliderInputs &inputs, const MPointArray &baseBellPoints,
                                  BellColliderOutputs &outputs)
{
    outputs.points.clear();
    outputs.meanDisplacement = MVector(0, 0, 0);
    if (!ColliderInput::validSubdivision(inputs.bellSubdivision) ||
        baseBellPoints.length() != static_cast<unsigned int>(2 * inputs.bellSubdivision + 1))
        return MS::kInvalidParameter;
    const MMatrix bellMatrix = inputs.bellMatrix;
    const int bellSubdivision = inputs.bellSubdivision;

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();
    const Plane bellPlane(bell_translate, bellNormal);
    const bool gate = inputs.smoothness > 0.0 || inputs.followGain > 0.0;

    vector<MPointArray> bellPointsList;
    deformPoints(inputs, inputs.rings, baseBellPoints, bellPlane, bellPointsList);

    MPointArray outBellPoints;
    averageDisplacements(bellSubdivision, baseBellPoints, bellPointsList, outBellPoints);

    if (gate)
    {
        const int startIndex = bellSubdivision + 1;
        vector<MVector> displacements(bellSubdivision);
        MVector weightedDisplacement(0, 0, 0);
        double displacementLengthSum = 0.0;

        for (int i = 0; i < bellSubdivision; i++)
        {
            const MVector displacement = outBellPoints[startIndex + i] - baseBellPoints[startIndex + i];
            const double displacementLength = displacement.length();
            displacements[i] = displacement;
            weightedDisplacement += displacement * displacementLength;
            displacementLengthSum += displacementLength;
        }

        const MVector meanDisplacement = displacementLengthSum < 1e-12
            ? MVector(0, 0, 0)
            : bellPlane.projectVector(weightedDisplacement / displacementLengthSum);
        outputs.meanDisplacement = meanDisplacement;

        for (int i = 0; i < bellSubdivision; i++)
            displacements[i] += meanDisplacement * inputs.followGain;

        smoothDisplacements(displacements, inputs.smoothness);

        for (int i = 0; i < bellSubdivision; i++)
            outBellPoints.set(baseBellPoints[startIndex + i] + displacements[i], startIndex + i);

        for (const auto &ring : inputs.rings)
            relaxTowardRingBoundary(outBellPoints, ring, inputs.collision, startIndex, bellSubdivision,
                                    inputs.capAtRingOrigin);
    }

    outputs.points = outBellPoints;
    return MS::kSuccess;
}

namespace
{
bool validRowTopology(const BellRowTopology &topology, size_t count)
{
    if (topology.vertices.size() != count)
        return false;
    std::vector<bool> mapped(topology.outputCount, false);
    std::vector<std::vector<int>> members(topology.components.size());
    for (size_t i = 0; i < count; ++i)
    {
        const auto &vertex = topology.vertices[i];
        if (!std::isfinite(vertex.materialU) || vertex.componentId < 0 ||
            static_cast<size_t>(vertex.componentId) >= members.size() || vertex.previous < -1 || vertex.next < -1 ||
            vertex.previous >= static_cast<int>(count) || vertex.next >= static_cast<int>(count) ||
            vertex.side.bank < -1 || vertex.side.bank > 1 || (vertex.side.bank != 0 && vertex.side.seamIndex < 0) ||
            vertex.outputDuplicates.empty())
            return false;
        for (unsigned int output : vertex.outputDuplicates)
        {
            if (output >= mapped.size() || mapped[output])
                return false;
            mapped[output] = true;
        }
        members[vertex.componentId].push_back(static_cast<int>(i));
    }
    for (bool present : mapped)
        if (!present)
            return false;
    for (size_t c = 0; c < members.size(); ++c)
    {
        const auto &component = topology.components[c];
        auto &indices = members[c];
        std::sort(indices.begin(), indices.end(),
                  [&](int a, int b) { return topology.vertices[a].materialU < topology.vertices[b].materialU; });
        if (!std::isfinite(component.startU) || !std::isfinite(component.endU) || component.endU <= component.startU ||
            indices.empty())
            return false;
        if (component.closed && (indices.size() < 2 || members.size() != 1 || component.endU != component.startU + 1.0))
            return false;
        for (size_t other = 0; other < c; ++other)
            if ((std::max)(component.startU, topology.components[other].startU) <
                (std::min)(component.endU, topology.components[other].endU))
                return false;
        for (size_t j = 0; j < indices.size(); ++j)
        {
            const auto &vertex = topology.vertices[indices[j]];
            if (vertex.materialU < component.startU - 1e-9 || vertex.materialU > component.endU + 1e-9 ||
                (j != 0 && vertex.materialU <= topology.vertices[indices[j - 1]].materialU))
                return false;
            const int previous = j != 0 ? indices[j - 1] : (component.closed ? indices.back() : -1);
            const int next = j + 1 != indices.size() ? indices[j + 1] : (component.closed ? indices.front() : -1);
            if (vertex.previous != previous || vertex.next != next || (component.closed && vertex.side.bank != 0) ||
                (vertex.side.bank == -1 && previous != -1) || (vertex.side.bank == 1 && next != -1))
                return false;
        }
    }
    return true;
}

bool finiteRowVectors(const std::vector<MVector> &values)
{
    for (const auto &value : values)
        if (!std::isfinite(value.x) || !std::isfinite(value.y) || !std::isfinite(value.z))
            return false;
    return true;
}

bool finiteRowPoints(const MPointArray &points)
{
    for (unsigned int i = 0; i < points.length(); ++i)
        if (!std::isfinite(points[i].x) || !std::isfinite(points[i].y) || !std::isfinite(points[i].z) ||
            !std::isfinite(points[i].w))
            return false;
    return true;
}

bool validRowControls(const BellRowInputs &inputs)
{
    return std::isfinite(inputs.smoothness) && inputs.smoothness >= 0.0 && inputs.smoothness <= 1.0 &&
           std::isfinite(inputs.followGain) && inputs.followGain >= 0.0 && std::isfinite(inputs.followRange) &&
           inputs.followRange >= 0.0 && std::isfinite(inputs.contactBlendWidth) && inputs.contactBlendWidth >= 0.0;
}

void rowDistances(const BellRowTopology &topology, const std::vector<unsigned int> &seeds, double range,
                  std::vector<double> &distances, std::vector<unsigned int> &touched)
{
    for (unsigned int index : touched)
        distances[index] = std::numeric_limits<double>::infinity();
    touched.clear();
    typedef std::pair<double, unsigned int> Entry;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pending;
    for (unsigned int seed : seeds)
    {
        if (distances[seed] == 0.0)
            continue;
        touched.push_back(seed);
        distances[seed] = 0.0;
        pending.push(Entry(0.0, seed));
    }
    while (!pending.empty())
    {
        const Entry current = pending.top();
        pending.pop();
        if (current.first > range)
            break;
        if (current.first != distances[current.second])
            continue;
        const auto &vertex = topology.vertices[current.second];
        for (int direction = 0; direction != 2; ++direction)
        {
            const int neighbor = direction == 0 ? vertex.previous : vertex.next;
            if (neighbor < 0)
                continue;
            const auto &other = topology.vertices[neighbor];
            double length = std::abs(other.materialU - vertex.materialU);
            if (topology.components[vertex.componentId].closed &&
                ((direction == 0 && other.materialU > vertex.materialU) ||
                 (direction == 1 && other.materialU < vertex.materialU)))
                length =
                    (std::min)(other.materialU, vertex.materialU) + 1.0 - (std::max)(other.materialU, vertex.materialU);
            const double candidate = current.first + length;
            if (candidate <= range && candidate < distances[neighbor])
            {
                if (!std::isfinite(distances[neighbor]))
                    touched.push_back(static_cast<unsigned int>(neighbor));
                distances[neighbor] = candidate;
                pending.push(Entry(candidate, static_cast<unsigned int>(neighbor)));
            }
        }
    }
}

// Per-component rotation weight. The ray from the bell centre through the
// bell contact point of the ring (the same sphere hit that drives the rotation)
// crosses the row polyline at a material position. A closed row keeps weight
// 1. An open component that does not contain that position has weight 0. A
// component that contains it fades with the material distance from the
// position to its nearer free end, reaching exactly 0 at the end and 1 at
// blendWidth: a free end near the contact has no hoop tension to lift, so the
// lift disappears towards a seam instead of switching between panels. The
// bell centre is used because every row surrounds it, whereas a ring close to
// or outside the row would see the CVs under magnified angles or miss the row.
// Without a forward crossing only the seed's component rotates.
std::vector<double> componentContactWeights(const BellRowTopology &topology, const std::vector<MVector> &offsets,
                                            const MVector &direction, unsigned int seed, double blendWidth)
{
    std::vector<double> weights(topology.components.size(), 1.0);
    if (topology.components.size() == 1 && topology.components[0].closed)
        return weights;
    const auto cross = [](const MVector &a, const MVector &b) { return a.z * b.x - a.x * b.z; };
    double nearest = std::numeric_limits<double>::infinity();
    double contactU = 0.0;
    // A crossing on the edge (i, j). Virtual edges join the two banks of a
    // seam so a ray through the gap between displaced banks still finds a
    // contact position; every point of such an edge has the seam's material U.
    const auto consider = [&](size_t i, size_t j, bool virtualEdge) {
        const auto &vertex = topology.vertices[i];
        const double ci = cross(offsets[i], direction);
        const double cj = cross(offsets[j], direction);
        if ((ci == 0.0 && cj == 0.0) || !((ci >= 0.0 && cj <= 0.0) || (ci <= 0.0 && cj >= 0.0)))
            return;
        MVector hit;
        double u;
        if (ci == 0.0)
        {
            hit = offsets[i];
            u = vertex.materialU;
        }
        else if (cj == 0.0)
        {
            hit = offsets[j];
            u = virtualEdge ? vertex.materialU : topology.vertices[j].materialU;
        }
        else
        {
            const double lambda = ci / (ci - cj);
            hit = offsets[i] + (offsets[j] - offsets[i]) * lambda;
            double span = virtualEdge ? 0.0 : topology.vertices[j].materialU - vertex.materialU;
            if (span < 0.0)
                span += 1.0;
            u = vertex.materialU + span * lambda;
        }
        const double t = hit * direction;
        if (t > 0.0 && t < nearest)
        {
            nearest = t;
            contactU = u;
        }
    };
    for (size_t i = 0; i < offsets.size(); ++i)
        if (topology.vertices[i].next >= 0)
            consider(i, static_cast<size_t>(topology.vertices[i].next), false);
    for (size_t i = 0; i < offsets.size(); ++i)
    {
        const auto &end = topology.vertices[i].side;
        if (end.bank != 1)
            continue;
        for (size_t j = 0; j < offsets.size(); ++j)
        {
            const auto &start = topology.vertices[j].side;
            if (start.bank == -1 && start.seamIndex == end.seamIndex)
                consider(i, j, true);
        }
    }
    if (!std::isfinite(nearest))
    {
        std::fill(weights.begin(), weights.end(), 0.0);
        weights[topology.vertices[seed].componentId] = 1.0;
        return weights;
    }
    for (size_t c = 0; c < weights.size(); ++c)
    {
        const auto &component = topology.components[c];
        if (component.closed)
            continue;
        const double span = component.endU - component.startU;
        double inside = contactU - component.startU;
        inside -= std::floor(inside);
        if (inside > span)
        {
            weights[c] = 0.0;
            continue;
        }
        const double edge = (std::min)(inside, span - inside);
        if (blendWidth <= 0.0 || edge >= blendWidth)
            continue;
        const double t = edge / blendWidth;
        weights[c] = t * t * (3.0 - 2.0 * t);
    }
    return weights;
}

double liftSmoothstep(double x)
{
    const double t = (std::min)(1.0, (std::max)(0.0, x));
    return (t * t) * (3.0 - 2.0 * t);
}

struct LiftContact
{
    bool found = false;
    bool extended = false;
    MPoint bellHit, ringHit, linePoint;
    double collisionDelta = 0.0;
    double R = 0.0, rOffset = 0.0, r = 0.0, H = 0.0;
    double byHeight = 0.0, byMargin = 0.0, a = 0.0;
};

LiftContact liftContact(const MMatrix &bellMatrix, const MMatrix &bellInverse, const Plane &bellPlane,
                        const MPointArray &baseRow, const PreparedBellRing &ring)
{
    LiftContact contact;
    const MVector axis = maxis(bellMatrix, 1);
    if (baseRow.length() == 0 || axis.length() <= 1e-5)
        return contact;
    contact.found = BellColliderSolver::collisionPoints(bellMatrix, bellInverse, bellPlane, ring, contact.bellHit,
                                                        contact.ringHit, contact.linePoint, &contact.extended);
    if (!contact.found)
        return contact;
    contact.collisionDelta =
        (bellPlane.distance(contact.ringHit) - bellPlane.distance(contact.bellHit)) / axis.length();
    if (!(contact.collisionDelta < 0.0))
        return contact;
    contact.R = (contact.bellHit - ring.translation).length();
    contact.rOffset = (contact.linePoint - ring.translation).length();
    const MMatrix matrix = ring.inverse.inverse();
    const MVector A = maxis(matrix, 0), B = maxis(matrix, 2);
    contact.r = (std::max)(A.length(), B.length());
    const double L = ring.direction.length();
    const MVector up = -bellPlane.normal;
    contact.H = std::hypot(A * up, B * up);
    for (unsigned int i = 0; i < baseRow.length(); ++i)
    {
        MPoint point = baseRow[i];
        if (point.w != 1.0)
            point.cartesianize();
        const MVector local = MVector(point - ring.translation) * ring.inverse;
        const double axial = local.y * L;
        const double fade = contact.r > 1e-12
                                ? liftSmoothstep((axial + 0.25 * contact.r) / (0.25 * contact.r))
                                : (axial >= 0.0 ? 1.0 : 0.0);
        const double pen = 1.0 - std::hypot(local.x, local.z);
        const MPoint centre = ring.translation + ring.direction * local.y;
        const double h = contact.H - (point - centre) * up;
        const double g = liftSmoothstep((pen + 0.10) / 0.10);
        const double b = contact.r > 1e-12
                             ? g * liftSmoothstep((h - 0.25 * contact.r) / (0.5 * contact.r))
                             : 0.0;
        if (fade <= 0.0 || g <= 0.0 || contact.r <= 1e-12)
            continue;
        const double contribution = b * fade;
        if (contribution > contact.byHeight)
            contact.byHeight = contribution;
    }
    contact.byMargin = !contact.extended && contact.r > 1e-12
                           ? liftSmoothstep((contact.R - contact.rOffset) / (0.5 * contact.r))
                           : 0.0;
    double amplitude = (std::max)(contact.byHeight, contact.byMargin);
    const MVector projected = bellPlane.projectVector(ring.direction);
    if (L > 1e-12)
        amplitude *= liftSmoothstep(projected.length() / (0.05 * L));
    contact.a = amplitude;
    return contact;
}

std::vector<unsigned int> rowContactSeeds(const BellRowInputs &inputs, const MPointArray &baseRow,
                                          const BellRowTopology &topology, const PreparedBellRing &ring)
{
    const MMatrix inverse = inputs.bellMatrix.inverse();
    const MVector axis = maxis(inputs.bellMatrix, 1);
    const Plane plane(taxis(inputs.bellMatrix), axis.normal());
    const LiftContact contact = liftContact(inputs.bellMatrix, inverse, plane, baseRow, ring);
    if (!contact.found || !std::isfinite(contact.a) || !(contact.a > 0.0))
        return {};
    const MPoint ringProjection = plane.projectPoint(ring.translation) * inverse;
    const MVector direction = (plane.projectVector(ring.direction) * inverse).normal();
    unsigned int seed = 0;
    double maximum = -std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < baseRow.length(); ++i)
    {
        const MPoint projection = plane.projectPoint(baseRow[i]);
        const MVector offset = projection * inverse - ringProjection;
        const double z = offset.normal() * direction;
        if (z > maximum || (z == maximum && topology.vertices[i].materialU < topology.vertices[seed].materialU))
        {
            maximum = z;
            seed = i;
        }
    }
    return {seed};
}
} // namespace

MStatus BellColliderSolver::smoothDisplacements(std::vector<MVector> &displacements, double smoothness,
                                                const BellRowTopology &topology)
{
    if (!validRowTopology(topology, displacements.size()) || !finiteRowVectors(displacements) ||
        !std::isfinite(smoothness) || smoothness < 0.0 || smoothness > 1.0)
        return MS::kInvalidParameter;
    if (smoothness == 0.0)
        return MS::kSuccess;
    const double alpha = smoothness / 2.0;
    std::vector<MVector> next(displacements.size());
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        for (size_t i = 0; i < displacements.size(); ++i)
        {
            const auto &vertex = topology.vertices[i];
            if (vertex.previous < 0 && vertex.next < 0)
                next[i] = displacements[i];
            else
            {
                const size_t previous = vertex.previous < 0 ? i : static_cast<size_t>(vertex.previous);
                const size_t following = vertex.next < 0 ? i : static_cast<size_t>(vertex.next);
                next[i] = displacements[i] * (1.0 - alpha) +
                          (displacements[previous] + displacements[following]) * (alpha / 2.0);
            }
        }
        displacements.swap(next);
    }
    return MS::kSuccess;
}

MStatus BellColliderSolver::computeLocalFollow(const std::vector<MVector> &directDisplacements,
                                               const BellRowTopology &topology, const MMatrix &bellMatrix,
                                               double followRange, std::vector<MVector> &follow)
{
    if (!validRowTopology(topology, directDisplacements.size()) || !finiteRowVectors(directDisplacements) ||
        !std::isfinite(followRange) || followRange < 0.0)
        return MS::kInvalidParameter;
    std::vector<MVector> result(directDisplacements.size(), MVector(0, 0, 0));
    std::vector<double> distances(result.size(), std::numeric_limits<double>::infinity());
    std::vector<unsigned int> touched;
    std::vector<size_t> order(result.size());
    for (size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const auto &left = topology.vertices[a];
        const auto &right = topology.vertices[b];
        if (left.materialU != right.materialU)
            return left.materialU < right.materialU;
        if (left.panelId != right.panelId)
            return left.panelId < right.panelId;
        return left.side.bank < right.side.bank;
    });
    const MVector axis = maxis(bellMatrix, 1);
    if (axis.length() > 1e-5)
    {
        const MVector normal = axis.normal();
        for (size_t i = 0; i < result.size(); ++i)
        {
            MVector value = directDisplacements[i];
            if (followRange > 0.0)
            {
                rowDistances(topology, {static_cast<unsigned int>(i)}, followRange, distances, touched);
                MVector sum(0, 0, 0);
                double weight = 0.0;
                for (size_t j : order)
                    if (distances[j] <= followRange)
                    {
                        const double length = directDisplacements[j].length();
                        sum += directDisplacements[j] * length;
                        weight += length;
                    }
                if (weight < 1e-12)
                    continue;
                value = sum / weight;
            }
            result[i] = value - normal * (value * normal);
        }
    }
    follow.swap(result);
    return MS::kSuccess;
}

MStatus BellColliderSolver::deformPoints(const BellRowInputs &inputs, const MPointArray &baseRow,
                                         const BellRowTopology &topology,
                                         const std::vector<std::vector<unsigned int>> &ringSeeds,
                                         std::vector<MPointArray> &ringPoints,
                                         std::vector<std::vector<double>> *componentWeights)
{
    if (!validRowTopology(topology, baseRow.length()) || ringSeeds.size() != inputs.rings.size() ||
        !finiteRowPoints(baseRow) || !validRowControls(inputs))
        return MS::kInvalidParameter;
    for (const auto &seeds : ringSeeds)
    {
        if (seeds.size() > 1)
            return MS::kInvalidParameter;
        for (unsigned int seed : seeds)
            if (seed >= baseRow.length())
                return MS::kInvalidParameter;
    }
    const MMatrix bellMatrix = inputs.bellMatrix;
    const MMatrix bellMatrixInverse = bellMatrix.inverse();
    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1);
    const MVector bellNormal = bellAxis.normal();
    const Plane bellPlane(bell_translate, bellNormal);
    std::vector<MPointArray> result;
    std::vector<std::vector<double>> weights(inputs.rings.size(),
                                             std::vector<double>(topology.components.size(), 1.0));
    for (size_t r = 0; r < inputs.rings.size(); ++r)
    {
        const auto &ring = inputs.rings[r];
        MPointArray points = baseRow;
        LiftContact contact;
        if (bellAxis.length() > 1e-5 && !ringSeeds[r].empty())
            contact = liftContact(bellMatrix, bellMatrixInverse, bellPlane, baseRow, ring);
        if (contact.found && std::isfinite(contact.a) && contact.a > 0.0)
        {
            const MPoint &ring_translate = ring.translation;
            const MPoint ring_translate_proj = bellPlane.projectPoint(ring_translate);
            const MVector ringDirection_proj = bellPlane.projectVector(ring.direction);
            MTransformationMatrix rotationMatrixFn;
            rotationMatrixFn.setTranslation(ring_translate, MSpace::kWorld);
            const MMatrix rotateMatrixInverse = rotationMatrixFn.asMatrixInverse();
            const MQuaternion fullQuat(contact.bellHit - ring_translate, contact.ringHit - ring_translate);
            const MQuaternion quat = contact.a >= 1.0
                                         ? fullQuat
                                         : slerp(MQuaternion::identity, fullQuat, contact.a);
            rotationMatrixFn.rotateBy(quat, MSpace::kTransform);
            const MMatrix rotateMatrix = rotationMatrixFn.asMatrix();
            const Plane upperBellPlane(bell_translate + bellAxis, bellNormal);
            const MPoint ring_translate_proj_bell = ring_translate_proj * bellMatrixInverse;
            const MVector ringDirection_proj_bell = (ringDirection_proj * bellMatrixInverse).normal();
            const MPointArray bellHits = findSphereLineIntersection(
                ring_translate_proj * bellMatrixInverse, ringDirection_proj * bellMatrixInverse, MPoint(0, 0, 0),
                1.001);
            // collisionPoints succeeded with the same arguments, so a hit exists;
            // an empty array would only mean the two call sites drifted apart.
            const MVector contactRay = bellHits.length() ? MVector(bellHits[0]) : MVector(0, 0, 0);
            std::vector<MVector> offsets(points.length());
            for (unsigned int i = 0; i < points.length(); ++i)
                offsets[i] = MVector(bellPlane.projectPoint(points[i]) * bellMatrixInverse);
            const std::vector<double> &componentWeights = weights[r] = componentContactWeights(
                topology, offsets, contactRay, ringSeeds[r][0], inputs.contactBlendWidth);
            bool partial = false;
            for (double componentWeight : componentWeights)
                partial = partial || (componentWeight > 0.0 && componentWeight < 1.0);
            for (unsigned int i = 0; i < points.length(); ++i)
            {
                if (componentWeights[topology.vertices[i].componentId] <= 0.0)
                    continue;
                const MPoint bellPoint_proj = bellPlane.projectPoint(points[i]);
                const MVector offset_proj = bellPoint_proj * bellMatrixInverse - ring_translate_proj_bell;
                double weight = offset_proj.normal() * ringDirection_proj_bell;
                if (weight > inputs.falloff)
                {
                    const double divisor = 1.0 - inputs.falloff;
                    weight = divisor > 1e-5 ? (weight - inputs.falloff) / divisor : 1.0;
                    if (inputs.smoothness > 0.0)
                        weight = weight * weight * (3.0 - 2.0 * weight);
                    const MPoint rp = points[i] * rotateMatrixInverse * rotateMatrix;
                    MPoint p = rp;
                    if (upperBellPlane.distance(rp) > 0)
                        p = upperBellPlane.projectPoint(rp);
                    points[i] = p * weight + points[i] * (1.0 - weight);
                }
            }
            relaxTowardRingBoundary(points, ring, inputs.collision, 0, static_cast<int>(points.length()),
                                    inputs.capAtRingOrigin, rowDirections(points, baseRow), {});
            if (partial)
            {
                // A partially weighted component blends the relaxed result of
                // the full rotation with the relaxed result of no rotation.
                // Blending before the relax would leave partially lifted points
                // inside the ring, where the relax pushes them out along a path
                // unrelated to the weight.
                MPointArray still = baseRow;
                relaxTowardRingBoundary(still, ring, inputs.collision, 0, static_cast<int>(still.length()),
                                        inputs.capAtRingOrigin,
                                        std::vector<MVector>(still.length(), MVector(0, 0, 0)), {});
                for (unsigned int i = 0; i < points.length(); ++i)
                {
                    const double componentWeight = componentWeights[topology.vertices[i].componentId];
                    if (componentWeight <= 0.0 || componentWeight >= 1.0)
                        continue;
                    const MPoint &lifted = points[i];
                    const MPoint &rest = still[i];
                    points[i] = MPoint(rest.x * (1.0 - componentWeight) + lifted.x * componentWeight,
                                       rest.y * (1.0 - componentWeight) + lifted.y * componentWeight,
                                       rest.z * (1.0 - componentWeight) + lifted.z * componentWeight, lifted.w);
                }
            }
            result.push_back(points);
            continue;
        }
        relaxTowardRingBoundary(points, ring, inputs.collision, 0, static_cast<int>(points.length()),
                                inputs.capAtRingOrigin, std::vector<MVector>(points.length(), MVector(0, 0, 0)), {});
        result.push_back(points);
    }
    ringPoints.swap(result);
    if (componentWeights)
        componentWeights->swap(weights);
    return MS::kSuccess;
}

// Final relax of a merged row for one ring. Points of a component that ring
// lifts with a partial weight sit on the chord between the relaxed still and
// relaxed lifted states, inside the ring; relaxing them fully would snap them
// back to the boundary and undo the blend, so the relax fades out towards the
// middle of the transition (f = 1 - 4m(1-m)) and is complete at both ends.
void BellColliderSolver::relaxRowFaded(MPointArray &points, const PreparedBellRing &ring, double collision,
                                       bool capAtRingOrigin, const BellRowTopology &topology,
                                       const std::vector<double> &componentWeights)
{
    bool partial = false;
    for (double weight : componentWeights)
        partial = partial || (weight > 0.0 && weight < 1.0);
    if (!partial || topology.vertices.size() != points.length())
    {
        relaxTowardRingBoundary(points, ring, collision, 0, static_cast<int>(points.length()), capAtRingOrigin);
        return;
    }
    MPointArray relaxed = points;
    relaxTowardRingBoundary(relaxed, ring, collision, 0, static_cast<int>(relaxed.length()), capAtRingOrigin);
    for (unsigned int i = 0; i < points.length(); ++i)
    {
        const double m = componentWeights[topology.vertices[i].componentId];
        if (m <= 0.0 || m >= 1.0)
        {
            points[i] = relaxed[i];
            continue;
        }
        const double f = 1.0 - 4.0 * m * (1.0 - m);
        const MPoint &before = points[i];
        const MPoint &after = relaxed[i];
        points[i] = MPoint(before.x * (1.0 - f) + after.x * f, before.y * (1.0 - f) + after.y * f,
                           before.z * (1.0 - f) + after.z * f, after.w);
    }
}

void BellColliderSolver::relaxRowFaded(MPointArray &points, const PreparedBellRing &ring, double collision,
                                       bool capAtRingOrigin, const BellRowTopology &topology,
                                       const std::vector<double> &componentWeights,
                                       const std::vector<MVector> &directions)
{
    bool partial = false;
    for (double weight : componentWeights)
        partial = partial || (weight > 0.0 && weight < 1.0);
    if (!partial || topology.vertices.size() != points.length())
    {
        relaxTowardRingBoundary(points, ring, collision, 0, static_cast<int>(points.length()), capAtRingOrigin,
                                directions, {});
        return;
    }
    // A component without a weight entry is treated as fully weighted (no partial blend).
    const auto weightAt = [&](unsigned int index) -> double {
        const int component = topology.vertices[index].componentId;
        return component >= 0 && static_cast<size_t>(component) < componentWeights.size() ? componentWeights[component]
                                                                                            : 1.0;
    };
    std::vector<double> betaScales(points.length(), 1.0);
    for (unsigned int i = 0; i < points.length(); ++i)
    {
        const double m = weightAt(i);
        if (m > 0.0 && m < 1.0)
            betaScales[i] = 1.0 - 4.0 * m * (1.0 - m);
    }
    MPointArray relaxed = points;
    relaxTowardRingBoundary(relaxed, ring, collision, 0, static_cast<int>(relaxed.length()), capAtRingOrigin,
                            directions, betaScales);
    for (unsigned int i = 0; i < points.length(); ++i)
    {
        const double m = weightAt(i);
        if (m <= 0.0 || m >= 1.0)
        {
            points[i] = relaxed[i];
            continue;
        }
        const double f = 1.0 - 4.0 * m * (1.0 - m);
        const MPoint &before = points[i];
        const MPoint &after = relaxed[i];
        points[i] = MPoint(before.x * (1.0 - f) + after.x * f, before.y * (1.0 - f) + after.y * f,
                           before.z * (1.0 - f) + after.z * f, after.w);
    }
}

MStatus BellColliderSolver::solveRow(const BellRowInputs &inputs, const MPointArray &baseRow,
                                     const BellRowTopology &topology, BellRowOutputs &outputs)
{
    if (!validRowTopology(topology, baseRow.length()) || !finiteRowPoints(baseRow) || !validRowControls(inputs))
        return MS::kInvalidParameter;
    BellRowOutputs result;
    result.points = baseRow;
    result.directDisplacements.assign(baseRow.length(), MVector(0, 0, 0));
    result.directField.topology = topology;
    result.directField.values.assign(baseRow.length(), MVector(0, 0, 0));
    if (baseRow.length() == 0 || inputs.rings.empty())
    {
        outputs = result;
        return MS::kSuccess;
    }
    std::vector<std::vector<unsigned int>> seeds;
    for (const auto &ring : inputs.rings)
        seeds.push_back(rowContactSeeds(inputs, baseRow, topology, ring));
    std::vector<MPointArray> ringPoints;
    std::vector<std::vector<double>> componentWeights;
    const MStatus status = deformPoints(inputs, baseRow, topology, seeds, ringPoints, &componentWeights);
    if (!status)
        return status;
    const bool gate = inputs.smoothness > 0.0 || inputs.followGain > 0.0;
    for (unsigned int i = 0; i < baseRow.length(); ++i)
        result.directDisplacements[i] = mergeDisplacement(baseRow, ringPoints, i);
    computeLocalFollow(result.directDisplacements, topology, inputs.bellMatrix, inputs.followRange,
                       result.directField.values);
    std::vector<MVector> displacements = result.directDisplacements;
    if (gate)
    {
        for (size_t i = 0; i < displacements.size(); ++i)
            displacements[i] += result.directField.values[i] * inputs.followGain;
        smoothDisplacements(displacements, inputs.smoothness, topology);
    }
    for (unsigned int i = 0; i < baseRow.length(); ++i)
        result.points[i] = baseRow[i] + displacements[i];
    if (gate)
        for (size_t r = 0; r < inputs.rings.size(); ++r)
            relaxRowFaded(result.points, inputs.rings[r], inputs.collision, inputs.capAtRingOrigin, topology,
                          componentWeights[r], rowDirections(result.points, baseRow));
    result.componentWeights.swap(componentWeights);
    outputs = result;
    return MS::kSuccess;
}

void BellColliderSolver::transferRowValuesImpl(const BellRowTopology &source, const BellRowTopology &destination,
                                               BellRowCorrespondence &correspondence)
{
    correspondence = BellRowCorrespondence();
    correspondence.sourceCount = source.vertices.size();
    correspondence.entries.resize(destination.vertices.size());
    std::vector<std::vector<size_t>> members(source.components.size());
    for (size_t i = 0; i < source.vertices.size(); ++i)
        members[source.vertices[i].componentId].push_back(i);
    for (auto &indices : members)
        std::sort(indices.begin(), indices.end(),
                  [&](size_t a, size_t b) { return source.vertices[a].materialU < source.vertices[b].materialU; });
    for (size_t i = 0; i < destination.vertices.size(); ++i)
    {
        auto &entry = correspondence.entries[i];
        const auto &target = destination.vertices[i];
        for (size_t c = 0; c < members.size(); ++c)
        {
            const auto &component = source.components[c];
            const auto &indices = members[c];
            const auto &first = source.vertices[indices.front()];
            const auto &last = source.vertices[indices.back()];
            double u = target.materialU;
            if (component.closed)
                u -= std::floor(u - first.materialU);
            else
            {
                constexpr double tolerance = 1e-9;
                const double startOffset = target.materialU - component.startU;
                const double endOffset = target.materialU - component.endU;
                const bool atStart = std::abs(startOffset - std::round(startOffset)) <= tolerance;
                const bool atEnd = std::abs(endOffset - std::round(endOffset)) <= tolerance;
                const auto sideMatch = [&](const BellRowVertex &endpoint) {
                    if (endpoint.side.bank == 0 || target.side.bank == 0)
                        return 0;
                    return endpoint.side.bank == target.side.bank && endpoint.side.seamIndex == target.side.seamIndex
                               ? 1
                               : -1;
                };
                if (atStart || atEnd)
                {
                    const int startMatch = atStart ? sideMatch(first) : -1;
                    const int endMatch = atEnd ? sideMatch(last) : -1;
                    if (startMatch < 0 && endMatch < 0)
                        continue;
                    const bool useEnd = endMatch > startMatch ||
                                        (endMatch == startMatch && std::abs(endOffset) < std::abs(startOffset));
                    u = useEnd ? component.endU : component.startU;
                }
                else
                {
                    if (u < component.startU || u > component.endU)
                        u += std::ceil(component.startU - u);
                    if (u < component.startU || u > component.endU)
                        continue;
                }
            }
            if (u <= first.materialU || indices.size() == 1)
                entry.left = indices.front();
            else if (u >= last.materialU && !component.closed)
                entry.left = indices.back();
            else
            {
                size_t right = 0;
                while (right < indices.size() && source.vertices[indices[right]].materialU < u)
                    ++right;
                if (right < indices.size() && source.vertices[indices[right]].materialU == u)
                    entry.left = indices[right];
                else
                {
                    const size_t leftIndex = indices[right - 1];
                    const size_t rightIndex = right == indices.size() ? indices.front() : indices[right];
                    const double leftU = source.vertices[leftIndex].materialU;
                    const double rightU = source.vertices[rightIndex].materialU + (right == indices.size() ? 1.0 : 0.0);
                    const double lambda = (u - leftU) / (rightU - leftU);
                    entry.left = leftIndex;
                    entry.right = rightIndex;
                    entry.lambda = lambda;
                    entry.interpolate = true;
                }
            }
            entry.matched = true;
            break;
        }
    }
    correspondence.valid = true;
}

void BellColliderSolver::applyRowCorrespondence(const std::vector<MVector> &values,
                                                const BellRowCorrespondence &correspondence, std::vector<MVector> &out,
                                                std::vector<bool> &matched)
{
    std::vector<MVector> result(correspondence.entries.size(), MVector(0, 0, 0));
    matched.assign(correspondence.entries.size(), false);
    for (size_t i = 0; i < correspondence.entries.size(); ++i)
    {
        const auto &entry = correspondence.entries[i];
        if (!entry.matched)
            continue;
        result[i] = entry.interpolate ? values[entry.left] * (1.0 - entry.lambda) + values[entry.right] * entry.lambda
                                      : values[entry.left];
        matched[i] = true;
    }
    out.swap(result);
}

MStatus BellColliderSolver::transferDirectField(const BellDirectField &source, const BellRowTopology &destination,
                                                BellRowTransfer &transfer)
{
    if (!validRowTopology(source.topology, source.values.size()) || !finiteRowVectors(source.values) ||
        !validRowTopology(destination, destination.vertices.size()))
        return MS::kInvalidParameter;
    BellRowCorrespondence correspondence;
    transferRowValuesImpl(source.topology, destination, correspondence);
    std::vector<bool> matched;
    applyRowCorrespondence(source.values, correspondence, transfer.values, matched);
    return MS::kSuccess;
}

MStatus BellColliderSolver::buildRowCorrespondence(const BellRowTopology &source, const BellRowTopology &destination,
                                                   BellRowCorrespondence &correspondence)
{
    if (!validRowTopology(source, source.vertices.size()) ||
        !validRowTopology(destination, destination.vertices.size()))
        return MS::kInvalidParameter;
    transferRowValuesImpl(source, destination, correspondence);
    return MS::kSuccess;
}

MStatus BellColliderSolver::transferRowValues(const std::vector<MVector> &values, const BellRowTopology &source,
                                              const BellRowTopology &destination, std::vector<MVector> &out,
                                              std::vector<bool> &matched)
{
    BellRowCorrespondence correspondence;
    const MStatus status = buildRowCorrespondence(source, destination, correspondence);
    if (!status)
        return status;
    return transferRowValues(values, correspondence, out, matched);
}

MStatus BellColliderSolver::transferRowValues(const std::vector<MVector> &values,
                                              const BellRowCorrespondence &correspondence, std::vector<MVector> &out,
                                              std::vector<bool> &matched)
{
    if (!correspondence.valid || values.size() != correspondence.sourceCount || !finiteRowVectors(values))
        return MS::kInvalidParameter;
    std::vector<MVector> result;
    std::vector<bool> resultMatched;
    applyRowCorrespondence(values, correspondence, result, resultMatched);
    if (!finiteRowVectors(result))
        return MS::kInvalidParameter;
    out.swap(result);
    matched.swap(resultMatched);
    return MS::kSuccess;
}

MStatus BellColliderSolver::projectSuspendedRow(const MPointArray &anchorFinal, const MPointArray &anchorBase,
                                                const BellRowTopology &anchorTopology, const MPointArray *aboveFinal,
                                                const BellRowTopology *aboveTopology, const MPointArray &base,
                                                const MPointArray &current, const BellRowTopology &topology,
                                                MPointArray &out)
{
    if ((aboveFinal == nullptr) != (aboveTopology == nullptr))
        return MS::kInvalidParameter;
    BellRowCorrespondence anchor, above;
    MStatus status = buildRowCorrespondence(anchorTopology, topology, anchor);
    if (status && aboveTopology)
        status = buildRowCorrespondence(*aboveTopology, topology, above);
    if (!status)
        return status;
    return projectSuspendedRow(anchorFinal, anchorBase, anchor, aboveFinal, aboveTopology ? &above : nullptr, base,
                               current, out);
}

MStatus BellColliderSolver::projectSuspendedRow(const MPointArray &anchorFinal, const MPointArray &anchorBase,
                                                const BellRowCorrespondence &anchor, const MPointArray *aboveFinal,
                                                const BellRowCorrespondence *above, const MPointArray &base,
                                                const MPointArray &current, MPointArray &out)
{
    if (!anchor.valid || anchorFinal.length() != anchor.sourceCount || anchorBase.length() != anchor.sourceCount ||
        base.length() != anchor.entries.size() || current.length() != base.length() || !finiteRowPoints(anchorFinal) ||
        !finiteRowPoints(anchorBase) || !finiteRowPoints(base) || !finiteRowPoints(current) ||
        (aboveFinal == nullptr) != (above == nullptr))
        return MS::kInvalidParameter;
    if (above && (!above->valid || aboveFinal->length() != above->sourceCount ||
                  above->entries.size() != base.length() || !finiteRowPoints(*aboveFinal)))
        return MS::kInvalidParameter;
    std::vector<MVector> finalValues, baseValues, aboveValues;
    for (unsigned int i = 0; i < anchorFinal.length(); ++i)
    {
        finalValues.push_back(MVector(anchorFinal[i]));
        baseValues.push_back(MVector(anchorBase[i]));
    }
    std::vector<MVector> a, ab, c;
    std::vector<bool> am, abm, cm;
    MStatus status = transferRowValues(finalValues, anchor, a, am);
    if (status)
        status = transferRowValues(baseValues, anchor, ab, abm);
    if (status && above)
    {
        for (unsigned int i = 0; i < aboveFinal->length(); ++i)
            aboveValues.push_back(MVector((*aboveFinal)[i]));
        status = transferRowValues(aboveValues, *above, c, cm);
    }
    if (!status)
        return status;
    MPointArray result = current;
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        if (!am[i] || !abm[i])
            continue;
        const double ell = (MVector(base[i]) - ab[i]).length();
        const MVector offset = MVector(current[i]) - a[i];
        const double d = offset.length();
        if (!std::isfinite(ell) || !std::isfinite(d))
            return MS::kInvalidParameter;
        if (d <= ell)
            continue;
        if (ell == 0.0)
        {
            result[i] = MPoint(a[i]);
            continue;
        }
        MVector e = offset / d;
        if (above && cm[i] && ell > 1e-12)
        {
            const MVector incoming = a[i] - c[i];
            const double length = incoming.length();
            if (!std::isfinite(length))
                return MS::kInvalidParameter;
            if (length > 1e-12)
            {
                const double x = (d - ell) / ell;
                const double t = (std::max)(0.0, (std::min)(1.0, x));
                const double k = 0.5 * t * t * (3.0 - 2.0 * t);
                const MVector q = e * (1.0 - k) + (incoming / length) * k;
                const double qLength = q.length();
                if (qLength > 1e-12)
                    e = q / qLength;
            }
        }
        result[i] = MPoint(a[i] + e * ell);
    }
    if (!finiteRowPoints(result))
        return MS::kInvalidParameter;
    out = result;
    return MS::kSuccess;
}

MStatus BellColliderSolver::smoothRowDisplacements(const MPointArray &upBase, const MPointArray &upCurrent,
                                                   const BellRowTopology &upTopology, const MPointArray &base,
                                                   const MPointArray &current, const BellRowTopology &topology,
                                                   const MPointArray &downBase, const MPointArray &downCurrent,
                                                   const BellRowTopology &downTopology, MPointArray &out)
{
    BellRowCorrespondence up, down;
    MStatus status = buildRowCorrespondence(upTopology, topology, up);
    if (status)
        status = buildRowCorrespondence(downTopology, topology, down);
    if (!status)
        return status;
    return smoothRowDisplacements(upBase, upCurrent, up, base, current, downBase, downCurrent, down, out);
}

MStatus BellColliderSolver::smoothRowDisplacements(const MPointArray &upBase, const MPointArray &upCurrent,
                                                   const BellRowCorrespondence &up, const MPointArray &base,
                                                   const MPointArray &current, const MPointArray &downBase,
                                                   const MPointArray &downCurrent, const BellRowCorrespondence &down,
                                                   MPointArray &out)
{
    if (!up.valid || !down.valid || upBase.length() != up.sourceCount || upCurrent.length() != up.sourceCount ||
        downBase.length() != down.sourceCount || downCurrent.length() != down.sourceCount ||
        base.length() != up.entries.size() || base.length() != down.entries.size() ||
        current.length() != base.length() || !finiteRowPoints(upBase) || !finiteRowPoints(upCurrent) ||
        !finiteRowPoints(base) || !finiteRowPoints(current) || !finiteRowPoints(downBase) ||
        !finiteRowPoints(downCurrent))
        return MS::kInvalidParameter;
    std::vector<MVector> upValues, downValues, u, d;
    for (unsigned int i = 0; i < upBase.length(); ++i)
        upValues.push_back(upCurrent[i] - upBase[i]);
    for (unsigned int i = 0; i < downBase.length(); ++i)
        downValues.push_back(downCurrent[i] - downBase[i]);
    std::vector<bool> um, dm;
    MStatus status = transferRowValues(upValues, up, u, um);
    if (status)
        status = transferRowValues(downValues, down, d, dm);
    if (!status)
        return status;
    MPointArray result = current;
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        if (!um[i] || !dm[i])
            continue;
        const MVector displacement = current[i] - base[i];
        const MVector next = displacement + ((u[i] + d[i]) * 0.5 - displacement) * 0.5;
        if (next.x != 0.0 || next.y != 0.0 || next.z != 0.0)
            result[i] = base[i] + next;
    }
    if (!finiteRowPoints(result))
        return MS::kInvalidParameter;
    out = result;
    return MS::kSuccess;
}
