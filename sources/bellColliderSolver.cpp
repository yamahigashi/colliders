#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <cmath>

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

namespace
{
template <class Ops>
BellColliderRelax::Vector<Ops> broadcast(double x, double y, double z)
{
    return {Ops::splat(x), Ops::splat(y), Ops::splat(z)};
}

template <class Ops>
BellColliderRelax::Ring<Ops> prepareRelaxRing(const PreparedBellRing &ring, double collision)
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

void BellColliderSolver::relaxTowardRingBoundary(MPointArray &points, const PreparedBellRing &ring, double collision,
                                                 int startIndex, int count)
{
    if (!(collision > 1e-5))
        return;

#ifdef YDD_RELAX_SSE2
    const auto prepared = prepareRelaxRing<BellColliderRelax::Pair>(ring, collision);
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
    const auto prepared = prepareRelaxRing<BellColliderRelax::Scalar>(ring, collision);
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

bool BellColliderSolver::collisionPoints(const MMatrix &bellMatrix, const MMatrix &bellMatrixInverse,
                                         const Plane &bellPlane, const PreparedBellRing &ring,
                                         MPoint &collisionPointBell, MPoint &collisionPointRing, MPoint &linePoint)
{
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
                                (int)bellPoints.length() - bellSubdivision - 1);

        bellPointsList.push_back(bellPoints);
    }
}

void BellColliderSolver::averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const vector<MPointArray>& bellPointsList, MPointArray& outBellPoints, bool useUnnormalized)
{
    outBellPoints = baseBellPoints;

    for (size_t i = bellSubdivision + 1; i < baseBellPoints.length(); i++)
    {
        double sum = 0;
        double maxDist = 0;
        for (const auto& bellPoints : bellPointsList)
        {
            const MVector vec = bellPoints[i] - baseBellPoints[i];
            const double d = vec.length();
            sum += pow(d, 2);

            if (d > maxDist)
                maxDist = d;
        }

        if (sum > 0)
        {
            MVector wp;
            for (const auto& bellPoints : bellPointsList)
            {
                const MVector vec = bellPoints[i] - baseBellPoints[i];
                const double d = pow(vec.length(), 2);
                const double w = d / sum;

                wp += vec * w;
            }

            if (useUnnormalized)
            {
                outBellPoints[i] += wp;
            }
            else
            {
                double wp_len = wp.length();
                if (wp_len > 1e-5)
                {
                    outBellPoints[i] += wp.normal() * maxDist;
                }
            }
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
    averageDisplacements(bellSubdivision, baseBellPoints, bellPointsList, outBellPoints, gate);

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

        const double alpha = inputs.smoothness * 0.5;
        if (alpha != 0.0)
        {
            vector<MVector> smoothedDisplacements(bellSubdivision);
            for (int iteration = 0; iteration < 3; iteration++)
            {
                for (int i = 0; i < bellSubdivision; i++)
                {
                    const int previous = (i + bellSubdivision - 1) % bellSubdivision;
                    const int next = (i + 1) % bellSubdivision;
                    smoothedDisplacements[i] = displacements[i] * (1.0 - alpha)
                        + (displacements[previous] + displacements[next]) * (alpha * 0.5);
                }
                displacements.swap(smoothedDisplacements);
            }
        }

        for (int i = 0; i < bellSubdivision; i++)
            outBellPoints.set(baseBellPoints[startIndex + i] + displacements[i], startIndex + i);

        for (const auto &ring : inputs.rings)
            relaxTowardRingBoundary(outBellPoints, ring, inputs.collision, startIndex, bellSubdivision);
    }

    outputs.points = outBellPoints;
    return MS::kSuccess;
}
