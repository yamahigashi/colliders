#include "bellColliderSolver.h"
#include "bellColliderRelaxKernel.h"

#include <maya/MPointArray.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>

static bool sameBits(double a, double b)
{
    std::uint64_t ua = 0, ub = 0;
    std::memcpy(&ua, &a, sizeof(ua));
    std::memcpy(&ub, &b, sizeof(ub));
    return ua == ub;
}

static void referenceRelax(MPointArray &points, const PreparedBellRing &ring, double collision, int startIndex, int count,
                           bool capAtRingOrigin)
{
    if (!(collision > 1e-5))
        return;
    const MMatrix &ringMatrixInverse = ring.inverse;
    const MPoint &ringTranslate = ring.translation;
    const Plane &ringPlane = ring.plane;
    for (int j = startIndex; j < startIndex + count; ++j)
    {
        if (capAtRingOrigin && !(ringPlane.distance(points[j]) > 0.0))
            continue;
        const MVector vec = ringPlane.projectPoint(points[j]) - ringTranslate;
        const double vecLength = vec.length();
        if (vecLength > 1e-5)
        {
            const double localLength = (vec * ringMatrixInverse).length();
            const double delta = localLength > 1e-5 ? vecLength / localLength : 1.0;
            const MVector vecProjectedScaled = vec.normal() * delta;
            if (vecProjectedScaled.length() > vecLength)
                points[j] += vec.normal() * (vecProjectedScaled.length() - vecLength) * collision;
        }
    }
}

static MPointArray makePoints(int count, int seed)
{
    MPointArray points;
    std::uint32_t state = static_cast<std::uint32_t>(seed);
    const auto next = [&state]() {
        state = state * 1664525u + 1013904223u;
        return static_cast<double>(state) / 4294967296.0;
    };
    for (int i = 0; i != count; ++i)
    {
        const double x = -4.0 + 8.0 * next();
        const double y = -3.0 + 6.0 * next();
        const double z = -5.0 + 10.0 * next();
        const double w = i % 29 == 0 ? 0.0 : (i % 31 == 0 ? 2.0 : (i % 37 == 0 ? -1.0 : 1.0));
        points.append(MPoint(x, y, z, w));
    }
    return points;
}

static bool compare(const MPointArray &expected, const MPointArray &actual, int start, int count)
{
    for (unsigned int i = 0; i != expected.length(); ++i)
    {
        const MPoint &a = expected[i];
        const MPoint &b = actual[i];
        if (!sameBits(a.x, b.x) || !sameBits(a.y, b.y) || !sameBits(a.z, b.z) || !sameBits(a.w, b.w))
        {
            std::cerr << "mismatch index=" << i << " range=" << start << ":" << count << "\n";
            return false;
        }
    }
    return true;
}

static void scalarRelax(MPointArray &points, const PreparedBellRing &ring, double collision, int start, int count,
                        bool capAtRingOrigin)
{
    if (!(collision > 1e-5))
        return;
    using Ops = BellColliderRelax::Scalar;
    BellColliderRelax::Ring<Ops> prepared;
    MPoint origin = ring.plane.orig;
    MPoint translation = ring.translation;
    if (origin.w != 1.0)
        origin.cartesianize();
    if (translation.w != 1.0)
        translation.cartesianize();
    prepared.origin = {origin.x, origin.y, origin.z};
    prepared.translation = {translation.x, translation.y, translation.z};
    prepared.normal = {ring.plane.normal.x, ring.plane.normal.y, ring.plane.normal.z};
    for (int column = 0; column < 3; ++column)
        prepared.inverseColumns[column] = {ring.inverse[0][column], ring.inverse[1][column], ring.inverse[2][column]};
    prepared.collision = collision;
    prepared.capAtRingOrigin = capAtRingOrigin;
    for (int j = start; j < start + count; ++j)
    {
        MPoint &point = points[j];
        MPoint cartesian = point;
        if (cartesian.w != 1.0)
            cartesian.cartesianize();
        const BellColliderRelax::Vector<Ops> raw = {point.x, point.y, point.z};
        const BellColliderRelax::Vector<Ops> input = {cartesian.x, cartesian.y, cartesian.z};
        const auto output = BellColliderRelax::relax(raw, input, prepared);
        point.x = output.x;
        point.y = output.y;
        point.z = output.z;
    }
}

static bool runCase(const PreparedBellRing &ring, const MPointArray &source, double collision, int start, int count)
{
    for (bool capAtRingOrigin : {false, true})
    {
        MPointArray expected = source;
        MPointArray actual = source;
        MPointArray scalar = source;
        referenceRelax(expected, ring, collision, start, count, capAtRingOrigin);
        BellColliderSolver::relaxTowardRingBoundary(actual, ring, collision, start, count, capAtRingOrigin);
        scalarRelax(scalar, ring, collision, start, count, capAtRingOrigin);
        if (!compare(expected, actual, start, count) || !compare(expected, scalar, start, count))
        {
            std::cerr << "capAtRingOrigin=" << capAtRingOrigin << " collision=" << collision << "\n";
            return false;
        }
    }
    return true;
}

static bool runSequentialCase(const PreparedBellRing &first, const PreparedBellRing &second, const MPointArray &source,
                              double collision)
{
    const int count = static_cast<int>(source.length());
    for (bool capAtRingOrigin : {false, true})
    {
        MPointArray expected = source;
        MPointArray actual = source;
        MPointArray scalar = source;
        referenceRelax(expected, first, collision, 0, count, capAtRingOrigin);
        referenceRelax(expected, second, collision, 0, count, capAtRingOrigin);
        BellColliderSolver::relaxTowardRingBoundary(actual, first, collision, 0, count, capAtRingOrigin);
        BellColliderSolver::relaxTowardRingBoundary(actual, second, collision, 0, count, capAtRingOrigin);
        scalarRelax(scalar, first, collision, 0, count, capAtRingOrigin);
        scalarRelax(scalar, second, collision, 0, count, capAtRingOrigin);
        if (!compare(expected, actual, 0, count) || !compare(expected, scalar, 0, count))
            return false;
    }
    return true;
}

static bool checkCapBoundary(const PreparedBellRing &ring)
{
    MPointArray source;
    source.append(MPoint(0.5, -1.0, 0.0));
    source.append(MPoint(0.5, 1.0, 0.0));
    source.append(MPoint(0.5, 0.0, 0.0));
    source.append(MPoint(0.5, -0.0, 0.0));
    source.append(MPoint(0.5, -1e-12, 0.0));
    source.append(MPoint(0.5, 1e-12, 0.0));
    if (!runCase(ring, source, 1.0, 0, static_cast<int>(source.length())))
        return false;
    for (bool capAtRingOrigin : {false, true})
    {
        MPointArray actual = source;
        MPointArray expected = source;
        for (unsigned int j = 0; j < expected.length(); ++j)
            if (!capAtRingOrigin || source[j].y > 0.0)
                // The kernel adds a signed-zero y displacement, so keep the source's zero sign bitwise.
                expected[j] += MVector(0.5, -0.0, 0.0);
        BellColliderSolver::relaxTowardRingBoundary(actual, ring, 1.0, 0,
                                                   static_cast<int>(actual.length()), capAtRingOrigin);
        if (!compare(expected, actual, 0, static_cast<int>(actual.length())))
            return false;
    }
    MPointArray defaults = source;
    MPointArray uncapped = source;
    BellColliderSolver::relaxTowardRingBoundary(defaults, ring, 1.0, 0, static_cast<int>(source.length()));
    BellColliderSolver::relaxTowardRingBoundary(uncapped, ring, 1.0, 0, static_cast<int>(source.length()), false);
    return !BellColliderInputs().capAtRingOrigin && compare(uncapped, defaults, 0, static_cast<int>(source.length()));
}

int main()
{
    const double rawA[4][4] = {
        {1.1, 0.2, -0.1, 0.0}, {-0.3, 0.9, 0.25, 0.0}, {0.15, -0.2, 1.3, 0.0}, {0.4, -0.7, 0.8, 1.0}};
    const double rawB[4][4] = {
        {0.8, -0.4, 0.3, 0.0}, {0.25, 1.4, -0.2, 0.0}, {-0.1, 0.35, 0.7, 0.0}, {-1.0, 0.3, -0.5, 1.0}};
    const PreparedBellRing rings[] = {PreparedBellRing(MMatrix(rawA)), PreparedBellRing(MMatrix(rawB))};
    const MPointArray points = makePoints(1207, 0x51a7);
    const double collisions[] = {0.0, 1e-5, 1.00001e-5, 0.15, 0.65, 1.0};
    const int starts[] = {0, 1, 7, 31, 600};
    const int counts[] = {1, 2, 17, 599, 607};

    for (const PreparedBellRing &ring : rings)
        for (double collision : collisions)
            for (int start : starts)
                for (int count : counts)
                    if (start + count <= static_cast<int>(points.length()) &&
                        !runCase(ring, points, collision, start, count))
                        return 1;

    const double identityRaw[4][4] = {
        {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
    const PreparedBellRing identity{MMatrix(identityRaw)};
    if (!checkCapBoundary(identity))
        return 1;

    // Tangential offsets exercise the radial threshold after plane projection.
    MPointArray threshold;
    threshold.append(MPoint(std::nextafter(1e-5, 0.0), 0.0, 0.0));
    threshold.append(MPoint(std::nextafter(1e-5, std::numeric_limits<double>::infinity()), 0.0, 0.0));
    if (!runCase(identity, threshold, 0.8, 0, 2))
        return 1;

    // Inside/outside, signed zero, and an odd count exercise both SIMD lanes
    // with one inactive lane and a final unpaired point.
    MPointArray radial;
    radial.append(MPoint(0.5, 0.0, 0.0));
    radial.append(MPoint(1.5, 0.0, 0.0));
    radial.append(MPoint(-0.0, 0.0, 0.0));
    if (!runCase(identity, radial, 0.8, 0, 3))
        return 1;

    // Ring application is ordered: the second ring must consume the first
    // ring's output, exactly as the solver's production loop does.
    if (!runSequentialCase(identity, rings[0], radial, 0.8))
        return 1;

    // A zero collision value must preserve every input bit, including w and
    // signed zero, before any coordinate conversion.
    MPointArray signedZero;
    signedZero.append(MPoint(-0.0, 1.0, 0.0, -0.0));
    if (!runCase(identity, signedZero, 0.0, 0, 1))
        return 1;
    // Keep vector length active while crossing the separate local-length gate.
    for (double value : {std::nextafter(2e-5, 0.0), std::nextafter(2e-5, std::numeric_limits<double>::infinity())})
    {
        PreparedBellRing scaled = identity;
        scaled.inverse[0][0] = value;
        if (!runCase(scaled, radial, 0.8, 0, 3))
            return 1;
    }

    MPointArray exceptional;
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    exceptional.append(MPoint(nan, 0.0, 0.0));
    exceptional.append(MPoint(0.5, 0.0, 0.0));
    exceptional.append(MPoint(infinity, 0.0, 0.0));
    exceptional.append(MPoint(0.5, 0.0, 0.0, nan));
    exceptional.append(MPoint(0.5, 0.0, 0.0, infinity));
    exceptional.append(MPoint(-0.0, 0.0, 0.0));
    for (double collision : {0.65, infinity, nan})
        if (!runCase(identity, exceptional, collision, 0, static_cast<int>(exceptional.length())))
            return 1;
    std::cout << "relax solver Maya reference comparison passed\n";
    return 0;
}
