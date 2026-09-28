#include "bellColliderSolver.h"
#include "bellColliderRelaxKernel.h"
#include "skirtLegProfile.h"

#include <maya/MPointArray.h>

#include <algorithm>
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

static std::vector<MVector> referenceDirections(const MPointArray &points, const MPointArray &base)
{
    std::vector<MVector> result(points.length());
    for (unsigned int i = 0; i < points.length(); ++i)
    {
        MPoint p = points[i], b = base[i];
        if (p.w != 1.0)
            p.cartesianize();
        if (b.w != 1.0)
            b.cartesianize();
        result[i] = p - b;
    }
    return result;
}

static void referenceDirectionalRelax(MPointArray &points, const PreparedBellRing &ring, double collision, int start,
                                      int count, bool cap, const std::vector<MVector> &directions,
                                      const std::vector<double> &scales)
{
    if (!(collision > 1e-5))
        return;
    MPoint origin = ring.plane.orig, translation = ring.translation;
    if (origin.w != 1.0)
        origin.cartesianize();
    if (translation.w != 1.0)
        translation.cartesianize();
    const MVector n = ring.plane.normal;
    const auto dot = [](const MVector &a, const MVector &b) { return (a.x * b.x + a.y * b.y) + a.z * b.z; };
    const auto length = [&](const MVector &v) { return std::sqrt(dot(v, v)); };
    const auto local = [&](const MVector &v) {
        return MVector((v.x * ring.inverse[0][0] + v.y * ring.inverse[1][0]) + v.z * ring.inverse[2][0],
                       (v.x * ring.inverse[0][1] + v.y * ring.inverse[1][1]) + v.z * ring.inverse[2][1],
                       (v.x * ring.inverse[0][2] + v.y * ring.inverse[1][2]) + v.z * ring.inverse[2][2]);
    };
    for (int i = start; i < start + count; ++i)
    {
        MPoint p = points[i];
        if (p.w != 1.0)
            p.cartesianize();
        const double distance = dot(p - origin, n);
        const MVector r(p.x - n.x * distance - translation.x, p.y - n.y * distance - translation.y,
                        p.z - n.z * distance - translation.z);
        const MVector pl = local(r);
        const double c2 = pl.x * pl.x + pl.z * pl.z - 1.0;
        if (!(c2 < 0.0) || (cap && !(distance > 0.0)))
            continue;
        double g_d = 1.0;
        if (ring.distalEnd)
        {
            const bool validWidth = ring.distalWidth > 1e-5;
            const double x = (validWidth ? distance - ring.distalLength : 0.0) /
                             (validWidth ? ring.distalWidth : 1.0);
            const double t_d = x > 1.0 ? 1.0 : (x > 0.0 ? x : 0.0);
            g_d = 1.0 - ((t_d * t_d) * (3.0 - 2.0 * t_d));
            if ((validWidth && distance > ring.distalLength + ring.distalWidth) || !(g_d > 0.0))
                continue;
        }
        const MVector perpendicular = directions[i] - n * dot(directions[i], n);
        const double D = length(perpendicular);
        MVector centre(0, 0, 0);
        double radius = 1.0, ramp = 0.0;
        bool shifted = false;
        if (D > 1e-5)
        {
            const MVector unit = perpendicular * (1.0 / D);
            const MVector m = local(unit * -1.0);
            const double norm = m.x * m.x + m.z * m.z;
            if (norm > 0.0 && std::isfinite(norm))
            {
                radius = 1.0 / std::sqrt(norm);
                shifted = true;
                const double t = (std::max)(0.0, (std::min)(1.0, D / radius));
                const double beta = (0.5 * (t * t) * (3.0 - 2.0 * t)) * (scales.empty() ? 1.0 : scales[i]);
                ramp = t * t * (3.0 - 2.0 * t);
                centre = unit * ((0.0 - beta) * radius);
            }
        }
        const MVector v = r - centre;
        const double vLength = length(v);
        if (!(vLength > 1e-5))
            continue;
        const MVector q = v * (1.0 / vLength);
        const MVector ql = local(q);
        const double a = ql.x * ql.x + ql.z * ql.z;
        const double b = 2.0 * (pl.x * ql.x + pl.z * ql.z);
        if (!(a > 1e-12))
            continue;
        const double disc = (std::max)(0.0, b * b - (4.0 * a) * c2);
        const double s = ((0.0 - b) + std::sqrt(disc)) / (2.0 * a);
        if (!(s > 0.0))
            continue;
        double g = 1.0;
        if (shifted)
        {
            const double u = (std::max)(0.0, (std::min)(1.0, vLength / (0.75 * radius)));
            g = 1.0 - ramp * (1.0 - u * u * (3.0 - 2.0 * u));
        }
        MVector displacement = ((q * s) * g) * collision;
        if (ring.distalEnd)
            displacement *= g_d;
        points[i].x += displacement.x;
        points[i].y += displacement.y;
        points[i].z += displacement.z;
    }
}

static bool directionalCompare(const MPointArray &expected, const MPointArray &actual, double tolerance = 1e-12)
{
    if (expected.length() != actual.length())
        return false;
    for (unsigned int i = 0; i < expected.length(); ++i)
    {
        const MPoint &a = expected[i], &b = actual[i];
        const auto close = [tolerance](double x, double y) {
            return sameBits(x, y) || (std::isfinite(x) && std::isfinite(y) && std::abs(x - y) <= tolerance);
        };
        if (!close(a.x, b.x) || !close(a.y, b.y) || !close(a.z, b.z) || !sameBits(a.w, b.w))
        {
            std::cerr << "directional mismatch index=" << i << "\n";
            return false;
        }
    }
    return true;
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

static void scalarDirectionalRelax(MPointArray &points, const PreparedBellRing &ring, double collision, int start,
                                   int count, bool capAtRingOrigin, const std::vector<MVector> &directions,
                                   const std::vector<double> &scales)
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
    prepared.distalEnd = ring.distalEnd;
    prepared.distalLength = Ops::splat(ring.distalLength);
    prepared.distalWidth = Ops::splat(ring.distalWidth);
    for (int j = start; j < start + count; ++j)
    {
        MPoint &point = points[j];
        MPoint cartesian = point;
        if (cartesian.w != 1.0)
            cartesian.cartesianize();
        const BellColliderRelax::Vector<Ops> raw = {point.x, point.y, point.z};
        const BellColliderRelax::Vector<Ops> input = {cartesian.x, cartesian.y, cartesian.z};
        const MVector &d = directions[j];
        const BellColliderRelax::Vector<Ops> direction = {d.x, d.y, d.z};
        const auto output = BellColliderRelax::relax(raw, input, prepared, direction, scales.empty() ? 1.0 : scales[j]);
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

static PreparedBellRing directionalCircle(double radius)
{
    MMatrix matrix;
    matrix[0][0] = matrix[2][2] = radius;
    return PreparedBellRing(matrix);
}

static bool directionalCase(const PreparedBellRing &ring, const MPointArray &source,
                            const std::vector<MVector> &directions, const std::vector<double> &scales, double collision,
                            int start, int count, bool cap)
{
    auto expected = source, actual = source, scalar = source;
    referenceDirectionalRelax(expected, ring, collision, start, count, cap, directions, scales);
    BellColliderSolver::relaxTowardRingBoundary(actual, ring, collision, start, count, cap, directions, scales);
    scalarDirectionalRelax(scalar, ring, collision, start, count, cap, directions, scales);
    if (!directionalCompare(expected, actual) || !compare(scalar, actual, start, count))
        return false;
    for (unsigned int i = 0; i < source.length(); ++i)
    {
        if (!sameBits(source[i].w, actual[i].w))
            return false;
        if (sameBits(source[i].x, expected[i].x) && sameBits(source[i].y, expected[i].y) &&
            sameBits(source[i].z, expected[i].z))
            if (!sameBits(source[i].x, actual[i].x) || !sameBits(source[i].y, actual[i].y) ||
                !sameBits(source[i].z, actual[i].z))
                return false;
    }
    return true;
}

static bool directionalReference()
{
    MMatrix ellipse;
    const double angle = 0.37, c = std::cos(angle), s = std::sin(angle);
    ellipse[0][0] = 9.61 * c;
    ellipse[0][1] = 9.61 * s;
    ellipse[1][0] = -s;
    ellipse[1][1] = c;
    ellipse[2][2] = 6.06;
    ellipse[3][0] = 0.4;
    ellipse[3][1] = -0.7;
    ellipse[3][2] = 0.8;
    MMatrix shear = ellipse;
    shear[0][1] += 0.5;
    const auto source = makePoints(127, 0x51a7);
    auto base = source;
    for (unsigned int i = 0; i < base.length(); ++i)
        base[i].x -= 3.0;
    const auto directions = referenceDirections(source, base);
    std::vector<double> scales(source.length());
    for (unsigned int i = 0; i < source.length(); ++i)
        scales[i] = (i % 5) * 0.25;
    for (const auto &ring : {PreparedBellRing(ellipse), PreparedBellRing(shear), directionalCircle(6.15)})
        for (double collision : {0.0, 1e-5, 1.00001e-5, 0.15, 0.65, 1.0})
            for (int start : {0, 1, 7})
                for (int count : {1, 2, 17, 119})
                    for (bool cap : {false, true})
                        if (!directionalCase(ring, source, directions, scales, collision, start, count, cap))
                            return false;
    return true;
}

static bool directionalZero()
{
    for (double radius : {1.0, 2.0, 6.15})
    {
        const auto ring = directionalCircle(radius);
        MPointArray source;
        for (double x : {0.1, 0.5, 0.9})
            source.append(MPoint(radius * x, 1, radius * 0.1));
        const std::vector<MVector> zero(source.length(), MVector(0, 0, 0));
        auto old = source, actual = source;
        BellColliderSolver::relaxTowardRingBoundary(old, ring, 1, 0, 3, true);
        BellColliderSolver::relaxTowardRingBoundary(actual, ring, 1, 0, 3, true, zero, {});
        if (!directionalCompare(old, actual))
            return false;
    }
    // The local cross-section equation also applies below the old local-length fallback.
    for (double radius : {2.0, 1e5})
    {
        MPointArray point, expected;
        point.append(MPoint(1.5e-5, 1, 0));
        expected.append(MPoint(radius, 1, 0));
        BellColliderSolver::relaxTowardRingBoundary(point, directionalCircle(radius), 1, 0, 1, false,
                                                    {MVector(0, 0, 0)}, {});
        if (!directionalCompare(expected, point, 1e-9))
            return false;
    }
    MMatrix shear;
    shear[0][1] = 2;
    PreparedBellRing ring(shear);
    MPointArray point, expected;
    point.append(MPoint(0.8, 1, 0));
    expected.append(MPoint(1, 1, 0));
    // Local y does not contribute to the elliptical cross-section.
    BellColliderSolver::relaxTowardRingBoundary(point, ring, 1, 0, 1, false, {MVector(0, 0, 0)}, {});
    return directionalCompare(expected, point);
}

static bool directionalAxisAndCentre()
{
    MPointArray source, expected;
    for (double x : {0.0, -3.075, 3.0})
        source.append(MPoint(x, 1, 0));
    const std::vector<MVector> directions(3, MVector(12.3, 0, 0));
    expected = source;
    referenceDirectionalRelax(expected, directionalCircle(6.15), 1, 0, 3, false, directions, {});
    auto actual = source;
    BellColliderSolver::relaxTowardRingBoundary(actual, directionalCircle(6.15), 1, 0, 3, false, directions, {});
    if (!directionalCompare(expected, actual, 1e-9))
        return false;
    for (const MVector d : {MVector(0, 0, 0), MVector(0, 12.3, 0)})
        if (!directionalCase(directionalCircle(6.15), source, std::vector<MVector>(3, d), {}, 1, 0, 3, false))
            return false;
    MPointArray ramp, rampExpected;
    ramp.append(MPoint(0, 1, 1));
    rampExpected = ramp;
    referenceDirectionalRelax(rampExpected, directionalCircle(2), 1, 0, 1, false, {MVector(1, 0, 0)}, {});
    BellColliderSolver::relaxTowardRingBoundary(ramp, directionalCircle(2), 1, 0, 1, false, {MVector(1, 0, 0)}, {});
    return directionalCompare(rampExpected, ramp);
}

static bool directionalEntryDepth()
{
    MPointArray source, expected;
    for (double h : {0.5, 0.1, 0.01})
        source.append(MPoint(-6.15 + h, 1, 0));
    expected = source;
    referenceDirectionalRelax(expected, directionalCircle(6.15), 1, 0, 3, false,
                              std::vector<MVector>(3, MVector(12.3, 0, 0)), {});
    auto actual = source;
    BellColliderSolver::relaxTowardRingBoundary(actual, directionalCircle(6.15), 1, 0, 3, false,
                                                std::vector<MVector>(3, MVector(12.3, 0, 0)), {});
    return directionalCompare(expected, actual, 1e-9);
}

static bool directionalBoundaryAndCap()
{
    for (double distance : {-1e-12, -0.0, 0.0, 1e-12})
        for (double collision : {0.0, 1e-5, 1.00001e-5, 0.15, 0.65, 1.0})
        {
            MMatrix matrix;
            matrix[0][0] = matrix[2][2] = 6.15;
            matrix[3][1] = 1 - distance;
            const PreparedBellRing ring(matrix);
            MPointArray source, expected;
            source.append(MPoint(0, 1, 0));
            source.append(MPoint(6.15, 1, 0));
            source.append(MPoint(7, 1, 0));
            expected = source;
            referenceDirectionalRelax(expected, ring, collision, 0, 3, true,
                                      std::vector<MVector>(3, MVector(12.3, 0, 0)), {});
            auto actual = source;
            const std::vector<MVector> directions(3, MVector(12.3, 0, 0));
            BellColliderSolver::relaxTowardRingBoundary(actual, ring, collision, 0, 3, true, directions, {});
            if (!directionalCompare(expected, actual, 1e-9) ||
                !directionalCase(ring, source, directions, {}, collision, 0, 3, true))
                return false;
        }
    MPointArray raw;
    raw.append(MPoint(-0.0, 1, -0.0, -0.0));
    const std::vector<MVector> directions(1, MVector(12.3, 0, 0));
    for (double collision : {0.0, 1e-5, std::numeric_limits<double>::quiet_NaN()})
        if (!directionalCase(directionalCircle(6.15), raw, directions, {}, collision, 0, 1, false))
            return false;
    return true;
}

static bool directionalFade()
{
    const auto ring = directionalCircle(2);
    const double above = std::nextafter(1.5, 2.0);
    for (double collision : {0.15, 0.65, 1.0})
    {
        MPointArray source, expected;
        source.append(MPoint(-1, 1, 0));
        source.append(MPoint(above - 1, 1, 0));
        source.append(MPoint(-0.25, 1, 0));
        source.append(MPoint(-0.25, 1, 0));
        source.append(MPoint(-0.25, 1, 0));
        const std::vector<MVector> directions = {MVector(4, 0, 0), MVector(4, 0, 0), MVector(4, 0, 0),
                                                  MVector(0, 0, 0), MVector(4, 0, 0)};
        expected.append(source[0]);
        const double q = above * (1.0 / std::sqrt(above * above));
        const double px = (above - 1) * 0.5, qx = q * 0.5;
        const double a = qx * qx, b = 2 * (px * qx), c2 = px * px - 1;
        const double s = ((0.0 - b) + std::sqrt(b * b - (4 * a) * c2)) / (2 * a);
        expected.append(MPoint(above - 1 + (q * s) * collision, 1, 0));
        expected.append(MPoint(-0.25 + 2.25 * 0.5 * collision, 1, 0));
        expected.append(MPoint(-0.25 - 1.75 * collision, 1, 0));
        expected.append(expected[2]);
        auto actual = source;
        BellColliderSolver::relaxTowardRingBoundary(actual, ring, collision, 0, 5, false, directions, {});
        if (!directionalCompare(expected, actual) ||
            !sameBits(expected[1].x, actual[1].x) ||
            !directionalCase(ring, source, directions, {}, collision, 0, 5, false))
            return false;
    }
    return true;
}

static bool directionalPair()
{
    const double infinity = std::numeric_limits<double>::infinity();
    const auto ring = directionalCircle(1);
    const auto paired = [&](const PreparedBellRing &current, const MPoint &stop, const MPoint &active,
                            const MVector &stopD, const MVector &activeD) {
        for (bool reverse : {false, true})
            for (bool cap : {false, true})
            {
                MPointArray source;
                source.append(MPoint(9, 9, 9));
                source.append(reverse ? active : stop);
                source.append(reverse ? stop : active);
                source.append(stop);
                const std::vector<MVector> directions = {MVector(0, 0, 0), reverse ? activeD : stopD,
                                                         reverse ? stopD : activeD, stopD};
                if (!directionalCase(current, source, directions, {1, 1, 1, 1}, 1, 1, 3, cap))
                    return false;
                auto original = source, changed = source;
                changed[1] = MPoint(0.4, 2, 0.3);
                BellColliderSolver::relaxTowardRingBoundary(original, current, 1, 1, 3, cap, directions, {});
                BellColliderSolver::relaxTowardRingBoundary(changed, current, 1, 1, 3, cap, directions, {});
                for (int i : {2, 3})
                    if (!sameBits(original[i].x, changed[i].x) || !sameBits(original[i].y, changed[i].y) ||
                        !sameBits(original[i].z, changed[i].z) || !sameBits(original[i].w, changed[i].w))
                        return false;
            }
        return true;
    };
    for (double threshold : {std::nextafter(1e-5, 0.0), std::nextafter(1e-5, infinity)})
    {
        if (!paired(ring, MPoint(0.5, 1, 0), MPoint(0.5, 1, 0), MVector(threshold, 0, 0), MVector(2, 0, 0)) ||
            !paired(ring, MPoint(threshold, 1, 0), MPoint(0.5, 1, 0), MVector(0, 0, 0), MVector(0, 0, 0)))
            return false;
    }
    for (double coefficient : {std::nextafter(1e-6, 0.0), std::nextafter(1e-6, infinity)})
    {
        auto large = ring;
        large.inverse[0][0] = coefficient;
        if (!paired(large, MPoint(0.5, 1, 0), MPoint(0, 1, 0.5), MVector(0, 0, 0), MVector(0, 0, 0)))
            return false;
    }
    double lo = 0, hi = 1;
    for (int i = 0; i < 80; ++i)
    {
        const double mid = (lo + hi) * 0.5;
        if (0.5 * (mid * mid) * (3 - 2 * mid) <= 1e-5)
            lo = mid;
        else
            hi = mid;
    }
    lo = std::nextafter(lo, 0.0);
    hi = std::nextafter(hi, infinity);
    for (double D : {lo, hi})
    {
        MPointArray actual, expected;
        actual.append(MPoint(0, 1, 0));
        expected = actual;
        referenceDirectionalRelax(expected, ring, 1, 0, 1, false, {MVector(D, 0, 0)}, {});
        BellColliderSolver::relaxTowardRingBoundary(actual, ring, 1, 0, 1, false, {MVector(D, 0, 0)}, {});
        if (!directionalCompare(expected, actual))
            return false;
    }
    if (!paired(ring, MPoint(0, 1, 0), MPoint(0, 1, 0), MVector(lo, 0, 0), MVector(hi, 0, 0)))
        return false;
    for (const MPoint p : {MPoint(2, 1, 0), MPoint(0.5, -1, 0), MPoint(-0.0, -0.0, 0), MPoint(infinity, 1, 0),
                           MPoint(std::numeric_limits<double>::quiet_NaN(), 1, 0)})
        if (!paired(ring, p, MPoint(0.5, 1, 0), MVector(2, 0, 0), MVector(2, 0, 0)))
            return false;
    for (double coefficient : {0.0, 1e308, infinity})
    {
        auto degenerate = ring;
        degenerate.inverse[0][0] = coefficient;
        if (!paired(degenerate, MPoint(0, 1, 0.5), MPoint(0, 1, 0), MVector(2, 0, 0), MVector(0, 0, 0)))
            return false;
    }
    return true;
}

static bool kneeEndPreparation()
{
    const SkirtLegProfile profile(1.2, 1.4, 1.5, 2.0, 0.8, 0.9, 0.5, 0.6, 0.5, 0.5, 4, 6);
    for (double scaleY : {0.75, 1.0, 1.5})
        for (auto kind : {SkirtLegProfile::Ring::Knee, SkirtLegProfile::Ring::Extended, SkirtLegProfile::Ring::Heel})
        {
            const double L = (kind == SkirtLegProfile::Ring::Knee ? 4.0 : 10.0) * scaleY;
            const auto radius = profile.forRing(0.7, kind, scaleY);
            MMatrix matrix;
            matrix[0][0] = 0.8 * radius.x;
            matrix[0][1] = 0.6 * radius.x;
            matrix[1][0] = -0.6 * L;
            matrix[1][1] = 0.8 * L;
            matrix[2][2] = 2.0 * radius.z;
            matrix[3][0] = 0.25;
            matrix[3][1] = -0.5;
            PreparedBellRing ring(matrix);
            if (ring.distalEnd || ring.distalLength != 0 || ring.distalWidth != 0)
                return false;
            const bool distal = kind == SkirtLegProfile::Ring::Knee || kind == SkirtLegProfile::Ring::Extended;
            BellColliderSolver::prepareDistalEnd(ring, matrix, distal);
            const double a = MVector(matrix[0][0], matrix[0][1], matrix[0][2]).length();
            const double b = MVector(matrix[2][0], matrix[2][1], matrix[2][2]).length();
            const double expectedWidth = (std::min)(L, (std::max)(a, b));
            if (ring.distalEnd != distal || std::abs(ring.distalLength - L) > 1e-12 ||
                std::abs(ring.distalWidth - expectedWidth) > 1e-12)
                return false;
            MPointArray source;
            const MVector radial(0.8, 0.6, 0);
            for (double distance : {-1.0, L * 0.5, L, L + ring.distalWidth * 0.5,
                                    L + ring.distalWidth, L + ring.distalWidth + 1, L + 2 * ring.distalWidth})
                source.append(ring.translation + ring.normal * distance + radial * 0.2);
            source.append(ring.translation + ring.normal * L + radial * 20);
            source.append(MPoint(0.4, 2 * L, -0.0, 2.0));
            const std::vector<MVector> directions(source.length(), radial * 4);
            const std::vector<double> scales = {1, 0.25, 0.5, 0.75, 1, 1, 1, 1, 0.25};
            for (bool cap : {false, true})
                for (double collision : {0.15, 0.65, 1.0})
                    for (int start : {0, 1, 2})
                        for (int count : {1, 2, 5, 7})
                            if (!directionalCase(ring, source, directions, scales, collision, start, count, cap))
                                return false;
            auto actual = source;
            BellColliderSolver::relaxTowardRingBoundary(actual, ring, 0.65, 1, 7, true, directions, scales);
            for (unsigned int i : {0u, 7u, 8u})
                if (!sameBits(actual[i].x, source[i].x) || !sameBits(actual[i].y, source[i].y) ||
                    !sameBits(actual[i].z, source[i].z) || !sameBits(actual[i].w, source[i].w))
                    return false;
            if (distal)
                for (unsigned int i : {5u, 6u})
                    if (!sameBits(actual[i].x, source[i].x) || !sameBits(actual[i].y, source[i].y) ||
                        !sameBits(actual[i].z, source[i].z))
                        return false;
            if (!distal)
            {
                PreparedBellRing disabled(matrix);
                auto expected = source;
                BellColliderSolver::relaxTowardRingBoundary(expected, disabled, 0.65, 1, 7, true, directions, scales);
                if (!compare(expected, actual, 1, 7))
                    return false;
            }
        }
    for (double L : {10.0, 4.0})
    {
        MMatrix matrix;
        matrix[0][0] = 5;
        matrix[1][1] = L;
        matrix[2][2] = 3;
        PreparedBellRing prepared(matrix);
        BellColliderSolver::prepareDistalEnd(prepared, matrix, true);
        const double a = MVector(matrix[0][0], matrix[0][1], matrix[0][2]).length();
        const double b = MVector(matrix[2][0], matrix[2][1], matrix[2][2]).length();
        const double width = (std::min)(L, (std::max)(a, b));
        if (!prepared.distalEnd || prepared.distalLength != L || prepared.distalWidth != width)
            return false;
    }
    {
        MMatrix matrix;
        matrix[0][0] = matrix[2][2] = 2;
        matrix[1][1] = 4;
        PreparedBellRing prepared(matrix);
        BellColliderSolver::prepareDistalEnd(prepared, matrix, true);
        if (prepared.distalLength != 4 || prepared.distalWidth != 2)
            return false;
        MPointArray source;
        for (double distance : {-1.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 7.0, 8.0})
            source.append(MPoint(-0.25, distance, -0.0));
        const std::vector<MVector> directions(source.length(), MVector(4, 0, 0));
        const int start = 1, count = 7;
        const double previousDistance = (source[start + count - 2] - prepared.plane.orig) * prepared.plane.normal;
        const double trailingDistance = (source[start + count - 1] - prepared.plane.orig) * prepared.plane.normal;
        if (!(previousDistance > prepared.distalLength &&
              previousDistance < prepared.distalLength + prepared.distalWidth &&
              trailingDistance > prepared.distalLength + prepared.distalWidth))
            return false;
        if (!directionalCase(prepared, source, directions, {}, 0.65, start, count, true))
            return false;
        auto actual = source, scalar = source;
        scalarDirectionalRelax(scalar, prepared, 0.65, start, count, true, directions, {});
        BellColliderSolver::relaxTowardRingBoundary(actual, prepared, 0.65, start, count, true, directions, {});
        for (unsigned int i = 0; i < source.length(); ++i)
        {
            const bool inRange = i >= static_cast<unsigned int>(start) && i < static_cast<unsigned int>(start + count);
            const MPoint &expected = inRange ? scalar[i] : source[i];
            if (!sameBits(actual[i].x, expected.x) || !sameBits(actual[i].y, expected.y) ||
                !sameBits(actual[i].z, expected.z) || !sameBits(actual[i].w, expected.w))
                return false;
        }
    }

    MMatrix wide;
    wide[0][0] = 3;
    wide[1][1] = 2;
    wide[2][2] = 5;
    PreparedBellRing ring(wide);
    BellColliderSolver::prepareDistalEnd(ring, wide, true);
    if (!ring.distalEnd || ring.distalWidth != 2 || ring.distalLength != 2)
        return false;
    for (double nonfinite : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        ring.direction = MVector(nonfinite, 0, 0);
        BellColliderSolver::prepareDistalEnd(ring, wide, true);
        if (ring.distalEnd || std::isfinite(ring.distalLength))
            return false;
    }
    wide[0][0] = std::numeric_limits<double>::infinity();
    ring.direction = MVector(0, std::numeric_limits<double>::infinity(), 0);
    BellColliderSolver::prepareDistalEnd(ring, wide, true);
    return !ring.distalEnd && !std::isfinite(ring.distalWidth);
}

int main()
{
    if (!kneeEndPreparation() || !directionalReference() || !directionalZero() || !directionalAxisAndCentre() ||
        !directionalEntryDepth() || !directionalBoundaryAndCap() || !directionalFade() || !directionalPair())
        return 1;
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
