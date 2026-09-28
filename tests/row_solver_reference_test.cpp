#include "bellColliderSolver.h"
#include "bellColliderRelaxKernel.h"
#include "skirtLegProfile.h"

#include <maya/MPointArray.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <utility>
#include <maya/MTransformationMatrix.h>

static bool sameBits(double a, double b)
{
    std::uint64_t ua = 0, ub = 0;
    std::memcpy(&ua, &a, sizeof(ua));
    std::memcpy(&ub, &b, sizeof(ub));
    return ua == ub;
}

static void referenceRelax(MPointArray &points, const PreparedBellRing &ring, double collision, int startIndex,
                           int count, bool capAtRingOrigin)
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
    if (expected.length() != actual.length())
        return false;
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

struct RowFixture
{
    std::vector<double> u;
    std::set<std::pair<int, int>> cuts;
    bool periodic = false;
};

static bool hasEdge(const RowFixture &fixture, int a, int b)
{
    if (fixture.cuts.count(std::make_pair(a, b)) || fixture.cuts.count(std::make_pair(b, a)))
        return false;
    const int count = static_cast<int>(fixture.u.size());
    return std::abs(a - b) == 1 || (fixture.periodic && count > 1 && std::abs(a - b) == count - 1);
}

static BellRowTopology makeTopology(const RowFixture &fixture)
{
    BellRowTopology topology;
    topology.outputCount = static_cast<unsigned int>(fixture.u.size());
    topology.vertices.resize(fixture.u.size());
    if (fixture.u.empty())
        return topology;
    int component = -1;
    for (size_t i = 0; i < fixture.u.size(); ++i)
    {
        if (i == 0 || !hasEdge(fixture, static_cast<int>(i - 1), static_cast<int>(i)))
        {
            ++component;
            BellRowComponent interval;
            interval.startU = fixture.u[i];
            interval.endU = fixture.u[i] + 1.0;
            topology.components.push_back(interval);
        }
        auto &vertex = topology.vertices[i];
        vertex.materialU = fixture.u[i];
        vertex.componentId = component;
        vertex.panelId = component;
        vertex.previous =
            i > 0 && hasEdge(fixture, static_cast<int>(i - 1), static_cast<int>(i)) ? static_cast<int>(i - 1) : -1;
        vertex.next = i + 1 < fixture.u.size() && hasEdge(fixture, static_cast<int>(i), static_cast<int>(i + 1))
                          ? static_cast<int>(i + 1)
                          : -1;
        vertex.outputDuplicates.push_back(static_cast<unsigned int>(i));
        topology.components[component].endU = fixture.u[i] + (vertex.previous < 0 ? 0.01 : 0.0);
    }
    // Each cut between consecutive vertices is a seam: the vertex before it
    // is the +1 bank, the vertex after it the -1 bank, as the node builds them.
    int seam = 0;
    for (const auto &cut : fixture.cuts)
    {
        const int a = (std::min)(cut.first, cut.second);
        const int b = (std::max)(cut.first, cut.second);
        if (b != a + 1 || b >= static_cast<int>(fixture.u.size()))
            continue;
        topology.vertices[a].side.seamIndex = seam;
        topology.vertices[a].side.bank = 1;
        topology.vertices[b].side.seamIndex = seam;
        topology.vertices[b].side.bank = -1;
        ++seam;
    }
    if (fixture.periodic && fixture.cuts.empty())
    {
        topology.components[0].closed = true;
        topology.components[0].startU = fixture.u.front();
        topology.components[0].endU = fixture.u.front() + 1.0;
        topology.vertices.front().previous = static_cast<int>(fixture.u.size() - 1);
        topology.vertices.back().next = 0;
    }
    return topology;
}

static RowFixture chain(int count, bool periodic = false)
{
    RowFixture result;
    result.periodic = periodic;
    for (int i = 0; i < count; ++i)
        result.u.push_back(static_cast<double>(i) / (periodic ? count : (std::max)(1, count - 1)));
    return result;
}

static std::vector<std::vector<double>> referenceDistances(const RowFixture &fixture)
{
    const size_t count = fixture.u.size();
    std::vector<std::vector<double>> distances(count,
                                               std::vector<double>(count, std::numeric_limits<double>::infinity()));
    for (size_t i = 0; i < count; ++i)
    {
        distances[i][i] = 0.0;
        for (size_t j = i + 1; j < count; ++j)
            if (hasEdge(fixture, static_cast<int>(i), static_cast<int>(j)))
            {
                double length = fixture.u[j] - fixture.u[i];
                if (fixture.periodic && i == 0 && j + 1 == count)
                {
                    const double wrap = fixture.u.front() + 1.0 - fixture.u.back();
                    length = count == 2 ? (std::min)(length, wrap) : wrap;
                }
                distances[i][j] = distances[j][i] = length;
            }
    }
    for (size_t k = 0; k < count; ++k)
        for (size_t i = 0; i < count; ++i)
            for (size_t j = 0; j < count; ++j)
                distances[i][j] = (std::min)(distances[i][j], distances[i][k] + distances[k][j]);
    return distances;
}

static void referenceSmooth(std::vector<MVector> &values, const RowFixture &fixture, double smoothness)
{
    if (smoothness == 0.0)
        return;
    const int count = static_cast<int>(values.size());
    const double alpha = smoothness / 2.0;
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        std::vector<MVector> next(values);
        for (int i = 0; i < count; ++i)
        {
            int previous = i - 1;
            int following = i + 1;
            if (fixture.periodic)
            {
                if (previous < 0)
                    previous = count - 1;
                if (following == count)
                    following = 0;
            }
            if (previous < 0 || !hasEdge(fixture, previous, i))
                previous = i;
            if (following >= count || !hasEdge(fixture, i, following))
                following = i;
            if (previous != i || following != i)
                next[i] = values[i] * (1.0 - alpha) + (values[previous] + values[following]) * (alpha / 2.0);
        }
        values.swap(next);
    }
}

static std::vector<MVector> referenceFollow(const std::vector<MVector> &values, const RowFixture &fixture,
                                            const MMatrix &bell, double range)
{
    std::vector<MVector> result(values.size(), MVector(0, 0, 0));
    const MVector axis = maxis(bell, 1);
    if (axis.length() <= 1e-5)
        return result;
    const MVector normal = axis.normal();
    const auto distances = referenceDistances(fixture);
    for (size_t i = 0; i < values.size(); ++i)
    {
        MVector value = values[i];
        if (range > 0.0)
        {
            MVector sum(0, 0, 0);
            double total = 0.0;
            for (size_t j = 0; j < values.size(); ++j)
                if (distances[i][j] <= range)
                {
                    const double weight = values[j].length();
                    sum += values[j] * weight;
                    total += weight;
                }
            if (total < 1e-12)
                continue;
            value = sum / total;
        }
        result[i] = value - normal * (value * normal);
    }
    return result;
}

static bool sameVectors(const std::vector<MVector> &a, const std::vector<MVector> &b)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!sameBits(a[i].x, b[i].x) || !sameBits(a[i].y, b[i].y) || !sameBits(a[i].z, b[i].z))
        {
            std::cerr << "vector mismatch index=" << i << "\n";
            return false;
        }
    return true;
}

static bool zero(const MVector &value)
{
    return value.x == 0.0 && value.y == 0.0 && value.z == 0.0;
}

static bool runRowSmoothing(const RowFixture &fixture, const std::vector<MVector> &source, double smoothness)
{
    auto expected = source;
    auto actual = source;
    referenceSmooth(expected, fixture, smoothness);
    return BellColliderSolver::smoothDisplacements(actual, smoothness, makeTopology(fixture)) &&
           sameVectors(expected, actual);
}

static bool runRowFollow(const RowFixture &fixture, const std::vector<MVector> &source, double range,
                         const MMatrix &bell = MMatrix())
{
    std::vector<MVector> actual;
    const auto expected = referenceFollow(source, fixture, bell, range);
    return BellColliderSolver::computeLocalFollow(source, makeTopology(fixture), bell, range, actual) &&
           sameVectors(expected, actual);
}

static bool smoothingCutEdge()
{
    auto fixture = chain(65);
    std::vector<MVector> source(65, MVector(0, 0, 0));
    source[0] = MVector(1, -2, 3);
    for (double smoothness : {0.0, 0.5, 1.0})
    {
        if (!runRowSmoothing(fixture, source, smoothness))
            return false;
        auto actual = source;
        if (!BellColliderSolver::smoothDisplacements(actual, smoothness, makeTopology(fixture)))
            return false;
        for (size_t i = 4; i < actual.size(); ++i)
            if (!zero(actual[i]))
                return false;
        fixture.periodic = true;
        if (!runRowSmoothing(fixture, source, smoothness))
            return false;
        actual = source;
        if (!BellColliderSolver::smoothDisplacements(actual, smoothness, makeTopology(fixture)) ||
            (smoothness > 0.0 && zero(actual.back())))
            return false;
        fixture.periodic = false;
    }
    return true;
}

static bool smoothingComponents()
{
    for (int cut = 0; cut < 16; ++cut)
    {
        auto fixture = chain(17);
        fixture.cuts.insert({cut, cut + 1});
        fixture.cuts.insert({12, 13});
        std::vector<MVector> source(17, MVector(0, 0, 0));
        const auto random = makePoints(17, 1241);
        for (int i = 0; i <= (std::min)(cut, 12); ++i)
            source[i] = MVector(random[i].x, random[i].y, random[i].z);
        if (!runRowSmoothing(fixture, source, 0.5))
            return false;
        auto actual = source;
        if (!BellColliderSolver::smoothDisplacements(actual, 0.5, makeTopology(fixture)))
            return false;
        for (int i = (std::min)(cut, 12) + 1; i < 17; ++i)
            if (!zero(actual[i]))
                return false;
    }
    for (int count : {0, 1, 2})
        for (bool periodic : {false, true})
        {
            const auto fixture = chain(count, periodic);
            std::vector<MVector> values(count, MVector(-0.0, 0.0, 2.0));
            if (count == 2)
                values[1] = MVector(3, 1, -4);
            if (count == 1 && periodic)
            {
                if (BellColliderSolver::smoothDisplacements(values, 1.0, makeTopology(fixture)) !=
                    MS::kInvalidParameter)
                    return false;
                continue;
            }
            if (!runRowSmoothing(fixture, values, 1.0) || !runRowFollow(fixture, values, 1.0))
                return false;
        }
    auto fixture = chain(4);
    fixture.cuts = {{0, 1}, {1, 2}, {2, 3}};
    std::vector<MVector> values(4, MVector(-0.0, 0.0, 3.0));
    if (!runRowSmoothing(fixture, values, 0.5))
        return false;
    auto topology = makeTopology(chain(4));
    topology.vertices[0].next = 3;
    if (BellColliderSolver::smoothDisplacements(values, 0.5, topology) != MS::kInvalidParameter)
        return false;
    topology = makeTopology(chain(4));
    topology.vertices[0].outputDuplicates = {4};
    return BellColliderSolver::smoothDisplacements(values, 0.5, topology) == MS::kInvalidParameter;
}

static bool localFollow()
{
    auto fixture = chain(65);
    std::vector<MVector> values(65, MVector(0, 0, 0));
    values[0] = MVector(2, 7, -1);
    for (bool periodic : {false, true})
    {
        fixture.periodic = periodic;
        if (!runRowFollow(fixture, values, 3.0 / 64.0))
            return false;
    }
    fixture.periodic = false;
    std::vector<MVector> actual;
    if (!BellColliderSolver::computeLocalFollow(values, makeTopology(fixture), MMatrix(), 3.0 / 64.0, actual) ||
        zero(actual[3]) || !zero(actual[4]) || !zero(actual.back()))
        return false;
    for (double length : {0.0, std::nextafter(1e-12, 0.0), 1e-12, std::nextafter(1e-12, 1.0)})
    {
        values[0] = MVector(length, 0, 0);
        if (!runRowFollow(fixture, values, 0.1) || !runRowFollow(fixture, values, 0.0))
            return false;
    }
    for (int count : {3, 16, 65})
    {
        const auto circle = chain(count, true);
        std::vector<MVector> input(count, MVector(0, 0, 0));
        input[0] = MVector(1, 2, 3);
        if (!runRowFollow(circle, input, 0.1875))
            return false;
        if (count == 3)
        {
            if (!BellColliderSolver::computeLocalFollow(input, makeTopology(circle), MMatrix(), 0.1875, actual) ||
                !zero(actual[1]) || !zero(actual[2]))
                return false;
        }
    }
    RowFixture uneven;
    uneven.u = {0.02, 0.3, 0.7, 0.98};
    uneven.periodic = true;
    values.assign(4, MVector(0, 0, 0));
    values[0] = MVector(1, 0, 0);
    const double wrap = uneven.u.front() + 1.0 - uneven.u.back();
    if (!runRowFollow(uneven, values, wrap) || !runRowFollow(uneven, values, std::nextafter(wrap, 0.0)))
        return false;
    uneven.u = {0, 0.3, 0.6, 0.61, 0.62};
    values = {MVector(1, 0, 1), MVector(3, 0, -2), MVector(-4, 0, 3), MVector(2, 0, 7), MVector(-5, 0, -6)};
    if (!runRowFollow(uneven, values, 0.5))
        return false;
    for (double axis : {std::nextafter(1e-5, 0.0), 1e-5, std::nextafter(1e-5, 1.0)})
    {
        MMatrix bell;
        bell[1][1] = axis;
        if (!runRowFollow(uneven, values, 0.0, bell) || !runRowFollow(uneven, values, 0.5, bell))
            return false;
    }
    return true;
}

static bool followDetour()
{
    for (int count : {5, 65})
    {
        const auto fixture = chain(count);
        std::vector<MVector> values(count, MVector(0, 0, 0));
        values[0] = MVector(1, 0, 0);
        const double range = count == 5 ? 0.25 : 3.0 / 64.0;
        auto expected = referenceFollow(values, fixture, MMatrix(), range);
        std::vector<MVector> actual;
        if (!BellColliderSolver::computeLocalFollow(values, makeTopology(fixture), MMatrix(), range, actual))
            return false;
        for (int i = 0; i < count; ++i)
        {
            expected[i] = values[i] + expected[i];
            actual[i] = values[i] + actual[i];
        }
        referenceSmooth(expected, fixture, 0.5);
        if (!BellColliderSolver::smoothDisplacements(actual, 0.5, makeTopology(fixture)) ||
            !sameVectors(expected, actual) || (count == 5 ? actual.back().x <= 0.0 : !zero(actual.back())))
            return false;
    }
    return true;
}

static BellRowTopology sampleTopology(const std::vector<double> &positions, const std::vector<int> &banks)
{
    RowFixture fixture;
    fixture.u = positions;
    for (size_t i = 1; i < positions.size(); ++i)
        fixture.cuts.insert({static_cast<int>(i - 1), static_cast<int>(i)});
    auto topology = makeTopology(fixture);
    for (size_t i = 0; i < positions.size(); ++i)
    {
        topology.components[i].startU = positions[i] - (banks[i] == 1 ? 1e-6 : 0.0);
        topology.components[i].endU = positions[i] + (banks[i] == 1 ? 0.0 : 1e-6);
        topology.vertices[i].side.bank = banks[i];
        topology.vertices[i].side.seamIndex = banks[i] == 0 ? -1 : 7;
    }
    return topology;
}

static std::vector<MVector> referenceTransfer(const BellDirectField &source, const BellRowTopology &destination)
{
    std::vector<MVector> result(destination.vertices.size(), MVector(0, 0, 0));
    for (size_t i = 0; i < destination.vertices.size(); ++i)
    {
        const auto &target = destination.vertices[i];
        bool matched = false;
        // Resolve identified banks before looking for an interpolation interval.
        if (target.side.bank != 0)
            for (size_t j = 0; j < source.values.size(); ++j)
            {
                const auto &sample = source.topology.vertices[j];
                const double offset = target.materialU - sample.materialU;
                if (sample.side.bank == target.side.bank && sample.side.seamIndex == target.side.seamIndex &&
                    std::abs(offset - std::round(offset)) <= 1e-9)
                {
                    result[i] = source.values[j];
                    matched = true;
                    break;
                }
            }
        if (matched)
            continue;
        for (size_t c = 0; c < source.topology.components.size(); ++c)
        {
            const auto &component = source.topology.components[c];
            std::vector<std::pair<double, size_t>> samples;
            for (size_t j = 0; j < source.values.size(); ++j)
                if (source.topology.vertices[j].componentId == static_cast<int>(c))
                    samples.emplace_back(source.topology.vertices[j].materialU, j);
            std::sort(samples.begin(), samples.end());
            if (component.closed)
                samples.emplace_back(samples.front().first + 1.0, samples.front().second);
            const double origin = component.closed ? samples.front().first : component.startU;
            double u = target.materialU - std::floor(target.materialU - origin);
            if (!component.closed)
            {
                if (u > component.endU + 1e-9)
                    continue;
                bool wrongSide = false;
                for (size_t endpoint : {samples.front().second, samples.back().second})
                {
                    const auto &sample = source.topology.vertices[endpoint];
                    if (std::abs(u - sample.materialU) <= 1e-9 && target.side.bank != 0 && sample.side.bank != 0 &&
                        (target.side.bank != sample.side.bank || target.side.seamIndex != sample.side.seamIndex))
                        wrongSide = true;
                }
                if (wrongSide)
                    continue;
            }
            if (u <= samples.front().first)
                result[i] = source.values[samples.front().second];
            else if (u >= samples.back().first)
                result[i] = source.values[samples.back().second];
            else
                for (size_t j = 1; j < samples.size(); ++j)
                    if (u <= samples[j].first)
                    {
                        if (u == samples[j].first)
                            result[i] = source.values[samples[j].second];
                        else
                        {
                            const double fraction =
                                (u - samples[j - 1].first) / (samples[j].first - samples[j - 1].first);
                            result[i] = source.values[samples[j - 1].second] * (1.0 - fraction) +
                                        source.values[samples[j].second] * fraction;
                        }
                        break;
                    }
            break;
        }
    }
    return result;
}

static void reverseStorage(BellRowTopology &topology)
{
    const int count = static_cast<int>(topology.vertices.size());
    std::reverse(topology.vertices.begin(), topology.vertices.end());
    std::reverse(topology.components.begin(), topology.components.end());
    for (auto &vertex : topology.vertices)
    {
        vertex.componentId = static_cast<int>(topology.components.size()) - 1 - vertex.componentId;
        if (vertex.previous >= 0)
            vertex.previous = count - 1 - vertex.previous;
        if (vertex.next >= 0)
            vertex.next = count - 1 - vertex.next;
    }
}

static bool runRowTransfer(const RowFixture &fixture, const std::vector<MVector> &values,
                           const std::vector<double> &positions, const std::vector<int> &banks)
{
    BellDirectField field;
    field.topology = makeTopology(fixture);
    field.values = values;
    if (!fixture.periodic)
        for (auto &vertex : field.topology.vertices)
        {
            if (vertex.previous < 0)
                vertex.side.bank = -1;
            if (vertex.next < 0)
                vertex.side.bank = 1;
            if (vertex.side.bank != 0)
                vertex.side.seamIndex = 7;
        }
    const auto destination = sampleTopology(positions, banks);
    const auto expected = referenceTransfer(field, destination);
    for (bool reversed : {false, true})
    {
        if (reversed)
        {
            reverseStorage(field.topology);
            std::reverse(field.values.begin(), field.values.end());
        }
        BellRowTransfer actual;
        if (!BellColliderSolver::transferDirectField(field, destination, actual) ||
            !sameVectors(expected, actual.values) || !sameVectors(expected, referenceTransfer(field, destination)))
            return false;
    }
    return true;
}

static bool directTransfer()
{
    RowFixture fixture;
    fixture.periodic = true;
    fixture.u = {0.125, 0.375, 0.625, 0.875};
    const std::vector<MVector> values = {MVector(1, 0, 2), MVector(3, 0, 4), MVector(5, 0, 6), MVector(7, 0, 8)};
    if (!runRowTransfer(fixture, values, {-0.01, 0.0, 0.01, 0.5, 0.5, 1.0}, {0, 0, 0, -1, 1, 0}))
        return false;
    fixture.periodic = false;
    fixture.u = {0.25, 0.5, 1.0, 1.25};
    const std::vector<double> positions = {0.25, 0.25, 0.375, 0.75, 1.125, 1.25, 1.25};
    const std::vector<int> banks = {-1, 1, 0, 0, 0, -1, 1};
    for (int seed : {0, 3})
    {
        std::vector<MVector> injected(4, MVector(0, 0, 0));
        injected[seed] = MVector(1, 0, 3);
        if (!runRowTransfer(fixture, injected, positions, banks))
            return false;
    }
    for (double offset : {-5e-10, 0.0, 5e-10})
        if (!runRowTransfer(fixture, values, {0.25 + offset, 0.25 + offset, 1.25 + offset, 1.25 + offset},
                            {-1, 1, -1, 1}))
            return false;
    fixture.u = {0.0, 0.25, 0.5, 0.5, 0.75, 1.0};
    fixture.cuts = {{2, 3}};
    std::vector<MVector> injected(6, MVector(0, 0, 0));
    injected[2] = MVector(9, 0, 3);
    if (!runRowTransfer(fixture, injected, {0.375, 0.5, 0.5, 0.625}, {0, -1, 1, 0}))
        return false;
    for (double offset : {-2e-9, -5e-10, 0.0, 5e-10, 2e-9})
        if (!runRowTransfer(fixture, injected, {0.5 + offset, 0.5 + offset}, {-1, 1}))
            return false;
    BellRowTransfer result;
    if (!BellColliderSolver::transferDirectField(BellDirectField(), sampleTopology({0.1, 0.5}, {0, 0}), result))
        return false;
    for (const auto &value : result.values)
        if (!zero(value))
            return false;
    // Component bounds can extend beyond the first and last supplied samples.
    BellDirectField held;
    RowFixture heldFixture;
    heldFixture.u = {0.25, 0.75};
    held.topology = makeTopology(heldFixture);
    held.topology.components[0].startU = 0.0;
    held.topology.components[0].endU = 1.0;
    held.values = {MVector(2, 0, 0), MVector(4, 0, 0)};
    return BellColliderSolver::transferDirectField(held, sampleTopology({0.125, 0.875}, {0, 0}), result) &&
           sameVectors(result.values, held.values);
}

static double referenceDirectionDot(const BellRowInputs &inputs, const PreparedBellRing &ring, const MPoint &point)
{
    const MMatrix inverse = inputs.bellMatrix.inverse();
    const Plane plane(taxis(inputs.bellMatrix), maxis(inputs.bellMatrix, 1).normal());
    const MPoint projection = plane.projectPoint(point);
    const MVector offset = projection * inverse - plane.projectPoint(ring.translation) * inverse;
    return offset.normal() * (plane.projectVector(ring.direction) * inverse).normal();
}

static int referenceSeed(const BellRowInputs &inputs, const PreparedBellRing &ring, const MPointArray &base,
                         const RowFixture &fixture)
{
    const MVector axis = maxis(inputs.bellMatrix, 1);
    const Plane plane(taxis(inputs.bellMatrix), axis.normal());
    MPoint bellHit, ringHit, lineHit;
    if (axis.length() <= 1e-5 ||
        !BellColliderSolver::collisionPoints(inputs.bellMatrix, inputs.bellMatrix.inverse(), plane, ring, bellHit,
                                             ringHit, lineHit) ||
        !((plane.distance(ringHit) - plane.distance(bellHit)) / axis.length() < 0.0))
        return -1;
    int seed = -1;
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        const double z = referenceDirectionDot(inputs, ring, base[i]);
        if (seed < 0 || z > referenceDirectionDot(inputs, ring, base[seed]) ||
            (z == referenceDirectionDot(inputs, ring, base[seed]) && fixture.u[i] < fixture.u[seed]))
            seed = static_cast<int>(i);
    }
    return seed;
}

static MVector referenceBellOffset(const BellRowInputs &inputs, const MPoint &point)
{
    const MMatrix inverse = inputs.bellMatrix.inverse();
    const Plane plane(taxis(inputs.bellMatrix), maxis(inputs.bellMatrix, 1).normal());
    return MVector(plane.projectPoint(point) * inverse);
}

static double wrapAngle(double angle)
{
    while (angle > 3.141592653589793)
        angle -= 6.283185307179586;
    while (angle <= -3.141592653589793)
        angle += 6.283185307179586;
    return angle;
}

static std::vector<double> referenceComponentWeights(const BellRowInputs &inputs, const PreparedBellRing &ring,
                                                     const MPointArray &base, const RowFixture &fixture)
{
    // The ray runs from the bell centre through the ring's bell contact point
    // (the forward hit of the ring line with the unit sphere, as in
    // collisionPoints). The crossing edge is chosen independently with polar
    // angles about the bell centre (the shorter arc of a connected edge
    // contains the ray; the nearest such crossing along the ray wins). The
    // material position of the crossing, the interval test and the fade repeat
    // the solver's arithmetic so partial weights match bitwise.
    const auto distances = referenceDistances(fixture);
    const MMatrix inverse = inputs.bellMatrix.inverse();
    const Plane plane(taxis(inputs.bellMatrix), maxis(inputs.bellMatrix, 1).normal());
    const MPointArray hits = findSphereLineIntersection(plane.projectPoint(ring.translation) * inverse,
                                                        plane.projectVector(ring.direction) * inverse,
                                                        MPoint(0, 0, 0), 1.001);
    const MVector direction = hits.length() ? MVector(hits[0]) : MVector(0, 0, 0);
    const double directionAngle = std::atan2(direction.z, direction.x);
    const auto cross = [](const MVector &a, const MVector &b) { return a.z * b.x - a.x * b.z; };
    std::vector<MVector> offsets(base.length());
    std::vector<double> polar(base.length());
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        offsets[i] = referenceBellOffset(inputs, base[i]);
        polar[i] = std::atan2(offsets[i].z, offsets[i].x);
    }
    bool single = true;
    for (unsigned int j = 0; j < base.length(); ++j)
        single = single && std::isfinite(distances[0][j]);
    std::vector<double> weights(base.length(), 1.0);
    if (single && fixture.periodic && fixture.cuts.empty())
        return weights;
    double nearest = std::numeric_limits<double>::infinity();
    double contactU = 0.0;
    // Cut pairs are the two banks of a seam; they form virtual edges whose
    // every point carries the seam's material U (the first bank's value).
    for (unsigned int j = 0; j + 1 < base.length(); ++j)
    {
        const unsigned int k = j + 1;
        const bool virtualEdge = !hasEdge(fixture, static_cast<int>(j), static_cast<int>(k));
        if (virtualEdge && !fixture.cuts.count(std::make_pair(static_cast<int>(j), static_cast<int>(k))) &&
            !fixture.cuts.count(std::make_pair(static_cast<int>(k), static_cast<int>(j))))
            continue;
        const double arc = wrapAngle(polar[k] - polar[j]);
        const double offset = wrapAngle(directionAngle - polar[j]);
        const bool contains = (arc >= 0.0 && offset >= 0.0 && offset <= arc) ||
                              (arc <= 0.0 && offset <= 0.0 && offset >= arc) ||
                              std::abs(wrapAngle(directionAngle - polar[k])) == 0.0;
        if (!contains)
            continue;
        const double cj = cross(offsets[j], direction);
        const double ck = cross(offsets[k], direction);
        MVector hit;
        double u;
        if (cj == 0.0)
        {
            hit = offsets[j];
            u = fixture.u[j];
        }
        else if (ck == 0.0)
        {
            hit = offsets[k];
            u = virtualEdge ? fixture.u[j] : fixture.u[k];
        }
        else
        {
            const double lambda = cj / (cj - ck);
            hit = offsets[j] + (offsets[k] - offsets[j]) * lambda;
            u = virtualEdge ? fixture.u[j] : fixture.u[j] + (fixture.u[k] - fixture.u[j]) * lambda;
        }
        const double t = hit * direction;
        if (t > 0.0 && t < nearest)
        {
            nearest = t;
            contactU = u;
        }
    }
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        double low = std::numeric_limits<double>::infinity();
        double high = -std::numeric_limits<double>::infinity();
        for (unsigned int j = 0; j < base.length(); ++j)
            if (std::isfinite(distances[i][j]))
            {
                low = (std::min)(low, fixture.u[j]);
                high = (std::max)(high, fixture.u[j]);
            }
        if (!std::isfinite(nearest))
        {
            const int seed = referenceSeed(inputs, ring, base, fixture);
            weights[i] = seed >= 0 && std::isfinite(distances[i][static_cast<unsigned int>(seed)]) ? 1.0 : 0.0;
            continue;
        }
        double inside = contactU - low;
        inside -= std::floor(inside);
        if (inside > high - low)
        {
            weights[i] = 0.0;
            continue;
        }
        const double edge = (std::min)(inside, high - low - inside);
        if (inputs.contactBlendWidth > 0.0 && edge < inputs.contactBlendWidth)
        {
            const double t = edge / inputs.contactBlendWidth;
            weights[i] = t * t * (3.0 - 2.0 * t);
        }
    }
    return weights;
}

static MPointArray referenceDeform(const BellRowInputs &inputs, const PreparedBellRing &ring, const MPointArray &base,
                                   const std::vector<double> &mask)
{
    MPointArray points = base;
    const MMatrix inverse = inputs.bellMatrix.inverse();
    const MPoint origin = taxis(inputs.bellMatrix);
    const MVector axis = maxis(inputs.bellMatrix, 1);
    const MVector normal = axis.normal();
    const Plane plane(origin, normal);
    MPoint bellHit, ringHit, lineHit;
    if (axis.length() > 1e-5 &&
        BellColliderSolver::collisionPoints(inputs.bellMatrix, inverse, plane, ring, bellHit, ringHit, lineHit) &&
        (plane.distance(ringHit) - plane.distance(bellHit)) / axis.length() < 0.0)
    {
        MTransformationMatrix rotation;
        rotation.setTranslation(ring.translation, MSpace::kWorld);
        const MMatrix back = rotation.asMatrixInverse();
        const MQuaternion quaternion(bellHit - ring.translation, ringHit - ring.translation);
        rotation.rotateBy(quaternion, MSpace::kTransform);
        const MMatrix forward = rotation.asMatrix();
        const Plane top(origin + axis, normal);
        const MPoint ringProjection = plane.projectPoint(ring.translation) * inverse;
        const MVector direction = (plane.projectVector(ring.direction) * inverse).normal();
        for (unsigned int i = 0; i < points.length(); ++i)
            if (mask[i] > 0.0)
            {
                const MPoint projection = plane.projectPoint(points[i]);
                const MVector offset = projection * inverse - ringProjection;
                double weight = offset.normal() * direction;
                if (weight > inputs.falloff)
                {
                    const double divisor = 1.0 - inputs.falloff;
                    weight = divisor > 1e-5 ? (weight - inputs.falloff) / divisor : 1.0;
                    if (inputs.smoothness > 0.0)
                        weight = weight * weight * (3.0 - 2.0 * weight);
                    const MPoint rotated = points[i] * back * forward;
                    MPoint constrained = rotated;
                    if (top.distance(rotated) > 0.0)
                        constrained = top.projectPoint(rotated);
                    points[i] = constrained * weight + points[i] * (1.0 - weight);
                }
            }
    }
    referenceDirectionalRelax(points, ring, inputs.collision, 0, static_cast<int>(points.length()),
                               inputs.capAtRingOrigin, referenceDirections(points, base), {});
    MPointArray still = base;
    referenceDirectionalRelax(still, ring, inputs.collision, 0, static_cast<int>(still.length()),
                               inputs.capAtRingOrigin, std::vector<MVector>(still.length(), MVector(0, 0, 0)), {});
    for (unsigned int i = 0; i < points.length(); ++i)
        if (mask[i] > 0.0 && mask[i] < 1.0)
            points[i] = MPoint(still[i].x * (1.0 - mask[i]) + points[i].x * mask[i],
                               still[i].y * (1.0 - mask[i]) + points[i].y * mask[i],
                               still[i].z * (1.0 - mask[i]) + points[i].z * mask[i], points[i].w);
    return points;
}

static BellRowOutputs referenceRow(const BellRowInputs &inputs, const MPointArray &base, const RowFixture &fixture)
{
    BellRowOutputs result;
    result.points = base;
    result.directDisplacements.assign(base.length(), MVector(0, 0, 0));
    result.directField.values.assign(base.length(), MVector(0, 0, 0));
    if (inputs.rings.empty())
        return result;
    std::vector<MPointArray> responses;
    std::vector<std::vector<double>> masks;
    for (const auto &ring : inputs.rings)
    {
        std::vector<double> mask(base.length(), 0.0);
        if (referenceSeed(inputs, ring, base, fixture) >= 0)
            mask = referenceComponentWeights(inputs, ring, base, fixture);
        responses.push_back(referenceDeform(inputs, ring, base, mask));
        masks.push_back(mask);
    }
    const bool gate = inputs.smoothness > 0.0 || inputs.followGain > 0.0;
    for (unsigned int i = 0; i < base.length(); ++i)
    {
        double total = 0.0, cubic = 0.0, longest = 0.0;
        for (const auto &response : responses)
        {
            const double length = (response[i] - base[i]).length();
            total += std::pow(length, 2);
            cubic += length * length * length;
            longest = (std::max)(longest, length);
        }
        if (total > 0.0)
        {
            MVector sum(0, 0, 0);
            for (const auto &response : responses)
            {
                const MVector vector = response[i] - base[i];
                const double squared = std::pow(vector.length(), 2);
                const double weight = squared / total;
                sum += vector * weight;
            }
            const double average = cubic / total;
            result.directDisplacements[i] = average > 0.0 ? sum * (longest / average) : sum;
        }
    }
    result.directField.values =
        referenceFollow(result.directDisplacements, fixture, inputs.bellMatrix, inputs.followRange);
    auto displacements = result.directDisplacements;
    if (gate)
    {
        for (size_t i = 0; i < displacements.size(); ++i)
            displacements[i] += result.directField.values[i] * inputs.followGain;
        referenceSmooth(displacements, fixture, inputs.smoothness);
    }
    for (unsigned int i = 0; i < base.length(); ++i)
        result.points[i] = base[i] + displacements[i];
    if (gate)
        for (size_t r = 0; r < inputs.rings.size(); ++r)
        {
            MPointArray relaxed = result.points;
            std::vector<double> scales(base.length(), 1.0);
            for (unsigned int i = 0; i < base.length(); ++i)
            {
                const double m = masks[r][i];
                if (m > 0.0 && m < 1.0)
                    scales[i] = 1.0 - 4.0 * m * (1.0 - m);
            }
            referenceDirectionalRelax(relaxed, inputs.rings[r], inputs.collision, 0, static_cast<int>(base.length()),
                                       inputs.capAtRingOrigin, referenceDirections(result.points, base), scales);
            for (unsigned int i = 0; i < base.length(); ++i)
            {
                const double m = masks[r][i];
                if (m <= 0.0 || m >= 1.0)
                {
                    result.points[i] = relaxed[i];
                    continue;
                }
                const double f = 1.0 - 4.0 * m * (1.0 - m);
                result.points[i] = MPoint(result.points[i].x * (1.0 - f) + relaxed[i].x * f,
                                          result.points[i].y * (1.0 - f) + relaxed[i].y * f,
                                          result.points[i].z * (1.0 - f) + relaxed[i].z * f, relaxed[i].w);
            }
        }
    return result;
}

static bool runRowCase(const BellRowInputs &inputs, const MPointArray &base, const RowFixture &fixture)
{
    const BellRowOutputs expected = referenceRow(inputs, base, fixture);
    BellRowOutputs actual;
    if (!(BellColliderSolver::solveRow(inputs, base, makeTopology(fixture), actual) &&
           compare(expected.points, actual.points, 0, static_cast<int>(base.length())) &&
           sameVectors(expected.directDisplacements, actual.directDisplacements) &&
           sameVectors(expected.directField.values, actual.directField.values)))
        return false;
    MPointArray displacedBase = base;
    for (unsigned int i = 0; i < displacedBase.length(); ++i)
        displacedBase[i] = MPoint(displacedBase[i].x + 0.125, displacedBase[i].y - 0.25,
                                  displacedBase[i].z + 0.375, displacedBase[i].w);
    const BellRowOutputs displacedExpected = referenceRow(inputs, displacedBase, fixture);
    BellRowOutputs displacedActual;
    return BellColliderSolver::solveRow(inputs, displacedBase, makeTopology(fixture), displacedActual) &&
           compare(displacedExpected.points, displacedActual.points, 0, static_cast<int>(displacedBase.length())) &&
           sameVectors(displacedExpected.directDisplacements, displacedActual.directDisplacements) &&
           sameVectors(displacedExpected.directField.values, displacedActual.directField.values);
}

static bool mergeBeforeFollow()
{
    const auto fixture = chain(5);
    MPointArray base;
    for (int i = 0; i < 5; ++i)
        base.append(MPoint(0.2 + i * 0.1, 1.0, 0.0));
    BellRowInputs inputs;
    inputs.collision = 1.0f;
    inputs.followRange = 0.5;
    MMatrix first;
    MMatrix second;
    second[3][0] = 0.8;
    inputs.rings.emplace_back(first);
    inputs.rings.emplace_back(second);
    for (double smoothness : {0.0, 0.5})
        for (double follow : {0.0, 0.75})
        {
            inputs.smoothness = smoothness;
            inputs.followGain = follow;
            if (!runRowCase(inputs, base, fixture))
                return false;
        }
    inputs.rings.clear();
    if (!runRowCase(inputs, base, fixture))
        return false;
    return runRowCase(inputs, MPointArray(), RowFixture());
}

static bool accumulateFields(const std::vector<double> &levels, const std::vector<BellDirectField> &fields, double t,
                             double follow, double tightness, const BellRowTopology &destination,
                             std::vector<MVector> &values)
{
    values.assign(destination.vertices.size(), MVector(0, 0, 0));
    if (levels.size() != fields.size() + 1)
        return false;
    for (size_t k = 0; k < fields.size(); ++k)
    {
        if (!(levels[k + 1] > levels[k]))
            return false;
        const double alpha = (std::max)(0.0, (std::min)(1.0, (t - levels[k]) / (levels[k + 1] - levels[k])));
        if (alpha == 0.0)
            continue;
        BellRowTransfer transferred;
        if (!BellColliderSolver::transferDirectField(fields[k], destination, transferred))
            return false;
        for (size_t i = 0; i < values.size(); ++i)
            values[i] += transferred.values[i] * alpha;
    }
    for (auto &value : values)
        value = value * (follow * (1.0 - tightness));
    return true;
}

static bool oneWayAndInsertions()
{
    const auto topology = makeTopology(chain(5));
    std::vector<BellDirectField> fields(2);
    for (size_t k = 0; k < fields.size(); ++k)
    {
        fields[k].topology = topology;
        fields[k].values.assign(5, MVector(k == 0 ? 2.0 : 4.0, 0.0, 0.0));
    }
    const std::vector<double> levels = {0.25, 0.5, 0.75};
    std::vector<MVector> baseline;
    if (!accumulateFields(levels, fields, 0.75, 1.0, 0.0, topology, baseline))
        return false;
    for (const auto &value : baseline)
        if (value.x != 6.0)
            return false;
    std::vector<MVector> actual;
    for (double t : {0.125, 0.375, 0.625, 0.75})
    {
        if (!accumulateFields(levels, fields, t, 1.0, 0.0, topology, actual))
            return false;
        const double expected = t == 0.125 ? 0.0 : t == 0.375 ? 1.0 : t == 0.625 ? 4.0 : 6.0;
        for (const auto &value : actual)
            if (value.x != expected)
                return false;
    }
    if (!sameVectors(baseline, actual))
        return false;
    std::vector<MVector> upperBefore;
    if (!accumulateFields(levels, fields, 0.375, 1.0, 0.0, topology, upperBefore))
        return false;
    fields[1].values.assign(5, MVector(1000, 0, -50));
    if (!accumulateFields(levels, fields, 0.375, 1.0, 0.0, topology, actual) || !sameVectors(upperBefore, actual) ||
        accumulateFields({0.25, 0.25, 0.75}, fields, 0.75, 1.0, 0.0, topology, actual))
        return false;
    const double bellScaleY = 2.0;
    const double safeDistance = (std::max)(0.00002, 0.0001);
    const MPoint bottom(0.0, 3.0, 0.0);
    const MPoint top(2.0, 3.0 + safeDistance * bellScaleY, 0.0);
    const double start = 0.25;
    const double end = 0.5;
    const double t = 0.375;
    const double alpha = (t - start) / (end - start);
    const MPoint inserted = bottom * (1.0 - alpha) + top * alpha;
    BellRowInputs inputs;
    inputs.bellMatrix[3][0] = inserted.x;
    inputs.bellMatrix[3][1] = inserted.y;
    inputs.bellMatrix[3][2] = inserted.z;
    inputs.bellMatrix[1][0] = top.x - bottom.x;
    inputs.bellMatrix[1][1] = top.y - bottom.y;
    inputs.followGain = 0.8 * (t / 1.0) * 0.5;
    inputs.physicalLevel = 1;
    MPointArray row;
    for (int i = 0; i < 5; ++i)
        row.append(inserted + MVector(0.1 * i, 0.0, 0.0));
    BellRowOutputs output;
    if (!sameBits(inserted.x, 1.0) || !sameBits(inserted.y, (bottom.y + top.y) * 0.5) ||
        inputs.followGain != 0.15000000000000002 || !BellColliderSolver::solveRow(inputs, row, topology, output) ||
        !compare(row, output.points, 0, 5))
        return false;
    const SkirtLegProfile profile(1.1, 0.9, 0.8, 0.7, 0.9, 0.8, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5);
    const double targetDistance = t;
    const double parameter = profile.levelParameter(targetDistance, 0.0);
    if (parameter != t)
        return false;
    const auto radius = profile.forRing(parameter, SkirtLegProfile::Ring::Knee, 1.0);
    MMatrix ringMatrix;
    ringMatrix[0][0] = radius.x;
    ringMatrix[2][2] = radius.z;
    ringMatrix[3][1] = targetDistance;
    inputs.rings.emplace_back(ringMatrix);
    inputs.collision = 1.0f;
    return runRowCase(inputs, row, chain(5));
}

static BellRowOutputs blendRows(const BellRowOutputs &with, const BellRowOutputs &without, double tightness)
{
    if (tightness == 0.0)
        return with;
    if (tightness == 1.0)
        return without;
    BellRowOutputs result = with;
    for (unsigned int i = 0; i < result.points.length(); ++i)
    {
        result.points[i] = with.points[i] * (1.0 - tightness) + without.points[i] * tightness;
        result.directField.values[i] =
            with.directField.values[i] * (1.0 - tightness) + without.directField.values[i] * tightness;
    }
    return result;
}

static bool directionalVectors(const std::vector<MVector> &a, const std::vector<MVector> &b)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!(std::abs(a[i].x - b[i].x) <= 1e-12 && std::abs(a[i].y - b[i].y) <= 1e-12 &&
              std::abs(a[i].z - b[i].z) <= 1e-12))
            return false;
    return true;
}

static bool longTightness()
{
    const auto fixture = chain(5);
    const auto topology = makeTopology(fixture);
    MPointArray base;
    for (int i = 0; i < 5; ++i)
        base.append(MPoint(0.2 + 0.03 * i, 1.0, 0.1));
    BellRowInputs withInputs;
    withInputs.collision = 1.0f;
    withInputs.followGain = 0.2;
    withInputs.smoothness = 0.5;
    withInputs.rings.emplace_back(MMatrix());
    withInputs.rings.emplace_back(MMatrix());
    BellRowInputs withoutInputs = withInputs;
    MMatrix extended;
    extended[3][0] = 0.7;
    withInputs.rings.emplace_back(extended);
    extended[3][2] = 0.2;
    withInputs.rings.emplace_back(extended);
    const auto expectedWith = referenceRow(withInputs, base, fixture);
    const auto expectedWithout = referenceRow(withoutInputs, base, fixture);
    for (double tightness : {0.0, 0.5, 1.0})
    {
        BellRowOutputs with;
        BellRowOutputs without;
        if (tightness != 1.0 && !BellColliderSolver::solveRow(withInputs, base, topology, with))
            return false;
        if (tightness != 0.0 && !BellColliderSolver::solveRow(withoutInputs, base, topology, without))
            return false;
        auto actual = blendRows(with, without, tightness);
        auto expected = blendRows(expectedWith, expectedWithout, tightness);
        BellColliderSolver::roundMeshPoints(actual.points);
        for (unsigned int i = 0; i < expected.points.length(); ++i)
            expected.points[i] =
                MPoint(static_cast<float>(expected.points[i].x), static_cast<float>(expected.points[i].y),
                       static_cast<float>(expected.points[i].z));
        if (!directionalCompare(expected.points, actual.points) ||
            !directionalVectors(expected.directField.values, actual.directField.values))
            return false;
        for (size_t r = 0; r < withInputs.rings.size(); ++r)
        {
            if (r >= 2 && tightness == 1.0)
                continue;
            const double collision = r < 2 ? 1.0 : 1.0 - tightness;
            // Post-follow relaxation uses the displacement after the preceding ring.
            referenceDirectionalRelax(expected.points, withInputs.rings[r], collision, 0, 5, true,
                                       referenceDirections(expected.points, base), {});
            BellColliderSolver::relaxRowFaded(actual.points, withInputs.rings[r], collision, true, topology, {},
                                               referenceDirections(actual.points, base));
        }
        if (!directionalCompare(expected.points, actual.points))
            return false;
        for (double materialT : {0.49, 0.5, 0.51, 0.75})
            for (double mappedT : {0.4, 0.8})
            {
                BellRowInputs membership = withoutInputs;
                membership.physicalLevel = materialT <= 0.5 ? 1 : 2;
                const double rowTightness = membership.physicalLevel == 2 ? tightness : 0.0;
                if (rowTightness != (materialT > 0.5 ? tightness : 0.0))
                    return false;
                MPointArray mappedBase = base;
                for (unsigned int i = 0; i < mappedBase.length(); ++i)
                    mappedBase[i].y = mappedT;
                if (!runRowCase(membership, mappedBase, fixture))
                    return false;
            }
    }
    BellRowOutputs a;
    BellRowOutputs b;
    a.points.append(MPoint(1.0 + std::ldexp(1.0, -24) - std::ldexp(1.0, -27), 0, 0));
    b.points.append(MPoint(1.0 + 3.0 * std::ldexp(1.0, -24) - std::ldexp(1.0, -27), 0, 0));
    a.directField.values = b.directField.values = {MVector(1, 0, 0)};
    auto after = blendRows(a, b, 0.5);
    BellColliderSolver::roundMeshPoints(after.points);
    BellColliderSolver::roundMeshPoints(a.points);
    BellColliderSolver::roundMeshPoints(b.points);
    auto before = blendRows(a, b, 0.5);
    BellColliderSolver::roundMeshPoints(before.points);
    return !sameBits(after.points[0].x, before.points[0].x);
}

static MPointArray waistComposition(const MPointArray &base, const BellRowTopology &topology,
                                    const std::vector<PreparedBellRing> &rings, double smoothness)
{
    MPointArray result = base;
    for (const auto &ring : rings)
        BellColliderSolver::relaxTowardRingBoundary(result, ring, 1.0, 0, static_cast<int>(result.length()), true,
                                                    std::vector<MVector>(result.length(), MVector(0, 0, 0)), {});
    if (smoothness > 0.0)
    {
        std::vector<MVector> displacement(base.length());
        for (unsigned int i = 0; i < base.length(); ++i)
            displacement[i] = result[i] - base[i];
        if (!BellColliderSolver::smoothDisplacements(displacement, smoothness, topology))
        {
            result.clear();
            return result;
        }
        for (unsigned int i = 0; i < base.length(); ++i)
            result[i] = base[i] + displacement[i];
        for (const auto &ring : rings)
            BellColliderSolver::relaxTowardRingBoundary(result, ring, 1.0, 0, static_cast<int>(result.length()), true,
                                                        referenceDirections(result, base), {});
    }
    return result;
}

static bool waistRow()
{
    RowFixture fixture;
    fixture.u = {0.3, 0.4, 0.6, 0.8, 1.0, 1.15, 1.25, 1.3};
    const auto topology = makeTopology(fixture);
    MPointArray base;
    base.append(MPoint(0.4, 1, 0));
    for (int i = 1; i < 8; ++i)
        base.append(MPoint(2, 1, 0));
    std::vector<PreparedBellRing> rings;
    rings.emplace_back(MMatrix());
    for (double smoothness : {0.0, 0.5})
    {
        auto expected = base;
        // Before smoothing every waist ring receives zero direction.
        referenceDirectionalRelax(expected, rings[0], 1.0, 0, 8, true,
                                   std::vector<MVector>(8, MVector(0, 0, 0)), {});
        if (smoothness > 0.0)
        {
            std::vector<MVector> displacement(8);
            for (int i = 0; i < 8; ++i)
                displacement[i] = expected[i] - base[i];
            referenceSmooth(displacement, fixture, smoothness);
            for (int i = 0; i < 8; ++i)
                expected[i] = base[i] + displacement[i];
            referenceDirectionalRelax(expected, rings[0], 1.0, 0, 8, true, referenceDirections(expected, base), {});
        }
        MPointArray baseline;
        for (double follow : {0.0, 1.0})
        {
            const auto actual = waistComposition(base, topology, rings, smoothness);
            if (actual.length() != 8 || !directionalCompare(expected, actual) || !sameBits(actual[7].x, base[7].x))
                return false;
            if (follow == 0.0)
                baseline = actual;
            else if (!compare(baseline, actual, 0, 8))
                return false;
        }
        base[7] = MPoint(0.5, 1, 0);
        const auto lastContact = waistComposition(base, topology, rings, smoothness);
        if (sameBits(lastContact[7].x, base[7].x))
            return false;
        base[7] = MPoint(2, 1, 0);
    }
    return true;
}

static bool contactGate()
{
    BellRowInputs inputs;
    inputs.falloff = -0.5f;
    // Ring 0 points at the shared seam point (-1.5, 1, 1.5) with exactly
    // opposite x and z so both banks register the crossing; ring 1 points at
    // +Z, the second point of the first component.
    const double tilt = std::sin(0.6) / std::sqrt(2.0);
    for (const MVector direction : {MVector(-tilt, std::cos(0.6), tilt), MVector(0, std::cos(0.6), std::sin(0.6))})
    {
        MTransformationMatrix transform;
        transform.rotateBy(MQuaternion(MVector(0, 1, 0), direction), MSpace::kTransform);
        MMatrix matrix = transform.asMatrix();
        matrix[1][0] = direction.x;
        matrix[1][1] = direction.y;
        matrix[1][2] = direction.z;
        inputs.rings.emplace_back(matrix);
    }
    auto fixture = chain(8);
    fixture.u = {0.0, 0.125, 0.25, 0.375, 0.375, 0.5, 0.625, 0.75};
    fixture.cuts = {{3, 4}};
    const auto topology = makeTopology(fixture);
    MPointArray base;
    for (const MPoint point : {MPoint(2, 1, 0), MPoint(1.5, 1, 1.5), MPoint(0, 1, 2), MPoint(-1.5, 1, 1.5),
                               MPoint(-1.5, 1, 1.5), MPoint(-2, 1, 0), MPoint(-1.5, 1, -1.5), MPoint(0, 1, -2)})
        base.append(point);
    if (referenceSeed(inputs, inputs.rings[0], base, fixture) != 3 ||
        referenceSeed(inputs, inputs.rings[1], base, fixture) != 2)
        return false;
    for (const auto &ring : inputs.rings)
    {
        auto relaxed = base;
        referenceRelax(relaxed, ring, 1.0, 0, 8, false);
        if (!compare(base, relaxed, 0, 8))
            return false;
    }
    const auto unchanged = [](const MPoint &a, const MPoint &b)
    { return sameBits(a.x, b.x) && sameBits(a.y, b.y) && sameBits(a.z, b.z) && sameBits(a.w, b.w); };
    // Ring 0 crosses the row at the shared seam point, material 0.375, the
    // free end of both components: both keep the full weight with a zero blend
    // width and neither rotates with a positive one. Ring 1 crosses at
    // material 0.25 inside the first component, 0.125 from its end: full
    // weight with a zero blend width and the exact weight 0.15625 with a blend
    // width of 0.5. The second component does not contain the crossing and
    // stays still.
    for (double blend : {0.0, 0.5})
        for (double range : {0.0, 0.01, 1.0})
            for (double smoothness : {0.0, 0.5})
                for (size_t ringCount : {size_t(1), size_t(2)})
                {
                    auto current = inputs;
                    current.rings.resize(ringCount, inputs.rings[0]);
                    current.followRange = range;
                    current.smoothness = smoothness;
                    current.contactBlendWidth = blend;
                    std::vector<std::vector<unsigned int>> seeds = {{3}};
                    if (ringCount == 2)
                        seeds.push_back({2});
                    std::vector<MPointArray> actual;
                    if (!BellColliderSolver::deformPoints(current, base, topology, seeds, actual) ||
                        actual.size() != ringCount || !runRowCase(current, base, fixture))
                        return false;
                    for (size_t r = 0; r < ringCount; ++r)
                    {
                        const auto weights = referenceComponentWeights(current, current.rings[r], base, fixture);
                        for (int i = 0; i < 8; ++i)
                        {
                            const double expectedWeight =
                                r == 0 ? (blend > 0.0 ? 0.0 : 1.0) : (i < 4 ? (blend > 0.0 ? 0.15625 : 1.0) : 0.0);
                            if (weights[i] != expectedWeight)
                                return false;
                            const bool moving = weights[i] > 0.0 &&
                                                referenceDirectionDot(current, current.rings[r], base[i]) >
                                                    current.falloff;
                            if (unchanged(base[i], actual[r][i]) == moving)
                                return false;
                        }
                        if (!compare(referenceDeform(current, current.rings[r], base, weights), actual[r], 0, 8))
                            return false;
                    }
                auto reversedTopology = topology;
                reverseStorage(reversedTopology);
                MPointArray reversedBase;
                for (unsigned int i = base.length(); i > 0; --i)
                    reversedBase.append(base[i - 1]);
                BellRowOutputs forward, reversed;
                if (!BellColliderSolver::solveRow(current, base, topology, forward) ||
                    !BellColliderSolver::solveRow(current, reversedBase, reversedTopology, reversed))
                    return false;
                for (unsigned int i = 0; i < base.length(); ++i)
                    if (!unchanged(forward.points[i], reversed.points[base.length() - 1 - i]))
                        return false;
            }
    for (float falloff : {-0.5f, 0.0f, 0.75f, 1.0f})
        for (double smoothness : {0.0, 0.5})
        {
            auto current = inputs;
            current.falloff = falloff;
            current.smoothness = smoothness;
            current.followRange = 0.0;
            const auto closed = chain(8, true);
            std::vector<MPointArray> actual;
            if (!BellColliderSolver::deformPoints(current, base, makeTopology(closed), {{0}, {6}}, actual) ||
                !runRowCase(current, base, closed))
                return false;
            BellColliderInputs generic;
            generic.bellMatrix = current.bellMatrix;
            generic.bellSubdivision = 8;
            generic.falloff = falloff;
            generic.smoothness = smoothness;
            MPointArray genericBase;
            for (int i = 0; i < 9; ++i)
                genericBase.append(MPoint(0, 0, 0));
            for (unsigned int i = 0; i < base.length(); ++i)
                genericBase.append(base[i]);
            const Plane plane(taxis(current.bellMatrix), maxis(current.bellMatrix, 1).normal());
            std::vector<MPointArray> expected;
            BellColliderSolver::deformPoints(generic, current.rings, genericBase, plane, expected);
            for (size_t r = 0; r < current.rings.size(); ++r)
                for (unsigned int i = 0; i < base.length(); ++i)
                {
                    const bool candidate = referenceDirectionDot(current, current.rings[r], base[i]) > falloff;
                    if (!unchanged(expected[r][9 + i], actual[r][i]) || unchanged(base[i], actual[r][i]) == candidate)
                        return false;
                }
        }
    for (double height : {1.0, 2.0})
    {
        auto current = inputs;
        MMatrix matrix;
        matrix[0][0] = 0.0;
        matrix[0][1] = -1.0;
        matrix[1][0] = 1.0;
        matrix[1][1] = 0.0;
        matrix[3][1] = height;
        current.rings = {PreparedBellRing(matrix)};
        const Plane plane(MPoint(0, 0, 0), MVector(0, 1, 0));
        MPoint bellHit, ringHit, lineHit;
        if (!BellColliderSolver::collisionPoints(current.bellMatrix, current.bellMatrix.inverse(), plane,
                                                 current.rings[0], bellHit, ringHit, lineHit))
            return false;
        const double delta = plane.distance(ringHit) - plane.distance(bellHit);
        if ((height == 1.0 && delta != 0.0) || (height == 2.0 && !(delta > 0.0)))
            return false;
        std::vector<MPointArray> actual;
        if (!BellColliderSolver::deformPoints(current, base, topology, {{0}}, actual) ||
            !compare(base, actual[0], 0, 8) || !runRowCase(current, base, fixture))
            return false;
    }
    for (double length : {std::nextafter(1e-5, 0.0), 1e-5})
    {
        auto current = inputs;
        current.bellMatrix[1][1] = length;
        current.collision = 1.0f;
        MPointArray intruding;
        for (int i = 0; i < 8; ++i)
            intruding.append(MPoint(0.2, 1, 0.1));
        std::vector<MPointArray> actual;
        if (!BellColliderSolver::deformPoints(current, intruding, topology, {{0}, {6}}, actual) ||
            !runRowCase(current, intruding, fixture))
            return false;
        for (size_t r = 0; r < current.rings.size(); ++r)
        {
            auto expected = intruding;
            // No-rotation rows use the cross-section quadratic with zero direction.
            referenceDirectionalRelax(expected, current.rings[r], 1.0, 0, 8, false,
                                       std::vector<MVector>(8, MVector(0, 0, 0)), {});
            if (!directionalCompare(expected, actual[r]) || unchanged(intruding[0], actual[r][0]))
                return false;
        }
    }
    auto failed = inputs;
    failed.rings = {PreparedBellRing(MMatrix())};
    std::vector<MPointArray> actual;
    if (!BellColliderSolver::deformPoints(failed, base, topology, {{0}}, actual) || !compare(base, actual[0], 0, 8) ||
        !runRowCase(failed, base, fixture))
        return false;
    if (BellColliderSolver::deformPoints(inputs, base, topology, {{0, 4}, {6}}, actual) != MS::kInvalidParameter)
        return false;
    return true;
}

// Sweep a ring direction around a circle cut into two open halves and compare
// the largest single-step point motion with and without the blend width.
static bool contactBlendContinuity()
{
    // Two half circles that share their end points, like the two banks of a
    // seam: material 0.5 and material 0 = 1 are each held by both components.
    RowFixture fixture;
    fixture.u = {0.0, 0.125, 0.25, 0.375, 0.5, 0.5, 0.625, 0.75, 0.875, 1.0};
    fixture.cuts = {{4, 5}};
    const auto topology = makeTopology(fixture);
    MPointArray base;
    for (int i = 0; i < 10; ++i)
    {
        const int step = i < 5 ? i : i - 1;
        const double angle = step * 6.283185307179586 / 8.0;
        base.append(MPoint(2.0 * std::cos(angle), 1.0, 2.0 * std::sin(angle)));
    }
    const auto sweep = [&](double blend, double &largestStep, double &largestLift) {
        largestStep = 0.0;
        largestLift = 0.0;
        MPointArray previous;
        const int steps = 1440;
        for (int step = 0; step <= steps; ++step)
        {
            // Half a step off the vertex angles: on a vertex the arc test and
            // the acos of a rounded dot may disagree about an exact zero.
            const double phi = (step + 0.5) * 6.283185307179586 / steps;
            BellRowInputs inputs;
            inputs.falloff = -0.5f;
            inputs.contactBlendWidth = blend;
            MTransformationMatrix transform;
            transform.rotateBy(MQuaternion(MVector(0, 1, 0), MVector(std::sin(0.6) * std::cos(phi), std::cos(0.6),
                                                                      std::sin(0.6) * std::sin(phi))),
                               MSpace::kTransform);
            inputs.rings.emplace_back(transform.asMatrix());
            BellRowOutputs output;
            if (!BellColliderSolver::solveRow(inputs, base, topology, output) || !runRowCase(inputs, base, fixture))
                return false;
            for (unsigned int i = 0; i < base.length(); ++i)
            {
                largestLift = (std::max)(largestLift, (output.points[i] - base[i]).length());
                if (step > 0)
                    largestStep = (std::max)(largestStep, (output.points[i] - previous[i]).length());
            }
            previous = output.points;
        }
        return true;
    };
    double binaryStep = 0.0, binaryLift = 0.0, blendedStep = 0.0, blendedLift = 0.0;
    if (!sweep(0.0, binaryStep, binaryLift) || !sweep(0.125, blendedStep, blendedLift))
        return false;
    // The binary mask moves a whole half circle in one step. With a blend width
    // of one material column the weight changes by at most 1.5 / 0.125 per unit
    // of material, i.e. about 0.0083 per step of this sweep, and the sweep itself
    // turns the lift by about 0.0044 rad per step, so a step stays below about
    // 1.5 percent of the lift; 3 percent leaves room without hiding a defect.
    return binaryLift > 0.1 && binaryStep > 0.5 * binaryLift && blendedStep < 0.03 * blendedLift;
}

static bool relaxOutsideGate()
{
    MPointArray source;
    for (double d : {-1e-12, 0.0, 1e-12})
        for (double radius : {0.0, std::nextafter(1e-5, 0.0), 1e-5, std::nextafter(1e-5, 1.0), 0.5, 1.0,
                              std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0)})
            source.append(MPoint(radius, d, -0.0));
    const auto fixture = chain(static_cast<int>(source.length()));
    BellRowInputs inputs;
    inputs.rings.emplace_back(MMatrix());
    inputs.collision = 1.0f;
    for (bool cap : {false, true})
    {
        inputs.capAtRingOrigin = cap;
        std::vector<MPointArray> actual;
        auto expected = source;
        // The row path uses the quadratic even at the old radial fallback threshold.
        referenceDirectionalRelax(expected, inputs.rings[0], 1.0, 0, static_cast<int>(source.length()), cap,
                                   std::vector<MVector>(source.length(), MVector(0, 0, 0)), {});
        if (!BellColliderSolver::deformPoints(inputs, source, makeTopology(fixture), {{}}, actual) ||
            !directionalCompare(expected, actual[0]) || !runRowCase(inputs, source, fixture))
            return false;
    }
    MPointArray base;
    base.append(MPoint(1, 1, 0));
    base.append(MPoint(1, 1, 0));
    base.append(MPoint(1, 1, 0));
    const auto small = chain(3);
    std::vector<MVector> displacement = {MVector(-0.5, 0, 0), MVector(0, 0, 0), MVector(0, 0, 0)};
    if (!BellColliderSolver::smoothDisplacements(displacement, 0.5, makeTopology(small)))
        return false;
    for (int i = 0; i < 3; ++i)
        base[i] += displacement[i];
    if (!(base[2].x < 1.0))
        return false;
    auto expected = base;
    referenceDirectionalRelax(expected, inputs.rings[0], 1.0, 0, 3, true,
                               std::vector<MVector>(3, MVector(0, 0, 0)), {});
    std::vector<MPointArray> actual;
    return BellColliderSolver::deformPoints(inputs, base, makeTopology(small), {{}}, actual) &&
           directionalCompare(expected, actual[0]) && actual[0][2].x > base[2].x;
}

static bool relaxLaneIndependence()
{
    MMatrix firstMatrix;
    firstMatrix[0][0] = 2.0;
    firstMatrix[2][2] = 0.7;
    MMatrix secondMatrix;
    secondMatrix[3][0] = 0.3;
    const PreparedBellRing first(firstMatrix);
    const PreparedBellRing second(secondMatrix);
    auto source = makePoints(129, 271828);
    source[3] = MPoint(-0.0, -0.0, 0.0, 2.0);
    source[4] = MPoint(std::numeric_limits<double>::infinity(), 1, 0, 1);
    source[5] = MPoint(std::numeric_limits<double>::quiet_NaN(), 1, 0, 1);
    for (double collision : {0.0, 1e-5, 0.5, 1.0})
        if (!runCase(first, source, collision, 0, 129) || !runCase(first, source, collision, 3, 123) ||
            !runSequentialCase(first, second, source, collision))
            return false;
    for (int changedLane : {0, 1})
        for (bool cap : {false, true})
        {
            MPointArray baseline;
            baseline.append(MPoint(-0.0, -0.0, 0.0, 2.0));
            baseline.append(MPoint(0.2, 1.0, 0.1, 3.0));
            MPointArray injected = baseline;
            injected[changedLane] = MPoint(0.4, 2.0, 0.3, -1.0);
            auto expected = baseline;
            referenceRelax(expected, first, 1.0, 0, 2, cap);
            BellColliderSolver::relaxTowardRingBoundary(baseline, first, 1.0, 0, 2, cap);
            BellColliderSolver::relaxTowardRingBoundary(injected, first, 1.0, 0, 2, cap);
            if (!compare(expected, baseline, 0, 2))
                return false;
            const int preserved = 1 - changedLane;
            const MPoint &a = baseline[preserved];
            const MPoint &b = injected[preserved];
            if (!sameBits(a.x, b.x) || !sameBits(a.y, b.y) || !sameBits(a.z, b.z) || !sameBits(a.w, b.w))
                return false;
        }
    return true;
}

static bool followSummationOrder()
{
    RowFixture fixture;
    fixture.u = {0.0, 0.25, 0.5, 0.5, 0.75, 1.0};
    fixture.cuts = {{2, 3}};
    for (bool samePanel : {false, true})
    {
        auto topology = makeTopology(fixture);
        topology.vertices[2].side = {7, 1};
        topology.vertices[3].side = {7, -1};
        if (samePanel)
            for (auto &vertex : topology.vertices)
                vertex.panelId = 0;
        std::vector<MVector> values = {MVector(1e8, 2, 0), MVector(-1e8, -2, 0), MVector(1, 0, 0),
                                       MVector(3, 0, 0),   MVector(-1e8, 0, 0),  MVector(1e8, 0, 0)};
        reverseStorage(topology);
        std::reverse(values.begin(), values.end());
        if (samePanel)
        {
            std::swap(topology.vertices[2], topology.vertices[3]);
            std::swap(values[2], values[3]);
            for (auto &vertex : topology.vertices)
            {
                if (vertex.previous == 2 || vertex.previous == 3)
                    vertex.previous = 5 - vertex.previous;
                if (vertex.next == 2 || vertex.next == 3)
                    vertex.next = 5 - vertex.next;
            }
        }
        std::vector<size_t> originalIndex = {5, 4, 3, 2, 1, 0};
        if (samePanel)
            std::swap(originalIndex[2], originalIndex[3]);
        std::vector<size_t> order;
        for (size_t i = 0; i < values.size(); ++i)
            order.push_back(i);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const auto &left = topology.vertices[a];
            const auto &right = topology.vertices[b];
            if (left.materialU != right.materialU)
                return left.materialU < right.materialU;
            if (left.panelId != right.panelId)
                return left.panelId < right.panelId;
            return left.side.bank == -1 && right.side.bank != -1;
        });
        if (topology.vertices[order[2]].panelId > topology.vertices[order[3]].panelId ||
            (samePanel && topology.vertices[order[2]].side.bank != -1))
            return false;
        const auto distances = referenceDistances(fixture);
        std::vector<MVector> expected(values.size(), MVector(0, 0, 0));
        for (size_t i = 0; i < values.size(); ++i)
        {
            MVector sum(0, 0, 0);
            double total = 0.0;
            for (size_t j : order)
                if (distances[originalIndex[i]][originalIndex[j]] <= 1.0)
                {
                    const double weight = values[j].length();
                    sum += values[j] * weight;
                    total += weight;
                }
            const MVector value = sum / total;
            const MVector normal(0, 1, 0);
            expected[i] = value - normal * (value * normal);
        }
        std::vector<MVector> actual;
        if (!BellColliderSolver::computeLocalFollow(values, topology, MMatrix(), 1.0, actual) ||
            !sameVectors(expected, actual))
            return false;
    }
    return true;
}

static bool invalidRowInputs()
{
    const auto topology = makeTopology(chain(4));
    MPointArray base;
    for (int i = 0; i < 4; ++i)
        base.append(MPoint(i, 1, 0));
    const std::vector<MVector> values(4, MVector(1, 2, 3));
    BellRowInputs inputs;
    inputs.rings.emplace_back(MMatrix());
    BellRowOutputs output;
    std::vector<MPointArray> points;
    std::vector<MVector> follow;
    BellRowTransfer transfer;
    const auto rejectedTopology = [&](const BellRowTopology &bad) {
        auto smooth = values;
        BellDirectField source;
        source.topology = bad;
        source.values = values;
        BellDirectField validSource;
        validSource.topology = topology;
        validSource.values = values;
        return BellColliderSolver::solveRow(inputs, base, bad, output) == MS::kInvalidParameter &&
               BellColliderSolver::deformPoints(inputs, base, bad, {{0}}, points) == MS::kInvalidParameter &&
               BellColliderSolver::smoothDisplacements(smooth, 0.5, bad) == MS::kInvalidParameter &&
               BellColliderSolver::computeLocalFollow(values, bad, MMatrix(), 0.5, follow) == MS::kInvalidParameter &&
               BellColliderSolver::transferDirectField(source, topology, transfer) == MS::kInvalidParameter &&
               BellColliderSolver::transferDirectField(validSource, bad, transfer) == MS::kInvalidParameter;
    };
    for (int kind = 0; kind < 7; ++kind)
    {
        auto bad = topology;
        if (kind == 0)
            bad.vertices.pop_back();
        else if (kind == 1)
            bad.vertices[1].previous = -1;
        else if (kind == 2)
            bad.vertices[1].next = -1;
        else if (kind == 3)
            bad.vertices[1].componentId = 1;
        else if (kind == 4)
            bad.vertices[1].outputDuplicates = {0};
        else if (kind == 5)
            bad.vertices[1].outputDuplicates = {1, 1};
        else
            bad.outputCount = 5;
        if (!rejectedTopology(bad))
            return false;
    }
    auto split = chain(4);
    split.cuts = {{1, 2}};
    auto overlapping = makeTopology(split);
    overlapping.components[0].endU = 0.8;
    if (!rejectedTopology(overlapping))
        return false;
    auto singleton = makeTopology(chain(1, true));
    std::vector<MVector> singleValue(1, MVector(0, 0, 0));
    if (BellColliderSolver::computeLocalFollow(singleValue, singleton, MMatrix(), 1.0, follow) != MS::kInvalidParameter)
        return false;
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()})
    {
        auto bad = topology;
        bad.vertices[1].materialU = invalid;
        if (!rejectedTopology(bad))
            return false;
        bad = topology;
        bad.components[0].startU = invalid;
        if (!rejectedTopology(bad))
            return false;
        bad = topology;
        bad.components[0].endU = invalid;
        if (!rejectedTopology(bad))
            return false;
        for (int coordinate = 0; coordinate < 4; ++coordinate)
        {
            auto badBase = base;
            if (coordinate == 0)
                badBase[1].x = invalid;
            if (coordinate == 1)
                badBase[1].y = invalid;
            if (coordinate == 2)
                badBase[1].z = invalid;
            if (coordinate == 3)
                badBase[1].w = invalid;
            if (BellColliderSolver::solveRow(inputs, badBase, topology, output) != MS::kInvalidParameter ||
                BellColliderSolver::deformPoints(inputs, badBase, topology, {{0}}, points) != MS::kInvalidParameter)
                return false;
            if (coordinate == 3)
                continue;
            auto badValues = values;
            if (coordinate == 0)
                badValues[1].x = invalid;
            if (coordinate == 1)
                badValues[1].y = invalid;
            if (coordinate == 2)
                badValues[1].z = invalid;
            BellDirectField source;
            source.topology = topology;
            source.values = badValues;
            if (BellColliderSolver::computeLocalFollow(badValues, topology, MMatrix(), 0.5, follow) !=
                    MS::kInvalidParameter ||
                BellColliderSolver::smoothDisplacements(badValues, 0.0, topology) != MS::kInvalidParameter ||
                BellColliderSolver::transferDirectField(source, topology, transfer) != MS::kInvalidParameter)
                return false;
        }
    }
    for (double invalid : {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()})
        for (int control = 0; control < 4; ++control)
        {
            auto bad = inputs;
            if (control == 0)
                bad.smoothness = invalid;
            if (control == 1)
                bad.followGain = invalid;
            if (control == 2)
                bad.followRange = invalid;
            if (control == 3)
                bad.contactBlendWidth = invalid;
            if (BellColliderSolver::solveRow(bad, base, topology, output) != MS::kInvalidParameter ||
                BellColliderSolver::deformPoints(bad, base, topology, {{0}}, points) != MS::kInvalidParameter)
                return false;
            auto smooth = values;
            if (control == 0 &&
                BellColliderSolver::smoothDisplacements(smooth, invalid, topology) != MS::kInvalidParameter)
                return false;
            if (control == 2 && BellColliderSolver::computeLocalFollow(values, topology, MMatrix(), invalid, follow) !=
                                    MS::kInvalidParameter)
                return false;
        }
    auto bad = inputs;
    bad.smoothness = std::nextafter(1.0, 2.0);
    auto smooth = values;
    if (BellColliderSolver::solveRow(bad, base, topology, output) != MS::kInvalidParameter ||
        BellColliderSolver::deformPoints(bad, base, topology, {{0}}, points) != MS::kInvalidParameter ||
        BellColliderSolver::smoothDisplacements(smooth, bad.smoothness, topology) != MS::kInvalidParameter ||
        BellColliderSolver::deformPoints(inputs, base, topology, {}, points) != MS::kInvalidParameter ||
        BellColliderSolver::deformPoints(inputs, base, topology, {{4}}, points) != MS::kInvalidParameter)
        return false;
    BellDirectField shortField;
    shortField.topology = topology;
    shortField.values.assign(3, MVector(0, 0, 0));
    smooth.resize(3);
    MPointArray shortBase;
    shortBase.append(base[0]);
    return BellColliderSolver::transferDirectField(shortField, topology, transfer) == MS::kInvalidParameter &&
           BellColliderSolver::smoothDisplacements(smooth, 0.5, topology) == MS::kInvalidParameter &&
           BellColliderSolver::computeLocalFollow(smooth, topology, MMatrix(), 0.5, follow) == MS::kInvalidParameter &&
           BellColliderSolver::solveRow(inputs, shortBase, topology, output) == MS::kInvalidParameter &&
           BellColliderSolver::deformPoints(inputs, shortBase, topology, {{0}}, points) == MS::kInvalidParameter;
}

static bool runDirectionalRowCase(const BellRowInputs &inputs, const MPointArray &base, const RowFixture &fixture)
{
    for (bool displaced : {false, true})
    {
        auto current = base;
        if (displaced)
            for (unsigned int i = 0; i < current.length(); ++i)
                current[i] = MPoint(current[i].x + 0.125, current[i].y - 0.25, current[i].z + 0.375, current[i].w);
        const auto expected = referenceRow(inputs, current, fixture);
        BellRowOutputs actual;
        if (!BellColliderSolver::solveRow(inputs, current, makeTopology(fixture), actual) ||
            !directionalCompare(expected.points, actual.points) ||
            !directionalVectors(expected.directDisplacements, actual.directDisplacements) ||
            !directionalVectors(expected.directField.values, actual.directField.values))
            return false;
    }
    return true;
}

static bool directionalRowOrder()
{
    auto fixture = chain(7);
    fixture.cuts = {{2, 3}};
    auto topology = makeTopology(fixture);
    topology.vertices[2].outputDuplicates.push_back(topology.outputCount++);
    MPointArray base, points;
    for (int i = 0; i < 7; ++i)
    {
        base.append(MPoint(-2 + i * 0.1, 1, -0.5));
        points.append(MPoint(-0.4 + i * 0.1, 1, 0.3));
    }
    BellRowInputs inputs;
    MMatrix first, second;
    first[0][0] = 2;
    first[2][2] = 1.5;
    second[3][0] = -0.6;
    inputs.rings = {PreparedBellRing(first), PreparedBellRing(second)};
    inputs.collision = 1;
    inputs.smoothness = 0.1;
    inputs.followGain = 0.2;
    if (!runDirectionalRowCase(inputs, base, fixture))
        return false;
    auto reversedTopology = topology;
    std::reverse(reversedTopology.vertices.begin(), reversedTopology.vertices.end());
    for (auto &vertex : reversedTopology.vertices)
    {
        if (vertex.previous >= 0)
            vertex.previous = 6 - vertex.previous;
        if (vertex.next >= 0)
            vertex.next = 6 - vertex.next;
    }
    MPointArray reversedBase, reversedPoints;
    for (unsigned int i = base.length(); i > 0; --i)
    {
        reversedBase.append(base[i - 1]);
        reversedPoints.append(points[i - 1]);
    }
    for (const auto &ring : inputs.rings)
    {
        const auto directions = referenceDirections(points, base);
        auto reversedDirections = directions;
        std::reverse(reversedDirections.begin(), reversedDirections.end());
        if (!directionalVectors(reversedDirections, referenceDirections(reversedPoints, reversedBase)))
            return false;
        BellColliderSolver::relaxRowFaded(points, ring, 1, true, topology, {0.25, 0.75}, directions);
        BellColliderSolver::relaxRowFaded(reversedPoints, ring, 1, true, reversedTopology, {0.25, 0.75},
                                          reversedDirections);
    }
    MPointArray restored;
    for (unsigned int i = reversedPoints.length(); i > 0; --i)
        restored.append(reversedPoints[i - 1]);
    if (!directionalCompare(points, restored))
        return false;
    BellRowOutputs forward, reversed;
    if (!BellColliderSolver::solveRow(inputs, base, topology, forward) ||
        !BellColliderSolver::solveRow(inputs, reversedBase, reversedTopology, reversed))
        return false;
    restored.clear();
    std::reverse(reversed.directDisplacements.begin(), reversed.directDisplacements.end());
    std::reverse(reversed.directField.values.begin(), reversed.directField.values.end());
    for (unsigned int i = reversed.points.length(); i > 0; --i)
        restored.append(reversed.points[i - 1]);
    return directionalCompare(forward.points, restored) &&
           directionalVectors(forward.directDisplacements, reversed.directDisplacements) &&
           directionalVectors(forward.directField.values, reversed.directField.values);
}

static bool directionalPartialBlend()
{
    RowFixture fixture = chain(5);
    fixture.cuts = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    auto topology = makeTopology(fixture);
    topology.vertices[2].outputDuplicates.push_back(topology.outputCount++);
    const std::vector<double> weights = {0, 0.25, 0.5, 0.75, 1};
    const std::vector<double> fades = {1, 0.25, 0, 0.25, 1};
    MPointArray base, points;
    for (int i = 0; i < 5; ++i)
    {
        base.append(MPoint(-2, 1, 0.2));
        points.append(MPoint(-0.2, 1, 0.3));
    }
    MMatrix matrix;
    matrix[0][0] = matrix[2][2] = 2;
    const PreparedBellRing ring(matrix);
    const auto directions = referenceDirections(points, base);
    auto after = points;
    referenceDirectionalRelax(after, ring, 1, 0, 5, true, directions, fades);
    auto expected = points;
    for (int i = 0; i < 5; ++i)
    {
        const double m = weights[i];
        const double f = m > 0 && m < 1 ? 1 - 4 * m * (1 - m) : 1;
        if (f != fades[i])
            return false;
        expected[i] = MPoint(points[i].x * (1 - f) + after[i].x * f, points[i].y * (1 - f) + after[i].y * f,
                             points[i].z * (1 - f) + after[i].z * f, after[i].w);
    }
    auto actual = points;
    BellColliderSolver::relaxRowFaded(actual, ring, 1, true, topology, weights, directions);
    if (!directionalCompare(expected, actual) || !sameBits(actual[2].x, points[2].x))
        return false;
    MPointArray expanded;
    expanded.setLength(topology.outputCount);
    for (size_t i = 0; i < topology.vertices.size(); ++i)
        for (unsigned int cv : topology.vertices[i].outputDuplicates)
            expanded[cv] = actual[static_cast<unsigned int>(i)];
    if (!sameBits(expanded[2].x, expanded[topology.outputCount - 1].x) ||
        topology.vertices[2].outputDuplicates.size() != 2)
        return false;
    for (bool mismatch : {false, true})
    {
        auto currentTopology = topology;
        if (mismatch)
            currentTopology.vertices.pop_back();
        auto full = points, oracle = points;
        referenceDirectionalRelax(oracle, ring, 1, 0, 5, true, directions, {});
        BellColliderSolver::relaxRowFaded(full, ring, 1, true, currentTopology,
                                          mismatch ? weights : std::vector<double>(), directions);
        if (!directionalCompare(oracle, full))
            return false;
    }
    auto still = points, stillExpected = points;
    const std::vector<MVector> zero(5, MVector(0, 0, 0));
    referenceDirectionalRelax(stillExpected, ring, 1, 0, 5, true, zero, {});
    BellColliderSolver::relaxTowardRingBoundary(still, ring, 1, 0, 5, true, zero, {});
    if (!directionalCompare(stillExpected, still))
        return false;

    const auto liftedFixture = chain(8, true);
    BellRowInputs input;
    input.collision = 1;
    input.falloff = -1;
    input.smoothness = 0.1;
    input.followGain = 0.2;
    input.contactBlendWidth = 0.125;
    MTransformationMatrix transform;
    transform.rotateBy(MQuaternion(MVector(0, 1, 0), MVector(0.6, 0.8, 0)), MSpace::kTransform);
    input.rings.emplace_back(transform.asMatrix());
    MPointArray liftedBase;
    for (int i = 0; i < 8; ++i)
    {
        const double angle = i * 2 * std::acos(-1.0) / 8;
        liftedBase.append(MPoint(0.8 * std::cos(angle), 1, 0.8 * std::sin(angle)));
    }
    if (!runDirectionalRowCase(input, liftedBase, liftedFixture))
        return false;
    auto cutFixture = chain(8);
    cutFixture.u = {0, 0.125, 0.25, 0.375, 0.375, 0.5, 0.625, 0.75};
    cutFixture.cuts = {{3, 4}};
    auto cutTopology = makeTopology(cutFixture);
    cutTopology.vertices[3].outputDuplicates.push_back(cutTopology.outputCount++);
    MPointArray cutBase;
    for (const MPoint point : {MPoint(2, 1, 0), MPoint(1.5, 1, 1.5), MPoint(0, 1, 2), MPoint(-1.5, 1, 1.5),
                               MPoint(-1.5, 1, 1.5), MPoint(-2, 1, 0), MPoint(-1.5, 1, -1.5), MPoint(0, 1, -2)})
        cutBase.append(point);
    MTransformationMatrix cutTransform;
    cutTransform.rotateBy(MQuaternion(MVector(0, 1, 0), MVector(0, std::cos(0.6), std::sin(0.6))), MSpace::kTransform);
    input.rings = {PreparedBellRing(cutTransform.asMatrix())};
    for (double m : {0.25, 0.5, 0.75, 1.0})
    {
        double low = 0, high = 1;
        for (int iteration = 0; iteration < 60; ++iteration)
        {
            const double t = (low + high) * 0.5;
            if (t * t * (3 - 2 * t) < m)
                low = t;
            else
                high = t;
        }
        input.contactBlendWidth = m == 1 ? 0 : 0.125 / ((low + high) * 0.5);
        std::vector<MPointArray> responses;
        std::vector<std::vector<double>> componentWeights;
        if (!BellColliderSolver::deformPoints(input, cutBase, cutTopology, {{2}}, responses, &componentWeights) ||
            componentWeights.size() != 1 || componentWeights[0].size() != 2 ||
            std::abs(componentWeights[0][0] - m) > 1e-12 || componentWeights[0][1] != 0)
            return false;
        const auto mask = referenceComponentWeights(input, input.rings[0], cutBase, cutFixture);
        if (!directionalCompare(referenceDeform(input, input.rings[0], cutBase, mask), responses[0]) ||
            !runDirectionalRowCase(input, cutBase, cutFixture))
            return false;
    }
    return true;
}

static bool directionalSweep(bool activeSecond)
{
    const double R = 6.06, H = 58, epsilon = 0.1;
    const double X = std::sqrt(R * R - epsilon * epsilon);
    const MPoint B(0, H, -12);
    MPoint previousP, previousOutput;
    double previousAngle = 0, maximumRatio = 0, maximumAngle = 0;
    const auto topology = makeTopology(chain(1));
    const int first = activeSecond ? 4 : -10;
    for (int degrees = first; degrees <= 10; ++degrees)
    {
        const double theta = degrees * std::acos(-1.0) / 180;
        const double cosine = std::cos(theta), sine = std::sin(theta);
        const auto rotate = [&](double x, double y, double z) {
            return MPoint(x * cosine - y * sine, x * sine + y * cosine, z);
        };
        MMatrix matrix;
        matrix[0][0] = R * cosine;
        matrix[0][1] = R * sine;
        matrix[1][0] = -sine;
        matrix[1][1] = cosine;
        matrix[2][2] = R;
        const PreparedBellRing ring1(matrix);
        const auto origin = rotate(activeSecond ? -X - 0.5 * R : -2 * X, 0, 0);
        matrix[3][0] = origin.x;
        matrix[3][1] = origin.y;
        const PreparedBellRing ring3(matrix);
        const auto A = rotate(X, H, epsilon), E = rotate(-X, H, epsilon);
        const MVector da = A - B, de = E - B;
        const double wa = da * da / (da * da + de * de);
        const MPoint P(wa * A.x + (1 - wa) * E.x, wa * A.y + (1 - wa) * E.y, wa * A.z + (1 - wa) * E.z);
        MPointArray base, actual;
        base.append(B);
        actual.append(P);
        auto expected = actual;
        const auto initialDirections = referenceDirections(actual, base);
        referenceDirectionalRelax(expected, ring1, 1, 0, 1, true, initialDirections, {});
        BellColliderSolver::relaxRowFaded(actual, ring1, 1, true, topology, {1}, initialDirections);
        if (!directionalCompare(expected, actual))
            return false;
        const MPoint firstOutput = actual[0];
        const MVector q = firstOutput - P;
        const double angle = std::atan2(q.z, q.x * cosine + q.y * sine);
        auto stale = actual;
        BellColliderSolver::relaxRowFaded(stale, ring3, 1, true, topology, {1}, initialDirections);
        referenceDirectionalRelax(expected, ring3, 1, 0, 1, true, referenceDirections(expected, base), {});
        BellColliderSolver::relaxRowFaded(actual, ring3, 1, true, topology, {1}, referenceDirections(actual, base));
        if (!directionalCompare(expected, actual))
            return false;
        if (activeSecond && (!((actual[0] - firstOutput).length() > 1e-3) || !((actual[0] - stale[0]).length() > 1e-6)))
            return false;
        if (!activeSecond && degrees != first)
        {
            const double ratio = (actual[0] - previousOutput).length() / (P - previousP).length();
            const double delta =
                std::abs(std::atan2(std::sin(angle - previousAngle), std::cos(angle - previousAngle))) * 180 /
                std::acos(-1.0);
            maximumRatio = (std::max)(maximumRatio, ratio);
            maximumAngle = (std::max)(maximumAngle, delta);
            if (!(ratio <= 2 && delta <= 30))
                return false;
        }
        previousP = P;
        previousOutput = actual[0];
        previousAngle = angle;
    }
    if (!activeSecond)
        std::cout << "directional sweep maximum ratio " << maximumRatio << " angle " << maximumAngle << "\n";
    return true;
}

static bool directionalTwoRingSweep()
{
    return directionalSweep(false);
}
static bool directionalTwoRingSweepActiveSecond()
{
    return directionalSweep(true);
}

static bool runDirectionalCases()
{
    struct DirectionalCase
    {
        const char *name;
        bool (*run)();
    };
    const DirectionalCase directionalCases[] = {
        {"directionalRowOrder", directionalRowOrder},
        {"directionalPartialBlend", directionalPartialBlend},
        {"directionalTwoRingSweep", directionalTwoRingSweep},
        {"directionalTwoRingSweepActiveSecond", directionalTwoRingSweepActiveSecond}};
    bool success = true;
    for (const auto &test : directionalCases)
    {
        const bool passed = test.run();
        std::cout << (passed ? "PASS " : "FAIL ") << test.name << "\n";
        success = passed && success;
    }
    return success;
}

static MVector mergeAtOrigin(const std::vector<MVector> &displacements)
{
    MPointArray base;
    base.append(MPoint(0, 0, 0));
    std::vector<MPointArray> rings;
    for (const auto &displacement : displacements)
    {
        MPointArray points;
        points.append(base[0] + displacement);
        rings.push_back(points);
    }
    return BellColliderSolver::mergeDisplacement(base, rings, 0);
}

static bool mergeClose(const MVector &actual, const MVector &expected, double tolerance = 1e-12)
{
    return std::isfinite(actual.x) && std::isfinite(actual.y) && std::isfinite(actual.z) &&
           (actual - expected).length() <= tolerance;
}

// Independent statement of the merge rule: squared-length weighted mean scaled by max|d| / (sum|d|^3 / sum|d|^2).
static MVector mergeClosedForm(const std::vector<MVector> &displacements)
{
    double total = 0.0, cubic = 0.0, longest = 0.0;
    MVector numerator(0, 0, 0);
    for (const auto &d : displacements)
    {
        const double length = d.length();
        total += length * length;
        cubic += length * length * length;
        longest = (std::max)(longest, length);
        numerator += d * (length * length);
    }
    return total > 0.0 ? numerator * (longest / cubic) : MVector(0, 0, 0);
}

static bool mergeExpectations()
{
    const MVector x(1, 0, 0), z(0, 0, 0);
    struct Example
    {
        std::vector<MVector> inputs;
        MVector expected;
    };
    const Example examples[] = {
        {{x, x}, x}, {{x, x * 0.5}, x}, {{x, x * -0.5}, x * (7.0 / 9.0)},
        {{x, -x}, z}, {{MVector(1, 2, 3)}, MVector(1, 2, 3)},
        {{z, z}, z}, {{z, x}, x}, {{x, -x, x * 2}, x * 1.6},
        {{x * 2, x}, x * 2}, {{x * 0.5, x * -0.25}, x * (7.0 / 18.0)},
        {{}, z}, {{x * 0.5e-5}, x * 0.5e-5}, {{x * 1e-5}, x * 1e-5}};
    for (const auto &example : examples)
        if (!mergeClose(mergeAtOrigin(example.inputs), example.expected))
            return false;
    const MVector single(1, 2, 3);
    return mergeClose(mergeAtOrigin({single}), single.normal() * single.length(),
                      1e-12 * (std::max)(1.0, single.length()));
}

static double cancellationAmplitude(double t)
{
    return 2.0 * (1.0 - 3.0 * t * t + 2.0 * t * t * t);
}

static double cancellationClosedForm(double amplitude)
{
    const double cube = amplitude * amplitude * amplitude;
    return (cube - 1.0) / (cube + 1.0) * (std::max)(amplitude, 1.0);
}

static bool mergeCancellation()
{
    const auto evaluate = [](double amplitude) {
        return mergeAtOrigin({MVector(amplitude, 0, 0), MVector(-1, 0, 0)});
    };
    for (double t : {0.4999, 0.5, 0.5001})
    {
        const double a = cancellationAmplitude(t);
        if (!mergeClose(evaluate(a), MVector(cancellationClosedForm(a), 0, 0), 1e-9))
            return false;
    }
    double previous = 0.0;
    double h = 0.001;
    for (int i = 0; i < 12; ++i, h *= 0.5)
    {
        const double a = cancellationAmplitude(0.5 + h);
        const double b = cancellationAmplitude(0.5 - h);
        const MVector right = evaluate(a), left = evaluate(b);
        if (!mergeClose(right, MVector(cancellationClosedForm(a), 0, 0), 1e-9) ||
            !mergeClose(left, MVector(cancellationClosedForm(b), 0, 0), 1e-9))
            return false;
        const double difference = (right - left).length();
        if (!(difference > 0.0) || (i > 0 && std::abs(difference / previous / 0.5 - 1.0) > 0.05))
            return false;
        previous = difference;
    }
    for (double epsilon : {-0.00001, -0.000005, 0.000005, 0.00001})
    {
        const double a = 1.0 + epsilon;
        const MVector actual = evaluate(a);
        if (actual.x == 0.0 || !mergeClose(actual, MVector(cancellationClosedForm(a), 0, 0), 1e-9))
            return false;
    }
    return true;
}

static bool mergeScalingRotationBound()
{
    const std::vector<MVector> original = {MVector(1, 2, 3), MVector(-2, 1, -1), MVector(0.5, -1, 2)};
    const MVector mean = (original[0] * 14.0 + original[1] * 6.0 + original[2] * 5.25) / 25.25;
    const double average = (14.0 * std::sqrt(14.0) + 6.0 * std::sqrt(6.0) + 5.25 * std::sqrt(5.25)) / 25.25;
    const MVector expected = mean * (std::sqrt(14.0) / average);
    if (!mergeClose(mergeClosedForm(original), expected, 1e-12))
        return false;
    double scale = 1.0;
    for (int i = 0; i < 32; ++i, scale *= 0.5)
    {
        std::vector<MVector> scaled, rotated;
        double maximum = 0.0;
        for (const auto &d : original)
        {
            scaled.push_back(d * scale);
            rotated.push_back(MVector(-d.y, d.x, d.z) * scale);
            maximum = (std::max)(maximum, d.length() * scale);
        }
        const MVector actual = mergeAtOrigin(scaled);
        if (!mergeClose(actual, expected * scale, 1e-12 * scale) ||
            !mergeClose(mergeAtOrigin(rotated), MVector(-expected.y, expected.x, expected.z) * scale,
                        1e-12 * scale) || actual.length() > maximum + 1e-12 * scale)
            return false;
    }
    for (double a : {-2.0, -1.0, 0.0, 0.5, 1.0, 2.0})
        for (double b : {-2.0, -1.0, 0.0, 0.5, 1.0, 2.0})
        {
            const MVector first(a, b, 0), second(-b, 0, a), third(0, -a, b);
            const double maximum = (std::max)(first.length(), (std::max)(second.length(), third.length()));
            if (mergeAtOrigin({first, second, third}).length() > maximum + 1e-12)
                return false;
        }
    return true;
}

static MMatrix parallelMergeRing(double radius, double centerX)
{
    MMatrix matrix;
    matrix[0][0] = radius;
    matrix[2][2] = radius;
    matrix[3][0] = centerX;
    return matrix;
}

static MVector radialMergeDisplacement(const MPoint &point, double radius, double centerX, double collision)
{
    const MVector radial(point.x - centerX, 0, point.z);
    const double distance = radial.length();
    return distance < radius ? radial * (collision * (radius - distance) / distance) : MVector(0, 0, 0);
}

static bool bellSolveMerge()
{
    BellColliderInputs inputs;
    inputs.bellSubdivision = 16;
    const MPointArray base = BellColliderSolver::makeBellPoints(inputs.bellMatrix, 1, inputs.bellSubdivision);
    const unsigned int start = inputs.bellSubdivision + 1;
    const Plane plane(MPoint(0, 0, 0), MVector(0, 1, 0));
    for (double collision : {1.0, 0.5})
        for (bool offset : {false, true})
        {
            inputs.collision = static_cast<float>(collision);
            const double radii[] = {2.0, offset ? 2.0 : 1.5};
            const double centers[] = {0.0, offset ? 2.0 : 0.0};
            std::vector<BellColliderOutputs> singles(2);
            for (int r = 0; r < 2; ++r)
            {
                inputs.rings.clear();
                inputs.rings.emplace_back(parallelMergeRing(radii[r], centers[r]));
                MPoint bellHit, ringHit, lineHit;
                if (BellColliderSolver::collisionPoints(inputs.bellMatrix, inputs.bellMatrix.inverse(), plane,
                                                        inputs.rings[0], bellHit, ringHit, lineHit) ||
                    !BellColliderSolver::solve(inputs, base, singles[r]))
                    return false;
                for (unsigned int i = 0; i < base.length(); ++i)
                {
                    const MVector expected = i < start ? MVector(0, 0, 0) :
                        radialMergeDisplacement(base[i], radii[r], centers[r], collision);
                    if (!mergeClose(singles[r].points[i] - base[i], expected) ||
                        (i < start && (!sameBits(singles[r].points[i].x, base[i].x) ||
                                       !sameBits(singles[r].points[i].y, base[i].y) ||
                                       !sameBits(singles[r].points[i].z, base[i].z))) ||
                        !sameBits(singles[r].points[i].w, base[i].w))
                        return false;
                }
                if (!zero(singles[r].meanDisplacement))
                    return false;
            }
            inputs.rings.clear();
            for (int r = 0; r < 2; ++r)
                inputs.rings.emplace_back(parallelMergeRing(radii[r], centers[r]));
            BellColliderOutputs actual;
            if (!BellColliderSolver::solve(inputs, base, actual) || !zero(actual.meanDisplacement))
                return false;
            for (unsigned int i = 0; i < base.length(); ++i)
            {
                std::vector<MVector> analytic, fromSingle;
                if (i >= start)
                    for (int r = 0; r < 2; ++r)
                    {
                        analytic.push_back(radialMergeDisplacement(base[i], radii[r], centers[r], collision));
                        fromSingle.push_back(singles[r].points[i] - base[i]);
                    }
                const MVector expected = mergeClosedForm(analytic);
                const MVector fromSingles = mergeClosedForm(fromSingle);
                if (!mergeClose(actual.points[i] - base[i], expected) ||
                    !mergeClose(actual.points[i] - base[i], fromSingles) ||
                    !sameBits(actual.points[i].w, base[i].w))
                    return false;
                if (i < start && (!sameBits(actual.points[i].x, base[i].x) ||
                                  !sameBits(actual.points[i].y, base[i].y) ||
                                  !sameBits(actual.points[i].z, base[i].z)))
                    return false;
            }
            // Same centre: parallel pushes keep the larger length. Offset: the pushes cancel exactly.
            if (!mergeClose(actual.points[start] - base[start], MVector(offset ? 0.0 : collision, 0, 0)))
                return false;
        }
    return true;
}

// Regression oracle for the gate path (smoothness and followGain positive), recorded from the accepted
// build with the bellSolveGateBits inputs. Rows contain point x/y/z/w in output order, followed by
// meanDisplacement x/y/z/0. On mismatch the test prints the actual rows in the same layout.
static const std::vector<std::array<double, 4>> bellGateExpected = {
    {0, 0, 0, 1},
    {1, 0, 0, 1},
    {6.123233995736766e-17, 0, 1, 1},
    {-1, 0, 1.2246467991473532e-16, 1},
    {-1.8369701987210297e-16, 0, -1, 1},
    {1.946616841943523, 1, 2.3606207941720405e-18, 1},
    {0.054056402839586126, 1, 1.6045772431830492, 1},
    {-1.5940710036129544, 1, 2.0081366929062648e-16, 1},
    {0.054056402839585731, 1, -1.6045772431830492, 1},
    {0.06410264417328905, 0, 1.3080123850768456e-17, 0.0},
};

static bool bellSolveGateBits()
{
    BellColliderInputs inputs;
    inputs.bellSubdivision = 4;
    inputs.collision = 0.5f;
    inputs.smoothness = 0.5;
    inputs.followGain = 0.25;
    inputs.rings.emplace_back(parallelMergeRing(2.0, 0.0));
    inputs.rings.emplace_back(parallelMergeRing(1.5, 0.75));
    const MPointArray base = BellColliderSolver::makeBellPoints(inputs.bellMatrix, 1, inputs.bellSubdivision);
    BellColliderOutputs actual;
    if (!BellColliderSolver::solve(inputs, base, actual) || bellGateExpected.size() != actual.points.length() + 1)
        return false;
    bool matched = true;
    for (unsigned int i = 0; i < actual.points.length(); ++i)
    {
        const auto &expected = bellGateExpected[i];
        const MPoint &point = actual.points[i];
        if (!sameBits(point.x, expected[0]) || !sameBits(point.y, expected[1]) ||
            !sameBits(point.z, expected[2]) || !sameBits(point.w, expected[3]))
            matched = false;
    }
    const auto &mean = bellGateExpected.back();
    if (!sameBits(actual.meanDisplacement.x, mean[0]) || !sameBits(actual.meanDisplacement.y, mean[1]) ||
        !sameBits(actual.meanDisplacement.z, mean[2]))
        matched = false;
    if (!matched)
    {
        std::cout << std::setprecision(17);
        for (unsigned int i = 0; i < actual.points.length(); ++i)
            std::cout << "    {" << actual.points[i].x << ", " << actual.points[i].y << ", " << actual.points[i].z
                      << ", " << actual.points[i].w << "},\n";
        std::cout << "    {" << actual.meanDisplacement.x << ", " << actual.meanDisplacement.y << ", "
                  << actual.meanDisplacement.z << ", 0.0},\n";
    }
    return matched;
}

static bool kneeEndRowPreparation()
{
    const SkirtLegProfile profile(1.2, 1.4, 1.5, 2, 0.8, 0.9, 0.5, 0.6, 0.5, 0.5, 4, 6);
    for (double scaleY : {0.75, 1.5})
        for (auto kind : {SkirtLegProfile::Ring::Knee, SkirtLegProfile::Ring::Extended, SkirtLegProfile::Ring::Heel})
        {
            const auto radius = profile.forRing(0.7, kind, scaleY);
            const double L = (kind == SkirtLegProfile::Ring::Knee ? 4.0 : 10.0) * scaleY;
            MMatrix matrix;
            matrix[0][0] = radius.x;
            matrix[1][1] = L;
            matrix[2][2] = 2 * radius.z;
            PreparedBellRing ring(matrix);
            const bool distal = kind == SkirtLegProfile::Ring::Knee || kind == SkirtLegProfile::Ring::Extended;
            BellColliderSolver::prepareDistalEnd(ring, matrix, distal);
            const double a = MVector(matrix[0][0], matrix[0][1], matrix[0][2]).length();
            const double b = MVector(matrix[2][0], matrix[2][1], matrix[2][2]).length();
            const double expectedWidth = (std::min)(L, (std::max)(a, b));
            if (ring.distalEnd != distal || ring.distalLength != L || ring.distalWidth != expectedWidth)
                return false;
            MPointArray source;
            for (double distance : {-1.0, L * 0.5, L, L + ring.distalWidth * 0.5,
                                    L + ring.distalWidth, L + ring.distalWidth + 1, L + 2 * ring.distalWidth, 1.0})
                source.append(MPoint(0.2, distance, -0.0));
            const MVector offset = source[4] - ring.plane.orig;
            const MVector normal = ring.plane.normal;
            const double measured = (offset.x * normal.x + offset.y * normal.y) + offset.z * normal.z;
            if (measured != L + expectedWidth)
                return false;
            const std::vector<MVector> directions(source.length(), MVector(4, 0, 0));
            auto actual = source, expected = source;
            referenceDirectionalRelax(expected, ring, 0.65, 1, 6, true, directions, {});
            BellColliderSolver::relaxTowardRingBoundary(actual, ring, 0.65, 1, 6, true, directions, {});
            if (!directionalCompare(expected, actual))
                return false;
            for (unsigned int i : {0u, 7u})
                if (!sameBits(source[i].x, actual[i].x) || !sameBits(source[i].y, actual[i].y) ||
                    !sameBits(source[i].z, actual[i].z))
                    return false;
            if (distal)
                for (unsigned int i : {4u, 5u, 6u})
                    if (!sameBits(source[i].x, actual[i].x) || !sameBits(source[i].y, actual[i].y) ||
                        !sameBits(source[i].z, actual[i].z))
                        return false;
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
        MPointArray source;
        for (double distance : {L, L + width * 0.5, L + width, L + width + 1})
            source.append(MPoint(1, distance, -0.0));
        const MVector offset = source[2] - prepared.plane.orig;
        const MVector normal = prepared.plane.normal;
        const double measured = (offset.x * normal.x + offset.y * normal.y) + offset.z * normal.z;
        if (measured != L + width)
            return false;
        const double t = (measured - L) / width;
        const double gain = 1 - ((t * t) * (3 - 2 * t));
        if (t != 1 || gain != 0)
            return false;
        const auto topology = makeTopology(chain(4));
        const std::vector<MVector> directions(source.length(), MVector(4, 0, 0));
        auto actual = source, expected = source;
        referenceDirectionalRelax(expected, prepared, 0.65, 0, 4, true, directions, {});
        BellColliderSolver::relaxRowFaded(actual, prepared, 0.65, true, topology,
                                         std::vector<double>(topology.components.size(), 1.0), directions);
        if (!directionalCompare(expected, actual))
            return false;
        for (unsigned int i : {2u, 3u})
            if (!sameBits(actual[i].x, source[i].x) || !sameBits(actual[i].y, source[i].y) ||
                !sameBits(actual[i].z, source[i].z) || !sameBits(actual[i].w, source[i].w))
                return false;
    }

    MMatrix wide;
    wide[0][0] = 3;
    wide[1][1] = 2;
    wide[2][2] = 5;
    PreparedBellRing ring(wide);
    BellColliderSolver::prepareDistalEnd(ring, wide, true);
    if (!ring.distalEnd || ring.distalWidth != 2)
        return false;
    for (double nonfinite : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        ring.direction = MVector(nonfinite, 0, 0);
        BellColliderSolver::prepareDistalEnd(ring, wide, true);
        if (ring.distalEnd || std::isfinite(ring.distalLength))
            return false;
    }
    return true;
}

static bool extended_keeps_leg_length_support()
{
    MMatrix kneeMatrix;
    kneeMatrix[0][0] = kneeMatrix[2][2] = 2;
    kneeMatrix[1][1] = 2;
    MMatrix extendedMatrix = kneeMatrix;
    extendedMatrix[1][1] = 8;
    PreparedBellRing knee(kneeMatrix), extended(extendedMatrix);
    BellColliderSolver::prepareDistalEnd(knee, kneeMatrix, true);
    BellColliderSolver::prepareDistalEnd(extended, extendedMatrix, true);
    PreparedBellRing previous = extended;
    previous.distalEnd = false;
    const std::vector<MVector> direction(1, MVector(4, 0, 0));
    for (double distance : {5.0, 8.0, 9.0, 10.0, 11.0})
    {
        MPointArray source;
        source.append(MPoint(-0.25, distance, -0.0));
        const double measured = (source[0] - knee.translation) * knee.normal;
        if (!(measured > knee.distalLength + knee.distalWidth))
            return false;
        auto kneeOut = source, extendedOut = source, oldOut = source, expected = source;
        BellColliderSolver::relaxTowardRingBoundary(kneeOut, knee, 0.65, 0, 1, true, direction, {});
        BellColliderSolver::relaxTowardRingBoundary(extendedOut, extended, 0.65, 0, 1, true, direction, {});
        BellColliderSolver::relaxTowardRingBoundary(oldOut, previous, 0.65, 0, 1, true, direction, {});
        referenceDirectionalRelax(expected, extended, 0.65, 0, 1, true, direction, {});
        if (!compare(source, kneeOut, 0, 1) || !directionalCompare(expected, extendedOut))
            return false;
        if (distance <= extended.distalLength)
        {
            if (!compare(oldOut, extendedOut, 0, 1) || !((extendedOut[0] - source[0]).length() > 0))
                return false;
        }
        const double gain = distance <= 8 ? 1 : (distance >= 10 ? 0 : 0.5);
        const double correction = ((2.25 * 0.5) * 0.65) * gain;
        if (std::abs(extendedOut[0].x - (source[0].x + correction)) > 1e-12)
            return false;
        if (distance > 10 && !compare(source, extendedOut, 0, 1))
            return false;
    }

    const auto fixture = chain(5);
    const auto topology = makeTopology(fixture);
    for (double distance : {5.0, 9.0, 11.0})
    {
        MPointArray base;
        for (int i = 0; i < 5; ++i)
            base.append(MPoint(-0.25 + 0.1 * i, distance, 0.1));
        for (double smoothness : {0.0, 0.1})
            for (double follow : {0.0, 0.2})
            {
                BellRowInputs inputs;
                inputs.rings = {knee, extended};
                inputs.collision = 0.65f;
                inputs.capAtRingOrigin = true;
                inputs.smoothness = smoothness;
                inputs.followGain = follow;
                if (!runDirectionalRowCase(inputs, base, fixture))
                    return false;
                auto expected = referenceRow(inputs, base, fixture);
                BellRowOutputs actual;
                if (!BellColliderSolver::solveRow(inputs, base, topology, actual))
                    return false;
                for (unsigned int i = 0; i < base.length(); ++i)
                {
                    expected.points[i] += expected.directField.values[i] * 0.2;
                    actual.points[i] += actual.directField.values[i] * 0.2;
                }
                for (size_t r = 0; r < inputs.rings.size(); ++r)
                {
                    const double collision = r == 0 ? 1.0 : 0.4;
                    referenceDirectionalRelax(expected.points, inputs.rings[r], collision, 0, 5, true,
                                               referenceDirections(expected.points, base), {});
                    BellColliderSolver::relaxRowFaded(actual.points, inputs.rings[r], collision, true, topology,
                                                      actual.componentWeights[r],
                                                      referenceDirections(actual.points, base));
                    if (!directionalCompare(expected.points, actual.points))
                        return false;
                }
            }
        for (double smoothness : {0.0, 0.5})
        {
            const std::vector<PreparedBellRing> rings = {knee, extended};
            auto expected = base;
            for (const auto &ring : rings)
                referenceDirectionalRelax(expected, ring, 1, 0, 5, true,
                                           std::vector<MVector>(5, MVector(0, 0, 0)), {});
            if (smoothness > 0)
            {
                std::vector<MVector> displacements(5);
                for (unsigned int i = 0; i < 5; ++i)
                    displacements[i] = expected[i] - base[i];
                referenceSmooth(displacements, fixture, smoothness);
                for (unsigned int i = 0; i < 5; ++i)
                    expected[i] = base[i] + displacements[i];
                for (const auto &ring : rings)
                    referenceDirectionalRelax(expected, ring, 1, 0, 5, true, referenceDirections(expected, base), {});
            }
            if (!directionalCompare(expected, waistComposition(base, topology, rings, smoothness)))
                return false;
        }
    }

    auto partialFixture = chain(5);
    partialFixture.cuts = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    const auto partialTopology = makeTopology(partialFixture);
    const std::vector<double> weights = {0, 0.25, 0.5, 0.75, 1};
    const std::vector<double> fades = {1, 0.25, 0, 0.25, 1};
    for (double collision : {1.0, 0.4})
    {
        MPointArray base, expected;
        for (int i = 0; i < 5; ++i)
        {
            base.append(MPoint(-2, 9, 0.2));
            expected.append(MPoint(-0.2, 9, 0.3));
        }
        auto actual = expected;
        for (const auto &ring : {knee, extended})
        {
            const auto before = expected;
            auto after = before;
            referenceDirectionalRelax(after, ring, collision, 0, 5, true, referenceDirections(before, base), fades);
            for (unsigned int i = 0; i < 5; ++i)
            {
                const double f = fades[i];
                expected[i] = MPoint(before[i].x * (1 - f) + after[i].x * f,
                                     before[i].y * (1 - f) + after[i].y * f,
                                     before[i].z * (1 - f) + after[i].z * f, after[i].w);
            }
            BellColliderSolver::relaxRowFaded(actual, ring, collision, true, partialTopology, weights,
                                              referenceDirections(actual, base));
            if (!directionalCompare(expected, actual))
                return false;
        }
    }

    auto cutFixture = chain(8);
    cutFixture.u = {0, 0.125, 0.25, 0.375, 0.375, 0.5, 0.625, 0.75};
    cutFixture.cuts = {{3, 4}};
    const auto cutTopology = makeTopology(cutFixture);
    MPointArray cutBase;
    for (const MPoint point : {MPoint(2, 1.25, 0), MPoint(1.5, 1.25, 1.5), MPoint(0, 1.25, 2), MPoint(-1.5, 1.25, 1.5),
                               MPoint(-1.5, 1.25, 1.5), MPoint(-2, 1.25, 0), MPoint(-1.5, 1.25, -1.5), MPoint(0, 1.25, -2)})
        cutBase.append(point);
    MTransformationMatrix transform;
    transform.rotateBy(MQuaternion(MVector(0, 1, 0), MVector(0, std::cos(0.6), std::sin(0.6))), MSpace::kTransform);
    for (double length : {0.5, 1.5})
    {
        MMatrix matrix = transform.asMatrix();
        for (int column = 0; column < 3; ++column)
            matrix[1][column] *= length;
        PreparedBellRing ring(matrix);
        BellColliderSolver::prepareDistalEnd(ring, matrix, true);
        BellRowInputs input;
        input.rings = {ring};
        input.collision = 1;
        input.capAtRingOrigin = true;
        input.falloff = -1;
        input.smoothness = 0.1;
        input.followGain = 0.2;
        for (double m : {0.25, 0.5, 0.75, 1.0})
        {
            double low = 0, high = 1;
            for (int iteration = 0; iteration < 60; ++iteration)
            {
                const double t = (low + high) * 0.5;
                if (t * t * (3 - 2 * t) < m)
                    low = t;
                else
                    high = t;
            }
            input.contactBlendWidth = m == 1 ? 0 : 0.125 / ((low + high) * 0.5);
            std::vector<MPointArray> responses;
            std::vector<std::vector<double>> componentWeights;
            if (!BellColliderSolver::deformPoints(input, cutBase, cutTopology, {{2}}, responses, &componentWeights) ||
                responses.size() != 1 || componentWeights.size() != 1 || componentWeights[0].size() != 2 ||
                std::abs(componentWeights[0][0] - m) > 1e-12 || componentWeights[0][1] != 0)
                return false;
            const auto mask = referenceComponentWeights(input, ring, cutBase, cutFixture);
            if (!directionalCompare(referenceDeform(input, ring, cutBase, mask), responses[0]) ||
                !runDirectionalRowCase(input, cutBase, cutFixture))
                return false;
            auto rotationOnly = input;
            rotationOnly.collision = 0;
            std::vector<MPointArray> lifted;
            if (!BellColliderSolver::deformPoints(rotationOnly, cutBase, cutTopology, {{2}}, lifted))
                return false;
            double motion = 0;
            for (unsigned int i = 0; i < cutBase.length(); ++i)
                motion += (lifted[0][i] - cutBase[i]).length();
            if (!(motion > 1e-6))
                return false;
            auto still = cutBase;
            referenceDirectionalRelax(still, ring, 1, 0, 8, true,
                                       std::vector<MVector>(8, MVector(0, 0, 0)), {});
            if (length > 1)
            {
                const double distance = (cutBase[2] - ring.translation) * ring.normal;
                if (!(distance > ring.distalLength && distance < ring.distalLength + ring.distalWidth) ||
                    !((still[2] - cutBase[2]).length() > 1e-6))
                    return false;
            }
            std::vector<MPointArray> noRotation;
            if (!BellColliderSolver::deformPoints(input, cutBase, cutTopology, {{}}, noRotation) ||
                !directionalCompare(still, noRotation[0]) ||
                !directionalCompare(referenceDeform(input, ring, cutBase, std::vector<double>(8, 0)), noRotation[0]))
                return false;
        }
    }
    return true;
}

int main()
{
    struct Case
    {
        const char *name;
        bool (*run)();
    };
    const Case cases[] = {{"knee_end_row_preparation", kneeEndRowPreparation},
                          {"extended_keeps_leg_length_support", extended_keeps_leg_length_support},
                          {"merge_expectations", mergeExpectations},
                          {"merge_cancellation", mergeCancellation},
                          {"merge_scaling_rotation_bound", mergeScalingRotationBound},
                          {"bell_solve_merge", bellSolveMerge},
                          {"bell_solve_gate_bits", bellSolveGateBits},
                          {"1 smoothing_cut_edge", smoothingCutEdge},
                          {"2 smoothing_components", smoothingComponents},
                          {"3 local_follow", localFollow},
                          {"3 follow_summation_order", followSummationOrder},
                          {"12 invalid_row_inputs", invalidRowInputs},
                          {"3 merge_before_follow", mergeBeforeFollow},
                          {"4 follow_detour", followDetour},
                          {"5 direct_transfer", directTransfer},
                          {"6 one_way_and_insertions", oneWayAndInsertions},
                          {"7 long_tightness", longTightness},
                          {"8 waist_row", waistRow},
                          {"9 contact_gate", contactGate},
                          {"9 contact_blend_continuity", contactBlendContinuity},
                          {"10 relax_outside_gate", relaxOutsideGate},
                          {"11 relax_lane_independence", relaxLaneIndependence}};
    bool success = runDirectionalCases();
    for (const auto &test : cases)
    {
        if (test.run == bellSolveGateBits && bellGateExpected.empty())
        {
            std::cout << "SKIP " << test.name << ": pre-change binary baseline table is empty\n";
            continue;
        }
        const bool passed = test.run();
        std::cout << (passed ? "PASS " : "FAIL ") << test.name << "\n";
        success = passed && success;
    }
    return success ? 0 : 1;
}
