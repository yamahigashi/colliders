#include "bellColliderSolver.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>

namespace
{
bool bits(double a, double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool bits(const MVector &a, const MVector &b)
{
    return bits(a.x, b.x) && bits(a.y, b.y) && bits(a.z, b.z);
}

bool bits(const MPoint &a, const MPoint &b)
{
    return bits(a.x, b.x) && bits(a.y, b.y) && bits(a.z, b.z) && bits(a.w, b.w);
}

bool same(const MPointArray &a, const MPointArray &b)
{
    if (a.length() != b.length())
        return false;
    for (unsigned int i = 0; i < a.length(); ++i)
        if (!bits(a[i], b[i]))
            return false;
    return true;
}

bool near(const MPoint &a, const MPoint &b)
{
    return (a - b).length() < 1e-12 && std::abs(a.w - b.w) < 1e-12;
}

BellRowTopology topology(const std::vector<double> &u, bool closed = false, double start = 0.0, double end = 1.0)
{
    BellRowTopology result;
    if (u.empty())
        return result;
    BellRowComponent component;
    component.closed = closed;
    component.startU = start;
    component.endU = end;
    result.components.push_back(component);
    result.outputCount = static_cast<unsigned int>(u.size());
    for (size_t i = 0; i < u.size(); ++i)
    {
        BellRowVertex vertex;
        vertex.materialU = u[i];
        vertex.previous = i ? static_cast<int>(i - 1) : closed ? static_cast<int>(u.size() - 1) : -1;
        vertex.next = i + 1 < u.size() ? static_cast<int>(i + 1) : closed ? 0 : -1;
        vertex.outputDuplicates.push_back(static_cast<unsigned int>(i));
        result.vertices.push_back(vertex);
    }
    return result;
}

MPointArray points(const MPoint &point, unsigned int count = 1)
{
    MPointArray result;
    for (unsigned int i = 0; i < count; ++i)
        result.append(point);
    return result;
}

bool projection()
{
    const auto t = topology({0.5});
    const auto anchor = points(MPoint(0, 0, 0));
    auto base = points(MPoint(0, 2, 0));
    MPointArray out;
    for (const MPoint &x : {MPoint(-0.0, 0, 0), MPoint(1, 0, 0), MPoint(0, 2, -0.0)})
    {
        const auto current = points(x);
        if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, nullptr, nullptr, base, current, t, out) ||
            !same(out, current))
            return false;
    }
    if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, nullptr, nullptr, base, points(MPoint(3, 4, 0)), t,
                                                 out) ||
        !near(out[0], MPoint(1.2, 1.6, 0)))
        return false;
    for (double d : {0.0, 1e-15, 3.0})
        if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, nullptr, nullptr, anchor,
                                                     points(MPoint(d, 0, 0)), t, out) ||
            !same(out, anchor))
            return false;
    base = points(MPoint(0, 1, 0));
    auto above = points(MPoint(0, -1, 0));
    const double distances[] = {1.0, 1.5, 2.0};
    const double shares[] = {0.0, 0.25, 0.5};
    for (unsigned int i = 0; i < 3; ++i)
    {
        const double k = shares[i];
        const MVector direction = MVector(1.0 - k, k, 0).normal();
        if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, &above, &t, base,
                                                     points(MPoint(distances[i], 0, 0)), t, out) ||
            !near(out[0], MPoint(direction)))
            return false;
    }
    for (const MPoint &c : {MPoint(0, 0, 0), MPoint(0, -1e-12, 0), MPoint(1, 0, 0)})
    {
        above = points(c);
        if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, &above, &t, base, points(MPoint(2, 0, 0)), t,
                                                     out) ||
            !near(out[0], MPoint(1, 0, 0)))
            return false;
    }
    base = points(MPoint(0, 1e-13, 0));
    above = points(MPoint(0, -1, 0));
    if (!BellColliderSolver::projectSuspendedRow(anchor, anchor, t, &above, &t, base, points(MPoint(2e-13, 0, 0)), t,
                                                 out) ||
        !bits(out[0].y, 0.0) || !bits(out[0].x, 1e-13))
        return false;
    BellRowCorrespondence cache;
    if (!BellColliderSolver::buildRowCorrespondence(t, t, cache))
        return false;
    MPointArray cached;
    return BellColliderSolver::projectSuspendedRow(anchor, anchor, cache, &above, &cache, base,
                                                   points(MPoint(2e-13, 0, 0)), cached) &&
           same(out, cached);
}

bool compareTransfers(const BellRowTopology &source, const BellRowTopology &destination,
                      const std::vector<MVector> &values, const std::vector<MVector> &expected,
                      const std::vector<bool> &expectedMask)
{
    BellRowCorrespondence cache;
    std::vector<MVector> a(31), b;
    std::vector<bool> am(19), bm;
    BellDirectField field;
    field.topology = source;
    field.values = values;
    BellRowTransfer direct;
    if (!BellColliderSolver::buildRowCorrespondence(source, destination, cache) ||
        !BellColliderSolver::transferRowValues(values, source, destination, a, am) ||
        !BellColliderSolver::transferRowValues(values, cache, b, bm) ||
        !BellColliderSolver::transferDirectField(field, destination, direct) || a.size() != expected.size() ||
        a.size() != b.size() || a.size() != direct.values.size() || am != expectedMask || am != bm)
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!bits(a[i], b[i]) || !bits(a[i], direct.values[i]) || !bits(a[i], expected[i]))
            return false;
    return true;
}

bool transfers()
{
    const auto closed = topology({0.125, 0.375, 0.625, 0.875}, true);
    const auto destination = topology({0.0, 0.125, 0.25, 0.5, 0.75, 0.875}, true);
    const std::vector<MVector> values = {MVector(-0.0, -0.0, -0.0), MVector(2, 0, 0), MVector(4, 0, 0),
                                         MVector(6, 0, 0)};
    if (!compareTransfers(
            closed, destination, values,
            {MVector(3, 0, 0), values[0], MVector(1, 0, 0), MVector(3, 0, 0), MVector(5, 0, 0), values[3]},
            std::vector<bool>(6, true)))
        return false;
    auto open = topology({0.25, 0.75});
    if (!compareTransfers(open, topology({0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0}), {values[0], values[1]},
                          {values[0], values[0], values[0], MVector(1, 0, 0), values[1], values[1], values[1]},
                          std::vector<bool>(7, true)))
        return false;
    auto seam = topology({0.0, 0.5}, false, 0, 0.5);
    BellRowComponent second;
    second.startU = 0.5;
    second.endU = 1.0;
    seam.components.push_back(second);
    auto tail = topology({0.5, 1.0}, false, 0.5, 1.0);
    for (auto v : tail.vertices)
    {
        v.componentId = 1;
        if (v.previous >= 0)
            v.previous += 2;
        if (v.next >= 0)
            v.next += 2;
        v.outputDuplicates[0] += 2;
        seam.vertices.push_back(v);
    }
    seam.outputCount = 4;
    for (size_t i = 0; i < seam.vertices.size(); ++i)
    {
        seam.vertices[i].side.bank = i % 2 ? 1 : -1;
        seam.vertices[i].side.seamIndex = 7;
    }
    const std::vector<MVector> seamValues = {MVector(1, 0, 0), MVector(2, 0, 0), MVector(3, 0, 0), MVector(4, 0, 0)};
    if (!compareTransfers(seam, seam, seamValues, seamValues, std::vector<bool>(4, true)))
        return false;
    auto single = topology({0.5}, false, 0.5, 0.6);
    single.vertices[0].side.bank = -1;
    single.vertices[0].side.seamIndex = 8;
    if (!compareTransfers(seam, single, seamValues, {MVector(0, 0, 0)}, {false}))
        return false;
    single.vertices[0].side.bank = 0;
    if (!compareTransfers(seam, single, seamValues, {seamValues[1]}, {true}))
        return false;
    std::swap(seam.components[0], seam.components[1]);
    for (auto &vertex : seam.vertices)
        vertex.componentId = 1 - vertex.componentId;
    if (!compareTransfers(seam, single, seamValues, {seamValues[2]}, {true}))
        return false;
    if (!compareTransfers(topology({0.1, 0.2}, false, 0.1, 0.2), topology({0.7}), {MVector(0, 0, 0), MVector(0, 0, 0)},
                          {MVector(0, 0, 0)}, {false}))
        return false;
    if (!compareTransfers(BellRowTopology(), destination, {}, std::vector<MVector>(6, MVector(0, 0, 0)),
                          std::vector<bool>(6, false)))
        return false;
    return compareTransfers(closed, BellRowTopology(), values, {}, {});
}

bool unmatched()
{
    const auto target = topology({0.75}, false, 0.7, 0.8);
    const auto source = topology({0.25}, false, 0.2, 0.3);
    const auto a = points(MPoint(0, 0, 0));
    const auto b = points(MPoint(0, 1, 0));
    const auto x = points(MPoint(3, 4, 0));
    MPointArray out;
    if (!BellColliderSolver::projectSuspendedRow(a, a, source, nullptr, nullptr, b, x, target, out) || !same(out, x))
        return false;
    if (!BellColliderSolver::projectSuspendedRow(a, a, target, &x, &source, b, x, target, out) ||
        !near(out[0], MPoint(0.6, 0.8, 0)))
        return false;
    return BellColliderSolver::smoothRowDisplacements(a, a, source, b, x, target, b, b, target, out) && same(out, x);
}

bool smoothing()
{
    const auto t = topology({0.5});
    MPointArray base[4], current[4];
    for (unsigned int i = 0; i < 4; ++i)
    {
        base[i] = points(MPoint(i * i, i * i * i, -0.0));
        current[i] = base[i];
    }
    current[1][0].x += 8;
    for (unsigned int i = 1; i < 3; ++i)
        if (!BellColliderSolver::smoothRowDisplacements(base[i - 1], current[i - 1], t, base[i], current[i], t,
                                                        base[i + 1], current[i + 1], t, current[i]))
            return false;
    if (current[0][0].x - base[0][0].x != 0 || current[1][0].x - base[1][0].x != 4 ||
        current[2][0].x - base[2][0].x != 1 || current[3][0].x - base[3][0].x != 0)
        return false;
    BellRowCorrespondence cache;
    if (!BellColliderSolver::buildRowCorrespondence(t, t, cache))
        return false;
    MPointArray uncached, cached;
    if (!BellColliderSolver::smoothRowDisplacements(base[0], current[0], t, base[1], current[1], t, base[2], current[2],
                                                    t, uncached) ||
        !BellColliderSolver::smoothRowDisplacements(base[0], current[0], cache, base[1], current[1], base[2],
                                                    current[2], cache, cached) ||
        !same(uncached, cached))
        return false;
    for (double zero : {-0.0, 0.0})
    {
        const auto b = points(MPoint(-zero, 3, zero));
        const auto x = points(MPoint(zero, 3, -zero));
        if (!BellColliderSolver::smoothRowDisplacements(b, x, t, b, x, t, b, x, t, cached) || !same(cached, x))
            return false;
    }
    // A cancelled displacement must copy the incoming point even when it differs from base.
    const auto b = points(MPoint(0, 0, 0));
    const auto x = points(MPoint(2, 0, 0));
    const auto opposite = points(MPoint(-2, 0, 0));
    return BellColliderSolver::smoothRowDisplacements(b, opposite, t, b, x, t, b, opposite, t, cached) &&
           same(cached, x);
}

bool invalidInputs()
{
    const auto t = topology({0.5});
    auto broken = t;
    broken.vertices[0].next = 0;
    std::vector<MVector> out;
    std::vector<bool> matched;
    BellRowCorrespondence cache;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (BellColliderSolver::buildRowCorrespondence(broken, t, cache) ||
        BellColliderSolver::buildRowCorrespondence(t, broken, cache) ||
        BellColliderSolver::transferRowValues({MVector(0, 0, 0)}, cache, out, matched) ||
        BellColliderSolver::transferRowValues({}, t, t, out, matched) ||
        BellColliderSolver::transferRowValues({MVector(nan, 0, 0)}, t, t, out, matched))
        return false;
    const auto p = points(MPoint(1, 0, 0));
    const auto bad = points(MPoint(1, 0, 0, nan));
    MPointArray result, empty;
    if (BellColliderSolver::projectSuspendedRow(p, p, t, &p, nullptr, p, p, t, result) ||
        BellColliderSolver::projectSuspendedRow(p, p, t, nullptr, &t, p, p, t, result) ||
        BellColliderSolver::projectSuspendedRow(p, p, t, &p, &broken, p, p, t, result) ||
        BellColliderSolver::projectSuspendedRow(p, p, t, &bad, &t, p, p, t, result) ||
        BellColliderSolver::projectSuspendedRow(p, empty, t, nullptr, nullptr, p, p, t, result) ||
        BellColliderSolver::projectSuspendedRow(p, p, t, nullptr, nullptr, bad, p, t, result) ||
        BellColliderSolver::smoothRowDisplacements(p, bad, t, p, p, t, p, p, t, result) ||
        BellColliderSolver::smoothRowDisplacements(p, p, t, p, p, broken, p, p, t, result) ||
        BellColliderSolver::smoothRowDisplacements(p, p, t, p, empty, t, p, p, t, result))
        return false;
    const double big = std::numeric_limits<double>::max();
    const auto positive = points(MPoint(big, 0, 0));
    const auto negative = points(MPoint(-big, 0, 0));
    return !BellColliderSolver::smoothRowDisplacements(negative, positive, t, p, p, t, p, p, t, result) &&
           !BellColliderSolver::projectSuspendedRow(negative, p, t, nullptr, nullptr, p, positive, t, result);
}

bool branchRelax()
{
    const auto t = topology({0.25, 0.5, 0.75});
    const auto base = points(MPoint(-2, 1, 0.2), 3);
    const auto chord = points(MPoint(-0.2, 1, 0.3), 3);
    MMatrix farMatrix, rightMatrix;
    farMatrix[3][0] = 100;
    rightMatrix[0][0] = rightMatrix[2][2] = 2;
    const PreparedBellRing heel(farMatrix), left(farMatrix), right(rightMatrix);
    const std::vector<PreparedBellRing> branchRings[] = {
        {heel, heel, left, right}, {heel, heel, left}, {heel, heel, right}, {heel, heel}};
    const std::vector<std::vector<double>> branchWeights[] = {
        {{1}, {1}, {1}, {0.5}}, {{1}, {1}, {1}}, {{1}, {1}, {0.5}}, {{1}, {1}}};
    const double bendLeft = 1.0, bendRight = 0.5, tightness = 0.5;
    const double pi[] = {(1 - bendLeft) * (1 - bendRight), (1 - bendLeft) * bendRight, bendLeft * (1 - bendRight),
                         bendLeft * bendRight};
    MPointArray mixed;
    bool first = true;
    for (unsigned int b = 0; b < 4; ++b)
    {
        if (!(pi[b] > 0))
            continue;
        auto current = chord;
        for (size_t r = 0; r < branchRings[b].size(); ++r)
            BellColliderSolver::relaxRowFaded(current, branchRings[b][r], r < 2 ? 1.0 : 1.0 - tightness, true, t,
                                              branchWeights[b][r], BellColliderSolver::rowDirections(current, base));
        if (!same(current, chord))
            return false;
        if (first)
            mixed.setLength(current.length());
        for (unsigned int i = 0; i < current.length(); ++i)
            mixed[i] = first ? current[i] * pi[b] : mixed[i] + current[i] * pi[b];
        first = false;
    }
    if (!same(mixed, chord))
        return false;
    auto wrong = chord;
    const auto &firstWeights = branchWeights[2];
    for (size_t r = 0; r < branchRings[0].size(); ++r)
        BellColliderSolver::relaxRowFaded(wrong, branchRings[0][r], r < 2 ? 1.0 : 1.0 - tightness, true, t,
                                          r < firstWeights.size() ? firstWeights[r] : std::vector<double>(),
                                          BellColliderSolver::rowDirections(wrong, base));
    return (wrong[0] - chord[0]).length() > 1e-6;
}

bool ringOrder()
{
    const auto t = topology({0.5});
    const double radius = 6.06, height = 58, epsilon = 0.1;
    const double x = std::sqrt(radius * radius - epsilon * epsilon);
    const double angle = 4.0 * std::acos(-1.0) / 180.0;
    const double c = std::cos(angle), s = std::sin(angle);
    const auto rotate = [&](double px, double py, double pz) { return MPoint(px * c - py * s, px * s + py * c, pz); };
    const MPoint b(0, height, -12);
    const MPoint left = rotate(x, height, epsilon), right = rotate(-x, height, epsilon);
    const MVector dl = left - b, dr = right - b;
    const double w = dl * dl / (dl * dl + dr * dr);
    const auto base = points(b);
    const auto start =
        points(MPoint(w * left.x + (1 - w) * right.x, w * left.y + (1 - w) * right.y, w * left.z + (1 - w) * right.z));
    MMatrix matrix;
    matrix[0][0] = radius * c;
    matrix[0][1] = radius * s;
    matrix[1][0] = -s;
    matrix[1][1] = c;
    matrix[2][2] = radius;
    const PreparedBellRing first(matrix);
    const MPoint origin = rotate(-x - 0.5 * radius, 0, 0);
    matrix[3][0] = origin.x;
    matrix[3][1] = origin.y;
    const PreparedBellRing second(matrix);
    auto ordered = start, reversed = start, stale = start, reference = start;
    const auto initialDirections = BellColliderSolver::rowDirections(start, base);
    for (const auto &ring : {first, second})
    {
        BellColliderSolver::relaxRowFaded(ordered, ring, 1, true, t, {1},
                                          BellColliderSolver::rowDirections(ordered, base));
        BellColliderSolver::relaxRowFaded(stale, ring, 1, true, t, {1}, initialDirections);
        BellColliderSolver::relaxTowardRingBoundary(reference, ring, 1, 0, 1, true,
                                                    BellColliderSolver::rowDirections(reference, base), {});
    }
    for (const auto &ring : {second, first})
        BellColliderSolver::relaxRowFaded(reversed, ring, 1, true, t, {1},
                                          BellColliderSolver::rowDirections(reversed, base));
    return same(ordered, reference) && (ordered[0] - stale[0]).length() > 1e-6 &&
           (ordered[0] - reversed[0]).length() > 1e-6;
}

bool pipeline()
{
    const auto t = topology({0.0, 0.25, 0.5, 0.75}, true);
    BellRowInputs input;
    MMatrix matrix;
    input.rings.push_back(PreparedBellRing(matrix));
    input.collision = 1;
    input.capAtRingOrigin = true;
    input.followGain = 0.5;
    input.smoothness = 0.1;
    MPointArray base[4], current[4];
    BellRowOutputs solutions[4];
    for (unsigned int r = 0; r < 4; ++r)
    {
        base[r] = points(MPoint(5 + r * r, r, 0), 4);
        BellColliderSolver::roundMeshPoints(base[r]);
        if (!BellColliderSolver::solveRow(input, base[r], t, solutions[r]))
            return false;
        current[r] = solutions[r].points;
        BellColliderSolver::roundMeshPoints(current[r]);
        BellRowTransfer follow;
        if (!BellColliderSolver::transferDirectField(solutions[r].directField, t, follow))
            return false;
        for (unsigned int i = 0; i < current[r].length(); ++i)
            current[r][i] += follow.values[i] * 0.5;
        BellColliderSolver::relaxRowFaded(current[r], input.rings[0], 1, true, t, solutions[r].componentWeights[0],
                                          BellColliderSolver::rowDirections(current[r], base[r]));
    }
    const auto relax = [&](unsigned int r) {
        BellColliderSolver::relaxRowFaded(current[r], input.rings[0], 1, true, t, solutions[r].componentWeights[0],
                                          BellColliderSolver::rowDirections(current[r], base[r]));
    };
    const auto suspend = [&]() {
        for (unsigned int r = 2; r < 4; ++r)
        {
            if (!BellColliderSolver::projectSuspendedRow(current[r - 1], base[r - 1], t, &current[r - 2], &t, base[r],
                                                         current[r], t, current[r]))
                return false;
            relax(r);
        }
        return true;
    };
    if (!suspend())
        return false;
    for (unsigned int iteration = 0; iteration < 2; ++iteration)
    {
        for (unsigned int r = 1; r < 3; ++r)
        {
            if (!BellColliderSolver::smoothRowDisplacements(base[r - 1], current[r - 1], t, base[r], current[r], t,
                                                            base[r + 1], current[r + 1], t, current[r]))
                return false;
            relax(r);
        }
        if (!suspend())
            return false;
    }
    for (unsigned int r = 0; r < 4; ++r)
    {
        if (!same(base[r], current[r]))
            return false;
        for (const auto &v : solutions[r].directField.values)
            if (!bits(v, MVector(0, 0, 0)))
                return false;
    }
    // Moving the finalised anchor must affect the next row immediately.
    const auto single = topology({0.5});
    auto a = points(MPoint(1, 0, 0));
    const auto ab = points(MPoint(0, 0, 0));
    const auto b = points(MPoint(0, 1, 0));
    auto x = points(MPoint(4, 0, 0));
    if (!BellColliderSolver::projectSuspendedRow(a, ab, single, nullptr, nullptr, b, x, single, x) ||
        !near(x[0], MPoint(2, 0, 0)))
        return false;
    a[0] = MPoint(0, 0, 0);
    return BellColliderSolver::projectSuspendedRow(a, ab, single, nullptr, nullptr, b, x, single, x) &&
           near(x[0], MPoint(1, 0, 0));
}
} // namespace

int main()
{
    const struct Test
    {
        const char *name;
        bool (*run)();
    } tests[] = {{"projection", projection},        {"transfers", transfers},
                 {"unmatched", unmatched},          {"smoothing", smoothing},
                 {"invalid inputs", invalidInputs}, {"branch relax", branchRelax},
                 {"ring order", ringOrder},         {"pipeline", pipeline}};
    for (const auto &test : tests)
        if (!test.run())
        {
            std::cerr << test.name << " failed\n";
            return 1;
        }
    std::cout << "row suspension passed\n";
    return 0;
}
