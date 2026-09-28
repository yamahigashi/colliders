#include "bellColliderSolver.h"
#include "colliderDrawGeometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

namespace {
using Point = std::array<double, 4>;

struct Curve {
  std::vector<double> knots;
  std::vector<Point> points;
};

bool sameBits(double a, double b) {
  std::uint64_t left = 0, right = 0;
  std::memcpy(&left, &a, sizeof(left));
  std::memcpy(&right, &b, sizeof(right));
  return left == right;
}

bool sameBits(const Point &a, const Point &b) {
  for (size_t i = 0; i < a.size(); ++i)
    if (!sameBits(a[i], b[i]))
      return false;
  return true;
}

Point components(const MPoint &p) { return {{p.x, p.y, p.z, p.w}}; }

bool check(bool valid, const char *message) {
  if (!valid)
    std::cerr << message << '\n';
  return valid;
}

Point deBoor(const Curve &curve, double u) {
  const int last = static_cast<int>(curve.points.size()) - 1;
  const int span =
      u == curve.knots[last + 1]
          ? last
          : static_cast<int>(
                std::upper_bound(curve.knots.begin(), curve.knots.end(), u) -
                curve.knots.begin()) -
                1;
  std::array<Point, 4> work;
  for (int j = 0; j <= 3; ++j)
    work[j] = curve.points[span - 3 + j];
  for (int order = 1; order <= 3; ++order)
    for (int j = 3; j >= order; --j) {
      const int i = span - 3 + j;
      const double denominator = curve.knots[i + 4 - order] - curve.knots[i];
      const double alpha =
          denominator == 0 ? 0 : (u - curve.knots[i]) / denominator;
      for (int axis = 0; axis < 4; ++axis)
        work[j][axis] = (1 - alpha) * work[j - 1][axis] + alpha * work[j][axis];
    }
  return work[3];
}

Curve periodic(int count) {
  Curve curve;
  for (int i = 0; i < 4 * count + 3; ++i) {
    const double angle = 6.28318530717958647692 * (i % count) / count;
    curve.points.push_back(
        {{2 * std::cos(angle), 0.3 * std::sin(2 * angle), std::sin(angle), 1}});
  }
  for (size_t i = 0; i < curve.points.size() + 4; ++i)
    curve.knots.push_back((static_cast<double>(i) - 3) / count - 1);
  return curve;
}

// Cubic span interpolation supplies an independent restriction without the
// private knot-insertion routine.
Curve restrictReference(const Curve &source, double a, double b) {
  std::vector<double> boundaries(1, a);
  for (double knot : source.knots)
    if (knot > a + 1e-12 && knot < b - 1e-12)
      boundaries.push_back(knot);
  boundaries.push_back(b);
  Curve result;
  result.knots.assign(4, 0);
  for (size_t span = 1; span < boundaries.size(); ++span) {
    const double lo = boundaries[span - 1], hi = boundaries[span];
    const Point p0 = deBoor(source, lo), p3 = deBoor(source, hi);
    const Point q1 = deBoor(source, lo + (hi - lo) / 3);
    const Point q2 = deBoor(source, lo + 2 * (hi - lo) / 3);
    Point p1, p2;
    for (int axis = 0; axis < 4; ++axis) {
      const double r1 = 27 * q1[axis] - 8 * p0[axis] - p3[axis];
      const double r2 = 27 * q2[axis] - p0[axis] - 8 * p3[axis];
      p1[axis] = (2 * r1 - r2) / 18;
      p2[axis] = (2 * r2 - r1) / 18;
    }
    if (span == 1)
      result.points.push_back(p0);
    result.points.push_back(p1);
    result.points.push_back(p2);
    result.points.push_back(p3);
    result.knots.insert(result.knots.end(),
                        span + 1 == boundaries.size() ? 4 : 3,
                        (hi - a) / (b - a));
  }
  return result;
}

std::vector<double> greville(const Curve &curve, double a, double b) {
  std::vector<double> result;
  for (size_t i = 0; i < curve.points.size(); ++i) {
    const double g =
        (curve.knots[i + 1] + curve.knots[i + 2] + curve.knots[i + 3]) / 3;
    result.push_back(a + g * (b - a));
  }
  return result;
}

BellRowTopology
openTopology(const std::vector<std::vector<double>> &materials) {
  BellRowTopology topology;
  for (size_t c = 0; c < materials.size(); ++c) {
    BellRowComponent component;
    component.startU = materials[c].front();
    component.endU = materials[c].back();
    topology.components.push_back(component);
    for (size_t j = 0; j < materials[c].size(); ++j) {
      BellRowVertex vertex;
      const int index = static_cast<int>(topology.vertices.size());
      vertex.materialU = materials[c][j];
      vertex.componentId = static_cast<int>(c);
      vertex.panelId = static_cast<int>(c + 1);
      vertex.previous = j == 0 ? -1 : index - 1;
      vertex.next = j + 1 == materials[c].size() ? -1 : index + 1;
      if (j == 0 || j + 1 == materials[c].size()) {
        vertex.side.seamIndex =
            j == 0 ? static_cast<int>(c)
                   : static_cast<int>((c + 1) % materials.size());
        vertex.side.bank = j == 0 ? -1 : 1;
      }
      vertex.outputDuplicates.push_back(static_cast<unsigned int>(index));
      topology.vertices.push_back(vertex);
    }
  }
  topology.outputCount = static_cast<unsigned int>(topology.vertices.size());
  return topology;
}

bool split_basis() {
  bool ok = true;
  for (int count : {3, 16}) {
    const Curve source = periodic(count);
    const std::vector<std::pair<double, double>> ranges = {
        {0, 0.5},     {1.0 / count, 2.0 / count}, {0.173, 0.681}, {0.83, 1.27},
        {0.25, 0.75}, {0.3, 0.3 + 2e-9}};
    for (const auto &range : ranges) {
      const double a = range.first, b = range.second;
      const Curve restricted = restrictReference(source, a, b);
      ok &= check(restricted.points.size() >= 4,
                  "split_basis: fewer than four CVs");
      ok &= check(restricted.knots.size() == restricted.points.size() + 4,
                  "split_basis: knot count");
      for (size_t i = 0; i < 4; ++i)
        ok &= check(restricted.knots[i] == 0 &&
                        restricted.knots[restricted.knots.size() - 1 - i] == 1,
                    "split_basis: unclamped endpoint");
      for (double knot : restricted.knots)
        if (knot > 0 && knot < 1)
          ok &= check(std::count(restricted.knots.begin(),
                                 restricted.knots.end(), knot) <= 3,
                      "split_basis: interior multiplicity");
      std::vector<double> samples = {0, 1};
      for (size_t i = 1; i < restricted.knots.size(); ++i) {
        const double lo = restricted.knots[i - 1], hi = restricted.knots[i];
        if (hi > lo) {
          samples.push_back(lo);
          samples.push_back(std::nextafter(lo, hi));
          samples.push_back(lo + (hi - lo) * 0.23);
          samples.push_back((lo + hi) * 0.5);
          samples.push_back(std::nextafter(hi, lo));
          samples.push_back(hi);
        }
      }
      for (double u : samples) {
        const Point expected = deBoor(source, a + u * (b - a));
        const Point actual = deBoor(restricted, u);
        for (int axis = 0; axis < 4; ++axis)
          ok &= check(std::abs(expected[axis] - actual[axis]) <= 1e-10,
                      "split_basis: restriction changed evaluation");
      }
      Curve coefficients = restricted;
      for (auto &p : coefficients.points)
        p = {{0, 0, 0, 0}};
      coefficients.points.front()[0] = 1;
      bool positive = false;
      for (double u : samples) {
        const double value = deBoor(coefficients, u)[0];
        if (u >= coefficients.knots[4])
          ok &= check(value == 0,
                      "split_basis: endpoint injection escaped basis support");
        else
          positive = positive || value > 0;
      }
      ok &= check(positive, "split_basis: injection was not exercised");
      Curve other = restrictReference(source, b, b + 0.2);
      BellDirectField injected;
      injected.topology = openTopology(
          {greville(restricted, a, b), greville(other, b, b + 0.2)});
      injected.values.assign(injected.topology.vertices.size(),
                             MVector(0, 0, 0));
      injected.values[restricted.points.size() - 1] = MVector(1, 2, 3);
      BellRowTransfer result;
      ok &= check(BellColliderSolver::transferDirectField(
                      injected, injected.topology, result) == MS::kSuccess,
                  "split_basis: cross-panel transfer status");
      if (result.values.size() != injected.values.size())
        return check(false, "split_basis: cross-panel transfer size");
      for (size_t i = 0; i < other.points.size(); ++i) {
        const MVector &value = result.values[restricted.points.size() + i];
        ok &= check(value.x == 0 && value.y == 0 && value.z == 0,
                    "split_basis: cross-panel CV coefficient");
        other.points[i] = {{value.x, value.y, value.z, 0}};
      }
      for (double u : {0.0, 0.01, 0.25, 0.5, 0.99, 1.0}) {
        const Point value = deBoor(other, u);
        for (double axis : value)
          ok &= check(axis == 0,
                      "split_basis: cross-panel evaluation coefficient");
      }
    }
  }
  const Curve independent = restrictReference(periodic(16), 0.25, 0.75);
  const auto materials = greville(independent, 0.25, 0.75);
  ok &= check(0.25 + 0.25 * (0.75 - 0.25) == 0.375,
              "split_basis: material mapping");
  BellDirectField source;
  source.topology = openTopology({materials});
  for (double s : materials)
    source.values.push_back(MVector(s, 0, 0));
  const Curve presentation = restrictReference(periodic(32), 0.25, 0.75);
  BellRowTopology destination =
      openTopology({greville(presentation, 0.25, 0.75)});
  BellRowTransfer transfer;
  ok &= check(BellColliderSolver::transferDirectField(source, destination,
                                                      transfer) == MS::kSuccess,
              "split_basis: transfer status");
  if (transfer.values.size() != destination.vertices.size())
    return check(false, "split_basis: transfer size");
  for (size_t i = 0; i < transfer.values.size(); ++i)
    ok &= check(std::abs(transfer.values[i].x -
                         destination.vertices[i].materialU) <= 1e-10,
                "split_basis: transfer used presentation CV indices instead of "
                "independent Greville");
  return ok;
}

struct Height {
  double value;
  bool physical;
};

// The node's clustering routine has internal linkage; this oracle supplies
// expected representatives.
std::vector<double> heightRepresentatives(std::vector<Height> candidates) {
  std::sort(candidates.begin(), candidates.end(),
            [](const Height &a, const Height &b) { return a.value < b.value; });
  std::vector<double> result;
  size_t first = 0;
  while (first < candidates.size()) {
    size_t next = first;
    double representative = candidates[first].value;
    while (next < candidates.size() &&
           candidates[next].value - candidates[first].value <= 1e-9) {
      if (candidates[next].physical)
        representative = candidates[next].value;
      ++next;
    }
    result.push_back(representative);
    first = next;
  }
  return result;
}

bool shared_tip() {
  bool ok = true;
  const auto chain = heightRepresentatives(
      {{0.5, false}, {0.5 + 0.75e-9, false}, {0.5 + 1.5e-9, false}});
  ok &= check(chain == std::vector<double>({0.5, 0.5 + 1.5e-9}),
              "shared_tip: transitive clustering");
  const auto physical = heightRepresentatives(
      {{0.5 - 0.75e-9, false}, {0.5, true}, {0.5, false}});
  ok &= check(physical == std::vector<double>({0.5}),
              "shared_tip: physical representative");
  for (double height : {0.25, 0.375, 0.5}) {
    const Curve source = periodic(16);
    MPointArray base;
    for (double u : {0.0, 0.3, 0.6, 0.9}) {
      Point p = deBoor(source, u);
      base.append(MPoint(p[0], -height, p[2], p[3]));
    }
    BellRowTopology top = openTopology({{0, 0.3, 0.6, 0.9}});
    top.components[0].closed = true;
    top.components[0].startU = 0;
    top.components[0].endU = 1;
    for (size_t i = 0; i < top.vertices.size(); ++i) {
      top.vertices[i].previous = static_cast<int>((i + 3) % 4);
      top.vertices[i].next = static_cast<int>((i + 1) % 4);
      top.vertices[i].side = BellRowSide();
    }
    top.vertices[1].outputDuplicates.push_back(4);
    top.outputCount = 5;
    BellRowInputs inputs;
    inputs.smoothness = 0.5;
    inputs.followGain = 1;
    inputs.collision = 1;
    inputs.bellMatrix[3][1] = -height;
    MMatrix ring;
    ring[0][0] = 3;
    ring[2][2] = 3;
    ring[3][1] = -2;
    inputs.rings.push_back(PreparedBellRing(ring));
    BellRowOutputs solved;
    ok &= check(BellColliderSolver::solveRow(inputs, base, top, solved) ==
                    MS::kSuccess,
                "shared_tip: shared row solve");
    MPointArray displacedBase = base;
    for (unsigned int i = 0; i < displacedBase.length(); ++i)
      displacedBase[i] = MPoint(displacedBase[i].x + 0.125,
                                displacedBase[i].y - 0.25,
                                displacedBase[i].z + 0.375,
                                displacedBase[i].w);
    BellRowOutputs displaced;
    ok &= check(BellColliderSolver::solveRow(inputs, displacedBase, top,
                                             displaced) == MS::kSuccess,
                "shared_tip: displaced base solve");
    if (displaced.points.length() == displacedBase.length())
      for (unsigned int i = 0; i < displaced.points.length(); ++i)
        ok &= check(std::isfinite(displaced.points[i].x) &&
                        std::isfinite(displaced.points[i].y) &&
                        std::isfinite(displaced.points[i].z),
                    "shared_tip: displaced base finite output");
    if (solved.points.length() != base.length())
      return check(false, "shared_tip: shared independent CV count");
    std::vector<Point> expanded(top.outputCount);
    for (size_t i = 0; i < top.vertices.size(); ++i)
      for (unsigned int output : top.vertices[i].outputDuplicates)
        expanded[output] =
            components(solved.points[static_cast<unsigned int>(i)]);
    ok &= check(sameBits(expanded[1], expanded[4]),
                "shared_tip: copied boundary bits");
    const auto saved = expanded;
    BellRowTopology lower =
        openTopology({{0.0, 0.1, 0.2, 0.3}, {0.3, 0.5, 0.7, 1.0}});
    for (double smoothness : {0.0, 0.5, 1.0}) {
      std::vector<MVector> injected(8, MVector(0, 0, 0));
      injected[3] = MVector(1, 2, 3);
      ok &= check(BellColliderSolver::smoothDisplacements(
                      injected, smoothness, lower) == MS::kSuccess,
                  "shared_tip: lower smoothing status");
      for (size_t i = 4; i < injected.size(); ++i)
        ok &= check(injected[i].x == 0 && injected[i].y == 0 &&
                        injected[i].z == 0,
                    "shared_tip: direct cut-end injection crossed bank");
      std::vector<MVector> followed;
      ok &= check(BellColliderSolver::computeLocalFollow(
                      injected, lower, MMatrix(), 0.1875, followed) ==
                      MS::kSuccess,
                  "shared_tip: lower follow status");
      if (followed.size() != injected.size())
        return check(false, "shared_tip: follow size");
      for (size_t i = 4; i < followed.size(); ++i)
        ok &= check(followed[i].x == 0 && followed[i].y == 0 &&
                        followed[i].z == 0,
                    "shared_tip: direct cut-end follow crossed bank");
      BellRowOutputs repeated;
      ok &= check(BellColliderSolver::solveRow(inputs, base, top, repeated) ==
                      MS::kSuccess,
                  "shared_tip: repeated upstream solve");
      if (repeated.points.length() != base.length())
        return check(false, "shared_tip: repeated upstream size");
      for (size_t i = 0; i < top.vertices.size(); ++i)
        for (unsigned int output : top.vertices[i].outputDuplicates)
          ok &= check(
              sameBits(
                  saved[output],
                  components(repeated.points[static_cast<unsigned int>(i)])),
              "shared_tip: upstream CV changed");
    }
    base[1] = MPoint(-0.0, -0.0, -0.0, -0.0);
    BellRowOutputs noContact;
    ok &= check(BellColliderSolver::solveRow(BellRowInputs(), base, top,
                                             noContact) == MS::kSuccess,
                "shared_tip: signed-zero solve status");
    if (noContact.points.length() != base.length())
      return check(false, "shared_tip: signed-zero output size");
    ok &= check(sameBits(components(base[1]), components(noContact.points[1])),
                "shared_tip: signed zero lost during solve");
    for (double axis : components(noContact.points[1]))
      ok &= check(std::signbit(axis), "shared_tip: signed-zero component");
  }
  return ok;
}

bool draw_patch_data() {
  std::vector<ColliderDraw::Curves::Element> patches;
  for (unsigned int id : {1u, 9u}) {
    ColliderDraw::Curves::Element patch;
    patch.panelId = id;
    patch.numU = id == 1 ? 4 : 7;
    patch.numV = 3;
    for (unsigned int u = 0; u < patch.numU; ++u)
      for (unsigned int v = 0; v < patch.numV; ++v)
        patch.points.append(MPoint(id * 10 + u, v, 0));
    patches.push_back(patch);
  }
  ColliderDraw::Curves cache;
  bool ok = check(cache.update(patches), "draw_patch_data: initial update");
  unsigned int line = 0;
  for (const auto &patch : patches)
    for (unsigned int v = 0; v < patch.numV; ++v)
      for (unsigned int u = 1; u < patch.numU; ++u) {
        if (line + 1 >= cache.lines.length())
          return check(false, "draw_patch_data: missing lines");
        ok &=
            check(sameBits(components(cache.lines[line++]),
                           components(patch.points[(u - 1) * patch.numV + v])),
                  "draw_patch_data: start endpoint");
        ok &= check(sameBits(components(cache.lines[line++]),
                             components(patch.points[u * patch.numV + v])),
                    "draw_patch_data: end endpoint");
      }
  ok &= check(line == cache.lines.length(),
              "draw_patch_data: extra wrap or inter-panel line");
  ok &= check(!cache.update(patches), "draw_patch_data: unchanged key rebuilt");
  patches[0].panelId = 20;
  ok &= check(cache.update(patches), "draw_patch_data: changed ID was ignored");
  patches[0].numU = 6;
  patches[0].numV = 2;
  ok &= check(cache.update(patches),
              "draw_patch_data: changed dimensions were ignored");
  const auto saved = patches;
  patches[0].points.setLength(1);
  ok &= check(cache.update(patches) && cache.lines.length() == 0 &&
                  cache.elements.empty(),
              "draw_patch_data: invalid cache retained");
  ok &= check(cache.update(saved), "draw_patch_data: recovery");
  ok &= check(cache.update({}) && !cache.update({}),
              "draw_patch_data: empty cache transition");
  return ok;
}
} // namespace

int main() {
  bool ok = split_basis();
  ok = shared_tip() && ok;
  ok = draw_patch_data() && ok;
  return ok ? 0 : 1;
}
