#include "colliderDrawGeometry.h"
#include "skirtRingFrames.h"
#include "skirtLegBuild.h"
#include <climits>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>

static void require(bool value) {
  if (!value)
    throw std::runtime_error("draw geometry assertion failed");
}
static MPointArray pointArray(std::initializer_list<MPoint> points) {
  MPointArray result;
  for (const auto &point : points)
    result.append(point);
  return result;
}

static void requireEndpoints(const MPointArray &actual,
                             std::initializer_list<MPoint> expected) {
  require(actual.length() == expected.size());
  unsigned int i = 0;
  for (const auto &point : expected) {
    require(actual[i].x == point.x && actual[i].y == point.y &&
            actual[i].z == point.z && actual[i].w == point.w);
    ++i;
  }
}

static void testPanelCurves() {
  using Element = ColliderDraw::Curves::Element;
  const MPoint a(0, 0, 0), b(0, 1, 0), c(1, 0, 0), d(1, 1, 0);
  const MPoint e(2, 0, 0), f(2, 1, 0), g(3, 0, 0), h(3, 1, 0);
  const MPoint tip(4, 0, 0), left(4, 1, -1), right(4, 1, 1);
  std::vector<Element> patches = {
      {3, pointArray({a, b, c, d, tip, left}), 3, 2, {}},
      {9, pointArray({tip, right, g, h}), 2, 2, {}}};
  ColliderDraw::Curves curves;
  require(curves.update(patches));
  require(!curves.update(patches));
  require(curves.elements.size() == 2);
  require(curves.elements[0].panelId == 3 && curves.elements[1].panelId == 9);
  requireEndpoints(curves.elements[0].lines, {a, c, c, tip, b, d, d, left});
  requireEndpoints(curves.elements[1].lines, {tip, g, right, h});
  requireEndpoints(curves.lines,
                   {a, c, c, tip, b, d, d, left, tip, g, right, h});

  patches[1].points[3] = f;
  require(curves.update(patches));
  require(!curves.update(patches));
  requireEndpoints(curves.elements[0].lines, {a, c, c, tip, b, d, d, left});
  requireEndpoints(curves.elements[1].lines, {tip, g, right, f});

  std::vector<Element> single = {
      {3, pointArray({tip, left, a, b, c, d, tip, right}), 4, 2, {}}};
  require(curves.update(single));
  require(!curves.update(single));
  requireEndpoints(curves.lines,
                   {tip, a, a, c, c, tip, left, b, b, d, d, right});

  const MPointArray samePoints = pointArray({a, b, c, d, e, f, g, h});
  single = {{3, samePoints, 4, 2, {}}};
  require(curves.update(single));
  requireEndpoints(curves.lines, {a, c, c, e, e, g, b, d, d, f, f, h});
  patches = {{3, pointArray({a, b, c, d}), 2, 2, {}},
             {9, pointArray({e, f, g, h}), 2, 2, {}}};
  require(curves.update(patches));
  require(!curves.update(patches));
  requireEndpoints(curves.lines, {a, c, b, d, e, g, f, h});
  patches[1].panelId = 12;
  require(curves.update(patches));
  require(!curves.update(patches));
  requireEndpoints(curves.lines, {a, c, b, d, e, g, f, h});
  std::swap(patches[0], patches[1]);
  require(curves.update(patches));
  requireEndpoints(curves.lines, {e, g, f, h, a, c, b, d});
  require(curves.update(single));
  requireEndpoints(curves.lines, {a, c, c, e, e, g, b, d, d, f, f, h});
  single[0].numU = 2;
  single[0].numV = 4;
  require(curves.update(single));
  require(!curves.update(single));
  requireEndpoints(curves.lines, {a, e, b, f, c, g, d, h});

  single[0].points[0].x = -0.0;
  require(curves.update(single));
  require(!curves.update(single));
  require(std::signbit(curves.elements[0].points[0].x));
  require(std::signbit(curves.lines[0].x));
  single[0].points[0].w = 2.0;
  require(curves.update(single));
  require(!curves.update(single));
  require(curves.lines[0].w == 2.0);

  for (int invalidCase = 0; invalidCase < 5; ++invalidCase) {
    require(curves.update(patches));
    require(curves.lines.length() == 8);
    require(!curves.update(patches));
    auto invalid = patches;
    if (invalidCase == 0)
      invalid[1].numU = 3;
    else if (invalidCase == 1)
      invalid[1].numV = 0;
    else if (invalidCase == 2)
      invalid[1].points.clear();
    else if (invalidCase == 3) {
      invalid[1].numU = std::numeric_limits<unsigned int>::max();
      invalid[1].numV = std::numeric_limits<unsigned int>::max();
    } else
      invalid[1].points[0].x = std::numeric_limits<double>::infinity();
    require(curves.update(invalid));
    require(curves.elements.empty() && curves.lines.length() == 0);
    require(!curves.update(invalid));
    require(!curves.update(std::vector<Element>{}));
    require(curves.elements.empty() && curves.lines.length() == 0);
    require(curves.update(patches));
    require(curves.lines.length() == 8);
    require(!curves.update(patches));
    require(curves.update(std::vector<Element>{}));
    require(curves.elements.empty() && curves.lines.length() == 0);
    require(!curves.update(std::vector<Element>{}));
  }
  require(curves.update(patches));
  require(curves.update(MPointArray(), 0, 0));
  require(curves.elements.empty() && curves.lines.length() == 0);
  require(!curves.update(MPointArray(), 0, 0));
}
int main() {
  testPanelCurves();
  ColliderDraw::Rings rings;
  std::vector<MMatrix> matrices(1);
  require(rings.update(matrices, 4));
  require(rings.geometry.triangles.length() == 36);
  require(rings.geometry.lines.length() == 24);
  require(rings.unitPoints.length() == 9);
  const int expected[] = {0, 1, 2, 1, 5, 6, 1, 6, 2, 0, 2, 3, 2, 6, 7, 2, 7, 3,
                          0, 3, 4, 3, 7, 8, 3, 8, 4, 0, 4, 1, 4, 8, 5, 4, 5, 1};
  for (unsigned int i = 0; i < 36; ++i) {
    require(rings.indices[i] == expected[i]);
    const MPoint actual(rings.geometry.triangles[i]);
    require(actual.isEquivalent(rings.unitPoints[expected[i]], 1e-6));
  }
  require(rings.geometry.lines[18].isEquivalent(rings.unitPoints[4], 1e-6));
  require(rings.geometry.lines[19].isEquivalent(rings.unitPoints[1], 1e-6));
  require(!rings.update(matrices, 4));
  matrices[0][3][0] = 3;
  require(rings.update(matrices, 4));
  require(rings.geometry.lines[0].isEquivalent(MPoint(4, 0, 0), 1e-6));
  matrices[0][3][0] = 0;
  require(rings.update(matrices, 4));
  matrices.push_back(MMatrix());
  require(rings.update(matrices, 4));
  require(rings.geometry.triangles.length() == 72);
  require(rings.update(matrices, 8));
  require(rings.geometry.lines.length() == 96);
  require(rings.update({}, 8));
  require(rings.geometry.lines.length() == 0);

  // Reject hostile subdivision values before any cache construction and make
  // sure a previously valid cache can recover on the next valid update.
  require(rings.update(matrices, 4));
  require(rings.subdivision == 4);
  for (int invalid : { -1, 0, 2, 4097, INT_MAX }) {
    require(rings.update(matrices, invalid));
    require(rings.subdivision == 0);
    require(rings.matrices.empty());
    require(rings.farMultipliers.empty());
    require(rings.unitPoints.length() == 0);
    require(rings.indices.length() == 0);
    require(rings.geometry.triangles.length() == 0);
    require(rings.geometry.lines.length() == 0);
    require(!rings.update(matrices, invalid));
    require(rings.update(matrices, 4));
  }
  std::vector<std::array<double, 2>> mismatchedFar(1, {{1.0, 1.0}});
  require(rings.update(matrices, 4, mismatchedFar));
  require(rings.subdivision == 0);
  require(rings.matrices.empty());
  require(rings.geometry.lines.length() == 0);
  require(rings.update(matrices, 4));
  require(rings.update({}, 4096));
  require(rings.subdivision == 4096);
  require(rings.unitPoints.length() == 8193);
  require(rings.indices.length() == 4096 * 9);

  MMatrix leftHip, leftKnee, leftHeel, rightHip, rightKnee, rightHeel;
  leftKnee[3][1] = -2;
  leftHeel[3][1] = -5;
  rightKnee[3][1] = -4;
  rightHeel[3][1] = -9;
  const SkirtRingFrames frames(leftHip, leftKnee, leftHeel, rightHip, rightKnee,
                               rightHeel, MVector(0.5, 2, 0.5), 0, 0, true);
  require(frames.thighLength == 3 && frames.calfLength == 4);
  require(std::abs(yaxis(frames.leftKnee).length() - 6) < 1e-6);
  require(std::abs(yaxis(frames.rightKnee).length() - 6) < 1e-6);
  require(std::abs(yaxis(frames.leftHeel).length() - 14) < 1e-6);
  require(std::abs(yaxis(frames.rightHeel).length() - 14) < 1e-6);
  const SkirtLegProfile profile(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, 0.25,
                                0.75, 3, 4);
  const MMatrix joints[2][3] = {{leftHip, leftKnee, leftHeel},
                                {rightHip, rightKnee, rightHeel}};
  const short axes[2] = {0, 0};
  const double radii[8] = {1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6};
  std::vector<std::array<double, 2>> farMultipliers;
  const auto drawLegs = [&](double scaleY, bool longSkirt, double thighPosition,
                           double calfPosition) {
    Leg legs[2];
    std::vector<LegSegment> rest;
    buildLegs(joints, joints, axes, MVector(0.5, scaleY, 0.5), radii,
              thighPosition, calfPosition, longSkirt, legs, rest);
    std::vector<std::array<LegVec3, 4>> rows;
    farMultipliers.clear();
    for (const auto &leg : legs)
      for (const auto &segment : leg)
        legRingMatrices(segment.current, rows, farMultipliers);
    std::vector<MMatrix> result;
    for (const auto &ring : rows)
      result.push_back(legRingMatrix(ring));
    return result;
  };
  const auto longMatrices = drawLegs(2, true, 0.25, 0.75);
  require(longMatrices.size() == 8 && farMultipliers.size() == 8);
  require(rings.update(longMatrices, 4, farMultipliers));
  require(rings.geometry.triangles.length() == 288);
  require(rings.geometry.lines.length() == 192);
  require(!rings.update(longMatrices, 4, farMultipliers));
  for (size_t leg = 0; leg < 2; ++leg) {
    const MPoint hip = taxis(joints[leg][0]), knee = taxis(joints[leg][1]),
                 heel = taxis(joints[leg][2]);
    const MPoint starts[] = {hip, hip + (knee - hip) * 0.5, knee,
                            knee + (heel - knee) * 1.5};
    const MPoint ends[] = {starts[1], hip + (knee - hip) * 2,
                          starts[3], knee + (heel - knee) * 2};
    for (size_t segment = 0; segment < 4; ++segment) {
      const size_t index = leg * 4 + segment;
      const auto &matrix = longMatrices[index];
      require(taxis(matrix).isEquivalent(starts[segment], 1e-9));
      require((MPoint(0, 1, 0) * matrix).isEquivalent(ends[segment], 1e-9));
      require(std::abs(xaxis(matrix).length() -
                       0.5 * profile.stations[segment].x) < 1e-9);
      require(std::abs(zaxis(matrix).length() -
                       0.5 * profile.stations[segment].z) < 1e-9);
      require(rings.geometry.lines[index * 24 + 5].isEquivalent(
          MPoint(farMultipliers[index][0], 1, 0) * matrix, 1e-6));
      require(rings.geometry.lines[index * 24 + 8].isEquivalent(
          MPoint(0, 1, farMultipliers[index][1]) * matrix, 1e-6));
      require(std::abs(xaxis(matrix).length() * farMultipliers[index][0] -
                       0.5 * profile.stations[segment + 1].x) < 1e-9);
      require(std::abs(zaxis(matrix).length() * farMultipliers[index][1] -
                       0.5 * profile.stations[segment + 1].z) < 1e-9);
    }
  }
  const auto savedFar = farMultipliers;
  farMultipliers[0][0] *= 0.5;
  require(rings.update(longMatrices, 4, farMultipliers));
  require(rings.geometry.lines[0].isEquivalent(
      MPoint(1, 0, 0) * longMatrices[0], 1e-6));
  require(rings.geometry.lines[5].isEquivalent(
      MPoint(farMultipliers[0][0], 1, 0) * longMatrices[0], 1e-6));
  require(!rings.update(longMatrices, 4, farMultipliers));
  const auto unitLengthMatrices = drawLegs(1, true, 0.25, 0.75);
  require(!ColliderDraw::sameMatrices(longMatrices, unitLengthMatrices));
  for (size_t i = 0; i < farMultipliers.size(); ++i)
    for (int axis = 0; axis < 2; ++axis)
      require(std::abs(farMultipliers[i][axis] - savedFar[i][axis]) < 1e-9);
  const auto shortMatrices = drawLegs(2, false, 0.25, 0.75);
  require(shortMatrices.size() == 4 && farMultipliers.size() == 4);
  require(rings.update(shortMatrices, 4, farMultipliers));
  require(rings.geometry.triangles.length() == 144);
  require(rings.geometry.lines.length() == 96);
  require(rings.update(shortMatrices, 8, farMultipliers));
  require(rings.geometry.triangles.length() == 288);
  require(rings.geometry.lines.length() == 192);

  const SkirtLegProfile defaults(1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 3, 4);
  for (double s : {-1.0, 0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.05}) {
    require(defaults.fX(s, 7) == 1.0 && defaults.fZ(s, 7) == 1.0);
  }
  require(profile.fX(-0.01, 7) == 1.0);
  require(profile.fZ(1.05, 7) == 0.6);
  for (double thighPosition : {0.0, 1.0}) {
    for (double calfPosition : {0.0, 1.0}) {
      const SkirtLegProfile ties(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6,
                                 thighPosition, calfPosition, 3, 4);
      require(std::abs(ties.fX(ties.knee * thighPosition, 7) -
                       (thighPosition == 0 ? 1.2 : 0.7)) < 1e-12);
      require(std::abs(ties.fX(ties.knee + (1 - ties.knee) * calfPosition, 7) -
                       (calfPosition == 0 ? 0.7 : 0.5)) < 1e-12);
      const auto thin = drawLegs(2, false, thighPosition, calfPosition);
      require(thin.size() == 2);
      rings.update(thin, 4, farMultipliers);
      for (unsigned int i = 0; i < rings.geometry.lines.length(); ++i)
        require(std::isfinite(rings.geometry.lines[i].x) &&
                std::isfinite(rings.geometry.lines[i].y) &&
                std::isfinite(rings.geometry.lines[i].z));
    }
  }
  const SkirtLegProfile noThigh(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, 0.5,
                                0.5, 0, 4);
  const SkirtLegProfile noCalf(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, 0.5, 0.5,
                               3, 0);
  const SkirtLegProfile noLeg(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, 0, 0, 0,
                              0);
  require(noThigh.fX(0, 4) == 0.7);
  require(noCalf.fX(1, 3) == 0.5);
  require(noLeg.fX(0, 0) == 1 && noLeg.fZ(1.1, 0) == 1);
  const SkirtLegProfile tinyLeg(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, 0, 0,
                                0.4e-6, 0.4e-6);
  require(tinyLeg.fX(0, 1) == 1 && tinyLeg.fZ(1.1, 1) == 1);
  require(noThigh.forRing(0, SkirtLegProfile::Ring::Knee, 1).x == 1);
  const SkirtLegProfile toleranceChain(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6,
                                       0.5, 0.5, 1 - 1.5e-9, 1.5e-9);
  require(std::abs(toleranceChain.fX(toleranceChain.knee, 1) - 0.7) < 1e-12);
  require(profile.fX(0.5, 0.999e-6) == 1);
  require(profile.fZ(0.5, 0.999e-6) == 1);
  require(profile.fX(0.5, 1e-6) != 1);
  require(profile.parameter(0.5, SkirtLegProfile::Ring::Knee, 2) == 3.0 / 7.0);
  require(profile.parameter(0.5, SkirtLegProfile::Ring::Heel, 2) == 1);
  require(profile.parameter(0.9, SkirtLegProfile::Ring::Extended, 2) ==
          profile.knee);
  require(std::abs(profile.forRing(0.9, SkirtLegProfile::Ring::Extended, 2).x -
                   0.7) < 1e-12);

  ColliderDraw::Curves curves;
  MPointArray points;
  for (int i = 0; i < 6; ++i)
    points.append(MPoint(i, 0, 0));
  require(curves.update(points, 3, 2));
  require(curves.lines.length() == 8);
  const int lineIndices[] = {0, 2, 2, 4, 1, 3, 3, 5};
  for (unsigned int i = 0; i < 8; ++i)
    require(curves.lines[i] == points[lineIndices[i]]);
  require(!curves.update(points, 3, 2));
  points[0].x = 8;
  require(curves.update(points, 3, 2));
  require(curves.lines[0].x == 8);
  require(curves.update(points, 2, 3));
  require(curves.lines.length() == 6);
  require(curves.update(MPointArray(), 0, 0));
  require(curves.lines.length() == 0);
  std::cout << "draw geometry and cache tests passed\n";
}
