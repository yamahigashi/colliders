#include "colliderDrawGeometry.h"
#include "skirtRingFrames.h"
#include <climits>
#include <iostream>
#include <stdexcept>

static void require(bool value) {
  if (!value)
    throw std::runtime_error("draw geometry assertion failed");
}
int main() {
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
  std::vector<std::array<double, 2>> farMultipliers;
  const auto longMatrices = frames.visibleMatrices(profile, farMultipliers);
  require(longMatrices.size() == 8 && farMultipliers.size() == 8);
  require(rings.update(longMatrices, 4, farMultipliers));
  require(rings.geometry.triangles.length() == 288);
  require(rings.geometry.lines.length() == 192);
  require(!rings.update(longMatrices, 4, farMultipliers));
  for (size_t leg = 0; leg < 2; ++leg) {
    const auto &joints = leg == 0 ? frames.leftJoints : frames.rightJoints;
    const MPoint positions[] = {
        joints[0], joints[0] + (joints[1] - joints[0]) * 0.25, joints[1],
        joints[1] + (joints[2] - joints[1]) * 0.75, joints[2]};
    for (size_t segment = 0; segment < 4; ++segment) {
      const size_t index = leg * 4 + segment;
      const auto &matrix = longMatrices[index];
      require(taxis(matrix).isEquivalent(positions[segment], 1e-6));
      require((MPoint(0, 1, 0) * matrix)
                  .isEquivalent(positions[segment + 1], 1e-6));
      require(std::abs(xaxis(matrix).length() -
                       0.5 * profile.stations[segment].x) < 1e-6);
      require(std::abs(zaxis(matrix).length() -
                       0.5 * profile.stations[segment].z) < 1e-6);
      const MPoint endpoint(farMultipliers[index][0], 1, 0);
      require(rings.geometry.lines[index * 24 + 5].isEquivalent(
          endpoint * matrix, 1e-6));
      require(rings.geometry.lines[index * 24 + 8].isEquivalent(
          MPoint(0, 1, farMultipliers[index][1]) * matrix, 1e-6));
      require(std::abs(xaxis(matrix).length() * farMultipliers[index][0] -
                       0.5 * profile.stations[segment + 1].x) < 1e-6);
      require(std::abs(zaxis(matrix).length() * farMultipliers[index][1] -
                       0.5 * profile.stations[segment + 1].z) < 1e-6);
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
  const SkirtRingFrames unitLengthFrames(leftHip, leftKnee, leftHeel, rightHip,
                                         rightKnee, rightHeel,
                                         MVector(0.5, 1, 0.5), 0, 0, true);
  require(ColliderDraw::sameMatrices(
      longMatrices, unitLengthFrames.visibleMatrices(profile, farMultipliers)));
  require(farMultipliers == savedFar);
  const SkirtRingFrames shortFrames(leftHip, leftKnee, leftHeel, rightHip,
                                    rightKnee, rightHeel, MVector(0.5, 2, 0.5),
                                    0, 0, false);
  const auto shortMatrices =
      shortFrames.visibleMatrices(profile, farMultipliers);
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
      const auto thin = shortFrames.visibleMatrices(ties, farMultipliers);
      require(thin.size() == 4);
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
