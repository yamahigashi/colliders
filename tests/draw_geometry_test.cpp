#include "colliderDrawGeometry.h"
#include "skirtRingFrames.h"
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
  require(frames.visibleMatrices().size() == 4);
  require(rings.update(frames.visibleMatrices(), 4));
  require(rings.geometry.lines[5].isEquivalent(
      rings.unitPoints[5] * frames.leftKnee, 1e-6));
  const SkirtRingFrames shortFrames(leftHip, leftKnee, leftHeel, rightHip,
                                    rightKnee, rightHeel, MVector(0.5, 2, 0.5),
                                    0, 0, false);
  require(shortFrames.visibleMatrices().size() == 2);

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
