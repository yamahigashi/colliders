#pragma once

#include "skirtLegProfile.h"
#include "utils.hpp"
#include <vector>

struct SkirtRingFrames {
  double thighLength;
  double calfLength;
  MMatrix leftKnee, rightKnee, leftHeel, rightHeel, leftExtended, rightExtended;
  bool longSkirt;
  std::array<MPoint, 3> leftJoints, rightJoints;

  SkirtRingFrames(const MMatrix &leftHip, const MMatrix &leftKneeJoint,
                  const MMatrix &leftHeelJoint, const MMatrix &rightHip,
                  const MMatrix &rightKneeJoint, const MMatrix &rightHeelJoint,
                  const MVector &scale, short leftAxis, short rightAxis,
                  bool isLong)
      : longSkirt(isLong) {
    const MPoint lh = taxis(leftHip), lk = taxis(leftKneeJoint),
                 le = taxis(leftHeelJoint);
    const MPoint rh = taxis(rightHip), rk = taxis(rightKneeJoint),
                 re = taxis(rightHeelJoint);
    leftJoints = {{lh, lk, le}};
    rightJoints = {{rh, rk, re}};
    thighLength = ((lk - lh).length() + (rk - rh).length()) * 0.5;
    calfLength = ((le - lk).length() + (re - rk).length()) * 0.5;
    const MVector thighScale(scale.x, thighLength * scale.y, scale.z);
    leftKnee = createRingMatrix(leftHip, thighScale, leftAxis, &lk);
    rightKnee = createRingMatrix(rightHip, thighScale, rightAxis, &rk);
    if (longSkirt) {
      const MVector legScale(scale.x, (thighLength + calfLength) * scale.y,
                             scale.z);
      leftHeel = createRingMatrix(leftHip, legScale, leftAxis, &le);
      rightHeel = createRingMatrix(rightHip, legScale, rightAxis, &re);
      leftExtended = createRingMatrix(leftHip, legScale, leftAxis, &lk);
      rightExtended = createRingMatrix(rightHip, legScale, rightAxis, &rk);
    }
  }

  std::vector<MMatrix>
  visibleMatrices(const SkirtLegProfile &profile,
                  std::vector<std::array<double, 2>> &farMultipliers) const {
    std::vector<MMatrix> result;
    farMultipliers.clear();
    auto appendLeg = [&](const std::array<MPoint, 3> &joints,
                         const MMatrix &kneeFrame, const MMatrix &heelFrame) {
      const std::array<MPoint, 5> positions{
          {joints[0],
           joints[0] + (joints[1] - joints[0]) * profile.thighPosition,
           joints[1],
           joints[1] + (joints[2] - joints[1]) * profile.calfPosition,
           joints[2]}};
      const int segments = longSkirt ? 4 : 2;
      for (int i = 0; i < segments; ++i) {
        MMatrix matrix = i < 2 ? kneeFrame : heelFrame;
        const auto &startRadius = profile.stations[i];
        const auto &endRadius = profile.stations[i + 1];
        for (unsigned int c = 0; c < 3; ++c) {
          matrix[0][c] *= startRadius.x;
          matrix[2][c] *= startRadius.z;
        }
        set_maxis(matrix, 1, positions[i + 1] - positions[i]);
        matrix[3][0] = positions[i].x;
        matrix[3][1] = positions[i].y;
        matrix[3][2] = positions[i].z;
        result.push_back(matrix);
        const double farX = startRadius.x > 1e-12 ? endRadius.x / startRadius.x : 1.0;
        const double farZ = startRadius.z > 1e-12 ? endRadius.z / startRadius.z : 1.0;
        farMultipliers.push_back({{farX, farZ}});
      }
    };
    appendLeg(leftJoints, leftKnee, leftHeel);
    appendLeg(rightJoints, rightKnee, rightHeel);
    return result;
  }
};
