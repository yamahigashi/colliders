#pragma once

#include "utils.hpp"
#include <vector>

struct SkirtRingFrames {
  double thighLength;
  double calfLength;
  MMatrix leftKnee, rightKnee, leftHeel, rightHeel, leftExtended, rightExtended;
  bool longSkirt;

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

  std::vector<MMatrix> visibleMatrices() const {
    std::vector<MMatrix> result = {leftKnee, rightKnee};
    if (longSkirt) {
      result.push_back(leftHeel);
      result.push_back(rightHeel);
    }
    return result;
  }
};
