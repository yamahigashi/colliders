#pragma once

#include "skirtLegModel.h"
#include <array>
#include "utils.hpp"

inline LegVec3 legVector(const MVector &value)
{
    return {value.x, value.y, value.z};
}

inline LegVec3 legPoint(const MPoint &value)
{
    return {value.x, value.y, value.z};
}

inline MVector mayaLegVector(const LegVec3 &value)
{
    return MVector(value.x, value.y, value.z);
}

inline LegJoint legJoint(const MMatrix &matrix)
{
    LegJoint joint;
    for (int row = 0; row < 3; ++row)
        joint.rows[row] = {matrix[row][0], matrix[row][1], matrix[row][2]};
    joint.origin = legPoint(taxis(matrix));
    return joint;
}

inline bool isFinitePoint(const MPoint &point)
{
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && std::isfinite(point.w);
}

inline bool validMatrix(const MMatrix &matrix)
{
    for (unsigned int row = 0; row < 4; row++)
        for (unsigned int column = 0; column < 4; column++)
            if (!std::isfinite(matrix[row][column]))
                return false;
    const double determinant = matrix.det4x4();
    return std::isfinite(determinant) && std::fabs(determinant) >= 1e-8;
}

inline MMatrix legRingMatrix(const std::array<LegVec3, 4> &rows)
{
    MMatrix matrix;
    for (int row = 0; row < 4; ++row)
        for (int column = 0; column < 3; ++column)
            matrix[row][column] = rows[row][column];
    return matrix;
}

inline void buildLegs(const MMatrix (&current)[2][3], const MMatrix (&rest)[2][3], const short (&axes)[2],
                      const MVector &ringScale, const double (&radii)[8], double thighPosition, double calfPosition,
                      bool longSkirt, Leg (&legs)[2], std::vector<LegSegment> &restSegments)
{
    legs[0].segments.clear();
    legs[1].segments.clear();
    restSegments.clear();
    for (int side = 0; side < 2; side++)
    {
        bool finiteCurrent = true, finiteRest = true;
        for (int joint = 0; joint < 3; joint++)
        {
            finiteCurrent = finiteCurrent && isFinitePoint(taxis(current[side][joint]));
            finiteRest = finiteRest && isFinitePoint(taxis(rest[side][joint]));
        }
        if (!finiteRest)
            continue;
        const auto makeProfile = [&](const MMatrix matrices[3]) {
            return SkirtLegProfile(radii[0], radii[1], radii[2], radii[3], radii[4], radii[5], radii[6], radii[7],
                                   thighPosition, calfPosition, (taxis(matrices[1]) - taxis(matrices[0])).length(),
                                   (taxis(matrices[2]) - taxis(matrices[1])).length());
        };
        const SkirtLegProfile restProfile = makeProfile(rest[side]);
        for (int segment = 0; segment < (longSkirt ? 2 : 1); segment++)
        {
            LegSegmentPair pair;
            pair.legId = side;
            pair.segmentId = segment;
            pair.currentJoint = legJoint(current[side][segment]);
            pair.restJoint = legJoint(rest[side][segment]);
            const bool restOk =
                validMatrix(rest[side][segment]) && validMatrix(rest[side][segment + 1]) &&
                makeLegSegment(pair.restJoint.rows, pair.restJoint.origin, legPoint(taxis(rest[side][segment + 1])),
                               axes[side], legVector(ringScale), restProfile, segment == 1, pair.rest);
            if (restOk)
                restSegments.push_back(pair.rest);
            if (finiteCurrent && restOk && validMatrix(current[side][segment]) &&
                validMatrix(current[side][segment + 1]) &&
                makeLegSegment(pair.currentJoint.rows, pair.currentJoint.origin,
                               legPoint(taxis(current[side][segment + 1])), axes[side], legVector(ringScale),
                               makeProfile(current[side]), segment == 1, pair.current))
                legs[side].segments.push_back(pair);
        }
    }
    for (auto &leg : legs)
        legPrepareLeg(leg);
}
