#include "skirtLegBuild.h"
#include <maya/MQuaternion.h>
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
int checks = 0;

void require(bool condition, const char *message)
{
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance = 1e-9)
{
    if (!std::isfinite(expected))
        require(!std::isfinite(actual), "nonfinite classification mismatch");
    else
    {
        if (!(std::isfinite(actual) && std::abs(actual - expected) <= tolerance))
            std::cerr << "scalar mismatch: actual " << actual << " expected " << expected << std::endl;
        require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, "scalar mismatch");
    }
}

void near(const LegVec3 &actual, const MVector &expected, double tolerance = 1e-9)
{
    if (!(std::abs(actual.x - expected.x) <= tolerance && std::abs(actual.y - expected.y) <= tolerance &&
          std::abs(actual.z - expected.z) <= tolerance))
        std::cerr << "vector mismatch: actual (" << actual.x << ", " << actual.y << ", " << actual.z << ") expected ("
                  << expected.x << ", " << expected.y << ", " << expected.z << ")" << std::endl;
    near(actual.x, expected.x, tolerance);
    near(actual.y, expected.y, tolerance);
    near(actual.z, expected.z, tolerance);
}

// Maya's half-turn quaternion keeps a rounded half-angle (w up to ~1e-8), so antiparallel
// frames agree with the exact half turn only to about 1e-7.
void checkFrame(const MMatrix &joint, short axis, const MPoint &target, double tolerance = 1e-9)
{
    const auto input = legJoint(joint);
    LegVec3 x, y, z;
    legRingFrame(input.rows, input.origin, axis, legPoint(target), x, y, z);
    const MMatrix expected = createRingMatrix(joint, MVector(1, 1, 1), axis, &target);
    near(x, xaxis(expected), tolerance);
    near(y, yaxis(expected), tolerance);
    near(z, zaxis(expected), tolerance);
}

bool mayaJointRotation(const MMatrix &joint, MMatrix &rotation)
{
    rotation.setToIdentity();
    for (short row = 0; row < 3; row++)
    {
        const MVector axis = getAxis(joint, row);
        const double norm = axis.length();
        if (!finiteVector(legVector(axis)) || !std::isfinite(norm) || norm < 1e-8)
            return false;
        set_maxis(rotation, row, axis / norm);
    }
    const double det = rotation.det3x3();
    return std::isfinite(det) && std::fabs(det) >= 0.5;
}

void checkJoint(const MMatrix &joint)
{
    MMatrix expected;
    LegMat3 actual;
    const bool valid = legJointRotation(legJoint(joint).rows, actual);
    require(valid == mayaJointRotation(joint, expected), "joint rotation validity");
    if (!valid)
        return;
    near(actual.det3x3(), expected.det3x3());
    const auto inverse = actual.inverse();
    const MMatrix expectedInverse = expected.inverse();
    for (int row = 0; row < 3; ++row)
    {
        near(actual[row], maxis(expected, row));
        near(inverse[row], maxis(expectedInverse, row));
    }
}

void testFrames()
{
    std::mt19937 random(72401);
    std::uniform_real_distribution<double> distribution(-1, 1);
    for (int i = 0; i < 1000; ++i)
    {
        MMatrix joint = MQuaternion(distribution(random) * 3.0,
                                    MVector(distribution(random), distribution(random), distribution(random)))
                            .asMatrix();
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                joint[row][column] += distribution(random) * 0.2;
        for (int column = 0; column < 3; ++column)
            joint[3][column] = distribution(random) * 5;
        const MPoint target = taxis(joint) + MVector(distribution(random), distribution(random), distribution(random));
        checkFrame(joint, static_cast<short>(i % 6), target);
        checkFrame(joint, static_cast<short>(i % 6), taxis(joint) + getAxis(joint, static_cast<short>(i % 6)) * 2);
        checkJoint(joint);
        const MVector from(distribution(random), distribution(random), distribution(random));
        const MVector to(distribution(random), distribution(random), distribution(random));
        const MVector vector(distribution(random), distribution(random), distribution(random));
        if (from.normal() * to.normal() >= -1.0 + 1e-6)
            near(legRotateVector(legVector(vector), legVector(from), legVector(to)),
                 vector.rotateBy(MQuaternion(from, to)));
    }
    for (const MVector &axis : {MVector(1, 0, 0), MVector(0, 1, 0), MVector(0, 0, 1), MVector(1, 2, 3), MVector(1, 1, 1),
                                MVector(0.5, 0.5, 0)})
        for (const MVector &vector : {MVector(1, 2, 3), MVector(-0.3, 0.7, 0.2)})
            near(legRotateVector(legVector(vector), legVector(axis), legVector(-axis)),
                 vector.rotateBy(MQuaternion(axis, -axis)), 1e-7);
    const MVector from(1, 0, 0), vector(1, 2, 3);
    for (const MVector &to : {MVector(1, 0, 0), MVector(0, 1, 0), MVector(-1, 0.002, 0), MVector(1, 1e-5, 0),
                             MVector(0, 0, 0)})
        near(legRotateVector(legVector(vector), legVector(from), legVector(to)),
             vector.rotateBy(MQuaternion(from, to)));
    MMatrix identity;
    for (short axis = 0; axis < 6; ++axis)
    {
        const MVector raw = getAxis(identity, axis);
        checkFrame(identity, axis, MPoint(-raw * 2), 1e-7);
        checkFrame(identity, axis, MPoint(raw * 2));
        checkFrame(identity, axis, MPoint(maxis(identity, (axis + 1) % 3) * 2));
        for (double length : {0.999999e-6, 1e-6, 1.000001e-6, 0.99e-5, 1.01e-5, 2e-5})
            checkFrame(identity, axis, MPoint(0, -length, 0));
    }
    for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
            for (int z = -1; z <= 1; ++z)
            {
                const MVector raw(x, y, z);
                if (raw.length() == 0)
                    continue;
                MMatrix joint;
                set_maxis(joint, 0, raw);
                const MVector basis = std::abs(x) < std::abs(y) ? MVector(1, 0, 0) : MVector(0, 1, 0);
                set_maxis(joint, 1, raw ^ basis);
                checkFrame(joint, 0, MPoint(-raw), 1e-7);
                checkFrame(joint, 0, MPoint(raw));
            }
    for (const MVector &raw : {MVector(1, 2, 3), MVector(0.5, 0.5, 0),
                              MVector(-0.38023345746295101, 0.89354277136842208, 0.41442833956952163)})
    {
        MMatrix joint;
        set_maxis(joint, 0, raw);
        checkFrame(joint, 0, MPoint(-raw), 1e-7);
        checkFrame(joint, 0, MPoint(raw));
    }
    for (double epsilon : {0.0, 5e-6, 2e-5})
        for (double sign : {-1.0, 1.0})
            checkFrame(identity, 0, MPoint(-std::cos(epsilon), sign * std::sin(epsilon), 0), 1e-7);
    MMatrix zero;
    set_maxis(zero, 0, MVector(0, 0, 0));
    checkFrame(zero, 0, MPoint(0, -2, 0));
    checkJoint(zero);
    for (double determinant : {0.5 - 1e-12, 0.5, 0.5 + 1e-12})
        for (double sign : {-1.0, 1.0})
        {
            MMatrix joint;
            joint[1][0] = std::sqrt(1 - determinant * determinant);
            joint[1][1] = sign * determinant;
            checkJoint(joint);
        }
    for (short axis : {-1, -2, -3})
        checkFrame(identity, axis, MPoint(1, 2, 3));
}

void testBuild()
{
    const short axes[2] = {0, 0};
    const double radii[8] = {1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6};
    const MVector scale(0.5, 1.3, 0.8);
    const MPoint restPoints[3] = {MPoint(1, 2, 3), MPoint(1, -2, 3), MPoint(1, -5, 7)};
    const MPoint currentPoints[3] = {MPoint(1, 2, 3), MPoint(2, -4, 4), MPoint(4, -6, 8)};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int fixture = 0; fixture < 13; ++fixture)
    {
        MMatrix current[2][3], rest[2][3];
        for (int side = 0; side < 2; ++side)
            for (int joint = 0; joint < 3; ++joint)
            {
                set_maxis(rest[side][joint], 3, MVector(restPoints[joint]));
                set_maxis(current[side][joint], 3, MVector(fixture == 3 ? currentPoints[joint] : restPoints[joint]));
            }
        if (fixture == 4 || fixture == 12)
            rest[0][2][3][0] = nan;
        else if (fixture == 5 || fixture == 8)
            current[0][2][3][0] = nan;
        else if (fixture == 6)
            current[0][2] = current[0][1];
        else if (fixture == 7)
            rest[0][2] = rest[0][1];
        else if (fixture == 9)
            current[0][2][0][0] = 0;
        else if (fixture == 10)
            rest[0][0][0][0] = 0;
        else if (fixture == 11)
            current[0][0][3][3] = nan;
        const bool longSkirt = fixture != 2 && fixture != 8 && fixture != 12;
        const double thighPosition = fixture == 1 ? 1.0 : 0.25;
        Leg legs[2];
        std::vector<LegSegment> restSegments;
        buildLegs(current, rest, axes, scale, radii, thighPosition, 0.75, longSkirt, legs, restSegments);
        const std::size_t expectedLegs[] = {2, 2, 1, 2, 0, 0, 1, 1, 0, 1, 1, 1, 0};
        const std::size_t expectedRest[] = {2, 2, 1, 2, 0, 2, 2, 1, 1, 2, 1, 2, 0};
        const std::size_t otherSide = longSkirt ? 2 : 1;
        require(legs[0].segments.size() == expectedLegs[fixture], "current branch");
        require(legs[1].segments.size() == otherSide, "other side retained");
        require(restSegments.size() == expectedRest[fixture] + otherSide, "rest branch");
        if (fixture == 3)
        {
            require(legs[0].segments[0].current.length != legs[0].segments[0].rest.length, "separate current profile");
            near(legs[0].segments[0].current.stations[1].z, legs[0].segments[0].current.length * 0.25);
            near(legs[0].segments[1].current.stations[1].z, legs[0].segments[1].current.length * 0.75);
        }
        if (fixture < 4)
            require(legs[0].segments[0].current.stations.size() == (fixture == 1 ? 2u : 3u), "merged stations");
        for (const auto &leg : legs)
            for (const auto &pair : leg)
            {
                const auto &joints = current[pair.legId];
                const int id = pair.segmentId;
                near(pair.current.origin, MVector(taxis(joints[id])));
                near(pair.current.length, (taxis(joints[id + 1]) - taxis(joints[id])).length() * 1.3);
            }
        if (fixture == 5 || fixture == 9)
            require(legLocalRadius(restSegments, {1.2, -4, 5}) > 0, "rest radius retained");
        if (fixture == 11)
            require(!validMatrix(current[0][0]), "matrix homogeneous component finite gate");
    }
    MPoint point(1, 2, 3);
    point.w = nan;
    require(!isFinitePoint(point), "point homogeneous component finite gate");
    point.w = 2;
    require(isFinitePoint(point), "finite point weight");
    near(legPoint(point), MVector(1, 2, 3));
}

}

int main()
{
    try
    {
        testFrames();
        testBuild();
        std::cout << "skirt leg frame: 2 test groups, " << checks << " checks passed\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "skirt leg frame: " << error.what() << " after " << checks << " checks" << std::endl;
        return 1;
    }
}
