#include "skirtLegModel.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

static double testCapsuleDistance(const LegVec3& point, const LegContactGeometry& geometry)
{
    LegCapsuleColumns columns;
    columns.prepare(geometry);
    return legCapsuleDistance(point, columns);
}

static double testCapsuleDistance(const Leg& leg, const LegVec3& point)
{
    return testCapsuleDistance(point, legPrepareContactGeometry(leg, 0.0));
}

static void require(bool value, const char *message)
{
    if (!value)
        throw std::runtime_error(message);
}

static void near(double actual, double expected)
{
    require(std::isfinite(actual) && std::abs(actual - expected) <= 1e-9, "scalar mismatch");
}

static void near(const LegVec3 &actual, const LegVec3 &expected)
{
    near(actual.x, expected.x);
    near(actual.y, expected.y);
    near(actual.z, expected.z);
}

static SkirtLegProfile profile(const LegVec3 (&points)[3], double thighPosition = 0.25)
{
    return SkirtLegProfile(1.2, 0.8, 0.7, 0.9, 1.1, 1.3, 0.5, 0.6, thighPosition, 0.75,
                           (points[1] - points[0]).length(), (points[2] - points[1]).length());
}

static void buildFixture(const LegVec3 (&current)[3], const LegVec3 (&rest)[3], bool longSkirt, double thighPosition,
                         Leg &leg, std::vector<LegSegment> &restSegments)
{
    leg.segments.clear();
    restSegments.clear();
    bool finiteCurrent = true, finiteRest = true;
    for (int i = 0; i < 3; ++i)
    {
        finiteCurrent = finiteCurrent && finiteVector(current[i]);
        finiteRest = finiteRest && finiteVector(rest[i]);
    }
    if (!finiteRest)
        return;
    const SkirtLegProfile restProfile = profile(rest, thighPosition);
    for (int i = 0; i < (longSkirt ? 2 : 1); ++i)
    {
        LegSegmentPair pair;
        pair.legId = 0;
        pair.segmentId = i;
        pair.restJoint.origin = rest[i];
        pair.currentJoint.origin = current[i];
        const LegVec3 scale(0.5, 1.3, 0.8);
        const bool restOk =
            makeLegSegment(pair.restJoint.rows, rest[i], rest[i + 1], 0, scale, restProfile, i == 1, pair.rest);
        if (restOk)
            restSegments.push_back(pair.rest);
        if (finiteCurrent && restOk &&
            makeLegSegment(pair.currentJoint.rows, current[i], current[i + 1], 0, scale,
                           profile(current, thighPosition), i == 1, pair.current))
            leg.segments.push_back(pair);
    }
    legPrepareLeg(leg);
}

static void checkSurface(const LegSegment &segment)
{
    std::vector<std::array<LegVec3, 4>> rows;
    std::vector<std::array<double, 2>> far;
    legRingMatrices(segment, rows, far);
    require(rows.size() + 1 == segment.stations.size() && far.size() == rows.size(), "ring count");
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        for (int end = 0; end < 2; ++end)
        {
            const double z = segment.stations[i + end].z;
            const LegVec3 center = rows[i][3] + rows[i][1] * end;
            const LegSectionFrame frame = legSectionFrame(segment, z);
            near(center, frame.center);
            for (int step = 0; step < 32; ++step)
            {
                const double theta = step * (6.28318530717958647692 / 32);
                const double x = std::cos(theta) * (end ? far[i][0] : 1.0);
                const double y = std::sin(theta) * (end ? far[i][1] : 1.0);
                const LegVec3 point = rows[i][0] * x + rows[i][1] * end + rows[i][2] * y + rows[i][3];
                near(legLocal(segment, point).eta, 0.0);
                near(point, legSurfacePoint(segment, z, theta));
                near(point, frame.center + frame.x * std::cos(theta) + frame.z * std::sin(theta));
            }
        }
    }
    near(legSectionFrame(segment, -1.0).center, segment.origin);
    near(legSectionFrame(segment, segment.length + 1.0).center, segment.origin + segment.axis * segment.length);
}

static void testDrawingFixtures()
{
    const LegVec3 rest[3] = {{1, 2, 3}, {1, -2, 3}, {1, -5, 7}};
    const LegVec3 current[3] = {{1, 2, 3}, {2, -4, 4}, {4, -6, 8}};
    for (int fixture = 0; fixture < 4; ++fixture)
    {
        Leg leg;
        std::vector<LegSegment> restSegments;
        const bool longSkirt = fixture != 2;
        buildFixture(fixture == 3 ? current : rest, rest, longSkirt, fixture == 1 ? 1.0 : 0.25, leg, restSegments);
        require(leg.segments.size() == (longSkirt ? 2u : 1u), "segment count");
        require(leg.segments[0].current.stations.size() == (fixture == 1 ? 2u : 3u), "merged stations");
        for (const auto &pair : leg)
        {
            checkSurface(pair.current);
            checkSurface(pair.rest);
            const auto &points = fixture == 3 ? current : rest;
            near(pair.current.length, (points[pair.segmentId + 1] - points[pair.segmentId]).length() * 1.3);
        }
        if (fixture == 3)
        {
            require(leg.segments[0].current.length != leg.segments[0].rest.length, "independent current profile");
            near(leg.segments[0].current.stations[1].z, leg.segments[0].current.length * 0.25);
            near(leg.segments[1].current.stations[1].z, leg.segments[1].current.length * 0.75);
        }
    }
}

static void testBuildBranches()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int fixture = 0; fixture < 6; ++fixture)
    {
        LegVec3 rest[3] = {{1, 2, 3}, {1, -2, 3}, {1, -5, 7}};
        LegVec3 current[3] = {{1, 2, 3}, {1, -2, 3}, {1, -5, 7}};
        if (fixture == 0 || fixture == 5)
            rest[2].x = nan;
        else if (fixture == 1 || fixture == 4)
            current[2].x = nan;
        else if (fixture == 2)
            current[2] = current[1];
        else
            rest[2] = rest[1];
        Leg leg;
        std::vector<LegSegment> restSegments;
        buildFixture(current, rest, fixture < 4, 0.25, leg, restSegments);
        const std::size_t expectedLegs[] = {0, 0, 1, 1, 0, 0};
        const std::size_t expectedRest[] = {0, 2, 2, 1, 1, 0};
        require(leg.segments.size() == expectedLegs[fixture], "current branch");
        require(restSegments.size() == expectedRest[fixture], "rest branch");
        for (const auto &pair : leg)
            checkSurface(pair.current);
    }
}

static void testDegenerateSegments()
{
    const LegMat3 rows;
    const LegVec3 origin;
    const SkirtLegProfile tiny(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0.4e-6, 0.4e-6);
    const auto parameters = legStationParameters(tiny);
    require(parameters.size() == 2 && parameters[0] == 0 && parameters[1] == 1, "tiny stations");
    LegSegment segment;
    require(makeLegSegment(rows, origin, LegVec3(0, 0.4e-6, 0), 0, LegVec3(1, 10, 1), tiny, true, segment),
            "scaled tiny segment");
    near(segment.axis, LegVec3(0, 1, 0));
    near(segment.x, LegVec3(0, 1, 0));
    near(segment.stations.back().z, 8e-6);
    require(LegVec3().normal().length() == 0.0, "zero normal stays zero");
    require(LegVec3(3e-11, 4e-11, 0).normal().x == 3e-11, "tiny normal stays unchanged");
    near(LegVec3(3e-9, 4e-9, 0).normal(), LegVec3(0.6, 0.8, 0));
    const LegVec3 points[3] = {{0, 0, 0}, {0, -4, 0}, {0, -7, 4}};
    for (int fixture = 0; fixture < 5; ++fixture)
    {
        LegSegment invalid;
        LegMat3 joint;
        LegVec3 target = points[1], scale(1, 1, 1);
        if (fixture == 0)
            target = origin;
        else if (fixture == 1)
            scale.x = 0;
        else if (fixture == 2)
            scale.z = -1;
        else if (fixture == 3)
            scale.y = std::numeric_limits<double>::infinity();
        else
            joint[0] = LegVec3();
        require(!makeLegSegment(joint, origin, target, 0, scale, profile(points), false, invalid),
                "invalid segment accepted");
    }
    const SkirtLegProfile chain(1, 1, 2, 2, 3, 3, 4, 4, 0.5, 0.5, 1 - 1.5e-9, 1.5e-9);
    const auto stations = legStationParameters(chain);
    require(stations.size() == 4 && stations[2] == chain.knee && stations[3] == 1.0, "station precedence");
}

static void testRotationBoundaries()
{
    const LegVec3 raw(-0.38023345746295101, 0.89354277136842208, 0.41442833956952163);
    LegMat3 rows;
    rows[0] = raw;
    LegVec3 x, y, z;
    legRingFrame(rows, LegVec3(), 0, raw, x, y, z);
    near(y, raw.normal());
    near(x, (rows[1] - (rows[1] * y) * y).normal());
    near(z, (y ^ x).normal());
    const LegVec3 vector(1, 2, 3);
    near(legRotateVector(vector, raw, LegVec3()), vector);
    for (double epsilon : {0.0, 5e-6, 2e-5})
    {
        const LegVec3 target(-std::cos(epsilon), std::sin(epsilon), 0);
        legRingFrame(LegMat3(), LegVec3(), 0, target, x, y, z);
        if (epsilon < 1e-5)
        {
            near(x, LegVec3(0, 1, 0));
            near(y, LegVec3(-std::cos(epsilon), 0, -std::sin(epsilon)));
        }
        else
        {
            near(x, LegVec3(-std::sin(epsilon), -std::cos(epsilon), 0));
            near(y, target);
        }
        near(z, (y ^ x).normal());
    }
}

static void testAlgebraAndContacts()
{
    LegMat3 matrix;
    matrix[0] = {2, 0.3, 0};
    matrix[1] = {0, 3, 0.4};
    matrix[2] = {0.1, 0, 4};
    const LegVec3 vector(0.2, -0.4, 0.7);
    near(vector * matrix * matrix.inverse(), vector);
    require((vector * matrix.transpose() - vector * matrix.inverse()).length() > 1.0, "true inverse");
    near(legRotateVector(LegVec3(1, 2, 3), LegVec3(1, 0, 0), LegVec3(0, 1, 0)), LegVec3(-2, 1, 3));
    LegSegmentPair pair;
    pair.current.origin = pair.rest.origin = LegVec3();
    pair.current.axis = pair.rest.axis = LegVec3(0, 1, 0);
    pair.current.x = pair.rest.x = LegVec3(1, 0, 0);
    pair.current.z = pair.rest.z = LegVec3(0, 0, 1);
    pair.current.length = pair.rest.length = 4;
    pair.current.stations = pair.rest.stations = {{0, 1, 1}, {4, 1, 1}};
    pair.legId = pair.segmentId = 0;
    pair.current.raySlopes = pair.rest.raySlopes = {{{0.0, 0.0}}};
    Leg leg;
    leg.segments.push_back(pair);
    legPrepareLeg(leg);
    for (double x : {0.5, 1.5})
    {
        LegContact contact;
        require(legConstraint(leg, LegVec3(x, 2, 0), LegVec3(x, 2, 0), LegVec3(1, 0, 0), 0.0, contact),
                "contact missing");
        near(contact.normal, LegVec3(1, 0, 0));
        near(contact.distance, 1.0 - x);
        near(contact.weight, 1.0);
    }
    near(legLocalRadius({pair.rest}, LegVec3(3, 2, 0)), 1.0);
    near(legLocalRadius({}, vector), 0.0);
    LegContact unused;
    require(!legConstraint({}, vector, vector, vector, 0, unused), "empty leg contact");
    pair.restJoint.rows[0] = pair.currentJoint.rows[0] = LegVec3();
    pair.current.axis = LegVec3(0, -1, 0);
    leg.segments[0] = pair;
    legPrepareLeg(leg);
    LegVec3 side;
    require(legTransportedSide(leg.segments[0], LegVec3(1, 2, 0), LegVec3(1, 0, 0), side), "antiparallel fallback");
    near(side, LegVec3(1, 0, 0));
    LegPolynomial polynomial;
    polynomial.degree = 2;
    polynomial.coefficients[0] = 1;
    polynomial.coefficients[1] = -2;
    polynomial.coefficients[2] = 1;
    LegScratch scratch;
    const auto &roots = legPolynomialRoots(polynomial, 0, 2, scratch);
    require(roots.size() == 1, "tangent root");
    near(roots[0], 1.0);
}

static void sameBits(double actual, double expected)
{
    require(std::memcmp(&actual, &expected, sizeof(double)) == 0, "contact bits differ");
}

static void sameBits(const LegVec3 &actual, const LegVec3 &expected)
{
    sameBits(actual.x, expected.x);
    sameBits(actual.y, expected.y);
    sameBits(actual.z, expected.z);
}

static void testPreparedContacts()
{
    const LegVec3 rest[3] = {{1, 2, 3}, {1, -2, 3}, {1, -5, 7}};
    const LegVec3 current[3] = {{1, 2, 3}, {2, -4, 4}, {4, -6, 8}};
    const LegVec3 radial(1, 0, 0);
    LegScratch scratch;
    for (int fixture = 0; fixture < 3; ++fixture)
    {
        Leg leg;
        std::vector<LegSegment> restSegments;
        const bool longSkirt = fixture != 0;
        buildFixture(current, rest, longSkirt, 0.25, leg, restSegments);
        if (fixture == 2)
            leg.segments.erase(leg.segments.begin());
        for (auto &pair : leg.segments)
        {
            pair.restJoint.rows[0] = {2, 0.3, 0};
            pair.restJoint.rows[1] = {0, 3, 0.4};
            pair.restJoint.rows[2] = {0.1, 0, 4};
            pair.currentJoint.rows[0] = {0.8, 0.6, 0};
            pair.currentJoint.rows[1] = {-0.6, 0.8, 0};
        }
        legPrepareLeg(leg);
        double bound = 0.0;
        for (const auto &pair : leg)
        {
            require(pair.rotationsValid, "prepared rotations invalid");
            double radius = 0.0;
            for (const auto &st : pair.current.stations)
                radius = (std::max)(radius, (std::max)(st.radiusX, st.radiusZ));
            bound = (std::max)(bound, std::hypot(pair.current.length, radius));
            for (std::size_t i = 0; i < pair.current.raySlopes.size(); ++i)
            {
                const auto &a = pair.current.stations[i], &b = pair.current.stations[i + 1];
                sameBits(pair.current.raySlopes[i][0], (b.radiusX - a.radiusX) / (b.z - a.z));
                sameBits(pair.current.raySlopes[i][1], (b.radiusZ - a.radiusZ) / (b.z - a.z));
            }
        }
        sameBits(leg.bound, bound);
        sameBits(leg.tmax,
                 2.0 * bound + (leg.segments.size() == 2
                                    ? 2.0 * (leg.segments[1].current.origin - leg.segments[0].current.origin).length()
                                    : 0.0));
        for (const auto &pair : leg)
        {
            const LegVec3 restPoint = pair.rest.origin + 0.4 * pair.rest.length * pair.rest.axis + 0.2 * pair.rest.x;
            LegSampleTransport transport;
            legSampleTransport(leg, restPoint, radial, transport);
            for (double kappa : {0.0, 0.2, 1.0})
            {
                const LegContactGeometry geometry = legPrepareContactGeometry(leg, kappa);
                for (double fraction : {0.2, 2.0})
                    for (double step : {0.0, 0.03})
                    {
                        const auto section = legSectionFrame(pair.current, 0.4 * pair.current.length);
                        const LegVec3 point = section.center + (fraction + step) * section.x;
                        double phiMin = std::numeric_limits<double>::infinity();
                        for (const auto &segment : leg)
                            phiMin = (std::min)(phiMin, legLocal(segment.current, point).phi);
                        require((phiMin <= 0.0) == (fraction < 1.0), "fixture classification");
                        LegContact wrapper, prepared, localScratch, sharedWrapper, fixedGeometry;
                        const bool expected = legConstraint(leg, point, restPoint, radial, kappa, wrapper);
                        require(expected, "fixture contact missing");
                        scratch.intervals.push_back({0.0, 1e10, 0, 0});
                        require(legConstraint(leg, point, transport, kappa, prepared, scratch) == expected,
                                "prepared contact validity");
                        require(legConstraint(leg, point, transport, kappa, localScratch) == expected,
                                "local scratch contact validity");
                        require(legConstraint(leg, point, restPoint, radial, kappa, sharedWrapper, scratch) == expected,
                                "shared wrapper contact validity");
                        require(legConstraint(leg, point, transport, kappa, fixedGeometry, scratch, geometry) ==
                                    expected,
                                "fixed geometry contact validity");
                        sameBits(testCapsuleDistance(point, geometry), testCapsuleDistance(leg, point));
                        for (const auto &actual : {prepared, localScratch, sharedWrapper, fixedGeometry})
                        {
                            sameBits(actual.normal, wrapper.normal);
                            sameBits(actual.distance, wrapper.distance);
                            sameBits(actual.weight, wrapper.weight);
                        }
                    }
            }
        }
    }
    LegContact contact;
    scratch.intervals.push_back({0.0, 1e10, 0, 0});
    require(!legConstraint(Leg(), LegVec3(), LegSampleTransport(), 0.2, contact, scratch), "empty prepared leg");
    require(scratch.intervals.empty() && scratch.raySites.empty(), "scratch not cleared on early return");
    for (const auto &level : scratch.polynomial)
        require(level.roots.empty() && level.sites.empty(), "polynomial scratch not cleared");
}

static void testPreparedTransportFallbacks()
{
    const LegVec3 points[3] = {{0, 0, 0}, {0, 4, 0}, {0, 8, 0}};
    Leg leg;
    std::vector<LegSegment> restSegments;
    buildFixture(points, points, false, 0.25, leg, restSegments);
    auto &pair = leg.segments[0];
    pair.currentJoint.rows[0] = {0, 1, 0};
    pair.currentJoint.rows[1] = {1, 0, 0};
    legPrepareLeg(leg);
    require(pair.rotationsValid, "fallback fixture rotations");
    LegSampleTransport transport;
    legSampleTransport(leg, {1, 2, 0}, {1, 0, 0}, transport);
    require(transport.valid[0], "projected side fallback missing");
    near(transport.sides[0], LegVec3(1, 0, 0));
    legSampleTransport(leg, {1, 2, 0}, {0, 1, 0}, transport);
    require(transport.valid[0], "rest point outward fallback missing");
    near(transport.sides[0], LegVec3(1, 0, 0));
    legSampleTransport(leg, {0, 2, 0}, {0, 1, 0}, transport);
    require(!transport.valid[0], "zero outward accepted");
    pair.restJoint.rows[0] = LegVec3();
    pair.current.axis = {0, -1, 0};
    legPrepareLeg(leg);
    require(!pair.rotationsValid && pair.antiparallel, "antiparallel preparation");
    legSampleTransport(leg, {1, 2, 0}, {1, 0, 0}, transport);
    require(transport.valid[0], "antiparallel transport missing");
    near(transport.sides[0], LegVec3(1, 0, 0));
}

static void testPolynomialScratch()
{
    LegScratch scratch;
    for (int repeat = 0; repeat < 3; ++repeat)
    {
        LegPolynomial quartic;
        const double coefficients[] = {24, -50, 35, -10, 1};
        std::copy(coefficients, coefficients + 5, quartic.coefficients);
        const auto &roots = legPolynomialRoots(quartic, 0, 5, scratch);
        require(roots.size() == 4, "quartic roots lost");
        for (std::size_t i = 0; i < roots.size(); ++i)
            near(roots[i], static_cast<double>(i + 1));
        const double tangent[] = {4, -12, 13, -6, 1};
        std::copy(tangent, tangent + 5, quartic.coefficients);
        const auto &touches = legPolynomialRoots(quartic, 0, 3, scratch);
        require(touches.size() == 2, "stationary roots lost");
        near(touches[0], 1.0);
        near(touches[1], 2.0);
        LegPolynomial linear;
        linear.coefficients[0] = -2;
        linear.coefficients[1] = 1;
        const auto &single = legPolynomialRoots(linear, 0, 3, scratch);
        require(single.size() == 1, "degree reduction roots");
        near(single[0], 2.0);
        require(legPolynomialRoots(LegPolynomial(), 0, 3, scratch).empty(), "constant retained stale roots");
    }
}

static Leg straightLeg()
{
    LegSegmentPair pair;
    pair.current.origin = pair.rest.origin = LegVec3();
    pair.current.axis = pair.rest.axis = LegVec3(0, 1, 0);
    pair.current.x = pair.rest.x = LegVec3(1, 0, 0);
    pair.current.z = pair.rest.z = LegVec3(0, 0, 1);
    pair.current.length = pair.rest.length = 4;
    pair.current.stations = pair.rest.stations = {{0, 1, 1}, {4, 1, 1}};
    pair.legId = pair.segmentId = 0;
    pair.current.raySlopes = pair.rest.raySlopes = {{{0.0, 0.0}}};
    Leg leg;
    leg.segments.push_back(pair);
    legPrepareLeg(leg);
    return leg;
}

static void testCapsuleGaps()
{
    const Leg leg = straightLeg();
    near(testCapsuleDistance(leg, LegVec3(1.5, 2, 0)), 0.5);
    near(testCapsuleDistance(leg, LegVec3(0.5, 2, 0)), 0.0);
    near(testCapsuleDistance(leg, LegVec3(0, 6, 0)), 1.0);
    near(testCapsuleDistance(leg, LegVec3(0, -2, 0)), 1.0);
    require(legLocal(leg.segments[0].current, LegVec3(0.5, 2, 0)).phi <= 0.0, "inside phi");
    require(!(legLocal(leg.segments[0].current, LegVec3(1.5, 2, 0)).phi <= 0.0), "outside phi");
    Leg tapered = straightLeg();
    tapered.segments[0].current.stations = {{0, 2, 1}, {4, 1, 2}};
    near(testCapsuleDistance(tapered, LegVec3(2.5, 2, 0)), 0.5);
    require(!std::isfinite(testCapsuleDistance(leg, LegVec3(std::numeric_limits<double>::quiet_NaN(), 0, 0))),
            "non-finite gap");
    require(testCapsuleDistance(Leg(), LegVec3(5, 5, 5)) == 0.0, "empty leg gap");
    Leg offset = straightLeg();
    offset.segments[0].current.stations = {{3, 1, 1}, {4, 1, 1}};
    near(testCapsuleDistance(offset, LegVec3(2, 0, 0)), 1.0);
    near(testCapsuleDistance(offset, LegVec3(2, 2, 0)), 1.0);
    near(testCapsuleDistance(offset, LegVec3(0, -2, 0)), 1.0);
    {
        const LegLocal l = legLocal(offset.segments[0].current, LegVec3(2, 0, 0));
        const LegVec3 anchor = legAnchor(offset.segments[0].current, l);
        require((anchor - LegVec3(2, 0, 0)).length() >= testCapsuleDistance(offset, LegVec3(2, 0, 0)) - 1e-12,
                "anchor distance covers the clamped end piece");
    }
    Leg huge = straightLeg();
    huge.segments[0].current.length = 1e155;
    huge.segments[0].current.stations = {{0, 1, 1}, {1e155, 1, 1}};
    require(!std::isfinite(testCapsuleDistance(huge, LegVec3(2, 1e153, 0))), "overflowing span is not finite");
    const LegContactGeometry hugeGeometry = legPrepareContactGeometry(huge, 0.2);
    require(!std::isfinite(testCapsuleDistance(LegVec3(2, 1e153, 0), hugeGeometry)),
            "prepared overflowing span is not finite");
    LegScratch scratch;
    for (double x : {1.5, 3.0, 6.0})
    {
        LegContact contact;
        require(legConstraint(leg, LegVec3(x, 2, 0), LegVec3(x, 2, 0), LegVec3(1, 0, 0), 0.0, contact, scratch),
                "outside contact missing");
        const double gap = testCapsuleDistance(leg, LegVec3(x, 2, 0));
        require(contact.distance <= -gap + 1e-9, "capsule must contain frustum");
    }
}

int main()
{
    testDrawingFixtures();
    testBuildBranches();
    testDegenerateSegments();
    testRotationBoundaries();
    testAlgebraAndContacts();
    testPreparedContacts();
    testPreparedTransportFallbacks();
    testPolynomialScratch();
    testCapsuleGaps();
    std::cout << "skirt leg model: 9 test groups passed\n";
}
