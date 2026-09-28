#pragma once

#include "skirtLegProfile.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#if defined(_M_X64) || defined(__SSE2__)
#include <emmintrin.h>
#endif

struct LegVec3
{
    double x, y, z;

    LegVec3() : x(0.0), y(0.0), z(0.0)
    {
    }
    LegVec3(double xValue, double yValue, double zValue) : x(xValue), y(yValue), z(zValue)
    {
    }
    double &operator[](unsigned int i)
    {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    double operator[](unsigned int i) const
    {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    LegVec3 operator+(const LegVec3 &v) const
    {
        return {x + v.x, y + v.y, z + v.z};
    }
    LegVec3 operator-(const LegVec3 &v) const
    {
        return {x - v.x, y - v.y, z - v.z};
    }
    LegVec3 operator-() const
    {
        return {-x, -y, -z};
    }
    LegVec3 operator*(double value) const
    {
        return {x * value, y * value, z * value};
    }
    LegVec3 operator/(double value) const
    {
        return {x / value, y / value, z / value};
    }
    double operator*(const LegVec3 &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    LegVec3 operator^(const LegVec3 &v) const
    {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }
    LegVec3 &operator+=(const LegVec3 &v)
    {
        return *this = *this + v;
    }
    LegVec3 &operator*=(double value)
    {
        return *this = *this * value;
    }
    double length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }
    LegVec3 normal() const
    {
        // Maya's MVector::normal leaves vectors shorter than its 1e-10 tolerance unchanged.
        const double l = length();
        return l < 1e-10 ? *this : *this / l;
    }
    void normalize()
    {
        *this = normal();
    }
};

inline LegVec3 operator*(double value, const LegVec3 &v)
{
    return v * value;
}

inline bool finiteVector(const LegVec3 &value)
{
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

struct LegMat3
{
    LegVec3 rows[3];

    LegMat3() : rows{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}
    {
    }
    LegVec3 &operator[](unsigned int row)
    {
        return rows[row];
    }
    const LegVec3 &operator[](unsigned int row) const
    {
        return rows[row];
    }
    double det3x3() const
    {
        return rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1]) -
               rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0]) +
               rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]);
    }
    LegMat3 transpose() const
    {
        LegMat3 result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result[i][j] = rows[j][i];
        return result;
    }
    LegMat3 inverse() const
    {
        const double determinant = det3x3();
        LegMat3 cofactors;
        cofactors[0] = (rows[1] ^ rows[2]) / determinant;
        cofactors[1] = (rows[2] ^ rows[0]) / determinant;
        cofactors[2] = (rows[0] ^ rows[1]) / determinant;
        return cofactors.transpose();
    }
};

inline LegVec3 operator*(const LegVec3 &v, const LegMat3 &matrix)
{
    return {v.x * matrix[0][0] + v.y * matrix[1][0] + v.z * matrix[2][0],
            v.x * matrix[0][1] + v.y * matrix[1][1] + v.z * matrix[2][1],
            v.x * matrix[0][2] + v.y * matrix[1][2] + v.z * matrix[2][2]};
}

struct LegJoint
{
    LegMat3 rows;
    LegVec3 origin;
};

inline double clampValue(double value, double low, double high)
{
    return (std::min)(high, (std::max)(low, value));
}

// Maya picks the 180-degree axis from the unit axis that follows the largest component by two
// (ties resolve to the later index): k=0 uses z x f, k=1 uses x x f, k=2 uses y x f.
inline LegVec3 legAntiparallelAxis(const LegVec3 &from)
{
    int largest = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(from[i]) >= std::abs(from[largest]))
            largest = i;
    LegVec3 basis;
    basis[(largest + 2) % 3] = 1.0;
    return basis ^ from;
}

inline LegVec3 legRotateVector(const LegVec3 &v, const LegVec3 &fromAxis, const LegVec3 &toAxis)
{
    // Maya's MQuaternion(from, to) is the identity when either vector is shorter than 1e-5 or the
    // unnormalised cross product is shorter than 1e-5 (measured), except for antiparallel inputs.
    if (fromAxis.length() < 1e-5 || toAxis.length() < 1e-5)
        return v;
    const LegVec3 from = fromAxis.normal(), to = toAxis.normal();
    const LegVec3 cross = from ^ to;
    const double cosine = from * to;
    if (1.0 + cosine >= 2.5e-11 && (fromAxis ^ toAxis).length() < 1e-5)
        return v;
    if (1.0 + cosine < 2.5e-11)
    {
        // Exact half turn; Maya's rounded half-angle can fall short of it by ~1e-8 rad.
        const LegVec3 axis = legAntiparallelAxis(from).normal();
        return axis * (2.0 * (axis * v)) - v;
    }
    return v + (cross ^ v) + (cross ^ (cross ^ v)) / (1.0 + cosine);
}

inline void legRingFrame(const LegMat3 &jointRows, const LegVec3 &jointOrigin, short axisIndex, const LegVec3 &target,
                         LegVec3 &x, LegVec3 &y, LegVec3 &z)
{
    const int yIndex = axisIndex % 3, xIndex = (axisIndex + 1) % 3;
    const LegVec3 rawY = yIndex < 0 ? LegVec3() : jointRows[yIndex];
    y = (axisIndex >= 3 ? -rawY : rawY).normal();
    LegVec3 uX = (xIndex < 0 ? LegVec3() : jointRows[xIndex]).normal();
    const LegVec3 v = target - jointOrigin;
    // Maya's MVector::rotateTo returns the identity when the target is shorter than 1e-5 or the
    // unnormalised cross product with the unit axis is shorter than 1e-5 (measured), except for
    // antiparallel inputs.
    if (v.length() >= 1e-5)
    {
        const LegVec3 from = y.normal(), to = v.normal();
        LegVec3 rotationAxis = from ^ to;
        // atan2 keeps parallel axes finite when their rounded dot product is below one.
        double angle = std::atan2(rotationAxis.length(), clampValue(from * to, -1.0, 1.0));
        if (1.0 + from * to < 2.5e-11)
            rotationAxis = legAntiparallelAxis(from);
        else if ((from ^ v).length() < 1e-5)
            angle = 0.0;
        LegMat3 rotation;
        if (angle != 0.0)
        {
            const LegVec3 q = rotationAxis.normal() * std::sin(angle / 2.0);
            const double w = std::cos(angle / 2.0);
            rotation[0] = {1.0 - 2.0 * (q.y * q.y + q.z * q.z), 2.0 * (q.x * q.y + q.z * w),
                           2.0 * (q.x * q.z - q.y * w)};
            rotation[1] = {2.0 * (q.x * q.y - q.z * w), 1.0 - 2.0 * (q.x * q.x + q.z * q.z),
                           2.0 * (q.y * q.z + q.x * w)};
            rotation[2] = {2.0 * (q.x * q.z + q.y * w), 2.0 * (q.y * q.z - q.x * w),
                           1.0 - 2.0 * (q.x * q.x + q.y * q.y)};
        }
        y = y * rotation;
        uX = uX * rotation;
    }
    x = (uX - (uX * y) * y).normal();
    z = (y ^ x).normal();
}

struct LegStation
{
    double z, radiusX, radiusZ;
};

struct LegSegment
{
    LegVec3 origin;
    LegVec3 axis, x, z;
    double length;
    std::vector<LegStation> stations;
    std::vector<std::array<double, 2>> raySlopes;
};

struct LegSegmentPair
{
    LegSegment current, rest;
    LegJoint currentJoint, restJoint;
    int legId, segmentId;
    bool rotationsValid = false, antiparallel = false;
    LegMat3 restInverse, currentRotation;
};

struct Leg
{
    std::vector<LegSegmentPair> segments;
    double bound = 0.0, tmax = 0.0;

    std::vector<LegSegmentPair>::const_iterator begin() const
    {
        return segments.begin();
    }
    std::vector<LegSegmentPair>::const_iterator end() const
    {
        return segments.end();
    }
};

struct LegContactGeometry
{
    struct Capsule
    {
        LegVec3 p0, span;
        double length2, radius;
    };
    struct Segment
    {
        std::vector<Capsule> capsules;
        std::vector<std::array<double, 2>> slopes;
        std::vector<double> halfWidths, widths;
    };
    std::vector<Segment> segments;
    bool capsulesValid = true;
};

inline LegContactGeometry::Segment legPrepareContactSegment(const LegSegment &c, double kappa,
                                                            bool &capsulesValid)
{
    LegContactGeometry::Segment segment;
    if (c.stations.empty())
    {
        capsulesValid = false;
        return segment;
    }
    const auto capsule = [&](double z0, double z1, double radius) {
        if (!std::isfinite(z0) || !std::isfinite(z1) || !std::isfinite(radius))
        {
            capsulesValid = false;
            return;
        }
        const LegVec3 p0 = c.origin + c.axis * z0, p1 = c.origin + c.axis * z1;
        const LegVec3 span = p1 - p0;
        const double length2 = span * span;
        if (!std::isfinite(length2))
        {
            capsulesValid = false;
            return;
        }
        segment.capsules.push_back({p0, span, length2, radius});
    };
    if (!std::isfinite(c.length))
        capsulesValid = false;
    else
    {
        segment.capsules.reserve(c.stations.size() + 1);
        const LegStation &front = c.stations.front(), &back = c.stations.back();
        if (front.z > 0.0)
            capsule(0.0, front.z, (std::max)(front.radiusX, front.radiusZ));
        for (std::size_t i = 0; i + 1 < c.stations.size(); ++i)
        {
            const double radius = (std::max)((std::max)(c.stations[i].radiusX, c.stations[i].radiusZ),
                                             (std::max)(c.stations[i + 1].radiusX, c.stations[i + 1].radiusZ));
            capsule(c.stations[i].z, c.stations[i + 1].z, radius);
        }
        if (back.z < c.length)
            capsule(back.z, c.length, (std::max)(back.radiusX, back.radiusZ));
    }
    segment.slopes.reserve(c.stations.size() - 1);
    segment.halfWidths.reserve(c.stations.size());
    segment.widths.reserve(c.stations.size());
    for (std::size_t i = 0; i + 1 < c.stations.size(); ++i)
    {
        const LegStation &a = c.stations[i], &b = c.stations[i + 1];
        segment.slopes.push_back({{(b.radiusX - a.radiusX) / (b.z - a.z),
                                   (b.radiusZ - a.radiusZ) / (b.z - a.z)}});
    }
    for (std::size_t j = 1; j + 1 < c.stations.size(); ++j)
    {
        const double f = kappa * (std::min)(c.stations[j].radiusX, c.stations[j].radiusZ);
        segment.halfWidths.push_back(f);
        segment.widths.push_back(2.0 * f);
    }
    return segment;
}

inline LegContactGeometry legPrepareContactGeometry(const Leg &leg, double kappa)
{
    LegContactGeometry geometry;
    geometry.segments.reserve(leg.segments.size());
    for (const auto &pair : leg.segments)
        geometry.segments.push_back(legPrepareContactSegment(pair.current, kappa, geometry.capsulesValid));
    return geometry;
}

struct LegContact
{
    LegVec3 normal;
    double distance, weight;
};

inline std::vector<double> legStationParameters(const SkirtLegProfile &profile)
{
    const double candidates[] = {0.0, profile.knee * profile.thighPosition, profile.knee,
                                 profile.knee + (1.0 - profile.knee) * profile.calfPosition, 1.0};
    const int precedence[] = {4, 2, 1, 3, 0};
    std::vector<double> parameters;
    for (int index : precedence)
    {
        bool coincident = false;
        for (double s : parameters)
            coincident = coincident || std::fabs(candidates[index] - s) <= 1e-9;
        if (!coincident)
            parameters.push_back(candidates[index]);
    }
    std::sort(parameters.begin(), parameters.end());
    return parameters;
}

inline bool makeLegSegment(const LegMat3 &jointRows, const LegVec3 &jointOrigin, const LegVec3 &target, short ringAxis,
                           const LegVec3 &scale, const SkirtLegProfile &profile, bool calf, LegSegment &cylinder)
{
    cylinder.stations.clear();
    cylinder.raySlopes.clear();
    cylinder.origin = jointOrigin;
    const LegVec3 end = target;
    const LegVec3 direction = end - cylinder.origin;
    cylinder.length = direction.length() * scale.y;
    if (!std::isfinite(cylinder.length) || cylinder.length < 1e-6)
        return false;
    cylinder.axis = direction.normal();
    LegMat3 frame;
    legRingFrame(jointRows, jointOrigin, ringAxis, end, frame[0], frame[1], frame[2]);
    for (int row = 0; row < 3; ++row)
        if (!finiteVector(frame[row]))
            return false;
    const double determinant = frame.det3x3();
    if (!std::isfinite(determinant) || std::fabs(determinant) < 1e-8)
        return false;
    cylinder.x = frame[0].normal();
    cylinder.z = frame[2].normal();
    const double start = calf ? profile.knee : 0.0;
    const double finish = calf ? 1.0 : profile.knee;
    for (double s : legStationParameters(profile))
    {
        if (s < start || s > finish)
            continue;
        const SkirtLegProfile::Radius radius = profile.evaluate(s, cylinder.length);
        const LegStation station = {(s - start) * profile.legLength * scale.y, scale.x * radius.x, scale.z * radius.z};
        if (!std::isfinite(station.z) || !std::isfinite(station.radiusX) || !std::isfinite(station.radiusZ) ||
            station.radiusX <= 0.0 || station.radiusZ <= 0.0)
            return false;
        cylinder.stations.push_back(station);
    }
    if (cylinder.stations.size() < 2)
        return false;
    for (std::size_t i = 0; i + 1 < cylinder.stations.size(); ++i)
    {
        const LegStation &first = cylinder.stations[i], &second = cylinder.stations[i + 1];
        cylinder.raySlopes.push_back({{(second.radiusX - first.radiusX) / (second.z - first.z),
                                       (second.radiusZ - first.radiusZ) / (second.z - first.z)}});
    }
    return true;
}

inline bool legJointRotation(const LegMat3 &joint, LegMat3 &rotation)
{
    rotation = LegMat3();
    for (short row = 0; row < 3; row++)
    {
        const LegVec3 axis = joint[row];
        const double norm = axis.length();
        if (!finiteVector(axis) || !std::isfinite(norm) || norm < 1e-8)
            return false;
        rotation[row] = axis / norm;
    }
    const double det = rotation.det3x3();
    return std::isfinite(det) && std::fabs(det) >= 0.5;
}

inline void legPrepareLeg(Leg &leg)
{
    leg.bound = 0.0;
    for (auto &pair : leg.segments)
    {
        LegMat3 restRotation;
        pair.rotationsValid = legJointRotation(pair.restJoint.rows, restRotation) &&
                              legJointRotation(pair.currentJoint.rows, pair.currentRotation);
        if (pair.rotationsValid)
            pair.restInverse = restRotation.inverse();
        pair.antiparallel = pair.rest.axis * pair.current.axis < -1.0 + 1e-6;
        double radius = 0.0;
        for (const auto &st : pair.current.stations)
            radius = (std::max)(radius, (std::max)(st.radiusX, st.radiusZ));
        leg.bound = (std::max)(leg.bound, std::hypot(pair.current.length, radius));
    }
    leg.tmax = 2.0 * leg.bound + (leg.segments.size() == 2
                                      ? 2.0 * (leg.segments[1].current.origin - leg.segments[0].current.origin).length()
                                      : 0.0);
}

inline bool legTransportedSide(const LegSegmentPair &pair, const LegVec3 &restPoint, const LegVec3 &waistRadial,
                               LegVec3 &side)
{
    const LegSegment &rest = pair.rest;
    const LegSegment &current = pair.current;
    LegVec3 outward = waistRadial - rest.axis * (waistRadial * rest.axis);
    if (outward.length() < 1e-6)
    {
        const LegVec3 offset = restPoint - rest.origin;
        outward = offset - rest.axis * (offset * rest.axis);
    }
    if (!finiteVector(outward) || outward.length() < 1e-8)
        return false;
    outward.normalize();
    if (pair.rotationsValid)
    {
        const LegVec3 transported = outward * pair.restInverse * pair.currentRotation;
        side = transported - current.axis * (transported * current.axis);
        if (finiteVector(side) && std::isfinite(side.length()) && side.length() >= 1e-6)
        {
            side.normalize();
            return true;
        }
    }
    const LegVec3 transported = pair.antiparallel ? outward : legRotateVector(outward, rest.axis, current.axis);
    side = transported - current.axis * (transported * current.axis);
    if (!finiteVector(side) || !std::isfinite(side.length()) || side.length() < 1e-6)
        return false;
    side.normalize();
    return true;
}

struct LegSampleTransport
{
    LegVec3 sides[2];
    bool valid[2] = {};
};

inline void legSampleTransport(const Leg &leg, const LegVec3 &restPoint, const LegVec3 &waistRadial,
                               LegSampleTransport &transport)
{
    transport = LegSampleTransport();
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
        transport.valid[k] = legTransportedSide(leg.segments[k], restPoint, waistRadial, transport.sides[k]);
}

inline double smoothRamp(double value)
{
    const double t = (std::min)(1.0, (std::max)(0.0, value));
    return t * t * (3.0 - 2.0 * t);
}

inline LegStation legSection(const LegSegment &cylinder, double z)
{
    z = clampValue(z, 0.0, cylinder.length);
    if (z <= cylinder.stations.front().z)
        return cylinder.stations.front();
    for (std::size_t j = 1; j < cylinder.stations.size(); ++j)
    {
        const LegStation &a = cylinder.stations[j - 1];
        const LegStation &b = cylinder.stations[j];
        if (z <= b.z)
        {
            const double t = (z - a.z) / (b.z - a.z);
            return {z, a.radiusX + t * (b.radiusX - a.radiusX), a.radiusZ + t * (b.radiusZ - a.radiusZ)};
        }
    }
    return cylinder.stations.back();
}

struct LegLocal
{
    double z, x, y, A, B, eta, r, rho, rhoMin, phi;
    LegVec3 e;
};

inline LegLocal legLocal(const LegSegment &c, const LegVec3 &point)
{
    const LegVec3 v = point - c.origin;
    LegLocal l;
    l.z = c.axis * v;
    l.x = c.x * v;
    l.y = c.z * v;
    const LegStation radii = legSection(c, l.z);
    l.A = radii.radiusX;
    l.B = radii.radiusZ;
    l.eta = 1.0 - std::hypot(l.x / l.A, l.y / l.B);
    l.r = std::hypot(l.x, l.y);
    l.e = l.r < 1e-8 ? c.x : (l.x * c.x + l.y * c.z) / l.r;
    const double ex = l.e * c.x, ez = l.e * c.z;
    l.rho = 1.0 / std::sqrt(ex * ex / (l.A * l.A) + ez * ez / (l.B * l.B));
    l.rhoMin = (std::min)(l.A, l.B);
    l.phi = (std::max)(-l.eta, (std::max)(-l.z / l.rhoMin, (l.z - c.length) / l.rhoMin));
    return l;
}

inline double legLocalPhi(const LegSegment &c, const LegVec3 &point)
{
    const LegVec3 v = point - c.origin;
    const double z = c.axis * v;
    const double x = c.x * v;
    const double y = c.z * v;
    const LegStation radii = legSection(c, z);
    const double eta = 1.0 - std::hypot(x / radii.radiusX, y / radii.radiusZ);
    const double rhoMin = (std::min)(radii.radiusX, radii.radiusZ);
    return (std::max)(-eta, (std::max)(-z / rhoMin, (z - c.length) / rhoMin));
}

struct LegCapsuleColumns
{
    std::vector<double> p0x, p0y, p0z, spanx, spany, spanz, length2, radius;
    bool valid = true;

    void prepare(const LegContactGeometry& geometry)
    {
        valid = geometry.capsulesValid;
        std::size_t count = 0;
        for (const auto& segment : geometry.segments)
            count += segment.capsules.size();
        p0x.resize(count); p0y.resize(count); p0z.resize(count);
        spanx.resize(count); spany.resize(count); spanz.resize(count);
        length2.resize(count); radius.resize(count);
        std::size_t i = 0;
        for (const auto& segment : geometry.segments)
            for (const auto& capsule : segment.capsules)
            {
                p0x[i] = capsule.p0.x; p0y[i] = capsule.p0.y; p0z[i] = capsule.p0.z;
                spanx[i] = capsule.span.x; spany[i] = capsule.span.y; spanz[i] = capsule.span.z;
                length2[i] = capsule.length2; radius[i] = capsule.radius;
                ++i;
            }
    }
};

inline double legCapsuleDistance(const LegVec3& point, const LegCapsuleColumns& columns)
{
    if (!finiteVector(point) || !columns.valid)
        return std::numeric_limits<double>::quiet_NaN();
    double minimum = std::numeric_limits<double>::infinity();
    std::size_t i = 0;
#if defined(_M_X64) || defined(__SSE2__)
    const __m128d px = _mm_set1_pd(point.x), py = _mm_set1_pd(point.y), pz = _mm_set1_pd(point.z);
    for (; i + 1 < columns.radius.size(); i += 2)
    {
        const __m128d dx = _mm_sub_pd(px, _mm_loadu_pd(columns.p0x.data() + i));
        const __m128d dy = _mm_sub_pd(py, _mm_loadu_pd(columns.p0y.data() + i));
        const __m128d dz = _mm_sub_pd(pz, _mm_loadu_pd(columns.p0z.data() + i));
        const __m128d projection = _mm_add_pd(
            _mm_add_pd(_mm_mul_pd(dx, _mm_loadu_pd(columns.spanx.data() + i)),
                       _mm_mul_pd(dy, _mm_loadu_pd(columns.spany.data() + i))),
            _mm_mul_pd(dz, _mm_loadu_pd(columns.spanz.data() + i)));
        double projected[2], t[2];
        _mm_storeu_pd(projected, projection);
        for (int lane = 0; lane < 2; ++lane)
        {
            const std::size_t index = i + static_cast<std::size_t>(lane);
            t[lane] = columns.length2[index] > 0.0
                ? clampValue(projected[lane] / columns.length2[index], 0.0, 1.0) : 0.0;
            if (!std::isfinite(t[lane]))
                return std::numeric_limits<double>::quiet_NaN();
        }
        const __m128d tv = _mm_loadu_pd(t);
        const __m128d cx = _mm_add_pd(_mm_loadu_pd(columns.p0x.data() + i),
                                       _mm_mul_pd(_mm_loadu_pd(columns.spanx.data() + i), tv));
        const __m128d cy = _mm_add_pd(_mm_loadu_pd(columns.p0y.data() + i),
                                       _mm_mul_pd(_mm_loadu_pd(columns.spany.data() + i), tv));
        const __m128d cz = _mm_add_pd(_mm_loadu_pd(columns.p0z.data() + i),
                                       _mm_mul_pd(_mm_loadu_pd(columns.spanz.data() + i), tv));
        const __m128d rx = _mm_sub_pd(px, cx), ry = _mm_sub_pd(py, cy), rz = _mm_sub_pd(pz, cz);
        const __m128d squared = _mm_add_pd(_mm_add_pd(_mm_mul_pd(rx, rx), _mm_mul_pd(ry, ry)),
                                            _mm_mul_pd(rz, rz));
        const __m128d distance = _mm_sub_pd(_mm_sqrt_pd(squared), _mm_loadu_pd(columns.radius.data() + i));
        double values[2];
        _mm_storeu_pd(values, distance);
        for (int lane = 0; lane < 2; ++lane)
        {
            const double value = values[lane];
            if (!std::isfinite(value))
                return std::numeric_limits<double>::quiet_NaN();
            minimum = (std::min)(minimum, value > 0.0 ? value : 0.0);
        }
    }
#endif
    for (; i < columns.radius.size(); ++i)
    {
        const LegVec3 p0(columns.p0x[i], columns.p0y[i], columns.p0z[i]);
        const LegVec3 span(columns.spanx[i], columns.spany[i], columns.spanz[i]);
        double t = 0.0;
        if (columns.length2[i] > 0.0)
            t = clampValue(((point - p0) * span) / columns.length2[i], 0.0, 1.0);
        if (!std::isfinite(t))
            return std::numeric_limits<double>::quiet_NaN();
        const double distance = (point - (p0 + span * t)).length() - columns.radius[i];
        if (!std::isfinite(distance))
            return std::numeric_limits<double>::quiet_NaN();
        minimum = (std::min)(minimum, distance > 0.0 ? distance : 0.0);
    }
    return std::isfinite(minimum) ? minimum : 0.0;
}

inline LegVec3 legSideNormal(const LegSegment &c, const LegLocal &l,
                            const LegContactGeometry::Segment &geometry)
{
    if (l.r < 1e-8)
        return l.e;
    const auto &st = c.stations;
    double previousA = geometry.slopes[0][0];
    double previousB = geometry.slopes[0][1];
    double da = previousA, db = previousB;
    const double z = clampValue(l.z, 0.0, c.length);
    for (std::size_t j = 1; j + 1 < st.size(); ++j)
    {
        const double nextA = geometry.slopes[j][0];
        const double nextB = geometry.slopes[j][1];
        const double width = geometry.widths[j - 1];
        const double w = geometry.halfWidths[j - 1] > 0.0 ? smoothRamp(0.5 + (z - st[j].z) / width)
                                                         : (z >= st[j].z ? 1.0 : 0.0);
        da += (nextA - previousA) * w;
        db += (nextB - previousB) * w;
        previousA = nextA;
        previousB = nextB;
    }
    const LegVec3 g = (2.0 * l.x / (l.A * l.A)) * c.x + (2.0 * l.y / (l.B * l.B)) * c.z +
                      (-2.0 * l.x * l.x * da / (l.A * l.A * l.A) - 2.0 * l.y * l.y * db / (l.B * l.B * l.B)) * c.axis;
    return g / g.length();
}

inline LegVec3 legSegmentNormal(const LegSegment &c, const LegLocal &l, double kappa, double exp0, double exp1,
                               const LegContactGeometry::Segment &geometry)
{
    const double f = kappa * l.rhoMin, side = l.r - l.rho;
    const auto corner = [&](double cap) {
        return kappa > 0.0 ? smoothRamp(0.5 + (cap - side) / (2.0 * f)) : (cap > side ? 1.0 : 0.0);
    };
    const double h = smoothRamp(l.z / c.length);
    const double a0 = (1.0 - h) * exp0 * corner(-l.z);
    const double a1 = h * exp1 * corner(l.z - c.length);
    const double t = kappa > 0.0 ? smoothRamp(l.r / l.rhoMin) : (l.r > 0.0 ? 1.0 : 0.0);
    return ((1.0 - a0 - a1) * t) * legSideNormal(c, l, geometry) - a0 * c.axis + a1 * c.axis;
}

inline LegVec3 legAnchor(const LegSegment &c, const LegLocal &l)
{
    double x = l.x, y = l.y;
    const double radius = std::hypot(x / l.A, y / l.B);
    if (radius > 1.0)
    {
        x /= radius;
        y /= radius;
    }
    return c.origin + clampValue(l.z, 0.0, c.length) * c.axis + (x * c.x + y * c.z);
}

inline void legExposures(const LegLocal (&values)[2], std::size_t count, double kappa, double shift, double exp0[2],
                         double exp1[2])
{
    exp0[0] = exp0[1] = exp1[0] = exp1[1] = 1.0;
    if (count == 2)
    {
        const double p0 = values[0].phi - shift, p1 = values[1].phi - shift;
        exp1[0] = kappa > 0.0 ? smoothRamp(p1 / kappa) : (p1 > 0.0 ? 1.0 : 0.0);
        exp0[1] = kappa > 0.0 ? smoothRamp(p0 / kappa) : (p0 > 0.0 ? 1.0 : 0.0);
    }
}

inline bool legNormal(const Leg &leg, const LegVec3 &point, const LegLocal (&values)[2], double kappa, LegVec3 &normal,
                      double &distance, double &weight, const LegContactGeometry &geometry)
{
    double phiMin = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
        phiMin = (std::min)(phiMin, values[k].phi);
    double exp0[2], exp1[2];
    legExposures(values, leg.segments.size(), kappa, 0.0, exp0, exp1);
    LegVec3 sum(0.0, 0.0, 0.0);
    double dmin = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
    {
        const LegLocal &l = values[k];
        const double rho = kappa > 0.0 ? smoothRamp(1.0 - (l.phi - phiMin) / kappa) : (l.phi == phiMin ? 1.0 : 0.0);
        sum += rho * legSegmentNormal(leg.segments[k].current, l, kappa, exp0[k], exp1[k], geometry.segments[k]);
        dmin = (std::min)(dmin, (legAnchor(leg.segments[k].current, l) - point).length());
    }
    const double magnitude = sum.length();
    distance = -dmin;
    weight = smoothRamp(magnitude / 0.25);
    if (magnitude <= 1e-15)
        return false;
    normal = sum / magnitude;
    return true;
}

inline bool legNormal(const Leg &leg, const LegVec3 &point, double kappa, LegVec3 &normal, double &distance,
                      double &weight, const LegContactGeometry &geometry)
{
    LegLocal values[2];
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
        values[k] = legLocal(leg.segments[k].current, point);
    return legNormal(leg, point, values, kappa, normal, distance, weight, geometry);
}

inline LegVec3 legDirection(const Leg &leg, const LegLocal (&values)[2], const LegSampleTransport &transport,
                            double kappa, double shift)
{
    double exp0[2], exp1[2], cont[2];
    LegVec3 dirs[2];
    legExposures(values, leg.segments.size(), kappa, shift, exp0, exp1);
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
    {
        const LegLocal &l = values[k];
        const LegSegment &c = leg.segments[k].current;
        LegVec3 b = transport.sides[k];
        const bool hasSide = transport.valid[k];
        const double alpha = hasSide ? smoothRamp(l.eta) : 0.0;
        if (!hasSide)
            b = LegVec3(0.0, 0.0, 0.0);
        const LegVec3 side = (1.0 - alpha) * l.e + alpha * b;
        const double gcap = l.rhoMin * (std::max)(l.eta, 0.0);
        const double zp = (std::max)(l.z, 0.0), hp = (std::max)(c.length - l.z, 0.0);
        const double g0 = gcap * exp0[k], g1 = gcap * exp1[k];
        const double t0 = kappa > 0.0 && g0 > 0.0 ? smoothRamp(zp / (zp + kappa * g0)) : 1.0;
        const double t1 = kappa > 0.0 && g1 > 0.0 ? smoothRamp(hp / (hp + kappa * g1)) : 1.0;
        dirs[k] = (t0 * t1) * side + ((1.0 - t1) - (1.0 - t0)) * c.axis;
        const double phi = l.phi - shift;
        cont[k] = kappa > 0.0 ? 1.0 - smoothRamp(phi / kappa) : (phi <= 0.0 ? 1.0 : 0.0);
    }
    if (leg.segments.size() == 1)
        return dirs[0];
    const LegStation knee = legSection(leg.segments[0].current, leg.segments[0].current.length);
    const double beta = kappa * (std::min)(knee.radiusX, knee.radiusZ);
    const double ht = leg.segments[0].current.length - values[0].z, hc = values[1].z;
    const double wc = beta > 0.0 ? smoothRamp(0.5 + (hc - ht) / (2.0 * beta)) : (hc > ht ? 1.0 : 0.0);
    return (cont[0] * (1.0 - wc)) * dirs[0] + (cont[1] * wc) * dirs[1];
}

struct LegInterval
{
    double start, end;
    int segmentId, stationId;
};

struct LegScratch
{
    struct PolynomialLevel
    {
        std::vector<double> roots, sites;
    };

    // A quartic has four derivative calls, including the constant leaf. Parents
    // retain their child's roots while constructing and evaluating their own sites.
    PolynomialLevel polynomial[5];
    std::vector<double> raySites;
    std::vector<LegInterval> intervals;

    void clear()
    {
        for (auto &level : polynomial)
        {
            level.roots.clear();
            level.sites.clear();
        }
        raySites.clear();
        intervals.clear();
    }
};

struct LegPolynomial
{
    double coefficients[5] = {};
    int degree = 4;
};

template <int Index> struct LegPolynomialHornerStep
{
    static void ordinary(const double *coefficients, double t, double &value, double &magnitude)
    {
        value = value * t + coefficients[Index];
        magnitude = magnitude * std::abs(t) + std::abs(coefficients[Index]);
        LegPolynomialHornerStep<Index - 1>::ordinary(coefficients, t, value, magnitude);
    }

    static void compensated(const double *coefficients, double t, double &value, double &correction)
    {
        const double coefficient = coefficients[Index];
        const double product = value * t;
        const double splitValue = 134217729.0 * value, splitT = 134217729.0 * t;
        const double highValue = splitValue - (splitValue - value), highT = splitT - (splitT - t);
        const double lowValue = value - highValue, lowT = t - highT;
        const double productError =
            ((highValue * highT - product) + highValue * lowT + lowValue * highT) + lowValue * lowT;
        const double total = product + coefficient;
        const double tail = total - product;
        const double sumError = (product - (total - tail)) + (coefficient - tail);
        correction = correction * t + (productError + sumError);
        value = total;
        LegPolynomialHornerStep<Index - 1>::compensated(coefficients, t, value, correction);
    }
};

template <> struct LegPolynomialHornerStep<-1>
{
    static void ordinary(const double *, double, double &, double &) {}
    static void compensated(const double *, double, double &, double &) {}
};

template <int Degree> inline double legPolynomialValueFixed(const LegPolynomial &polynomial, double t)
{
    double value = polynomial.coefficients[Degree];
    double magnitude = std::abs(value);
    LegPolynomialHornerStep<Degree - 1>::ordinary(polynomial.coefficients, t, value, magnitude);
    if (std::abs(value) > 8 * std::numeric_limits<double>::epsilon() * magnitude)
        return value;
    // Compensated Horner preserves signs near cancellation and multiple roots.
    value = polynomial.coefficients[Degree];
    double correction = 0.0;
    LegPolynomialHornerStep<Degree - 1>::compensated(polynomial.coefficients, t, value, correction);
    return value + correction;
}

inline double legPolynomialValue(const LegPolynomial &polynomial, double t)
{
    switch (polynomial.degree)
    {
    case 0: return legPolynomialValueFixed<0>(polynomial, t);
    case 1: return legPolynomialValueFixed<1>(polynomial, t);
    case 2: return legPolynomialValueFixed<2>(polynomial, t);
    case 3: return legPolynomialValueFixed<3>(polynomial, t);
    default: return legPolynomialValueFixed<4>(polynomial, t);
    }
}

inline const std::vector<double> &legPolynomialRoots(LegPolynomial polynomial, double lo, double hi,
                                                     LegScratch &scratch, std::size_t depth = 0)
{
    auto &roots = scratch.polynomial[depth].roots;
    auto &sites = scratch.polynomial[depth].sites;
    roots.clear();
    sites.clear();
    while (polynomial.degree > 0 && polynomial.coefficients[polynomial.degree] == 0.0)
        --polynomial.degree;
    if (polynomial.degree == 0)
        return roots;
    LegPolynomial derivative;
    derivative.degree = polynomial.degree - 1;
    for (int i = 1; i <= polynomial.degree; ++i)
        derivative.coefficients[i - 1] = i * polynomial.coefficients[i];
    const auto &critical = legPolynomialRoots(derivative, lo, hi, scratch, depth + 1);
    sites.assign(critical.begin(), critical.end());
    sites.push_back(lo);
    sites.push_back(hi);
    std::sort(sites.begin(), sites.end());
    sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
    double values[6];
    for (std::size_t i = 0; i < sites.size(); ++i)
    {
        values[i] = legPolynomialValue(polynomial, sites[i]);
        if (std::binary_search(critical.begin(), critical.end(), sites[i]))
        {
            double magnitude = std::abs(polynomial.coefficients[polynomial.degree]);
            for (int j = polynomial.degree - 1; j >= 0; --j)
                magnitude = magnitude * std::abs(sites[i]) + std::abs(polynomial.coefficients[j]);
            // A stationary root need not be representable; account for the squared
            // roundoff of its adjacent-double location and compensated evaluation.
            const double epsilon = std::numeric_limits<double>::epsilon();
            if (std::abs(values[i]) <= 32 * (epsilon * epsilon) * magnitude)
                values[i] = 0.0;
        }
        if (values[i] == 0.0)
            roots.push_back(sites[i]);
    }
    for (std::size_t i = 1; i < sites.size(); ++i)
    {
        double a = sites[i - 1], b = sites[i], fa = values[i - 1], fb = values[i];
        if (fa == 0.0 || fb == 0.0 || (fa < 0.0) == (fb < 0.0))
            continue;
        while (true)
        {
            const double mid = a + (b - a) * 0.5;
            if (mid == a || mid == b)
            {
                roots.push_back(b);
                break;
            }
            const double fm = legPolynomialValue(polynomial, mid);
            if (fm == 0.0)
            {
                roots.push_back(mid);
                break;
            }
            if ((fm < 0.0) == (fa < 0.0))
            {
                a = mid;
                fa = fm;
            }
            else
                b = mid;
        }
    }
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    return roots;
}

inline void legRayIntervals(const LegSegmentPair &segment, const LegVec3 &p, const LegVec3 &n, double tmax,
                            std::vector<LegInterval> &intervals, LegScratch &scratch)
{
    const LegSegment &c = segment.current;
    const LegVec3 offset = p - c.origin;
    const double z = offset * c.axis, dz = n * c.axis;
    const double x = offset * c.x, dx = n * c.x, y = offset * c.z, dy = n * c.z;
    for (std::size_t station = 0; station + 1 < c.stations.size(); ++station)
    {
        const LegStation &first = c.stations[station];
        const LegStation &second = c.stations[station + 1];
        const double lower = (std::max)(0.0, first.z), upper = (std::min)(c.length, second.z);
        double lo = 0.0, hi = tmax;
        if (lower > upper)
            continue;
        if (dz == 0.0)
        {
            if (z < lower || z > upper)
                continue;
        }
        else
        {
            const double a = (lower - z) / dz, b = (upper - z) / dz;
            lo = (std::max)(lo, (std::min)(a, b));
            hi = (std::min)(hi, (std::max)(a, b));
            if (lo > hi)
                continue;
        }
        const double slopeA = c.raySlopes[station][0];
        const double slopeB = c.raySlopes[station][1];
        const double a = first.radiusX + slopeA * (z - first.z), da = slopeA * dz;
        const double b = first.radiusZ + slopeB * (z - first.z), db = slopeB * dz;
        const double xx[3] = {x * x, 2.0 * x * dx, dx * dx};
        const double yy[3] = {y * y, 2.0 * y * dy, dy * dy};
        const double aa[3] = {a * a, 2.0 * a * da, da * da};
        const double bb[3] = {b * b, 2.0 * b * db, db * db};
        LegPolynomial polynomial;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                polynomial.coefficients[i + j] += xx[i] * bb[j] + yy[i] * aa[j] - aa[i] * bb[j];
        const auto &roots = legPolynomialRoots(polynomial, lo, hi, scratch);
        auto &sites = scratch.raySites;
        sites.assign(roots.begin(), roots.end());
        sites.push_back(lo);
        sites.push_back(hi);
        std::sort(sites.begin(), sites.end());
        sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
        for (double t : sites)
            if (std::binary_search(roots.begin(), roots.end(), t) || legPolynomialValue(polynomial, t) <= 0.0)
                intervals.push_back({t, t, segment.segmentId, static_cast<int>(station)});
        for (std::size_t i = 1; i < sites.size(); ++i)
        {
            const double a = sites[i - 1], b = sites[i];
            if (legPolynomialValue(polynomial, a + (b - a) * 0.5) <= 0.0)
                intervals.push_back({a, b, segment.segmentId, static_cast<int>(station)});
        }
    }
}

inline bool legUnionExit(const Leg &leg, const LegVec3 &p, const LegVec3 &n, double &distance, LegScratch &scratch)
{
    const double tmax = leg.tmax;
    auto &intervals = scratch.intervals;
    intervals.clear();
    for (const auto &segment : leg)
        legRayIntervals(segment, p, n, tmax, intervals, scratch);
    std::sort(intervals.begin(), intervals.end(), [](const LegInterval &a, const LegInterval &b) {
        return a.start < b.start || (a.start == b.start && (a.segmentId < b.segmentId ||
                                                            (a.segmentId == b.segmentId && a.stationId < b.stationId)));
    });
    distance = 0.0;
    bool found = std::any_of(intervals.begin(), intervals.end(), [=](const LegInterval &interval) {
        return interval.start <= 1e-7 * tmax && interval.end >= 0.0;
    });
    // Radial classification and polynomial expansion can round to opposite sides
    // at the origin; an inside classification still owns the closed point [0, 0].
    if (!found)
        found = std::any_of(leg.begin(), leg.end(),
                            [&](const LegSegmentPair &segment) { return legLocal(segment.current, p).phi <= 0.0; });
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (const auto &interval : intervals)
        {
            if (interval.start <= distance + 1e-7 * tmax && interval.end > distance)
            {
                distance = interval.end;
                changed = true;
            }
        }
    }
    return found;
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegSampleTransport &transport, double kappa,
                          LegContact &constraint, LegScratch &scratch, const LegContactGeometry &geometry)
{
    scratch.clear();
    if (leg.segments.empty())
        return false;
    LegLocal values[2];
    double phiMin = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < leg.segments.size(); ++k)
    {
        values[k] = legLocal(leg.segments[k].current, point);
        phiMin = (std::min)(phiMin, values[k].phi);
    }
    const bool inside = phiMin <= 0.0;
    const LegVec3 raw = legDirection(leg, values, transport, kappa, inside ? 0.0 : phiMin);
    const double magnitude = raw.length();
    double normalWeight;
    if (!inside)
    {
        if (!legNormal(leg, point, values, kappa, constraint.normal, constraint.distance, normalWeight, geometry))
            return false;
    }
    else
    {
        if (magnitude <= 1e-12)
            return false;
        const LegVec3 n = raw / magnitude;
        double d;
        if (!legUnionExit(leg, point, n, d, scratch))
            return false;
        const LegVec3 q = point + d * n;
        double unused;
        if (!legNormal(leg, q, kappa, constraint.normal, unused, normalWeight, geometry) ||
            constraint.normal * n <= 1e-6)
            return false;
        constraint.distance = constraint.normal * (q - point);
    }
    constraint.weight = smoothRamp(magnitude / 0.25) * normalWeight;
    return true;
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegSampleTransport &transport, double kappa,
                          LegContact &constraint, LegScratch &scratch)
{
    return legConstraint(leg, point, transport, kappa, constraint, scratch, legPrepareContactGeometry(leg, kappa));
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegSampleTransport &transport, double kappa,
                          LegContact &constraint)
{
    LegScratch scratch;
    return legConstraint(leg, point, transport, kappa, constraint, scratch);
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegVec3 &restPoint, const LegVec3 &waistRadial,
                          double kappa, LegContact &constraint, LegScratch &scratch,
                          const LegContactGeometry &geometry)
{
    LegSampleTransport transport;
    legSampleTransport(leg, restPoint, waistRadial, transport);
    return legConstraint(leg, point, transport, kappa, constraint, scratch, geometry);
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegVec3 &restPoint, const LegVec3 &waistRadial,
                          double kappa, LegContact &constraint, LegScratch &scratch)
{
    return legConstraint(leg, point, restPoint, waistRadial, kappa, constraint, scratch,
                         legPrepareContactGeometry(leg, kappa));
}

inline bool legConstraint(const Leg &leg, const LegVec3 &point, const LegVec3 &restPoint, const LegVec3 &waistRadial,
                          double kappa, LegContact &constraint)
{
    LegScratch scratch;
    return legConstraint(leg, point, restPoint, waistRadial, kappa, constraint, scratch);
}

inline double legLocalRadius(const std::vector<LegSegment> &rest, const LegVec3 &point)
{
    if (rest.empty())
        return 0.0;
    // Two sides with at most two segments each; larger inputs fall back to the heap.
    double fixedExponents[4], fixedRadii[4];
    std::vector<double> heapExponents, heapRadii;
    double *exponents = fixedExponents, *radii = fixedRadii;
    if (rest.size() > 4)
    {
        heapExponents.resize(rest.size());
        heapRadii.resize(rest.size());
        exponents = heapExponents.data();
        radii = heapRadii.data();
    }
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < rest.size(); ++k)
    {
        const auto &cylinder = rest[k];
        const LegLocal l = legLocal(cylinder, point);
        const double zc = clampValue(l.z, 0.0, cylinder.length);
        const double dist = std::hypot(l.r, l.z - zc);
        const double ratio = dist / l.rhoMin;
        const double exponent = -ratio * ratio;
        exponents[k] = exponent;
        radii[k] = l.rhoMin;
        maximum = (std::max)(maximum, exponent);
    }
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t k = 0; k < rest.size(); ++k)
    {
        const double w = std::exp(exponents[k] - maximum);
        numerator += w * radii[k];
        denominator += w;
    }
    return numerator / denominator;
}

struct LegSectionFrame
{
    LegVec3 center, x, z;
};

inline LegSectionFrame legSectionFrame(const LegSegment &segment, double z)
{
    z = clampValue(z, 0.0, segment.length);
    const LegStation radius = legSection(segment, z);
    return {segment.origin + segment.axis * z, segment.x * radius.radiusX, segment.z * radius.radiusZ};
}

inline LegVec3 legSurfacePoint(const LegSegment &segment, double z, double theta)
{
    const LegStation radius = legSection(segment, z);
    return segment.origin + radius.z * segment.axis + radius.radiusX * std::cos(theta) * segment.x +
           radius.radiusZ * std::sin(theta) * segment.z;
}

inline void legRingMatrices(const LegSegment &segment, std::vector<std::array<LegVec3, 4>> &rows,
                            std::vector<std::array<double, 2>> &farMultipliers)
{
    for (std::size_t i = 0; i + 1 < segment.stations.size(); ++i)
    {
        const double z0 = segment.stations[i].z, z1 = segment.stations[i + 1].z;
        const LegStation near = legSection(segment, z0), far = legSection(segment, z1);
        rows.push_back({{segment.x * near.radiusX, segment.axis * (z1 - z0), segment.z * near.radiusZ,
                         segment.origin + segment.axis * z0}});
        farMultipliers.push_back({{near.radiusX > 1e-12 ? far.radiusX / near.radiusX : 1.0,
                                   near.radiusZ > 1e-12 ? far.radiusZ / near.radiusZ : 1.0}});
    }
}
