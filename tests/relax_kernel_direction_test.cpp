#include "bellColliderRelaxKernel.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using Ops = BellColliderRelax::Scalar;
using Vec = BellColliderRelax::Vector<Ops>;
using Ring = BellColliderRelax::Ring<Ops>;
static int checks = 0;
static int failures = 0;

static bool sameBits(double a, double b)
{
    std::uint64_t aa, bb;
    std::memcpy(&aa, &a, sizeof(aa));
    std::memcpy(&bb, &b, sizeof(bb));
    return aa == bb;
}

static void check(bool result, const char *name)
{
    ++checks;
    if (!result)
    {
        ++failures;
        std::cerr << "FAIL " << name << '\n';
    }
}

static void near(const Vec &a, const Vec &b, double tolerance, const char *name)
{
    check(std::abs(a.x - b.x) <= tolerance && std::abs(a.y - b.y) <= tolerance && std::abs(a.z - b.z) <= tolerance,
          name);
}

static double smooth(double x)
{
    x = (std::max)(0.0, (std::min)(1.0, x));
    return x * x * (3.0 - 2.0 * x);
}

static Ring circle(double radius)
{
    Ring ring;
    ring.origin = ring.translation = {0, 0, 0};
    ring.normal = {0, 1, 0};
    ring.inverseColumns[0] = {1 / radius, 0, 0};
    ring.inverseColumns[1] = {0, 1, 0};
    ring.inverseColumns[2] = {0, 0, 1 / radius};
    ring.collision = 1;
    return ring;
}

static Vec apply(const Vec &p, const Ring &ring, const Vec &d, double scale = 1)
{
    return BellColliderRelax::relax(p, p, ring, d, scale);
}

#ifdef YDD_RELAX_SSE2
static void pairCheck(const Ring &ring, const Vec &p, const Vec &other, const Vec &d, const Vec &otherD,
                      double scale = 1, double otherScale = 1)
{
    using Pair = BellColliderRelax::Pair;
    using PV = BellColliderRelax::Vector<Pair>;
    const auto broadcast = [](const Vec &v) -> PV { return {Pair::splat(v.x), Pair::splat(v.y), Pair::splat(v.z)}; };
    const auto pack = [](const Vec &a, const Vec &b) -> PV {
        return {_mm_set_pd(b.x, a.x), _mm_set_pd(b.y, a.y), _mm_set_pd(b.z, a.z)};
    };
    BellColliderRelax::Ring<Pair> pr;
    pr.origin = broadcast(ring.origin);
    pr.translation = broadcast(ring.translation);
    pr.normal = broadcast(ring.normal);
    for (int i = 0; i < 3; ++i)
        pr.inverseColumns[i] = broadcast(ring.inverseColumns[i]);
    pr.collision = Pair::splat(ring.collision);
    pr.capAtRingOrigin = ring.capAtRingOrigin;
    pr.distalEnd = ring.distalEnd;
    pr.distalLength = Pair::splat(ring.distalLength);
    pr.distalWidth = Pair::splat(ring.distalWidth);
    for (int reverse = 0; reverse < 2; ++reverse)
    {
        const Vec points[] = {reverse ? other : p, reverse ? p : other};
        const Vec directions[] = {reverse ? otherD : d, reverse ? d : otherD};
        const double scales[] = {reverse ? otherScale : scale, reverse ? scale : otherScale};
        const PV raw = pack(points[0], points[1]);
        const PV output = BellColliderRelax::relax(raw, raw, pr, pack(directions[0], directions[1]),
                                                   _mm_set_pd(scales[1], scales[0]));
        double x[2], y[2], z[2];
        _mm_storeu_pd(x, output.x);
        _mm_storeu_pd(y, output.y);
        _mm_storeu_pd(z, output.z);
        for (int lane = 0; lane < 2; ++lane)
        {
            const Vec expected = apply(points[lane], ring, directions[lane], scales[lane]);
            check(sameBits(expected.x, x[lane]) && sameBits(expected.y, y[lane]) && sameBits(expected.z, z[lane]),
                  "Scalar/Pair bits and lane independence");
        }
    }
}
#endif

static void bits(const Vec &a, const Vec &b, const char *name)
{
    check(sameBits(a.x, b.x) && sameBits(a.y, b.y) && sameBits(a.z, b.z), name);
}

static void knee_end_gain_continuity()
{
    Ring ring = circle(2);
    ring.distalEnd = true;
    ring.distalLength = 2;
    ring.distalWidth = 2;
    Ring previous = ring;
    previous.distalEnd = false;
    const Vec d = {4, 0, 0};
    const auto correction = [&](double distance) {
        const Vec p = {-0.25, distance, 0};
        check((p - ring.origin).dot(ring.normal) == distance, "continuity exact axial distance");
        return apply(p, ring, d) - p;
    };
    for (double boundary : {2.0, 4.0})
    {
        const Vec at = correction(boundary);
        double previousJump = std::numeric_limits<double>::infinity();
        double previousLeft = previousJump, previousRight = previousJump;
        for (int exponent = 3; exponent <= 20; ++exponent)
        {
            const double h = std::ldexp(1.0, -exponent);
            const Vec left = correction(boundary - h), right = correction(boundary + h);
            const double jump = (right - left).length();
            const double leftSlope = (at - left).length() / h;
            const double rightSlope = (right - at).length() / h;
            check(jump <= previousJump * 0.3, "end correction difference converges quadratically");
            check(leftSlope <= previousLeft * 0.6 && rightSlope <= previousRight * 0.6,
                  "both correction slopes converge to zero");
            previousJump = jump;
            previousLeft = leftSlope;
            previousRight = rightSlope;
        }
        check(previousLeft < 1e-6 && previousRight < 1e-6, "end limiting correction slopes");
    }
    for (double distance : {1.0, 2.0})
    {
        const Vec p = {-0.25, distance, -0.0};
        check((p - ring.origin).dot(ring.normal) <= ring.distalLength, "before band distance");
        bits(apply(p, ring, d), apply(p, previous, d), "before band unchanged bits");
    }
    check(correction(4).length() <= 1e-12, "end correction bound");
    for (double distance : {4.001, 6.0})
    {
        const Vec p = {-0.25, distance, -0.0};
        check((p - ring.origin).dot(ring.normal) > ring.distalLength + ring.distalWidth, "beyond band distance");
        bits(apply(p, ring, d), p, "beyond band raw bits");
    }
    ring.distalLength = ring.distalWidth = 0.5;
    const Vec rounded = {-0.0, std::nextafter(1.0, 0.0), -0.0};
    const double t = ((rounded - ring.origin).dot(ring.normal) - ring.distalLength) / ring.distalWidth;
    check(t < 1 && 1 - ((t * t) * (3 - 2 * t)) == 0, "subunit distal parameter rounds gain to zero");
    bits(apply(rounded, ring, d), rounded, "rounded zero gain raw bits");
#ifdef YDD_RELAX_SSE2
    const Vec band = {-0.25, 0.75, 0};
    const double bandT = ((band - ring.origin).dot(ring.normal) - ring.distalLength) / ring.distalWidth;
    const double bandGain = 1 - ((bandT * bandT) * (3 - 2 * bandT));
    check(bandGain > 0 && bandGain < 1, "mixed lane has partial distal gain");
    const Vec excluded = {-0.0, 1.5, -0.0};
    check((excluded - ring.origin).dot(ring.normal) > ring.distalLength + ring.distalWidth,
          "mixed lane beyond distal support");
    pairCheck(ring, rounded, excluded, d, d);
    pairCheck(ring, band, rounded, d, d);
#endif
}

static void knee_end_scalar_pair_agreement()
{
    const Vec recordedPoints[] = {{-6.15 + 0.5, 1, 0}, {-6.15 + 0.01, 1, 0}, {3, 1, 0},
                                  {7, 1, 0}, {0.5, -1, 0}, {7, 1, -0.0}};
    const Vec recordedOutputs[] = {{-5.9435013824400338, 1, 0}, {-6.1473784536522995, 1, 0},
                                   {6.1500000000000004, 1, 0}, {7, 1, 0}, {0.5, -1, 0}, {7, 1, -0.0}};
    Ring recordedRing = circle(6.15);
    recordedRing.capAtRingOrigin = true;
    recordedRing.distalLength = 2;
    recordedRing.distalWidth = 2;
    const Vec recordedDirection = {12.3, 0, 0};
    for (bool distal : {false, true})
    {
        recordedRing.distalEnd = distal;
        for (unsigned int i = 0; i < sizeof(recordedPoints) / sizeof(recordedPoints[0]); ++i)
        {
            const Vec &p = recordedPoints[i];
            check((p - recordedRing.origin).dot(recordedRing.normal) <= recordedRing.distalLength,
                  "stored output fixture within full support");
            bits(apply(p, recordedRing, recordedDirection), recordedOutputs[i], "stored kernel output bits");
#ifdef YDD_RELAX_SSE2
            pairCheck(recordedRing, p, recordedPoints[(i + 1) % 6], recordedDirection, recordedDirection);
#endif
        }
    }
    const Vec d = {4, 0, 0};
    const Vec points[] = {{-0.25, 1, -0.0}, {-0.25, 2, 0.2}, {-0.25, 3, 0.2},
                          {-0.0, 4, -0.0}, {-0.0, 5, -0.0}, {3, 3, -0.0},
                          {-0.0, -1, -0.0}, {-0.0, -0.0, -0.0}};
    for (bool cap : {false, true})
    {
        Ring baseline = circle(2);
        baseline.capAtRingOrigin = cap;
        Ring ring = baseline;
        ring.distalLength = 2;
        ring.distalWidth = 2;
        for (const Vec &p : points)
        {
            bits(apply(p, ring, d), apply(p, baseline, d), "disabled distal fields ignored");
#ifdef YDD_RELAX_SSE2
            pairCheck(ring, p, {-0.25, 3, 0}, d, d);
#endif
        }
        ring.distalEnd = true;
        for (const Vec &p : points)
        {
            if ((p - ring.origin).dot(ring.normal) <= ring.distalLength)
                bits(apply(p, ring, d), apply(p, baseline, d), "full support old path bits");
#ifdef YDD_RELAX_SSE2
            for (const Vec &other : points)
                pairCheck(ring, p, other, d, {0, 0, 0}, 1, 0.25);
#endif
        }
        const Vec interior = {-0.25, 3, 0};
        bits(BellColliderRelax::relax(interior, interior, ring),
             BellColliderRelax::relax(interior, interior, baseline), "old overload ignores distal band");
        check(!sameBits(apply(interior, ring, d).x, apply(interior, baseline, d).x),
              "old overload fixture has nontrivial distal gain");
        Ring nonfinite = baseline;
        nonfinite.distalEnd = true;
        nonfinite.distalLength = nonfinite.distalWidth = std::numeric_limits<double>::quiet_NaN();
        bits(BellColliderRelax::relax(interior, interior, nonfinite),
             BellColliderRelax::relax(interior, interior, baseline), "old overload ignores nonfinite distal fields");
        for (double width : {0.0, 0.5e-5, 1e-5, std::numeric_limits<double>::quiet_NaN()})
        {
            ring.distalWidth = width;
            for (const Vec &p : points)
            {
                bits(apply(p, ring, d), apply(p, baseline, d), "degenerate width keeps old correction");
#ifdef YDD_RELAX_SSE2
                pairCheck(ring, p, interior, d, d);
#endif
            }
        }
    }
}

static void knee_end_directional_inward_ray()
{
    Ring ring = circle(6.15);
    ring.distalEnd = true;
    ring.distalLength = 0.75;
    ring.distalWidth = 0.5;
    const Vec d = {12.3, 0, 0};
    for (double collision : {0.15, 0.65, 1.0})
        for (double depth : {0.5, 0.1, 0.01})
        {
            ring.collision = collision;
            const Vec p = {-6.15 + depth, 1, 0};
            const double distance = (p - ring.origin).dot(ring.normal);
            const double t = (distance - ring.distalLength) / ring.distalWidth;
            const double g_d = 1 - (t * t) * (3 - 2 * t);
            const double g = smooth((3.075 - depth) / (0.75 * 6.15));
            check(g > 0 && g < 1 && g_d > 0 && g_d < 1, "inward ray both gains active");
            const double push = ((depth * g) * collision) * g_d;
            near(apply(p, ring, d, 1), {p.x - push, p.y, p.z}, 1e-12, "inward ray q s g collision distal gain");
#ifdef YDD_RELAX_SSE2
            pairCheck(ring, p, {0, 2, -0.0}, d, d);
#endif
        }
}

int main()
{
    std::cout << std::setprecision(17);
    knee_end_gain_continuity();
    knee_end_scalar_pair_agreement();
    knee_end_directional_inward_ray();
    const Vec zero = {0, 0, 0};
    const Vec d = {12.3, 0, 0};
    Ring ring = circle(6.15);
    const Vec axis = {0, 1, 0};
    const Vec exit = apply(axis, ring, d);
    near(exit, {6.15 * smooth(0.5 / 0.75), 1, 0}, 1e-9, "axis exit");
    std::cout << "T2 axis exit: " << exit.x << ',' << exit.y << ',' << exit.z << '\n';
    for (double depth : {0.5, 0.1, 0.01})
    {
        const Vec p = {-6.15 + depth, 1, 0};
        const Vec out = apply(p, ring, d);
        const double push = depth * smooth((3.075 - depth) / (0.75 * 6.15));
        near(out, {p.x - push, 1, 0}, 1e-9, "entry exit");
        check(std::abs((out - p).length() - push) <= 1e-9, "entry depth");
        std::cout << "T2 depth " << depth << " displacement " << (out - p).length() << '\n';
    }
    const Vec centre = {-3.075, 1, 0};
    near(apply(centre, ring, d), centre, 0, "shifted centre stop");
    std::cout << "T2 shifted centre displacement " << (apply(centre, ring, d) - centre).length() << '\n';
    near(apply({3, 1, 0}, ring, d), {6.15, 1, 0}, 1e-9, "far side");
    for (double x : {6.15, 7.0})
        near(apply({x, 1, 0}, ring, d), {x, 1, 0}, 1e-9, "boundary/outside");
    for (const Vec direction : {zero, Vec{0, 12.3, 0}})
    {
        near(apply(axis, ring, direction), axis, 0, "axis degenerate direction");
        near(apply({3, 1, 0}, ring, direction), {6.15, 1, 0}, 1e-12, "radial fallback");
    }
    for (double distance : {-1e-12, -0.0, 0.0, 1e-12})
    {
        Ring cap = ring;
        cap.capAtRingOrigin = true;
        cap.origin.y = 1 - distance;
        cap.translation = cap.origin;
        near(apply(axis, cap, d), distance > 0 ? Vec{6.15 * smooth(0.5 / 0.75), 1, 0} : axis, 1e-9, "cap distance");
    }
    ring = circle(2);
    const Vec rampPoint = {0, 1, 1};
    const Vec rampDirection = {1, 0, 0};
    const double root = (-1 / std::sqrt(5.0) + std::sqrt(0.95)) / 0.5;
    const double rampFade = 1 - smooth(0.5) * (1 - smooth(std::sqrt(5.0) / 3));
    const Vec rampExpected = {root / std::sqrt(5.0) * rampFade, 1, 1 + 2 * root / std::sqrt(5.0) * rampFade};
    const Vec ramp = apply(rampPoint, ring, rampDirection);
    near(ramp, rampExpected, 1e-12, "ramp beta .25 centre -.5");
    std::cout << "T3 ramp beta 0.25 centre -0.5 output " << ramp.x << ',' << ramp.y << ',' << ramp.z << '\n';
    for (double collision : {0.15, 0.65, 1.0})
    {
        Ring fadeRing = circle(2);
        fadeRing.collision = collision;
        const Vec direction = {4, 0, 0};
        const Vec singular = {-1, 1, 0};
        const double above = std::nextafter(1.5, 2.0);
        const Vec outsideFade = {above - 1, 1, 0};
        const Vec halfway = {-0.25, 1, 0};
        near(apply(singular, fadeRing, direction), singular, 0, "fade zero vector");
        const double q = above * (1.0 / std::sqrt(above * above));
        const double px = (above - 1) * 0.5, qx = q * 0.5;
        const double a = qx * qx, b = 2 * (px * qx), c2 = px * px - 1;
        const double s = ((0.0 - b) + std::sqrt(b * b - (4 * a) * c2)) / (2 * a);
        const Vec previous = {outsideFade.x + (q * s) * collision, 1, 0};
        const Vec full = apply(outsideFade, fadeRing, direction);
        near(full, previous, 1e-12, "fade above radius unchanged directional exit");
        check(sameBits(full.x, previous.x) && sameBits(full.y, previous.y) && sameBits(full.z, previous.z),
              "fade above radius directional bits");
        near(apply(halfway, fadeRing, direction), {-0.25 + 2.25 * smooth(0.5) * collision, 1, 0},
             1e-12, "fade half radius");
        near(apply(halfway, fadeRing, zero), {-0.25 - 1.75 * collision, 1, 0}, 1e-12,
             "unshifted fade is one");
#ifdef YDD_RELAX_SSE2
        pairCheck(fadeRing, singular, halfway, direction, direction);
        pairCheck(fadeRing, outsideFade, halfway, direction, direction);
        pairCheck(fadeRing, halfway, halfway, zero, direction);
        pairCheck(fadeRing, halfway, halfway, direction, direction);
#endif
    }
    const double infinity = std::numeric_limits<double>::infinity();
    for (double threshold : {std::nextafter(1e-5, 0.0), std::nextafter(1e-5, infinity)})
    {
        const Vec p = {threshold, 1, 0};
        const Vec result = apply(p, circle(1), zero);
        near(result, threshold <= 1e-5 ? p : Vec{1, 1, 0}, 1e-12, "v threshold");
        near(apply(axis, circle(1), {threshold, 0, 0}), axis, 0, "small D axis stop");
#ifdef YDD_RELAX_SSE2
        pairCheck(circle(1), p, {0.5, 1, 0}, zero, d);
        pairCheck(circle(1), {0.5, 1, 0}, {0.5, 1, 0}, {threshold, 0, 0}, d);
        pairCheck(circle(1), p, p, zero, zero);
#endif
    }
    for (double localQ : {std::nextafter(1e-6, 0.0), std::nextafter(1e-6, infinity)})
    {
        Ring large = circle(1);
        large.inverseColumns[0].x = localQ;
        const Vec p = {0.5, 1, 0};
        near(apply(p, large, zero), localQ * localQ <= 1e-12 ? p : Vec{1 / localQ, 1, 0}, 1e-9, "a threshold");
#ifdef YDD_RELAX_SSE2
        pairCheck(large, p, {0, 1, 0.5}, zero, zero);
        pairCheck(large, p, p, zero, zero);
#endif
    }
    double lo = 0, hi = 1;
    for (int i = 0; i < 80; ++i)
    {
        const double mid = (lo + hi) * 0.5;
        if (0.5 * (mid * mid) * (3 - 2 * mid) <= 1e-5)
            lo = mid;
        else
            hi = mid;
    }
    lo = std::nextafter(lo, 0.0);
    hi = std::nextafter(hi, infinity);
    near(apply(axis, circle(1), {lo, 0, 0}), axis, 0, "axis shift below threshold");
    const double axisFade = 1 - smooth(hi) * (1 - smooth(0.5 * smooth(hi) / 0.75));
    near(apply(axis, circle(1), {hi, 0, 0}), {axisFade, 1, 0}, 1e-12, "axis shift above threshold");
    std::cout << "T3 axis D threshold below " << lo << " above " << hi << " output x " << axisFade << '\n';
    for (double radius : {1.0, 2.0, 6.15})
        for (double x : {0.1, 0.5, 0.9})
        {
            const Vec p = {radius * x, 1, radius * 0.1};
            const Ring current = circle(radius);
            near(apply(p, current, zero), BellColliderRelax::relax(p, p, current), 1e-12, "zero old equality");
        }
    for (double norm : {0.0, infinity, std::numeric_limits<double>::quiet_NaN()})
    {
        Ring degenerate = circle(1);
        // A finite inverse column can overflow its squared directional norm.
        degenerate.inverseColumns[0] = {norm == infinity ? 1e308 : norm, 0, 0};
        const Vec p = {0, 1, 0.5};
        const Vec result = apply(p, degenerate, d);
        check(std::isfinite(result.x) && std::isfinite(result.y) && std::isfinite(result.z),
              "degenerate finite result");
        near(result, std::isnan(norm) ? p : Vec{0, 1, 1}, 1e-12, "degenerate zero centre");
#ifdef YDD_RELAX_SSE2
        pairCheck(degenerate, p, axis, d, zero);
#endif
    }
#ifdef YDD_RELAX_SSE2
    for (const Vec p : {axis, centre, Vec{7, 1, 0}, Vec{0.1, -1, 0}, Vec{-0.0, 0, -0.0}, Vec{infinity, 1, 0},
                        Vec{std::numeric_limits<double>::quiet_NaN(), 1, 0}})
        for (bool cap : {false, true})
        {
            Ring current = circle(6.15);
            current.capAtRingOrigin = cap;
            pairCheck(current, p, {0.5, 1, 0.25}, d, d, 1, 0.25);
            pairCheck(current, p, p, d, d);
        }
    for (double collision : {0.15, 0.65, 1.0, infinity, std::numeric_limits<double>::quiet_NaN()})
    {
        Ring current = circle(6.15);
        current.collision = collision;
        pairCheck(current, {0.5, 1, 0.25}, {7, 1, -0.0}, d, zero);
    }
    Ring ellipse = circle(1);
    const double c = std::cos(0.37), sn = std::sin(0.37);
    ellipse.normal = {-sn, c, 0};
    ellipse.origin = ellipse.translation = {0.4, -0.7, 0.8};
    ellipse.inverseColumns[0] = {c / 9.61, sn / 9.61, 0};
    ellipse.inverseColumns[1] = {-sn, c, 0};
    ellipse.inverseColumns[2] = {0, 0, 1 / 6.06};
    for (int i = 0; i < 32; ++i)
    {
        const double angle = i * 0.19;
        const Vec p = {0.4 + c * 2 * std::cos(angle) - sn, -0.7 + sn * 2 * std::cos(angle) + c,
                       0.8 + 3 * std::sin(angle)};
        near(apply(p, ellipse, zero), BellColliderRelax::relax(p, p, ellipse), 1e-12, "ellipse zero equality");
        pairCheck(ellipse, p, {8, 3, -4}, {2, 0.5, -3}, {-1, 2, 3}, 0.25, 0.75);
    }
    pairCheck(circle(1), axis, axis, {lo, 0, 0}, {hi, 0, 0});
#else
    check(false, "SSE2 Pair unavailable");
#endif
    std::cout << checks - failures << " checks passed; " << failures << " failed\n";
    return failures ? 1 : 0;
}
