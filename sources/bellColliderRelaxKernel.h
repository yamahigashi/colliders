#pragma once

#include <cmath>
#include <limits>

#if defined(_M_X64) || defined(__SSE2__)
#include <emmintrin.h>
#define YDD_RELAX_SSE2 1
#endif

namespace BellColliderRelax
{
struct Scalar
{
    using Value = double;
    using Mask = bool;
    static Value splat(double value)
    {
        return value;
    }
    static Value add(Value a, Value b)
    {
        return a + b;
    }
    static Value sub(Value a, Value b)
    {
        return a - b;
    }
    static Value mul(Value a, Value b)
    {
        return a * b;
    }
    static Value div(Value a, Value b)
    {
        return a / b;
    }
    static Value sqrt(Value value)
    {
        return std::sqrt(value);
    }
    static Mask greater(Value a, Value b)
    {
        return a > b;
    }
    static Mask both(Mask a, Mask b)
    {
        return a && b;
    }
    static Mask andNot(Mask a, Mask b)
    {
        return a && !b;
    }
    static bool any(Mask mask)
    {
        return mask;
    }
    static Value select(Mask mask, Value yes, Value no)
    {
        return mask ? yes : no;
    }
};

#ifdef YDD_RELAX_SSE2
struct Pair
{
    using Value = __m128d;
    using Mask = __m128d;
    static Value splat(double value)
    {
        return _mm_set1_pd(value);
    }
    static Value add(Value a, Value b)
    {
        return _mm_add_pd(a, b);
    }
    static Value sub(Value a, Value b)
    {
        return _mm_sub_pd(a, b);
    }
    static Value mul(Value a, Value b)
    {
        return _mm_mul_pd(a, b);
    }
    static Value div(Value a, Value b)
    {
        return _mm_div_pd(a, b);
    }
    static Value sqrt(Value value)
    {
        return _mm_sqrt_pd(value);
    }
    static Mask greater(Value a, Value b)
    {
        return _mm_cmpgt_pd(a, b);
    }
    static Mask both(Mask a, Mask b)
    {
        return _mm_and_pd(a, b);
    }
    static Mask andNot(Mask a, Mask b)
    {
        return _mm_andnot_pd(b, a);
    }
    static bool any(Mask mask)
    {
        return _mm_movemask_pd(mask) != 0;
    }
    static Value select(Mask mask, Value yes, Value no)
    {
        return _mm_or_pd(_mm_and_pd(mask, yes), _mm_andnot_pd(mask, no));
    }
};
#endif

template <class Ops> struct Vector
{
    using Value = typename Ops::Value;
    Value x, y, z;

    Vector operator-(const Vector &other) const
    {
        return {Ops::sub(x, other.x), Ops::sub(y, other.y), Ops::sub(z, other.z)};
    }
    Vector operator*(Value scale) const
    {
        return {Ops::mul(x, scale), Ops::mul(y, scale), Ops::mul(z, scale)};
    }
    Value dot(const Vector &other) const
    {
        // Preserve the reference x/y/z sum order; do not contract to a fused multiply-add.
        return Ops::add(Ops::add(Ops::mul(x, other.x), Ops::mul(y, other.y)), Ops::mul(z, other.z));
    }
    Value length() const
    {
        return Ops::sqrt(dot(*this));
    }
};

template <class Ops> struct Ring
{
    using Value = typename Ops::Value;
    Vector<Ops> origin, normal, translation, inverseColumns[3];
    Value collision;
    bool capAtRingOrigin = false;
    bool distalEnd = false;
    Value distalLength = Ops::splat(0.0);
    Value distalWidth = Ops::splat(0.0);
};

// Cartesian input is separate from raw output: MPoint += MVector preserves w.
// A lane represents one vertex; reductions never cross vertices or rings.
template <class Ops>
inline Vector<Ops> relax(const Vector<Ops> &raw, const Vector<Ops> &cartesian, const Ring<Ops> &ring)
{
    using Value = typename Ops::Value;
    const Value one = Ops::splat(1.0);
    const Value threshold = Ops::splat(1e-5);
    const Value distance = (cartesian - ring.origin).dot(ring.normal);
    const Vector<Ops> vector = (cartesian - ring.normal * distance) - ring.translation;
    const Value length = vector.length();
    const auto nonzero = Ops::greater(length, threshold);
    if (!Ops::any(nonzero))
        return raw;

    const Vector<Ops> local = {vector.dot(ring.inverseColumns[0]), vector.dot(ring.inverseColumns[1]),
                               vector.dot(ring.inverseColumns[2])};
    const Value localLength = local.length();
    const auto validLocal = Ops::greater(localLength, threshold);
    const Value delta = Ops::select(validLocal, Ops::div(length, Ops::select(validLocal, localLength, one)), one);
    // length > 1e-5 implies Maya's normalize() lensq > 1e-20 condition.
    // Keep reciprocal followed by multiplication, rather than x / length.
    const Vector<Ops> normal = vector * Ops::div(one, Ops::select(nonzero, length, one));
    const Value projectedLength = (normal * delta).length();
    auto push = Ops::both(nonzero, Ops::greater(projectedLength, length));
    if (ring.capAtRingOrigin)
        push = Ops::both(push, Ops::greater(distance, Ops::splat(0.0)));
    if (!Ops::any(push))
        return raw;

    const Vector<Ops> displacement = (normal * Ops::sub(projectedLength, length)) * ring.collision;
    // Selecting the original bits also preserves signed zeros on inactive lanes.
    return {Ops::select(push, Ops::add(raw.x, displacement.x), raw.x),
            Ops::select(push, Ops::add(raw.y, displacement.y), raw.y),
            Ops::select(push, Ops::add(raw.z, displacement.z), raw.z)};
}
template <class Ops>
inline Vector<Ops> relax(const Vector<Ops> &raw, const Vector<Ops> &cartesian, const Ring<Ops> &ring,
                         const Vector<Ops> &d, typename Ops::Value betaScale = Ops::splat(1.0))
{
    using Value = typename Ops::Value;
    const Value zero = Ops::splat(0.0);
    const Value one = Ops::splat(1.0);
    const Value two = Ops::splat(2.0);
    const Value threshold = Ops::splat(1e-5);
    const Value distance = (cartesian - ring.origin).dot(ring.normal);
    const Vector<Ops> vector = (cartesian - ring.normal * distance) - ring.translation;
    const Vector<Ops> local = {vector.dot(ring.inverseColumns[0]), vector.dot(ring.inverseColumns[1]),
                               vector.dot(ring.inverseColumns[2])};
    const Value c2 = Ops::sub(Ops::add(Ops::mul(local.x, local.x), Ops::mul(local.z, local.z)), one);
    auto active = Ops::greater(zero, c2);
    if (ring.capAtRingOrigin)
        active = Ops::both(active, Ops::greater(distance, zero));
    Value g_d = one;
    if (ring.distalEnd)
    {
        const Value L = ring.distalLength;
        const Value w = ring.distalWidth;
        const auto distalMask = Ops::greater(Ops::splat(ring.distalEnd ? 1.0 : 0.0), zero);
        const auto validWidth = Ops::both(distalMask, Ops::greater(w, threshold));
        const Value x = Ops::div(Ops::select(validWidth, Ops::sub(distance, L), zero),
                                 Ops::select(validWidth, w, one));
        const Value t_d = Ops::select(Ops::greater(x, one), one, Ops::select(Ops::greater(x, zero), x, zero));
        g_d = Ops::sub(one, Ops::mul(Ops::mul(t_d, t_d), Ops::sub(Ops::splat(3.0), Ops::mul(two, t_d))));
        const auto exclude = Ops::both(validWidth, Ops::greater(distance, Ops::add(L, w)));
        active = Ops::andNot(active, exclude);
    }
    if (!Ops::any(active))
        return raw;

    const Vector<Ops> perpendicular = d - ring.normal * d.dot(ring.normal);
    const Value D = Ops::sqrt(Ops::select(active, perpendicular.dot(perpendicular), zero));
    auto shifted = Ops::both(active, Ops::greater(D, threshold));
    const Vector<Ops> unit = perpendicular * Ops::div(one, Ops::select(shifted, D, one));
    const Vector<Ops> negative = unit * Ops::splat(-1.0);
    const Value mx = negative.dot(ring.inverseColumns[0]);
    const Value mz = negative.dot(ring.inverseColumns[2]);
    const Value norm = Ops::add(Ops::mul(mx, mx), Ops::mul(mz, mz));
    shifted = Ops::both(shifted, Ops::both(Ops::greater(norm, zero),
                                           Ops::greater(Ops::splat(std::numeric_limits<double>::infinity()), norm)));
    const Value radius = Ops::div(one, Ops::sqrt(Ops::select(shifted, norm, one)));
    const Value ratio = Ops::div(Ops::select(shifted, D, zero), radius);
    const Value positive = Ops::select(Ops::greater(ratio, zero), ratio, zero);
    const Value t = Ops::select(Ops::greater(positive, one), one, positive);
    const Value beta = Ops::select(
        shifted,
        Ops::mul(Ops::mul(Ops::mul(Ops::splat(0.5), Ops::mul(t, t)), Ops::sub(Ops::splat(3.0), Ops::mul(two, t))),
                 betaScale),
        zero);
    const Value shift = Ops::mul(Ops::sub(zero, beta), radius);
    const Vector<Ops> c = {Ops::select(shifted, Ops::mul(unit.x, shift), zero),
                           Ops::select(shifted, Ops::mul(unit.y, shift), zero),
                           Ops::select(shifted, Ops::mul(unit.z, shift), zero)};
    const Vector<Ops> v = vector - c;
    const Value vLength = Ops::sqrt(Ops::select(active, v.dot(v), zero));
    active = Ops::both(active, Ops::greater(vLength, threshold));
    if (!Ops::any(active))
        return raw;
    const Vector<Ops> q = v * Ops::div(one, Ops::select(active, vLength, one));
    const Vector<Ops> qLocal = {q.dot(ring.inverseColumns[0]), q.dot(ring.inverseColumns[1]),
                                q.dot(ring.inverseColumns[2])};
    const Value a = Ops::add(Ops::mul(qLocal.x, qLocal.x), Ops::mul(qLocal.z, qLocal.z));
    const Value b = Ops::mul(two, Ops::add(Ops::mul(local.x, qLocal.x), Ops::mul(local.z, qLocal.z)));
    active = Ops::both(active, Ops::greater(a, Ops::splat(1e-12)));
    if (!Ops::any(active))
        return raw;
    const Value discriminant = Ops::sub(Ops::mul(b, b), Ops::mul(Ops::mul(Ops::splat(4.0), a), c2));
    const Value disc = Ops::select(Ops::both(active, Ops::greater(discriminant, zero)), discriminant, zero);
    const Value s = Ops::div(Ops::add(Ops::sub(zero, b), Ops::sqrt(disc)), Ops::mul(two, Ops::select(active, a, one)));
    active = Ops::both(active, Ops::greater(s, zero));
    constexpr double rho = 0.75;
    const Value fadeRatio = Ops::div(Ops::select(shifted, vLength, zero), Ops::mul(Ops::splat(rho), radius));
    const Value fadePositive = Ops::select(Ops::greater(fadeRatio, zero), fadeRatio, zero);
    const Value u = Ops::select(Ops::greater(fadePositive, one), one, fadePositive);
    const Value smoothU = Ops::mul(Ops::mul(u, u), Ops::sub(Ops::splat(3.0), Ops::mul(two, u)));
    const Value smoothT = Ops::mul(Ops::mul(t, t), Ops::sub(Ops::splat(3.0), Ops::mul(two, t)));
    const Value g = Ops::select(shifted, Ops::sub(one, Ops::mul(smoothT, Ops::sub(one, smoothU))), one);
    Vector<Ops> displacement = ((q * Ops::select(active, s, zero)) * g) * ring.collision;
    auto keep = active;
    if (ring.distalEnd)
    {
        displacement = displacement * g_d;
        keep = Ops::both(active, Ops::greater(g_d, zero));
    }
    // Returning raw bits avoids adding zero to a signed zero when the gain vanishes.
    return {Ops::select(keep, Ops::add(raw.x, displacement.x), raw.x),
            Ops::select(keep, Ops::add(raw.y, displacement.y), raw.y),
            Ops::select(keep, Ops::add(raw.z, displacement.z), raw.z)};
}

} // namespace BellColliderRelax
