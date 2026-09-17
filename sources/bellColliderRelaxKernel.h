#pragma once

#include <cmath>

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
} // namespace BellColliderRelax
