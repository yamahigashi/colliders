"""Pure Python per-leg geometry fields and weighted, coupled collision corrections."""

import math
import warnings


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, value):
    return tuple(x * value for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def length(a):
    return math.sqrt(dot(a, a))


def unit(a):
    return mul(a, 1.0 / length(a))


def radial(offset, axis):
    return sub(offset, mul(axis, dot(offset, axis)))


def rotate_minimum(vector, source, target):
    cosine = dot(source, target)
    if cosine < -1 + 1e-6:
        return vector
    pivot = cross(source, target)
    return add(add(vector, cross(pivot, vector)), mul(cross(pivot, cross(pivot, vector)), 1 / (1 + cosine)))


def profile_stations(thigh=5.0, calf=5.0, **profile):
    knee = thigh / (thigh + calf)
    positions = dict(
        hip=0.0,
        thigh=knee * profile.get("thighPosition", 0.5),
        knee=knee,
        calf=knee + (1 - knee) * profile.get("calfPosition", 0.5),
        ankle=1.0,
    )
    selected = []
    for name in ("ankle", "knee", "thigh", "calf", "hip"):
        s = positions[name]
        if all(abs(s - entry[0]) > 1e-9 for entry in selected):
            selected.append((s, profile.get(name + "RadiusX", 1.0), profile.get(name + "RadiusZ", 1.0)))
    return sorted(selected)


IDENTITY = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
ZERO = (0.0, 0.0, 0.0)


def finite(vector):
    return all(math.isfinite(value) for value in vector)


def determinant(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def row_vector(vector, matrix):
    return tuple(sum(vector[i] * matrix[i][j] for i in range(3)) for j in range(3))


def inverse(matrix):
    volume = dot(matrix[0], cross(matrix[1], matrix[2]))
    columns = (cross(matrix[1], matrix[2]), cross(matrix[2], matrix[0]), cross(matrix[0], matrix[1]))
    return tuple(tuple(columns[j][i] / volume for j in range(3)) for i in range(3))


def joint_rotation(matrix):
    norms = [length(row) for row in matrix]
    if any(not finite(row) for row in matrix) or any(not math.isfinite(c) or c < 1e-8 for c in norms):
        return None
    rotation = tuple(mul(row, 1 / c) for row, c in zip(matrix, norms))
    det = determinant(rotation)
    return rotation if math.isfinite(det) and abs(det) >= 0.5 else None


class Cylinder:
    def __init__(self, origin, axis, x_axis, z_axis, extent, stations, joint=IDENTITY):
        self.origin = tuple(origin)
        self.axis = unit(axis)
        self.x_axis = unit(x_axis)
        self.z_axis = unit(z_axis)
        self.extent = extent
        self.stations = tuple(tuple(station) for station in stations)
        self.joint = tuple(tuple(row) for row in joint)
        if (
            not math.isfinite(extent)
            or extent < 1e-6
            or len(stations) < 2
            or any(not finite(station) or min(station[1:]) <= 0 for station in stations)
        ):
            raise ValueError("Invalid cylinder")

    def ellipse_support(self, normal, station):
        _, rx, rz = station
        x, z = rx * dot(normal, self.x_axis), rz * dot(normal, self.z_axis)
        return math.sqrt(x * x + z * z)

    def support(self, normal):
        return dot(normal, self.origin) + max(
            station[0] * dot(normal, self.axis) + self.ellipse_support(normal, station) for station in self.stations
        )

    def section_radii(self, z):
        if z <= self.stations[0][0]:
            return self.stations[0][1:]
        for first, second in zip(self.stations, self.stations[1:]):
            if z <= second[0]:
                t = (z - first[0]) / (second[0] - first[0])
                return tuple(first[i] + t * (second[i] - first[i]) for i in (1, 2))
        return self.stations[-1][1:]

    def inside(self, point, tolerance=0.0):
        offset = sub(point, self.origin)
        z = dot(offset, self.axis)
        a, b = self.section_radii(z)
        return (
            0 < z < self.extent
            and math.hypot(dot(offset, self.x_axis) / a, dot(offset, self.z_axis) / b) < 1 - tolerance
        )


def rest_contact(point, cylinder, waist=(0, 11, 0), waist_axis=(0, -1, 0)):
    waist_radial = radial(sub(point, waist), unit(waist_axis))
    if length(waist_radial) < 1e-8:
        return None
    outward = radial(unit(waist_radial), cylinder.axis)
    if length(outward) < 1e-6:
        outward = radial(sub(point, cylinder.origin), cylinder.axis)
    return unit(outward) if finite(outward) and length(outward) >= 1e-8 else None


def transported_side(b0, current, rest, quaternion_only=False):
    if b0 is None:
        return None, "degenerate"
    r0, r1 = joint_rotation(rest.joint), joint_rotation(current.joint)
    if not quaternion_only and r0 is not None and r1 is not None:
        side = radial(row_vector(row_vector(b0, inverse(r0)), r1), current.axis)
        if finite(side) and math.isfinite(length(side)) and length(side) >= 1e-6:
            return unit(side), "joint"
    side = radial(rotate_minimum(b0, rest.axis, current.axis), current.axis)
    if finite(side) and math.isfinite(length(side)) and length(side) >= 1e-6:
        return unit(side), "quaternion"
    return None, "degenerate"


def smooth_ramp(value):
    t = min(1.0, max(0.0, value))
    return t * t * (3.0 - 2.0 * t)


S = smooth_ramp


def C(t):
    return min(1.0, max(0.0, t))


def det3(M):
    return determinant(M)


def project_w(cons):
    planes = [(n, 1.25 * h, 4 * s) for n, h, s in cons if s != 0]
    best, best_obj = (ZERO, math.inf)
    best_any, best_any_obj = (ZERO, math.inf)
    for mask in range(1 << len(planes)):
        act = [i for i in range(len(planes)) if mask & 1 << i]
        A = [[1.0 if r == c else 0.0 for c in range(3)] for r in range(3)]
        rhs = [0.0] * 3
        for i in act:
            n, hh, w = planes[i]
            for r in range(3):
                rhs[r] += w * hh * n[r]
                for c in range(3):
                    A[r][c] += w * n[r] * n[c]
        cand = ZERO
        if mask:
            dt = det3(A)
            cand = tuple(
                (det3([[rhs[r] if j == c else A[r][j] for j in range(3)] for r in range(3)]) / dt for c in range(3))
            )
        ok = True
        obj = dot(cand, cand)
        for i, (n, hh, w) in enumerate(planes):
            pr = dot(n, cand)
            if pr > hh + 1.25e-09 if i in act else pr < hh - 1.25e-09:
                ok = False
            res = max(0.0, hh - pr)
            obj += w * res * res
        if obj < best_any_obj:
            best_any, best_any_obj = (cand, obj)
        if ok and obj < best_obj:
            best, best_obj = (cand, obj)
    return best if math.isfinite(best_obj) else best_any


class Seg:
    def __init__(self, c, c0, ident, free_start, free_end):
        self.c, self.c0, self.id, self.free_start, self.free_end = (c, c0, ident, free_start, free_end)


def sec(c, z):
    return c.section_radii(min(max(z, 0.0), c.extent))


def slopes(c, z, kappa):
    z = min(max(z, 0.0), c.extent)
    st = c.stations
    if len(st) < 2:
        return (0.0, 0.0)
    ms = [
        (
            (st[j + 1][1] - st[j][1]) / (st[j + 1][0] - st[j][0]),
            (st[j + 1][2] - st[j][2]) / (st[j + 1][0] - st[j][0]),
        )
        for j in range(len(st) - 1)
    ]
    a, b = ms[0]
    for j in range(1, len(st) - 1):
        f = kappa * min(st[j][1], st[j][2])
        w = S(C(0.5 + (z - st[j][0]) / (2 * f))) if f > 0 else 1.0 if z >= st[j][0] else 0.0
        a += (ms[j][0] - ms[j - 1][0]) * w
        b += (ms[j][1] - ms[j - 1][1]) * w
    return (a, b)


def local(c, X):
    v = sub(X, c.origin)
    z = dot(v, c.axis)
    x, y = (dot(v, c.x_axis), dot(v, c.z_axis))
    A, B = sec(c, z)
    eta = 1 - math.hypot(x / A, y / B)
    r = math.hypot(x, y)
    e = c.x_axis if r < 1e-08 else mul(add(mul(c.x_axis, x), mul(c.z_axis, y)), 1.0 / r)
    ex, ez = (dot(e, c.x_axis), dot(e, c.z_axis))
    rho = 1 / math.sqrt(ex * ex / (A * A) + ez * ez / (B * B))
    rmin = min(A, B)
    phi = max(-eta, -z / rmin, (z - c.extent) / rmin)
    return dict(z=z, x=x, y=y, A=A, B=B, eta=eta, r=r, e=e, rho=rho, rho_min=rmin, phi=phi)


def side_normal(c, values, kappa):
    if values["r"] < 1e-08:
        return values["e"]
    A, B, x, y = (values["A"], values["B"], values["x"], values["y"])
    dA, dB = slopes(c, values["z"], kappa)
    g = add(
        add(mul(c.x_axis, 2 * x / (A * A)), mul(c.z_axis, 2 * y / (B * B))),
        mul(c.axis, -2 * x * x * dA / A**3 - 2 * y * y * dB / B**3),
    )
    return unit(g)


def seg_normal(seg, X, kappa, exp0, exp1):
    """normal field of one segment with two-cap corner rounding; exposures scale the cap contributions."""
    c = seg.c
    values = local(c, X)
    f = kappa * values["rho_min"]
    ns = side_normal(c, values, kappa)
    dside = values["r"] - values["rho"]

    def xc(dcap):
        return S(C(0.5 + (dcap - dside) / (2 * f))) if f > 0 else 1.0 if dcap > dside else 0.0

    H = S(C(values["z"] / c.extent))
    a0 = (1 - H) * exp0 * xc(-values["z"])
    a1 = H * exp1 * xc(values["z"] - c.extent)
    T = S(C(values["r"] / values["rho_min"])) if kappa > 0 else 1.0 if values["r"] > 0 else 0.0
    n = add(add(mul(ns, (1 - a0 - a1) * T), mul(c.axis, -a0)), mul(c.axis, a1))
    return n


def anchor(seg, X):
    c = seg.c
    values = local(c, X)
    zc = min(max(values["z"], 0.0), c.extent)
    A, B = sec(c, zc)
    x, y = (values["x"], values["y"])
    k = math.hypot(x / A, y / B)
    if k > 1:
        x, y = (x / k, y / k)
    return add(add(c.origin, mul(c.axis, zc)), add(mul(c.x_axis, x), mul(c.z_axis, y)))


def inside_closed(c, X):
    values = local(c, X)
    return values["phi"] <= 0


def polynomial_value(coefficients, t):
    value = coefficients[-1]
    magnitude = abs(value)
    for coefficient in reversed(coefficients[:-1]):
        value = value * t + coefficient
        magnitude = magnitude * abs(t) + abs(coefficient)
    if abs(value) > 8 * 2.0**-52 * magnitude:
        return value
    # Compensated Horner preserves signs near cancellation and multiple roots.
    value, correction = coefficients[-1], 0.0
    for coefficient in reversed(coefficients[:-1]):
        product = value * t
        split_value, split_t = 134217729.0 * value, 134217729.0 * t
        high_value = split_value - (split_value - value)
        high_t = split_t - (split_t - t)
        low_value, low_t = value - high_value, t - high_t
        product_error = ((high_value * high_t - product) + high_value * low_t + low_value * high_t) + low_value * low_t
        total = product + coefficient
        tail = total - product
        sum_error = (product - (total - tail)) + (coefficient - tail)
        correction = correction * t + (product_error + sum_error)
        value = total
    return value + correction


def polynomial_roots(coefficients, lo, hi):
    coefficients = list(coefficients)
    while len(coefficients) > 1 and coefficients[-1] == 0.0:
        coefficients.pop()
    if len(coefficients) == 1:
        return []
    critical = polynomial_roots([i * coefficients[i] for i in range(1, len(coefficients))], lo, hi)
    sites = sorted(set([lo, *critical, hi]))
    values = [polynomial_value(coefficients, t) for t in sites]
    for i, t in enumerate(sites):
        if t in critical:
            magnitude = abs(coefficients[-1])
            for coefficient in reversed(coefficients[:-1]):
                magnitude = magnitude * abs(t) + abs(coefficient)
            # A stationary root need not be representable; account for the squared
            # roundoff of its adjacent-double location and compensated evaluation.
            if abs(values[i]) <= 32 * (2.0**-52) ** 2 * magnitude:
                values[i] = 0.0
    roots = [t for t, value in zip(sites, values) if value == 0.0]
    for a, b, fa, fb in zip(sites, sites[1:], values, values[1:]):
        if fa == 0.0 or fb == 0.0 or (fa < 0.0) == (fb < 0.0):
            continue
        while True:
            mid = a + (b - a) * 0.5
            if mid == a or mid == b:
                roots.append(b)
                break
            fm = polynomial_value(coefficients, mid)
            if fm == 0.0:
                roots.append(mid)
                break
            if (fm < 0.0) == (fa < 0.0):
                a, fa = mid, fm
            else:
                b = mid
    return sorted(set(roots))


def polynomial_intervals(coefficients, lo, hi):
    roots = polynomial_roots(coefficients, lo, hi)
    sites = sorted(set([lo, *roots, hi]))
    out = [(t, t) for t in sites if t in roots or polynomial_value(coefficients, t) <= 0.0]
    for a, b in zip(sites, sites[1:]):
        if polynomial_value(coefficients, a + (b - a) * 0.5) <= 0.0:
            out.append((a, b))
    return sorted(out)


def ray_intervals(c, p, n, tmax):
    offset = sub(p, c.origin)
    z, dz = dot(offset, c.axis), dot(n, c.axis)
    x, dx = dot(offset, c.x_axis), dot(n, c.x_axis)
    y, dy = dot(offset, c.z_axis), dot(n, c.z_axis)
    out = []
    for station, (first, second) in enumerate(zip(c.stations, c.stations[1:])):
        lower, upper = max(0.0, first[0]), min(c.extent, second[0])
        lo, hi = 0.0, tmax
        if lower > upper:
            continue
        if dz == 0.0:
            if z < lower or z > upper:
                continue
        else:
            a, b = (lower - z) / dz, (upper - z) / dz
            lo, hi = max(lo, min(a, b)), min(hi, max(a, b))
            if lo > hi:
                continue
        slopes = [(second[i] - first[i]) / (second[0] - first[0]) for i in (1, 2)]
        a, b = [first[i + 1] + slopes[i] * (z - first[0]) for i in range(2)]
        da, db = slopes[0] * dz, slopes[1] * dz

        def square(v, dv):
            return [v * v, 2.0 * v * dv, dv * dv]

        xx, yy, aa, bb = square(x, dx), square(y, dy), square(a, da), square(b, db)
        coefficients = [0.0] * 5
        for i in range(3):
            for j in range(3):
                coefficients[i + j] += xx[i] * bb[j] + yy[i] * aa[j] - aa[i] * bb[j]
        for start, end in polynomial_intervals(coefficients, lo, hi):
            out.append((start, end, station))
    return out


def union_exit(segs, p, n):
    tmax = 2 * max((math.hypot(s.c.extent, max((max(st[1], st[2]) for st in s.c.stations))) for s in segs)) + (
        2 * length(sub(segs[1].c.origin, segs[0].c.origin)) if len(segs) > 1 else 0.0
    )
    ivs = []
    for s in segs:
        for a, b, station in ray_intervals(s.c, p, n, tmax):
            ivs.append((a, s.id, station, b, s))
    ivs.sort(key=lambda x: (x[0], x[1], x[2]))
    tol = 1e-07 * tmax
    t = 0.0
    who = next((s for a, _, _, b, s in ivs if a <= tol and b >= 0.0), None)
    # Radial classification and polynomial expansion can round to opposite sides
    # at the origin; an inside classification still owns the closed point [0, 0].
    if who is None:
        who = next((s for s in segs if inside_closed(s.c, p)), None)
    changed = True
    while changed:
        changed = False
        for a, _, _, b, s in ivs:
            if a <= t + tol and b > t:
                t = b
                who = s
                changed = True
    return (t, who)


def exposures(segs, L, kappa, shift=0.0):
    """exp_end[0] (thigh knee face), exp_start[1] (calf knee face); free ends = 1. shift = phi_min for the outside evaluation."""
    e0 = [1.0] * len(segs)
    e1 = [1.0] * len(segs)
    if len(segs) == 2:
        p1 = L[1]["phi"] - shift
        p0 = L[0]["phi"] - shift
        e1[0] = S(C(p1 / kappa)) if kappa > 0 else 1.0 if p1 > 0 else 0.0
        e0[1] = S(C(p0 / kappa)) if kappa > 0 else 1.0 if p0 > 0 else 0.0
    return (e0, e1)


def leg_normal_and_anchor(segs, X, kappa):
    L = [local(s.c, X) for s in segs]
    e0, e1 = exposures(segs, L, kappa)
    phimin = min((values["phi"] for values in L))
    rho_w = [
        S(C(1 - (values["phi"] - phimin) / kappa)) if kappa > 0 else 1.0 if values["phi"] == phimin else 0.0
        for values in L
    ]
    nsum = ZERO
    dmin = math.inf
    for k, (s, values) in enumerate(zip(segs, L)):
        nk = seg_normal(s, X, kappa, e0[k], e1[k])
        ck = anchor(s, X)
        nsum = add(nsum, mul(nk, rho_w[k]))
        dmin = min(dmin, length(sub(ck, X)))
    mag = length(nsum)
    return (mul(nsum, 1 / mag) if mag > 1e-15 else (0.0, 0.0, 0.0), -dmin, S(C(mag / 0.25)))


def direction(segs, L, e0, e1, p, rp, waist, wa, kappa, use_cont, shift=0.0):
    dirs = []
    hk = []
    cont = []
    for k, (s, values) in enumerate(zip(segs, L)):
        c = s.c
        b0 = rest_contact(rp, s.c0, waist, wa)
        b, _ = transported_side(b0, c, s.c0)
        alpha = S(C(values["eta"])) if b is not None else 0.0
        side = add(mul(values["e"], 1 - alpha), mul(b or ZERO, alpha))
        gcap = values["rho_min"] * max(values["eta"], 0.0)
        zp = max(values["z"], 0.0)
        hp = max(c.extent - values["z"], 0.0)
        g0 = gcap * e0[k]
        g1 = gcap * e1[k]
        t0 = S(C(zp / (zp + kappa * g0))) if kappa > 0 and g0 > 0 else 1.0
        t1 = S(C(hp / (hp + kappa * g1))) if kappa > 0 and g1 > 0 else 1.0
        dirs.append(add(mul(side, t0 * t1), mul(c.axis, 1 - t1 - (1 - t0))))
        hk.append(c.extent - values["z"] if k == 0 else values["z"])
        ph = values["phi"] - shift
        cont.append(1 - (S(C(ph / kappa)) if kappa > 0 else 1.0 if ph > 0 else 0.0) if use_cont else 1.0)
    if len(segs) == 2:
        beta = kappa * min(sec(segs[0].c, segs[0].c.extent))
        wc = S(C(0.5 + (hk[1] - hk[0]) / (2 * beta))) if beta > 0 else 1.0 if hk[1] > hk[0] else 0.0
        return add(mul(dirs[0], cont[0] * (1 - wc)), mul(dirs[1], cont[1] * wc))
    return dirs[0]


def leg_constraints(segs, p, rp, waist, wa, kappa):
    L = [local(s.c, p) for s in segs]
    e0, e1 = exposures(segs, L, kappa)
    inside = any((values["phi"] <= 0 for values in L))
    if not inside:
        nu, d, sn = leg_normal_and_anchor(segs, p, kappa)
        if nu == ZERO:
            return []
        phimin = min((values["phi"] for values in L))
        e0s, e1s = exposures(segs, L, kappa, phimin)
        m_out = length(direction(segs, L, e0s, e1s, p, rp, waist, wa, kappa, True, phimin))
        return [(nu, d, S(C(m_out / 0.25)) * sn)]
    n_raw = direction(segs, L, e0, e1, p, rp, waist, wa, kappa, True)
    m = length(n_raw)
    if m <= 1e-12:
        return []
    n = mul(n_raw, 1 / m)
    d, who = union_exit(segs, p, n)
    if who is None:
        return []
    q = add(p, mul(n, d))
    nu, _, sn = leg_normal_and_anchor(segs, q, kappa)
    if dot(nu, n) <= 1e-06:
        return []
    return [(nu, dot(nu, sub(q, p)), S(C(m / 0.25)) * sn)]


def make_segs(current, rest, long_leg=True):
    """Cylinder maps use (leg, segment) keys: left/right and thigh/calf are 0/1."""
    return [
        [
            Seg(current[(leg, segment)], rest[(leg, segment)], segment, segment == 0, segment == 1)
            for segment in range(2 if long_leg else 1)
            if (leg, segment) in current and (leg, segment) in rest
        ]
        for leg in range(2)
    ]


def all_cons(legs, p, rp, waist, wa, kappa):
    return [
        constraint for segments in legs if segments for constraint in leg_constraints(segments, p, rp, waist, wa, kappa)
    ]


def correction(cons, x):
    shifted = [(n, d - dot(n, x), s) for n, d, s in cons if s != 0]
    if all(d <= 0 for _, d, _ in shifted):
        return ZERO
    return project_w(shifted)


def local_radius(rest_flat, rp):
    exponents, radii = [], []
    for cylinder in rest_flat:
        values = local(cylinder, rp)
        zc = min(max(values["z"], 0.0), cylinder.extent)
        radius = min(cylinder.section_radii(zc))
        dist = math.hypot(values["r"], values["z"] - zc)
        ratio = dist / radius
        exponents.append(-ratio * ratio)
        radii.append(radius)
    if not radii:
        return 0.0
    maximum = max(exponents)
    weights = [math.exp(a - maximum) for a in exponents]
    return sum(w * r for w, r in zip(weights, radii)) / sum(weights)


def point_weights(points, rest_points, q, waist=ZERO, waist_axis=(0, 0, 1)):
    axis = unit(waist_axis)
    result = []
    for p, rp, weight in zip(points, rest_points, q):
        valid = finite(p) and finite(rp)
        weight = min(1.0, max(0.0, weight)) if math.isfinite(weight) else 0.0
        if not valid or length(radial(sub(rp, waist), axis)) < 1e-8:
            weight = 0.0
        result.append(weight)
    return result


def couple(points, rest_points, cons, adj, kappa, rest_flat, q, iters=30):
    count = len(points)
    valid = [finite(p) and finite(rp) for p, rp in zip(points, rest_points)]
    weights = [min(1.0, max(0.0, w)) if math.isfinite(w) and valid[i] else 0.0 for i, w in enumerate(q)]
    radii = [local_radius(rest_flat, rp) if valid[i] else 0.0 for i, rp in enumerate(rest_points)]
    neighbours = [set() for _ in points]
    for i in range(count):
        for j in adj[i]:
            if 0 <= j < count and valid[i] and valid[j]:
                neighbours[i].add(j)
                neighbours[j].add(i)
    edges = [[] for _ in points]
    for i in range(count):
        for j in sorted(neighbours[i]):
            scale = kappa * (radii[i] + radii[j]) / 2.0
            offset = sub(rest_points[i], rest_points[j])
            value = scale * scale / max(dot(offset, offset), 1e-12) if kappa > 0 and rest_flat else 0.0
            edges[i].append((j, value))
    u, delta = [ZERO] * count, [ZERO] * count
    for _ in range(iters):
        for i in range(count):
            if weights[i] == 0:
                continue
            total, numerator = 0.0, ZERO
            for j, value in edges[i]:
                total += value
                numerator = add(numerator, mul(u[j], value))
            x = mul(numerator, 1.0 / (1.0 + total))
            delta[i] = add(x, correction(cons[i], x))
            u[i] = mul(delta[i], weights[i])
    for i in range(count):
        if weights[i] > 0:
            delta[i] = add(delta[i], correction(cons[i], delta[i]))
            u[i] = mul(delta[i], weights[i])
    return u


def surface_adjacency(nu, nv, du, dv, periodic_u, periodic_v):
    unique_u, unique_v = nu - du if periodic_u else nu, nv - dv if periodic_v else nv
    result = {}
    for u in range(unique_u):
        for v in range(unique_v):
            neighbours = set()
            for a, b in ((u - 1, v), (u + 1, v), (u, v - 1), (u, v + 1)):
                if periodic_u:
                    a %= unique_u
                if periodic_v:
                    b %= unique_v
                if 0 <= a < unique_u and 0 <= b < unique_v:
                    neighbours.add(a * nv + b)
            result[u * nv + v] = sorted(neighbours)
    return result


def spline_basis(knots, degree, count, parameter):
    """Evaluate Maya's knot convention, including the exact upper endpoint."""
    full = [knots[0]] + list(knots) + [knots[-1]]
    span = count - 1
    if parameter < full[count]:
        span = next(i for i in range(degree, count) if full[i] <= parameter < full[i + 1])
    values = [1.0] + [0.0] * degree
    for j in range(1, degree + 1):
        saved = 0.0
        for r in range(j):
            right = full[span + r + 1] - parameter
            left = parameter - full[span + 1 - j + r]
            temp = values[r] / (right + left) if right + left else 0.0
            values[r] = saved + right * temp
            saved = left * temp
        values[j] = saved
    return {span - degree + j: value for j, value in enumerate(values) if value != 0.0}


def surface_basis(surface, with_quadrature=False):
    nu, nv, du, dv = (surface[k] for k in ("nu", "nv", "du", "dv"))
    pu, pv = surface["formU"] == 3, surface["formV"] == 3
    uu, vv = nu - du if pu else nu, nv - dv if pv else nv
    weights = surface.get("weights", [1.0] * (nu * nv))
    if len(weights) != nu * nv or any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Surface weights must be finite and positive")

    def sites(knots, degree, count, divisions):
        full = [knots[0]] + list(knots) + [knots[-1]]
        return [
            full[i] + (full[i + 1] - full[i]) * (k - 0.5) / divisions
            for i in range(degree, count)
            if full[i + 1] > full[i]
            for k in range(1, divisions + 1)
        ]

    us = sites(surface["ku"], du, nu, 2)
    full_v = [surface["kv"][0]] + list(surface["kv"]) + [surface["kv"][-1]]
    offset = 0.5 * math.sqrt(3.0 / 5.0)
    vs = [
        (full_v[i] + (full_v[i + 1] - full_v[i]) * fraction, weight)
        for i in range(dv, nv)
        if full_v[i + 1] > full_v[i]
        for fraction, weight in ((0.5 - offset, 10.0 / 9.0), (0.5, 16.0 / 9.0), (0.5 + offset, 10.0 / 9.0))
    ]
    parameters = [(u, v) for u in us for v, _ in vs]
    if not pu:
        parameters += [
            (u, v)
            for u in (surface["ku"][du - 1], surface["ku"][nu - 1])
            for v, _ in vs
        ]
    quadrature = [weight for _ in us for _, weight in vs]
    if not pu:
        quadrature += [weight for _ in range(2) for _, weight in vs]
    rows = []
    for u, v in parameters:
        bu = spline_basis(surface["ku"], du, nu, u)
        bv = spline_basis(surface["kv"], dv, nv, v)
        row = {}
        total = 0.0
        for a, x in bu.items():
            for b, y in bv.items():
                value = x * y * weights[a * nv + b]
                index = (a % uu) * nv + b % vv
                row[index] = row.get(index, 0.0) + value
                total += value
        if not math.isfinite(total) or total <= 0:
            raise ValueError("Invalid surface basis")
        row = {i: w / total for i, w in sorted(row.items())}
        if any(not math.isfinite(w) or w < 0 for w in row.values()) or abs(sum(row.values()) - 1) > 1e-9:
            raise ValueError("Invalid surface basis")
        rows.append(row)
    return (parameters, rows, quadrature) if with_quadrature else (parameters, rows)


def surface_groups(surface):
    nu, nv, du, dv = (surface[k] for k in ("nu", "nv", "du", "dv"))
    adjacency = surface_adjacency(nu, nv, du, dv, surface["formU"] == 3, surface["formV"] == 3)
    parent = {i: i for i in adjacency}

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def join(i, j):
        a, b = sorted((root(i), root(j)))
        parent[b] = a

    for i in adjacency:
        u, v = divmod(i, nv)
        if surface.get("closedU", False) and surface["formU"] != 3 and u == 0:
            join(i, (nu - 1) * nv + v)
    groups = {}
    for i in adjacency:
        groups.setdefault(root(i), []).append(i)
    return list(groups.values()), adjacency


def solve_group(mass, neighbours, terms, initial=ZERO, trace=None):
    """Terms are (normal, target, weight, basis); neighbours are (weight, displacement)."""
    diagonal = mass + sum(w for w, _ in neighbours)
    rhs0 = tuple(sum(w * v[k] for w, v in neighbours) for k in range(3))

    def objective(z):
        return (
            mass * dot(z, z)
            + sum(w * dot(sub(z, v), sub(z, v)) for w, v in neighbours)
            + sum(w * max(0.0, t - b * dot(n, z)) ** 2 for n, t, w, b in terms)
        )

    z = initial
    for _ in range(50):
        active = [t - b * dot(n, z) > 0 for n, t, w, b in terms]
        matrix = [[diagonal if r == c else 0.0 for c in range(3)] for r in range(3)]
        rhs = list(rhs0)
        for enabled, (n, t, w, b) in zip(active, terms):
            if enabled:
                for r in range(3):
                    rhs[r] += w * b * t * n[r]
                    for c in range(3):
                        matrix[r][c] += w * b * b * n[r] * n[c]
        det = determinant(matrix)
        candidate = tuple(
            determinant([[rhs[r] if c == k else matrix[r][c] for c in range(3)] for r in range(3)]) / det
            for k in range(3)
        )
        consistent = all(enabled == (t - b * dot(n, candidate) > 0) for enabled, (n, t, w, b) in zip(active, terms))
        old_value = objective(z)
        if objective(candidate) >= old_value:
            step = sub(candidate, z)
            candidate = z
            for k in range(1, 31):
                trial = add(z, mul(step, 2.0**-k))
                if objective(trial) < old_value:
                    candidate = trial
                    break
            consistent = all(enabled == (t - b * dot(n, candidate) > 0) for enabled, (n, t, w, b) in zip(active, terms))
        small = length(sub(candidate, z)) <= 1e-12 * (1 + length(z))
        z = candidate
        if trace is not None:
            trace.append((z, objective(z)))
        if consistent or small:
            break
    return z


def _resolve_schedule(iters, regenerations):
    if iters is None and regenerations is None:
        iters = 3
        regeneration_sweeps = ()
    else:
        iters = 10 if iters is None else iters
        regeneration_count = 1 if regenerations is None else regenerations
        regeneration_sweeps = tuple(
            j * iters // (regeneration_count + 1) for j in range(1, regeneration_count + 1)
        )
    return iters, regeneration_sweeps


def sample_couple(points, rest_points, rows, groups, edges, q, generate, iters=None, events=None, regenerations=None, sample_weights=None):
    iters, schedule = _resolve_schedule(iters, regenerations)
    if sample_weights is None:
        sample_weights = [1.0] * len(rows)
    group_of = {i: g for g, members in enumerate(groups) for i in members}
    beta = []
    attachments = [[] for _ in groups]
    for s, row in enumerate(rows):
        grouped = {}
        for i, b in row.items():
            g = group_of[i]
            grouped[g] = grouped.get(g, 0.0) + b
        beta.append(sorted(grouped.items()))
        for g, b in beta[-1]:
            attachments[g].append((s, b))
    denominators = [sum(b * b for b in row.values()) for row in rows]

    def evaluate(row, values):
        return tuple(sum(b * values[i][k] for i, b in row.items()) for k in range(3))

    ps, rs = [evaluate(row, points) for row in rows], [evaluate(row, rest_points) for row in rows]
    paint = [min(q[i] for i in members) for members in groups]
    v = [ZERO] * len(groups)
    schedule = set(schedule)
    constraints = []
    for sweep in range(iters):
        if sweep == 0 or sweep in schedule:
            constraints = []
            for s in range(len(rows)):
                displacement = tuple(sum(b * v[g][k] for g, b in beta[s]) for k in range(3))
                position = add(ps[s], displacement)
                if events is not None:
                    events.append((sweep, s, position))
                constraints.append(
                    [
                        [n, d + dot(n, displacement), 4 * sample_weights[s] * strength / denominators[s], dot(n, displacement)]
                        for n, d, strength in generate(position, rs[s])
                        if strength != 0
                    ]
                )
        for g, members in enumerate(groups):
            if paint[g] == 0:
                continue
            terms = [
                (n, 1.25 * d - (cache - b * dot(n, v[g])), w, b)
                for s, b in attachments[g]
                for n, d, w, cache in constraints[s]
            ]
            z = solve_group(len(members), [(w, v[h]) for h, w in edges[g]], terms, v[g])
            updated = mul(z, paint[g])
            change = sub(updated, v[g])
            v[g] = updated
            for s, b in attachments[g]:
                for constraint in constraints[s]:
                    constraint[3] += b * dot(constraint[0], change)
    return {i: v[g] for i, g in group_of.items()}


def surface_displacements(
    surface,
    points,
    rest_points,
    legs,
    rest_flat,
    waist,
    waist_axis,
    kappa,
    q=None,
    members=None,
    iters=None,
    events=None,
    regenerations=None,
):
    groups, adjacency = surface_groups(surface)
    try:
        _, rows, sample_weights = surface_basis(surface, with_quadrature=True)
    except ValueError as error:
        warnings.warn(str(error), RuntimeWarning, stacklevel=2)
        return {i: ZERO for group in groups for i in group}
    q = point_weights(points, rest_points, q if q is not None else [1.0] * len(points), waist, waist_axis)
    if members is not None:
        q = [w if i in members else 0.0 for i, w in enumerate(q)]
    group_of = {i: g for g, group in enumerate(groups) for i in group}
    radii = {
        i: local_radius(rest_flat, rest_points[i]) if finite(rest_points[i]) and finite(points[i]) else 0.0
        for i in adjacency
    }
    edges = [{} for _ in groups]
    for i, neighbours in adjacency.items():
        for j in neighbours:
            if not all(finite(values[k]) for values in (points, rest_points) for k in (i, j)):
                continue
            g, h = group_of[i], group_of[j]
            if g == h:
                continue
            scale = kappa * (radii[i] + radii[j]) / 2
            offset = sub(rest_points[i], rest_points[j])
            w = scale * scale / max(dot(offset, offset), 1e-12) if kappa > 0 and rest_flat else 0.0
            edges[g][h] = edges[g].get(h, 0.0) + w

    def generate(p, r):
        if not finite(p) or not finite(r) or length(radial(sub(r, waist), unit(waist_axis))) < 1e-8:
            return []
        return all_cons(legs, p, r, waist, waist_axis, kappa)

    return sample_couple(
        points,
        rest_points,
        rows,
        groups,
        [sorted(e.items()) for e in edges],
        q,
        generate,
        iters,
        events,
        regenerations,
        sample_weights=sample_weights,
    )


def fixture_cylinders(data):
    """Construct object-space cylinders from the serialized node attributes."""
    result = []
    scale = data["ringScale"]
    names = [station + "Radius" + axis for station in ("thigh", "knee", "calf", "ankle") for axis in ("X", "Z")]
    profile = dict(zip(names, data["radii"]))
    profile.update({k: data[k] for k in ("thighPosition", "calfPosition")})
    for prefix in ("", "rest"):
        cylinders = {}
        for leg, side in enumerate(("left", "right")):
            stem = prefix + side.capitalize() if prefix else side
            matrices = [data["matrices"][stem + joint + "Matrix"] for joint in ("Hip", "Knee", "Heel")]
            positions = [tuple(m[12:15]) for m in matrices]
            lengths = [length(sub(b, a)) for a, b in zip(positions, positions[1:])]
            total = sum(lengths)
            knee = lengths[0] / total
            stations = profile_stations(*lengths, **profile)
            axis_index = data[side + "RingAxis"]
            for segment in range(1 + data["skirtType"]):
                joint = matrices[segment]
                rows = [tuple(joint[i : i + 3]) for i in (0, 4, 8)]
                axis = unit(sub(positions[segment + 1], positions[segment]))
                y = mul(unit(rows[axis_index % 3]), -1 if axis_index >= 3 else 1)
                x = rotate_minimum(unit(rows[(axis_index + 1) % 3]), y, axis)
                x = unit(radial(x, axis))
                z = unit(cross(axis, x))
                start, end = (0, knee) if segment == 0 else (knee, 1)
                values = [
                    ((t - start) * total * scale[1], a * scale[0], b * scale[2])
                    for t, a, b in stations
                    if start <= t <= end
                ]
                cylinders[leg, segment] = Cylinder(
                    positions[segment], axis, x, z, lengths[segment] * scale[1], values, rows
                )
        result.append(cylinders)
    return result


def fixture_output(data, closed=True, iters=None, regenerations=None):
    surface = dict(data["input"])
    surface["closedU"] = closed
    current, rest = fixture_cylinders(data)
    bell = data["matrices"]["restBellMatrix"]
    points = [tuple(p[:3]) for p in surface["cvs"]]
    rps = [tuple(p[:3]) for p in data["rest"]["cvs"]]
    u = surface_displacements(
        surface,
        points,
        rps,
        make_segs(current, rest, bool(data["skirtType"])),
        [rest[k] for k in sorted(rest)],
        tuple(bell[12:15]),
        tuple(bell[4:7]),
        data["falloff"],
        iters=iters,
        regenerations=regenerations,
    )
    return [add(p, mul(u[i], data["envelope"])) for i, p in enumerate(points)]
