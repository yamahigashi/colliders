#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace NurbsRefit
{
inline int findSpan(int nCV, int p, double u, const std::vector<double>& T)
{
    if (u >= T[nCV])
    {
        int span = nCV - 1;
        while (span > p && T[span] == T[span + 1])
            --span;
        return span;
    }

    int low = p;
    int high = nCV;
    while (low + 1 < high)
    {
        const int middle = low + (high - low) / 2;
        if (u < T[middle])
            high = middle;
        else
            low = middle;
    }
    return low;
}

inline void basisFunctions(int span, double u, int p, const std::vector<double>& T,
    std::vector<double>& out)
{
    out.assign(p + 1, 0.0);
    std::vector<double> left(p + 1), right(p + 1);
    out[0] = 1.0;
    for (int j = 1; j <= p; ++j)
    {
        left[j] = u - T[span + 1 - j];
        right[j] = T[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r)
        {
            const double term = out[r] / (right[r + 1] + left[j - r]);
            out[r] = saved + right[r + 1] * term;
            saved = left[j - r] * term;
        }
        out[j] = saved;
    }
}

inline std::vector<double> fullKnots(const std::vector<double>& mayaKnots)
{
    std::vector<double> result;
    result.reserve(mayaKnots.size() + 2);
    result.push_back(mayaKnots.front());
    result.insert(result.end(), mayaKnots.begin(), mayaKnots.end());
    result.push_back(mayaKnots.back());
    return result;
}

inline std::vector<double> clampedUniformKnotsMaya(int spans, int degree)
{
    std::vector<double> result(degree, 0.0);
    for (int j = 1; j < spans; ++j)
        result.push_back(static_cast<double>(j) / spans);
    result.insert(result.end(), degree, 1.0);
    return result;
}

inline std::vector<double> grevilleAbscissae(const std::vector<double>& Tfull,
    int degree, int nCV)
{
    std::vector<double> result(nCV, 0.0);
    for (int k = 0; k < nCV; ++k)
    {
        for (int j = 1; j <= degree; ++j)
            result[k] += Tfull[k + j];
        result[k] /= degree;
    }
    return result;
}

inline bool sampleMatrix(int nv, int q, const std::vector<double>& KvMaya,
    int m, std::vector<double>& S)
{
    const std::vector<double> inputKnots = fullKnots(KvMaya);
    const double vmin = inputKnots[q];
    const double length = inputKnots[nv] - vmin;
    if (length <= 1e-12)
        return false;

    const std::vector<double> outputKnots = fullKnots(clampedUniformKnotsMaya(m - 3, 3));
    const std::vector<double> greville = grevilleAbscissae(outputKnots, 3, m);
    S.assign(m * nv, 0.0);
    std::vector<double> basis;
    for (int k = 0; k < m; ++k)
    {
        const double v = vmin + greville[k] * length;
        const int span = findSpan(nv, q, v, inputKnots);
        basisFunctions(span, v, q, inputKnots, basis);
        for (int j = 0; j <= q; ++j)
            S[k * nv + span - q + j] = basis[j];
    }
    return true;
}

inline void refitColumns(int nu, int nv, const std::vector<double>& P, int m,
    const std::vector<double>& S, std::vector<double>& C)
{
    C.assign(nu * m * 3, 0.0);
    for (int u = 0; u < nu; ++u)
        for (int k = 0; k < m; ++k)
            for (int j = 0; j < nv; ++j)
                for (int c = 0; c < 3; ++c)
                    C[(u * m + k) * 3 + c] += S[k * nv + j] * P[(u * nv + j) * 3 + c];
}

inline std::vector<double> protectedBoundaries(const std::vector<double>& protectedV)
{
    std::vector<double> bounds;
    bounds.reserve(protectedV.size() + 2);
    bounds.push_back(0.0);
    bounds.insert(bounds.end(), protectedV.begin(), protectedV.end());
    bounds.push_back(1.0);
    return bounds;
}

inline int protectedBandCount(int protectedCount)
{
    return protectedCount + 1;
}

inline int protectedCvCount(int spansV, int protectedCount)
{
    return spansV + 2 * protectedBandCount(protectedCount) + 1;
}

inline std::vector<int> allocateBandSpans(int spansV, const std::vector<double>& bounds)
{
    const int bands = static_cast<int>(bounds.size()) - 1;
    std::vector<int> spans(bands, 1);
    if (bands <= 0)
        return spans;

    double sumLength = 0.0;
    std::vector<double> lengths(bands, 0.0);
    for (int band = 0; band < bands; ++band)
    {
        lengths[band] = bounds[band + 1] - bounds[band];
        sumLength += lengths[band];
    }

    const int remaining = spansV - bands;
    if (remaining <= 0 || sumLength <= 0.0)
        return spans;

    std::vector<double> fractional(bands, 0.0);
    int used = 0;
    for (int band = 0; band < bands; ++band)
    {
        const double share = remaining * lengths[band] / sumLength;
        const double whole = std::floor(share);
        spans[band] = 1 + static_cast<int>(whole);
        fractional[band] = share - whole;
        used += static_cast<int>(whole);
    }

    const int leftover = remaining - used;
    std::vector<int> order(bands);
    for (int band = 0; band < bands; ++band)
        order[band] = band;
    std::sort(order.begin(), order.end(), [&](int left, int right) {
        return fractional[left] > fractional[right];
    });
    // Anchor each tolerance group to its maximum: approximate equality is not transitive.
    for (int first = 0; first < bands;)
    {
        int end = first + 1;
        const double maximum = fractional[order[first]];
        while (end < bands && maximum - fractional[order[end]] <= 1e-12 * maximum)
            ++end;
        std::sort(order.begin() + first, order.begin() + end);
        first = end;
    }
    for (int i = 0; i < leftover; ++i)
        ++spans[order[i]];
    return spans;
}

inline std::vector<double> protectedKnotsMaya(int spansV, const std::vector<double>& protectedV)
{
    const std::vector<double> bounds = protectedBoundaries(protectedV);
    const std::vector<int> spans = allocateBandSpans(spansV, bounds);
    const int bands = static_cast<int>(spans.size());
    std::vector<double> knots(3, 0.0);
    for (int band = 0; band < bands; ++band)
    {
        const double start = bounds[band];
        const double length = bounds[band + 1] - start;
        for (int j = 1; j < spans[band]; ++j)
            knots.push_back(start + length * static_cast<double>(j) / spans[band]);
        if (band + 1 < bands)
            knots.insert(knots.end(), 3, bounds[band + 1]);
    }
    knots.insert(knots.end(), 3, 1.0);
    return knots;
}

inline std::vector<double> normalizedInputRows(int nv, int q, const std::vector<double>& KvMaya)
{
    const std::vector<double> inputKnots = fullKnots(KvMaya);
    const double vmin = inputKnots[q];
    const double length = inputKnots[nv] - vmin;
    std::vector<double> rows = grevilleAbscissae(inputKnots, q, nv);
    for (int j = 0; j < nv; ++j)
        rows[j] = (rows[j] - vmin) / length;
    return rows;
}

inline int protectedFitCondition(int nv, int q, const std::vector<double>& KvMaya,
    int spansV, const std::vector<double>& protectedV)
{
    for (size_t i = 0; i < protectedV.size(); ++i)
    {
        if (!std::isfinite(protectedV[i]) || protectedV[i] <= 0.0 || protectedV[i] >= 1.0)
            return 10;
        if (i > 0 && !(protectedV[i] > protectedV[i - 1]))
            return 10;
    }

    const std::vector<double> inputKnots = fullKnots(KvMaya);
    const double vmin = inputKnots[q];
    const double domainLength = inputKnots[nv] - vmin;
    double previous = 0.0;
    for (double boundary : protectedV)
    {
        if ((boundary - previous) * domainLength <= 1e-12)
            return 11;
        previous = boundary;
    }
    if ((1.0 - previous) * domainLength <= 1e-12)
        return 11;

    if (spansV <= 0 || protectedV.size() >= static_cast<size_t>(spansV))
        return 12;
    if (q != 1)
        return 13;

    const std::vector<double> rows = normalizedInputRows(nv, q, KvMaya);
    std::vector<bool> matched(nv, false);
    for (size_t i = 0; i < protectedV.size(); ++i)
    {
        int matches = 0;
        for (int j = 0; j < nv; ++j)
        {
            if (std::fabs(rows[j] - protectedV[i]) <= 1e-9)
            {
                if (matched[j])
                    return 14;
                matched[j] = true;
                ++matches;
            }
        }
        if (matches != 1)
            return 14;
    }
    return 0;
}

inline bool sampleMatrix(int nv, int q, const std::vector<double>& KvMaya,
    int spansV, const std::vector<double>& protectedV, std::vector<double>& S)
{
    if (protectedFitCondition(nv, q, KvMaya, spansV, protectedV) != 0)
        return false;

    const std::vector<double> bounds = protectedBoundaries(protectedV);
    const std::vector<int> spans = allocateBandSpans(spansV, bounds);
    const int bands = static_cast<int>(spans.size());
    const int m = protectedCvCount(spansV, static_cast<int>(protectedV.size()));
    S.assign(static_cast<size_t>(m) * nv, 0.0);

    std::vector<double> rowV = normalizedInputRows(nv, q, KvMaya);
    for (size_t i = 0; i < protectedV.size(); ++i)
    {
        for (int j = 0; j < nv; ++j)
        {
            if (std::fabs(rowV[j] - protectedV[i]) <= 1e-9)
            {
                rowV[j] = protectedV[i];
                break;
            }
        }
    }

    int outRow = 0;
    for (int band = 0; band < bands; ++band)
    {
        const int localCount = spans[band] + 3;
        const std::vector<double> localFull = fullKnots(clampedUniformKnotsMaya(spans[band], 3));
        const std::vector<double> greville = grevilleAbscissae(localFull, 3, localCount);
        const int localStart = band == 0 ? 0 : 1;
        for (int local = localStart; local < localCount; ++local, ++outRow)
        {
            const double sample = local == localCount - 1 ? bounds[band + 1]
                : bounds[band] + (bounds[band + 1] - bounds[band]) * greville[local];
            double* coefficients = &S[static_cast<size_t>(outRow) * nv];
            if (outRow == 0)
            {
                coefficients[0] = 1.0;
                continue;
            }
            if (outRow == m - 1)
            {
                coefficients[nv - 1] = 1.0;
                continue;
            }

            std::vector<int> inBand;
            inBand.reserve(nv);
            for (int j = 0; j < nv; ++j)
            {
                if (rowV[j] >= bounds[band] && rowV[j] <= bounds[band + 1])
                    inBand.push_back(j);
            }

            bool exact = false;
            for (int index = 0; index < static_cast<int>(inBand.size()); ++index)
            {
                const int row = inBand[index];
                if (std::fabs(rowV[row] - sample) <= 1e-12)
                {
                    coefficients[row] = 1.0;
                    exact = true;
                    break;
                }
            }
            if (exact)
                continue;

            bool bracketed = false;
            for (int index = 0; index + 1 < static_cast<int>(inBand.size()); ++index)
            {
                const int left = inBand[index];
                const int right = inBand[index + 1];
                if (rowV[left] <= sample && sample <= rowV[right])
                {
                    const double width = rowV[right] - rowV[left];
                    if (width == 0.0)
                        return false;
                    const double lambda = (sample - rowV[left]) / width;
                    coefficients[left] = 1.0 - lambda;
                    coefficients[right] = lambda;
                    bracketed = true;
                    break;
                }
            }
            if (!bracketed)
                return false;
        }
    }
    return outRow == m;
}
}
