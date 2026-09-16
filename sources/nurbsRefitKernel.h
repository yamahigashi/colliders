#pragma once

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
}
