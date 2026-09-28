#include "nurbsRefitKernel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const char* message)
{
    require(std::isfinite(actual) && std::fabs(actual - expected) <= tolerance, message);
}

bool sameBits(double left, double right)
{
    return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void partition(const std::vector<double>& basis)
{
    double sum = 0.0;
    for (double value : basis)
    {
        require(std::isfinite(value) && value >= 0.0, "basis must be finite and nonnegative");
        sum += value;
    }
    close(sum, 1.0, 1e-12, "basis must sum to one");
}

double deBoor(const std::vector<double>& knots, int degree, int nCV,
    const std::vector<double>& points, int column, int component, double parameter)
{
    int span = nCV - 1;
    for (int i = degree; i < nCV; ++i)
    {
        if (knots[i] <= parameter && parameter < knots[i + 1])
        {
            span = i;
            break;
        }
    }
    std::vector<double> values(degree + 1);
    for (int j = 0; j <= degree; ++j)
        values[j] = points[(column * nCV + span - degree + j) * 3 + component];
    for (int r = 1; r <= degree; ++r)
    {
        for (int j = degree; j >= r; --j)
        {
            const int i = span - degree + j;
            const double alpha = (parameter - knots[i]) / (knots[i + degree - r + 1] - knots[i]);
            values[j] = (1.0 - alpha) * values[j - 1] + alpha * values[j];
        }
    }
    return values[degree];
}

void testBasisFunctions()
{
    const std::vector<double> linearMaya = {0.0, 2.5, 4.0, 7.5};
    const std::vector<double> linear = NurbsRefit::fullKnots(linearMaya);
    const std::vector<double> cubic = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(4, 3));
    require(linear == std::vector<double>({0.0, 0.0, 2.5, 4.0, 7.5, 7.5}), "full knot endpoint duplication");
    for (int sample = 0; sample < 50; ++sample)
    {
        const double fraction = sample / 49.0;
        std::vector<double> basis;
        NurbsRefit::basisFunctions(NurbsRefit::findSpan(7, 3, fraction, cubic), fraction, 3, cubic, basis);
        require(basis.size() == 4, "cubic basis size");
        partition(basis);

        const double v = 7.5 * fraction;
        const int segment = v < 2.5 ? 0 : (v < 4.0 ? 1 : 2);
        const double alpha = (v - linearMaya[segment]) / (linearMaya[segment + 1] - linearMaya[segment]);
        const int span = NurbsRefit::findSpan(4, 1, v, linear);
        require(span == segment + 1, "nonuniform linear span");
        NurbsRefit::basisFunctions(span, v, 1, linear, basis);
        require(basis.size() == 2, "linear basis size");
        partition(basis);
        close(basis[0], 1.0 - alpha, 1e-12, "left linear interpolation coefficient");
        close(basis[1], alpha, 1e-12, "right linear interpolation coefficient");
    }
}

void testKnotsAndGreville()
{
    for (int spans : {1, 4, 9})
    {
        const int m = spans + 3;
        const std::vector<double> knots = NurbsRefit::clampedUniformKnotsMaya(spans, 3);
        std::vector<double> expected = {0.0, 0.0, 0.0};
        for (int i = 1; i < spans; ++i)
            expected.push_back(static_cast<double>(i) / spans);
        expected.insert(expected.end(), 3, 1.0);
        require(knots == expected, "uniform Maya knot values");
        require(static_cast<int>(knots.size()) == m + 2, "Maya knot count");
        const std::vector<double> full = NurbsRefit::fullKnots(knots);
        require(static_cast<int>(full.size()) == m + 4, "full knot count");
        const std::vector<double> greville = NurbsRefit::grevilleAbscissae(full, 3, m);
        require(static_cast<int>(greville.size()) == m, "Greville count");
        for (int k = 0; k < m; ++k)
        {
            double sum = 0.0;
            for (int j = 1; j <= 3; ++j)
                sum += std::max(0, std::min(spans, k + j - 3)) / static_cast<double>(spans);
            close(greville[k], sum / 3.0, 1e-15, "Greville values");
        }
        require(greville.front() == 0.0 && greville.back() == 1.0, "Greville endpoints");
    }
    require(NurbsRefit::clampedUniformKnotsMaya(3, 1) ==
        std::vector<double>({0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0}), "linear uniform knots");
}

void testLinearPrecision()
{
    for (int degree : {1, 3})
    {
        const std::vector<double> maya = degree == 1
            ? std::vector<double>({2.0, 4.5, 6.0, 9.5})
            : std::vector<double>({2.0, 2.0, 2.0, 4.5, 6.0, 9.5, 9.5, 9.5});
        const int nv = static_cast<int>(maya.size()) - degree + 1;
        const std::vector<double> knots = NurbsRefit::fullKnots(maya);
        const std::vector<double> parameters = NurbsRefit::grevilleAbscissae(knots, degree, nv);
        const double origins[2][3] = {{1.0, -3.0, 8.0}, {-5.0, 2.0, 0.5}};
        const double directions[2][3] = {{2.0, -0.5, 0.0}, {-1.0, 3.0, 0.25}};
        std::vector<double> points(2 * nv * 3);
        for (int u = 0; u < 2; ++u)
            for (int v = 0; v < nv; ++v)
                for (int c = 0; c < 3; ++c)
                    points[(u * nv + v) * 3 + c] = origins[u][c] + parameters[v] * directions[u][c];
        for (int spans : {1, 4, 9})
        {
            const int m = spans + 3;
            std::vector<double> samples, output;
            require(NurbsRefit::sampleMatrix(nv, degree, maya, m, samples), "linear sample matrix");
            NurbsRefit::refitColumns(2, nv, points, m, samples, output);
            require(static_cast<int>(output.size()) == 2 * m * 3, "refit output count");
            const std::vector<double> outputKnots = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(spans, 3));
            for (int sample = 0; sample < 50; ++sample)
            {
                const double fraction = sample / 49.0;
                const double inputParameter = 2.0 + fraction * 7.5;
                for (int u = 0; u < 2; ++u)
                    for (int c = 0; c < 3; ++c)
                        close(deBoor(outputKnots, 3, m, output, u, c, fraction),
                            origins[u][c] + inputParameter * directions[u][c], 1e-10, "de Boor linear precision");
            }
        }
    }
}

void testConvexCombinations()
{
    const std::vector<double> maya = {0.0, 2.5, 4.0, 7.5};
    const int nu = 3;
    const int nv = 4;
    const std::vector<double> points = {
        -9, 2, 1, 3, -5, 8, 7, 6, -4, -2, 1, 3,
        20, -8, 12, 13, 2, -3, 25, 9, 1, 17, -6, 9,
        -12, 11, -5, -3, -4, 7, -8, 5, 16, -15, 2, 0
    };
    for (int spans : {1, 4, 9, 256})
    {
        const int m = spans + 3;
        std::vector<double> samples, output;
        require(NurbsRefit::sampleMatrix(nv, 1, maya, m, samples), "convex sample matrix");
        require(static_cast<int>(samples.size()) == m * nv, "sample matrix dimensions");
        for (int k = 0; k < m; ++k)
            partition(std::vector<double>(samples.begin() + k * nv, samples.begin() + (k + 1) * nv));
        NurbsRefit::refitColumns(nu, nv, points, m, samples, output);
        for (int u = 0; u < nu; ++u)
        {
            for (int c = 0; c < 3; ++c)
            {
                double low = points[u * nv * 3 + c];
                double high = low;
                for (int j = 1; j < nv; ++j)
                {
                    low = std::min(low, points[(u * nv + j) * 3 + c]);
                    high = std::max(high, points[(u * nv + j) * 3 + c]);
                }
                for (int k = 0; k < m; ++k)
                {
                    const double value = output[(u * m + k) * 3 + c];
                    const double tolerance = 1e-12 * std::max(1.0, std::max(std::fabs(low), std::fabs(high)));
                    require(std::isfinite(value) && value >= low - tolerance && value <= high + tolerance,
                        "component convex bounds");
                }
            }
        }
    }
}

void testMinimalRowsMaximalSpans()
{
    const std::vector<double> maya = {0.0, 3.0};
    const int nu = 2;
    const int nv = 2;
    const std::vector<double> points = {1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6};
    const int m = 256 + 3;
    std::vector<double> samples, output;
    require(NurbsRefit::sampleMatrix(nv, 1, maya, m, samples), "minimal rows sample matrix");
    require(static_cast<int>(samples.size()) == m * nv, "minimal rows sample dimensions");
    NurbsRefit::refitColumns(nu, nv, points, m, samples, output);
    require(static_cast<int>(output.size()) == nu * m * 3, "minimal rows output size");
    for (int u = 0; u < nu; ++u)
        for (int k = 0; k < m; ++k)
            for (int c = 0; c < 3; ++c)
            {
                const double low = std::min(points[u * nv * 3 + c], points[(u * nv + 1) * 3 + c]);
                const double high = std::max(points[u * nv * 3 + c], points[(u * nv + 1) * 3 + c]);
                const double value = output[(u * m + k) * 3 + c];
                require(std::isfinite(value) && value >= low - 1e-12 && value <= high + 1e-12, "minimal rows bounds");
            }
}

void testRepeatedKnots()
{
    const std::vector<double> full = {0.0, 0.0, 1.0, 1.0, 2.0, 3.0};
    require(NurbsRefit::findSpan(4, 1, 2.0, full) == 3, "right domain endpoint span");
    require(NurbsRefit::findSpan(4, 1, 1.0, full) == 3, "skip internal empty span");
    for (double parameter : {1.0, 2.0})
    {
        std::vector<double> basis;
        NurbsRefit::basisFunctions(NurbsRefit::findSpan(4, 1, parameter, full), parameter, 1, full, basis);
        partition(basis);
    }

    const std::vector<double> maya = {0.0, 2.5, 5.0, 5.0};
    const std::vector<double> repeatedEnd = NurbsRefit::fullKnots(maya);
    require(NurbsRefit::findSpan(4, 1, 5.0, repeatedEnd) == 2, "skip trailing empty span");
    require(NurbsRefit::findSpan(4, 1, 6.0, repeatedEnd) == 2, "right of domain span");
    std::vector<double> basis;
    NurbsRefit::basisFunctions(2, 5.0, 1, repeatedEnd, basis);
    partition(basis);
    require(basis == std::vector<double>({0.0, 1.0}), "repeated endpoint weights");
    std::vector<double> samples;
    require(NurbsRefit::sampleMatrix(4, 1, maya, 7, samples), "repeated endpoint sample matrix");
    for (int k = 0; k < 7; ++k)
        partition(std::vector<double>(samples.begin() + k * 4, samples.begin() + (k + 1) * 4));
}

void testParameterDomainValidation()
{
    for (int spans : {1, 4, 9})
    {
        for (double knot : {0.0, 2.5})
        {
            std::vector<double> samples;
            require(!NurbsRefit::sampleMatrix(3, 1, {knot, knot, knot}, spans + 3, samples),
                "reject all-equal V knots");
        }
    }
    std::vector<double> samples;
    require(!NurbsRefit::sampleMatrix(2, 1, {0.0, 0.5e-12}, 7, samples), "reject subthreshold domain");
    require(!NurbsRefit::sampleMatrix(2, 1, {0.0, 1e-12}, 7, samples), "reject threshold domain");
    require(NurbsRefit::sampleMatrix(2, 1, {0.0, 2e-12}, 7, samples), "accept domain above threshold");
}

void testClampedRightEndpoint()
{
    const std::vector<double> maya = NurbsRefit::clampedUniformKnotsMaya(1, 3);
    require(maya == std::vector<double>({0.0, 0.0, 0.0, 1.0, 1.0, 1.0}), "single span knots");
    const std::vector<double> full = NurbsRefit::fullKnots(maya);
    const int span = NurbsRefit::findSpan(4, 3, 1.0, full);
    require(span == 3, "clamped right span");
    std::vector<double> basis;
    NurbsRefit::basisFunctions(span, 1.0, 3, full, basis);
    require(basis == std::vector<double>({0.0, 0.0, 0.0, 1.0}), "clamped last CV basis");
}

std::vector<double> slitPoints()
{
    return {
        0.0, 0.0, 0.0,
        0.0, 0.5, 0.0,
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.5, 0.0,
        1.0, 1.0, 0.5
    };
}

double columnDifference(const std::vector<double>& knots, int nCV, const std::vector<double>& points,
    double parameter)
{
    return deBoor(knots, 3, nCV, points, 1, 2, parameter)
        - deBoor(knots, 3, nCV, points, 0, 2, parameter);
}

void requireJunctions(int nu, int nv, const std::vector<double>& rows, int spansV,
    const std::vector<double>& protectedV, const std::vector<double>& samples,
    const std::vector<double>& points, const std::vector<double>& output)
{
    const std::vector<double> bounds = NurbsRefit::protectedBoundaries(protectedV);
    const std::vector<int> spans = NurbsRefit::allocateBandSpans(spansV, bounds);
    const int m = NurbsRefit::protectedCvCount(spansV, static_cast<int>(protectedV.size()));
    std::vector<double> referenceOutput(nu * m * 3, 0.0);
    int boundary = 0;
    for (size_t band = 1; band < spans.size(); ++band)
    {
        boundary += spans[band - 1] + 2;
        const std::vector<double> upper = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(spans[band - 1], 3));
        const std::vector<double> lower = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(spans[band], 3));
        const int last = spans[band - 1] + 2;
        const double upperWidth = (bounds[band] - bounds[band - 1]) * (upper[last + 3] - upper[last]);
        const double lowerWidth = (bounds[band + 1] - bounds[band]) * (lower[4] - lower[1]);
        const double rho = lowerWidth / upperWidth;
        const double effective = std::min(rho, 4.0);
        int inputBoundary = -1;
        for (int j = 0; j < nv; ++j)
        {
            if (rows[j] == bounds[band])
                inputBoundary = j;
            close(samples[(boundary + 1) * nv + j],
                (1.0 + effective) * samples[boundary * nv + j] - effective * samples[(boundary - 1) * nv + j],
                1e-12, "junction coefficient uses preceding endpoint and penultimate row");
            if (rows[j] < bounds[band - 1] || rows[j] > bounds[band])
                require(samples[(boundary + 1) * nv + j] == 0.0, "junction forbidden dependency");
        }
        require(inputBoundary >= 0, "boundary input exists");
        for (int j = 0; j < nv; ++j)
            require(samples[boundary * nv + j] == (j == inputBoundary ? 1.0 : 0.0), "boundary one-hot");
        for (int u = 0; u < nu; ++u)
        {
            double up[3], down[3];
            for (int c = 0; c < 3; ++c)
            {
                const int index = (u * m + boundary) * 3 + c;
                require(sameBits(output[index], points[(u * nv + inputBoundary) * 3 + c]), "boundary CV input bits");
                up[c] = output[index] - output[index - 3];
                down[c] = output[index + 3] - output[index];
                if (rho <= 4.0)
                    close(3.0 * down[c] / lowerWidth, 3.0 * up[c] / upperWidth, 1e-12, "one-sided derivatives match");
                else
                    close(down[c], 4.0 * up[c], 1e-12, "capped tangent displacement");
            }
            if (rho > 4.0)
            {
                double crossSquared = 0.0, upSquared = 0.0, downSquared = 0.0, dot = 0.0;
                for (int c = 0; c < 3; ++c)
                {
                    const double cross = up[(c + 1) % 3] * down[(c + 2) % 3] - up[(c + 2) % 3] * down[(c + 1) % 3];
                    crossSquared += cross * cross;
                    upSquared += up[c] * up[c];
                    downSquared += down[c] * down[c];
                    dot += up[c] * down[c];
                }
                close(std::sqrt(crossSquared), 0.0, 1e-12, "capped tangents parallel");
                close(std::sqrt(downSquared), 4.0 * std::sqrt(upSquared), 1e-12, "capped tangent length");
                require(dot > 0.0, "capped tangents same direction");
            }
        }
        std::vector<double> changed = points, changedSamples, changedOutput;
        for (int u = 0; u < nu; ++u)
            for (int j = 0; j < nv; ++j)
                if (rows[j] > bounds[band])
                    for (int c = 0; c < 3; ++c)
                        changed[(u * nv + j) * 3 + c] += (u + 1) * (j + 1) * (c + 1);
        require(NurbsRefit::sampleMatrix(nv, 1, rows, spansV, protectedV, changedSamples), "changed lower input matrix");
        NurbsRefit::refitColumns(nu, nv, changed, m, changedSamples, changedOutput);
        for (int k = 0; k <= boundary; ++k)
        {
            for (int j = 0; j < nv; ++j)
            {
                require(sameBits(samples[k * nv + j], changedSamples[k * nv + j]), "upper matrix bits independent");
                if (rows[j] > bounds[band])
                    require(samples[k * nv + j] == 0.0, "upper matrix has no lower dependency");
            }
            for (int u = 0; u < nu; ++u)
                for (int c = 0; c < 3; ++c)
                    require(sameBits(output[(u * m + k) * 3 + c], changedOutput[(u * m + k) * 3 + c]),
                        "upper CV bits independent of lower input");
        }
    }
    int start = 0;
    for (size_t band = 0; band < spans.size(); ++band)
    {
        const int count = spans[band] + 3;
        const std::vector<double> full = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(spans[band], 3));
        const std::vector<double> greville = NurbsRefit::grevilleAbscissae(full, 3, count);
        for (int local = band == 0 ? 0 : 1; local < count; ++local)
        {
            const int k = start + local;
            const bool junction = band > 0 && local == 1;
            double sum = 0.0;
            for (int j = 0; j < nv; ++j)
            {
                require(std::isfinite(samples[k * nv + j]), "finite coefficients");
                sum += samples[k * nv + j];
            }
            double rho = 0.0;
            if (junction)
            {
                const std::vector<double> upper = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(spans[band - 1], 3));
                const int last = spans[band - 1] + 2;
                rho = ((bounds[band + 1] - bounds[band]) * (full[4] - full[1]))
                    / ((bounds[band] - bounds[band - 1]) * (upper[last + 3] - upper[last]));
            }
            close(sum, 1.0, 4.0 * std::numeric_limits<double>::epsilon() * (1.0 + rho), "coefficient row sum");
            if (junction)
            {
                for (int u = 0; u < nu; ++u)
                    for (int c = 0; c < 3; ++c)
                    {
                        const int index = (u * m + k) * 3 + c;
                        referenceOutput[index] = (1.0 + std::min(rho, 4.0)) * referenceOutput[index - 3]
                            - std::min(rho, 4.0) * referenceOutput[index - 6];
                        close(output[index], referenceOutput[index], 1e-12, "junction CV extrapolation");
                    }
                continue;
            }
            const double v = local == count - 1 ? bounds[band + 1]
                : bounds[band] + (bounds[band + 1] - bounds[band]) * greville[local];
            std::vector<double> oldRow(nv, 0.0);
            bool exact = false;
            for (int j = 0; j < nv; ++j)
                if (rows[j] >= bounds[band] && rows[j] <= bounds[band + 1] && std::fabs(rows[j] - v) <= 1e-12)
                {
                    oldRow[j] = 1.0;
                    exact = true;
                    break;
                }
            if (!exact)
                for (int j = 0; j + 1 < nv; ++j)
                    if (rows[j] <= v && v <= rows[j + 1])
                    {
                        const double weight = (v - rows[j]) / (rows[j + 1] - rows[j]);
                        oldRow[j] = 1.0 - weight;
                        oldRow[j + 1] = weight;
                        break;
                    }
            partition(oldRow);
            for (int j = 0; j < nv; ++j)
                require(sameBits(samples[k * nv + j], oldRow[j]), "old sample row bits preserved");
            for (int u = 0; u < nu; ++u)
                for (int c = 0; c < 3; ++c)
                {
                    double expected = 0.0;
                    for (int j = 0; j < nv; ++j)
                        expected += oldRow[j] * points[(u * nv + j) * 3 + c];
                    require(sameBits(output[(u * m + k) * 3 + c], expected), "old sampled CV bits preserved");
                    referenceOutput[(u * m + k) * 3 + c] = expected;
                }
        }
        const std::vector<double> knots = NurbsRefit::fullKnots(NurbsRefit::protectedKnotsMaya(spansV, protectedV));
        for (double t : {0.0, 0.2, 0.5, 0.8, 1.0})
        {
            const double v = bounds[band] + t * (bounds[band + 1] - bounds[band]);
            for (int u = 0; u < nu; ++u)
                for (int c = 0; c < 3; ++c)
                    close(deBoor(knots, 3, m, output, u, c, v),
                        deBoor(knots, 3, m, referenceOutput, u, c, v), 1e-10, "general band de Boor from expected CVs");
        }
        start += spans[band] + 2;
    }
}

void testUnprotectedBaseline()
{
    const std::vector<double> expectedSamples = {
        1, 0, 0, 0.33333333333333337, 0.66666666666666663, 0,
        0, 0.66666666666666674, 0.33333333333333326, 0, 0, 1
    };
    const std::vector<double> expectedKnots = {0, 0, 0, 1, 1, 1};
    const std::vector<double> expectedOutput = {
        0, 0, 0, 0, 0.33333333333333331, 0, 0, 0.66666666666666663, 0, 0, 1, 0,
        1, 0, 0, 1, 0.33333333333333331, 0, 1, 0.66666666666666663, 0.16666666666666663, 1, 1, 0.5
    };
    std::vector<double> samples, output;
    require(NurbsRefit::sampleMatrix(3, 1, {0, 0.5, 1}, 4, samples), "unprotected baseline matrix");
    NurbsRefit::refitColumns(2, 3, slitPoints(), 4, samples, output);
    const std::vector<double> knots = NurbsRefit::clampedUniformKnotsMaya(1, 3);
    require(samples.size() == expectedSamples.size() &&
        std::memcmp(samples.data(), expectedSamples.data(), samples.size() * sizeof(double)) == 0, "saved unprotected matrix bytes");
    require(knots.size() == expectedKnots.size() &&
        std::memcmp(knots.data(), expectedKnots.data(), knots.size() * sizeof(double)) == 0, "saved unprotected knot bytes");
    require(output.size() == expectedOutput.size() &&
        std::memcmp(output.data(), expectedOutput.data(), output.size() * sizeof(double)) == 0, "saved unprotected CV bytes");
}

void testProtectedFit()
{
    const std::vector<double> maya = {0.0, 0.5, 1.0};
    const int nv = 3;
    const std::vector<double> points = slitPoints();

    std::vector<double> plainSamples, plainOutput;
    require(NurbsRefit::sampleMatrix(nv, 1, maya, 4, plainSamples), "plain sample matrix");
    NurbsRefit::refitColumns(2, nv, points, 4, plainSamples, plainOutput);
    const std::vector<double> plainKnots = NurbsRefit::fullKnots(NurbsRefit::clampedUniformKnotsMaya(1, 3));
    close(columnDifference(plainKnots, 4, plainOutput, 0.25), 1.0 / 32.0, 1e-12, "plain V=0.25 difference");

    const std::vector<double> protectedV = {0.5};
    const int spansV = 2;
    const int m = NurbsRefit::protectedCvCount(spansV, 1);
    require(m == 7, "protected CV count");
    std::vector<double> samples, output;
    require(NurbsRefit::sampleMatrix(nv, 1, maya, spansV, protectedV, samples), "protected sample matrix");
    require(static_cast<int>(samples.size()) == m * nv, "protected sample dimensions");
    NurbsRefit::refitColumns(2, nv, points, m, samples, output);
    require(static_cast<int>(output.size()) == 2 * m * 3, "protected output count");

    for (int k = 0; k < 4; ++k)
    {
        const double difference = output[(1 * m + k) * 3 + 2] - output[(0 * m + k) * 3 + 2];
        require(difference == 0.0, "upper band CV difference is zero");
    }
    for (int u = 0; u < 2; ++u)
        for (int c = 0; c < 3; ++c)
            require(sameBits(output[(u * m + 3) * 3 + c], points[(u * nv + 1) * 3 + c]),
                "protected tip row bit match");

    const std::vector<double> outputKnots = NurbsRefit::fullKnots(NurbsRefit::protectedKnotsMaya(spansV, protectedV));
    require(NurbsRefit::protectedKnotsMaya(spansV, protectedV)
            == std::vector<double>({0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}),
        "protected Maya knots");
    require(columnDifference(outputKnots, m, output, 0.25) == 0.0, "upper interior difference is zero");
    require(columnDifference(outputKnots, m, output, 0.5 - 1e-8) == 0.0, "just above protected difference is zero");
    require(columnDifference(outputKnots, m, output, 0.5) == 0.0, "protected position difference is zero");
    for (int u = 0; u < 2; ++u)
        for (int c = 0; c < 3; ++c)
            close(deBoor(outputKnots, 3, m, output, u, c, 0.5), points[(u * nv + 1) * 3 + c], 1e-10,
                "de Boor at protected position");
    require(columnDifference(outputKnots, m, output, 0.5 + 1e-8) != 0.0, "just below protected can differ");
    close(columnDifference(outputKnots, m, output, 1.0), 0.5, 1e-10, "hem difference");
    const double expectedRow[] = {-1.0 / 3.0, 4.0 / 3.0, 0.0};
    for (int j = 0; j < nv; ++j)
        close(samples[4 * nv + j], expectedRow[j], 1e-12, "junction extrapolation coefficients");
    const double expectedY[] = {0, 1.0 / 6.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1};
    const double expectedZ[] = {0, 0, 0, 0, 0, 1.0 / 3.0, 0.5};
    for (int u = 0; u < 2; ++u)
        for (int k = 0; k < m; ++k)
        {
            close(output[(u * m + k) * 3], u, 1e-12, "slit CV x");
            close(output[(u * m + k) * 3 + 1], expectedY[k], 1e-12, "slit CV y");
            close(output[(u * m + k) * 3 + 2], u * expectedZ[k], 1e-12, "slit CV z");
        }
    for (double t : {0.0, 0.1, 0.5, 0.75, 1.0})
        close(columnDifference(outputKnots, m, output, 0.5 + t * 0.5), t * t - t * t * t / 2.0,
            1e-12, "lower-band cubic difference");
    close(columnDifference(outputKnots, m, output, 0.75), 3.0 / 16.0, 1e-12, "lower midpoint difference");
    requireJunctions(2, nv, maya, spansV, protectedV, samples, points, output);
    std::vector<double> outside;
    NurbsRefit::refitColumns(1, 3, {0, 0, 0, 1, 0, 0, 1, 0, 0}, m, samples, outside);
    close(outside[4 * 3], 4.0 / 3.0, 1e-12, "junction may exceed input convex hull");
    require(outside[4 * 3] > 1.0, "extrapolation above input maximum");
}

void requireConvexOutput(int nu, int nv, int m, const std::vector<double>& points,
    const std::vector<double>& output, const std::vector<int>& junctionRows)
{
    for (int u = 0; u < nu; ++u)
    {
        for (int c = 0; c < 3; ++c)
        {
            double low = points[(u * nv) * 3 + c];
            double high = low;
            double bound = std::fabs(low);
            for (int j = 1; j < nv; ++j)
            {
                const double value = points[(u * nv + j) * 3 + c];
                low = std::min(low, value);
                high = std::max(high, value);
                bound = std::max(bound, std::fabs(value));
            }
            const double tolerance = 1e-12 * std::max(1.0, bound);
            for (int k = 0; k < m; ++k)
            {
                if (std::find(junctionRows.begin(), junctionRows.end(), k) != junctionRows.end())
                    continue;
                const double value = output[(u * m + k) * 3 + c];
                require(std::isfinite(value) && value >= low - tolerance && value <= high + tolerance,
                    "protected convex bounds");
            }
        }
    }
}

void testFitMatrixAndValidation()
{
    require(NurbsRefit::allocateBandSpans(5, {0.0, 0.5, 1.0}) == std::vector<int>({3, 2}),
        "equal bands leftover to waist");
    require(NurbsRefit::allocateBandSpans(5, {0.0, 0.25, 1.0}) == std::vector<int>({2, 3}),
        "unequal band span allocation");

    require(NurbsRefit::allocateBandSpans(5, {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0})
        == std::vector<int>({2, 2, 1}), "thirds prefer waist on rounded ties");
    require(NurbsRefit::allocateBandSpans(4, {0.0, 0.1, 0.4, 0.7, 1.0})
        == std::vector<int>({1, 1, 1, 1}), "minimum band allocation");
    require(NurbsRefit::allocateBandSpans(5, {0.0, 0.1, 0.4, 0.7, 1.0})
        == std::vector<int>({1, 2, 1, 1}), "rounded equal remainders prefer waist");
    require(NurbsRefit::allocateBandSpans(5, {0.0, 0.1, 0.4, 0.7 - 1e-10, 1.0})
        == std::vector<int>({1, 1, 1, 2}), "distinct remainders retain order");
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 3,
        {0.5 - 0.25e-9, 0.5 + 0.25e-9}) == 14, "one row matching two boundaries is 14");
    std::vector<double> excessive;
    for (int i = 1; i <= 257; ++i)
        excessive.push_back(i / 258.0);
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 256, excessive) == 12,
        "oversized protection array is 12");
    for (const std::vector<double>& rows : {
        std::vector<double>{0.0, 0.1, 0.4, 0.7, 1.0},
        std::vector<double>{0.0, 0.1, 0.2, 0.2, 0.4, 0.7, 1.0}})
    {
        const int count = static_cast<int>(rows.size());
        std::vector<double> matrix;
        require(NurbsRefit::sampleMatrix(count, 1, rows, 4, {0.1, 0.4, 0.7}, matrix),
            "rounded boundary and repeated interior rows accepted");
        const int boundaryRow = rows.size() == 5 ? 2 : 4;
        for (int j = 0; j < count; ++j)
            require(matrix[6 * count + j] == (j == boundaryRow ? 1.0 : 0.0),
                "boundary coefficients are exactly one-hot");
        if (count == 7)
            for (int j = 0; j < count; ++j)
                close(matrix[4 * count + j], j == 0 ? -1.0 : (j == 1 ? 2.0 : 0.0),
                    j < 2 ? 1e-12 : 0.0, "duplicate Greville junction extrapolation");
    }

    std::vector<double> snappedSamples;
    require(NurbsRefit::sampleMatrix(4, 1, {0.0, 0.1, 0.3 + 5e-10, 1.0}, 3,
        {0.1, 0.3}, snappedSamples), "non-dyadic snapped boundary accepted");
    for (int j = 0; j < 4; ++j)
        require(snappedSamples[6 * 4 + j] == (j == 2 ? 1.0 : 0.0),
            "snapped non-dyadic boundary is one-hot");

    const std::vector<double> unitMaya = {0.0, 0.25, 0.5, 0.75, 1.0};
    const int nv = 5;
    const std::vector<double> protectedV = {0.25, 0.5, 0.75};
    const int spansV = 5;
    const int m = NurbsRefit::protectedCvCount(spansV, 3);
    require(m == 14, "four-band CV count");

    std::vector<double> points(2 * nv * 3, 0.0);
    std::mt19937 random(173);
    for (double& value : points)
        value = static_cast<double>(random()) / random.max();

    std::vector<double> samples, output;
    require(NurbsRefit::sampleMatrix(nv, 1, unitMaya, spansV, protectedV, samples), "multi-band sample matrix");
    require(static_cast<int>(samples.size()) == m * nv, "multi-band sample dimensions");
    const std::vector<int> bandSpans = NurbsRefit::allocateBandSpans(spansV, {0.0, 0.25, 0.5, 0.75, 1.0});
    require(bandSpans == std::vector<int>({2, 1, 1, 1}), "four equal bands leftover to waist");
    int rowIndex = 0;
    for (int band = 0; band < 4; ++band)
    {
        const int localCount = bandSpans[band] + 3;
        const int localStart = band == 0 ? 0 : 1;
        for (int local = localStart; local < localCount; ++local, ++rowIndex)
        {
            std::vector<double> row(samples.begin() + rowIndex * nv, samples.begin() + (rowIndex + 1) * nv);
            const bool junction = band > 0 && local == 1;
            if (!junction)
                partition(row);
            for (int j = 0; j < nv; ++j)
            {
                if (junction ? (j < band - 1 || j > band) : (j < band || j > band + 1))
                    require(row[j] == 0.0, "out-of-band coefficient is zero");
            }
        }
    }
    require(rowIndex == m, "band row concatenation");
    NurbsRefit::refitColumns(2, nv, points, m, samples, output);
    requireConvexOutput(2, nv, m, points, output, {5, 8, 11});
    requireJunctions(2, nv, unitMaya, spansV, protectedV, samples, points, output);
    for (int r : {4, 7, 10})
        for (int j = 0; j < nv; ++j)
            close(samples[(r + 1) * nv + j], (r == 4 ? 3.0 : 2.0) * samples[r * nv + j]
                - (r == 4 ? 2.0 : 1.0) * samples[(r - 1) * nv + j], 1e-12, "chained junction reads Q2");
    for (int sourceColumn = 0; sourceColumn < 2; ++sourceColumn)
        for (int sourceRow = 0; sourceRow < nv; ++sourceRow)
            for (int sourceComponent = 0; sourceComponent < 3; ++sourceComponent)
            {
                std::vector<double> impulse(2 * nv * 3, 0.0), response;
                const int sourceIndex = (sourceColumn * nv + sourceRow) * 3 + sourceComponent;
                impulse[sourceIndex] = 1.0;
                NurbsRefit::refitColumns(2, nv, impulse, m, samples, response);
                for (int u = 0; u < 2; ++u)
                    for (int k = 0; k < m; ++k)
                        for (int c = 0; c < 3; ++c)
                            require(response[(u * m + k) * 3 + c]
                                == (u == sourceColumn && c == sourceComponent
                                    ? samples[k * nv + sourceRow] : 0.0),
                                "U column impulse coefficients");
                std::vector<double> changed = points, changedOutput;
                changed[sourceIndex] += 1.0;
                NurbsRefit::refitColumns(2, nv, changed, m, samples, changedOutput);
                for (int k = 0; k < m; ++k)
                    for (int c = 0; c < 3; ++c)
                    {
                        const int index = ((1 - sourceColumn) * m + k) * 3 + c;
                        require(sameBits(output[index], changedOutput[index]), "other U column unchanged");
                    }
            }

    const std::vector<double> scaledMaya = {2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<double> scaledSamples;
    require(NurbsRefit::sampleMatrix(nv, 1, scaledMaya, spansV, protectedV, scaledSamples),
        "non-unit domain sample matrix");
    require(scaledSamples == samples, "normalized non-unit domain matches unit domain");

    const std::vector<double> outputKnots = NurbsRefit::fullKnots(NurbsRefit::protectedKnotsMaya(spansV, protectedV));
    const double breaks[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (int band = 0; band < 4; ++band)
    {
        const double start = breaks[band];
        const double end = breaks[band + 1];
        const double interior = 0.5 * (start + end);
        for (int u = 0; u < 2; ++u)
            for (int c = 0; c < 3; ++c)
            {
                const double b = points[(u * nv + band) * 3 + c];
                const double next = points[(u * nv + band + 1) * 3 + c];
                const double expected = band == 0 ? 0.5 * b + 0.5 * next
                    : -points[(u * nv + band - 1) * 3 + c] / 8.0 + 3.0 * b / 4.0 + 3.0 * next / 8.0;
                close(deBoor(outputKnots, 3, m, output, u, c, interior), expected, 1e-10,
                    "independent de Boor inside band");
            }
        for (int u = 0; u < 2; ++u)
            for (int c = 0; c < 3; ++c)
                close(deBoor(outputKnots, 3, m, output, u, c, end),
                    points[(u * nv + band + 1) * 3 + c], 1e-10, "independent de Boor at band end");
        if (band + 1 < 4)
        {
            const double above = end - 1e-8;
            const double below = end + 1e-8;
            for (int u = 0; u < 2; ++u)
                for (int c = 0; c < 3; ++c)
                {
                    const double startValue = points[(u * nv + band) * 3 + c];
                    const double midValue = points[(u * nv + band + 1) * 3 + c];
                    const double nextValue = points[(u * nv + band + 2) * 3 + c];
                    const double tAbove = (above - start) / (end - start);
                    const double tBelow = (below - end) / (breaks[band + 2] - end);
                    const double aboveExpected = (1.0 - tAbove) * startValue + tAbove * midValue
                        + (band == 0 ? 0.0 : tAbove * (1.0 - tAbove) * (1.0 - tAbove)
                            * (2.0 * startValue - points[(u * nv + band - 1) * 3 + c] - midValue));
                    const double belowExpected = (1.0 - tBelow) * midValue + tBelow * nextValue
                        + tBelow * (1.0 - tBelow) * (1.0 - tBelow) * (2.0 * midValue - startValue - nextValue);
                    close(deBoor(outputKnots, 3, m, output, u, c, above), aboveExpected, 1e-10,
                        "independent de Boor just above protected");
                    close(deBoor(outputKnots, 3, m, output, u, c, below), belowExpected, 1e-10,
                        "independent de Boor just below protected");
                }
        }
    }

    for (const std::vector<double>& bounds : {
        std::vector<double>{0, 0.25, 1}, std::vector<double>{0, 0.5, 1},
        std::vector<double>{0, 0.25, 0.5, 1}, std::vector<double>{0, 0.01, 1}})
    {
        const int count = static_cast<int>(bounds.size());
        const int total = count == 4 ? 3 : (bounds[1] == 0.01 ? 2 : 5);
        const std::vector<double> protectedRows(bounds.begin() + 1, bounds.end() - 1);
        const int outCount = NurbsRefit::protectedCvCount(total, count - 2);
        std::vector<double> input(2 * count * 3), matrix, cvs;
        for (int u = 0; u < 2; ++u)
            for (int j = 0; j < count; ++j)
                for (int c = 0; c < 3; ++c)
                    input[(u * count + j) * 3 + c] = (u + 1) * (c + 1) * j * j / 16.0;
        require(NurbsRefit::sampleMatrix(count, 1, bounds, total, protectedRows, matrix), "junction ratio fixture");
        NurbsRefit::refitColumns(2, count, input, outCount, matrix, cvs);
        requireJunctions(2, count, bounds, total, protectedRows, matrix, input, cvs);
    }

    const std::vector<double> cubicMaya = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    require(NurbsRefit::protectedFitCondition(5, 3, cubicMaya, 4, {0.5}) == 13, "higher V degree is 13");
    std::vector<double> rejected;
    require(!NurbsRefit::sampleMatrix(5, 3, cubicMaya, 4, std::vector<double>({0.5}), rejected),
        "reject higher V degree");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 3, protectedV) == 12, "spansV < B is 12");
    require(!NurbsRefit::sampleMatrix(nv, 1, unitMaya, 3, protectedV, rejected), "reject spansV < B");

    const double justBelow = std::nextafter(1e-12, 0.0);
    const double justAbove = std::nextafter(1e-12, 1.0);
    const auto domainKnots = [](double length) {
        return std::vector<double>({0.0, 0.5 * length, length});
    };
    require(NurbsRefit::protectedFitCondition(3, 1, domainKnots(2.0 * justBelow), 2, {0.5}) == 11,
        "band length just below 1e-12 is 11");
    require(NurbsRefit::protectedFitCondition(3, 1, domainKnots(2e-12), 2, {0.5}) == 11,
        "band length at 1e-12 is 11");
    require(NurbsRefit::protectedFitCondition(3, 1, domainKnots(2.0 * justAbove), 2, {0.5}) == 0,
        "band length just above 1e-12 is accepted");

    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {std::numeric_limits<double>::quiet_NaN()}) == 10,
        "non-finite protected is 10");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {std::numeric_limits<double>::infinity()}) == 10,
        "infinite protected is 10");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {0.5, 0.25}) == 10, "reversed protected is 10");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {0.0}) == 10, "protected 0 is 10");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {1.0}) == 10, "protected 1 is 10");
    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, 4, {0.5, 0.5}) == 10, "duplicate protected is 10");
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 2, {0.3}) == 14, "missing row is 14");
    require(NurbsRefit::protectedFitCondition(4, 1, {0.0, 0.5, 0.5 + 0.5e-9, 1.0}, 2, {0.5}) == 14,
        "multiple matching rows is 14");

    const double matchAt = 0.5 + 1e-9;
    const double matchBefore = std::nextafter(matchAt, 0.5);
    const double matchAfter = std::nextafter(matchAt, 1.0);
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 2, {matchBefore}) == 0,
        "match just inside 1e-9");
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 2, {matchAt}) == 0, "match at 1e-9");
    require(NurbsRefit::protectedFitCondition(3, 1, {0.0, 0.5, 1.0}, 2, {matchAfter}) == 14,
        "no match just outside 1e-9");

    const std::vector<double> offsetMaya = {0.0, 0.5, 1.0};
    std::vector<double> offsetPoints = {
        0.0, 0.0, 0.0,
        9.0, 9.0, 9.0,
        1.0, 1.0, 1.0,
        2.0, 2.0, 2.0,
        8.0, 8.0, 8.0,
        3.0, 3.0, 3.0
    };
    std::vector<double> offsetSamples, offsetOutput;
    require(NurbsRefit::sampleMatrix(3, 1, offsetMaya, 2, std::vector<double>({matchAt}), offsetSamples),
        "tolerance match sample matrix");
    NurbsRefit::refitColumns(2, 3, offsetPoints, 7, offsetSamples, offsetOutput);
    for (int u = 0; u < 2; ++u)
        for (int c = 0; c < 3; ++c)
            require(sameBits(offsetOutput[(u * 7 + 3) * 3 + c], offsetPoints[(u * 3 + 1) * 3 + c]),
                "matched row curve not averaged");

    require(NurbsRefit::protectedFitCondition(nv, 1, unitMaya, spansV, protectedV) == 0, "valid recovery");
    require(NurbsRefit::sampleMatrix(nv, 1, unitMaya, spansV, protectedV, samples), "valid recovery matrix");
}

}

int main()
{
    try
    {
        testBasisFunctions();
        testKnotsAndGreville();
        testLinearPrecision();
        testConvexCombinations();
        testMinimalRowsMaximalSpans();
        testRepeatedKnots();
        testParameterDomainValidation();
        testClampedRightEndpoint();
        testUnprotectedBaseline();
        testProtectedFit();
        testFitMatrixAndValidation();
        std::cout << "nurbs_refit_kernel: 11 tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "nurbs_refit_kernel: " << error.what() << '\n';
        return 1;
    }
}
