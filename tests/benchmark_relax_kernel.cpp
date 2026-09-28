#include "bellColliderRelaxKernel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

struct Point
{
    double x, y, z;
};

static bool sameBits(double a, double b)
{
    std::uint64_t ua = 0, ub = 0;
    std::memcpy(&ua, &a, sizeof(ua));
    std::memcpy(&ub, &b, sizeof(ub));
    return ua == ub;
}

template <class Ops> static BellColliderRelax::Ring<Ops> makeRing()
{
    BellColliderRelax::Ring<Ops> ring{};
    ring.origin = {Ops::splat(0.35), Ops::splat(-0.8), Ops::splat(1.1)};
    ring.normal = {Ops::splat(0.0), Ops::splat(1.0), Ops::splat(0.0)};
    ring.translation = ring.origin;
    ring.inverseColumns[0] = {Ops::splat(1.0 / 1.5), Ops::splat(0.0), Ops::splat(0.0)};
    ring.inverseColumns[1] = {Ops::splat(0.0), Ops::splat(1.0), Ops::splat(0.0)};
    ring.inverseColumns[2] = {Ops::splat(0.0), Ops::splat(0.0), Ops::splat(1.0 / 0.8)};
    ring.collision = Ops::splat(0.65);
    return ring;
}

static std::vector<Point> makePoints(int count, int mode)
{
    std::vector<Point> points;
    points.reserve(count);
    std::uint32_t state = 0x4f1bbcdcU + static_cast<std::uint32_t>(mode);
    for (int i = 0; i != count; ++i)
    {
        state = state * 1664525U + 1013904223U;
        const double u = static_cast<double>(state) / 4294967296.0;
        state = state * 1664525U + 1013904223U;
        const double v = static_cast<double>(state) / 4294967296.0;
        state = state * 1664525U + 1013904223U;
        const double w = static_cast<double>(state) / 4294967296.0;
        const bool inside = mode == 0 || (mode == 2 && i % 2 == 0);
        const double radius = inside ? 0.25 + u * 0.5 : 1.25 + u;
        const double angle = v * 6.283185307179586;
        points.push_back({0.35 + 1.5 * radius * std::cos(angle), -0.8 + w, 1.1 + 0.8 * radius * std::sin(angle)});
    }
    return points;
}

static void scalarRun(const std::vector<Point> &input, std::vector<Point> &output)
{
    using Ops = BellColliderRelax::Scalar;
    const auto ring = makeRing<Ops>();
    for (std::size_t i = 0; i < input.size(); ++i)
    {
        const Point &p = input[i];
        const BellColliderRelax::Vector<Ops> raw{p.x, p.y, p.z};
        const auto result = BellColliderRelax::relax(raw, raw, ring);
        output[i] = {result.x, result.y, result.z};
    }
}

#ifdef YDD_RELAX_SSE2
static void pairRun(const std::vector<Point> &input, std::vector<Point> &output)
{
    using Ops = BellColliderRelax::Pair;
    const auto ring = makeRing<Ops>();
    for (std::size_t i = 0; i < input.size(); i += 2)
    {
        const std::size_t j = std::min(i + 1, input.size() - 1);
        const auto raw =
            BellColliderRelax::Vector<Ops>{_mm_setr_pd(input[i].x, input[j].x), _mm_setr_pd(input[i].y, input[j].y),
                                           _mm_setr_pd(input[i].z, input[j].z)};
        const auto result = BellColliderRelax::relax(raw, raw, ring);
        alignas(16) double x[2], y[2], z[2];
        _mm_store_pd(x, result.x);
        _mm_store_pd(y, result.y);
        _mm_store_pd(z, result.z);
        output[i] = {x[0], y[0], z[0]};
        if (i + 1 < input.size())
            output[i + 1] = {x[1], y[1], z[1]};
    }
}
#endif

static double checksum(const std::vector<Point> &points)
{
    volatile double result = 0.0;
    for (const Point &p : points)
        result += p.x * 0.5 + p.y * 0.25 + p.z * 0.125;
    return result;
}

using Run = void (*)(const std::vector<Point> &, std::vector<Point> &);

static double measure(Run function, const std::vector<Point> &input, std::vector<Point> &output)
{
    // The volatile function pointer prevents calls with identical inputs being hoisted.
    Run volatile run = function;
    const int iterations = 5000;
    const auto start = std::chrono::steady_clock::now();
    for (int iteration = 0; iteration < iterations; ++iteration)
        run(input, output);
    const auto stop = std::chrono::steady_clock::now();
    (void)checksum(output);
    return std::chrono::duration<double, std::nano>(stop - start).count() / (iterations * input.size());
}

static void printSamples(const std::vector<double> &samples)
{
    std::cout << "[";
    for (std::size_t i = 0; i < samples.size(); ++i)
        std::cout << (i ? "," : "") << samples[i];
    std::cout << "]";
}

int main()
{
#ifndef YDD_RELAX_SSE2
    std::cerr << "This comparison requires SSE2. The production scalar path remains available.\n";
    return 1;
#else
    const int counts[] = {17, 65, 257};
    const char *modes[] = {"inside", "outside", "mixed"};
    std::cout << std::setprecision(12) << "{\"iterations\":5000,\"batches\":15,\"cases\":[";
    bool first = true;
    for (int mode = 0; mode < 3; ++mode)
        for (int count : counts)
        {
            const auto input = makePoints(count, mode);
            std::vector<Point> scalar(input.size()), pair(input.size());
            scalarRun(input, scalar);
            pairRun(input, pair);
            for (std::size_t i = 0; i != scalar.size(); ++i)
                if (!sameBits(scalar[i].x, pair[i].x) || !sameBits(scalar[i].y, pair[i].y) ||
                    !sameBits(scalar[i].z, pair[i].z))
                {
                    std::cerr << "bitwise mismatch mode=" << mode << " count=" << count << " index=" << i << "\n";
                    return 1;
                }
            for (int warm = 0; warm != 100; ++warm)
            {
                scalarRun(input, scalar);
                pairRun(input, pair);
            }
            std::vector<double> scalarNs, pairNs;
            for (int batch = 0; batch != 15; ++batch)
            {
                if (batch % 2 == 0)
                {
                    scalarNs.push_back(measure(scalarRun, input, scalar));
                    pairNs.push_back(measure(pairRun, input, pair));
                }
                else
                {
                    pairNs.push_back(measure(pairRun, input, pair));
                    scalarNs.push_back(measure(scalarRun, input, scalar));
                }
            }
            if (!first)
                std::cout << ",";
            first = false;
            std::cout << "{\"mode\":\"" << modes[mode] << "\",\"count\":" << count
                      << ",\"all_exact\":true,\"scalar_batch_ns_per_point\":";
            printSamples(scalarNs);
            std::cout << ",\"pair_batch_ns_per_point\":";
            printSamples(pairNs);
            std::sort(scalarNs.begin(), scalarNs.end());
            std::sort(pairNs.begin(), pairNs.end());
            std::cout << ",\"scalar_ns_per_point\":" << scalarNs[scalarNs.size() / 2]
                      << ",\"pair_ns_per_point\":" << pairNs[pairNs.size() / 2]
                      << ",\"scalar_checksum\":" << checksum(scalar) << ",\"pair_checksum\":" << checksum(pair) << "}";
        }
    std::cout << "]}\n";
    return 0;
#endif
}
