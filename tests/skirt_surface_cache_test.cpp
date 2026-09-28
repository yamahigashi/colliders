#include "skirtSurfaceCache.h"

#include <atomic>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

struct Topology
{
    int count;
    std::vector<double> knots;
};

struct Preparation
{
    std::shared_ptr<const Topology> topology;
    std::vector<double> weights;
    std::vector<double> rest;
};

using Cache = SkirtSurfaceCache<Topology, Preparation>;

static bool sameTopology(const Topology& a, const Topology& b)
{
    return a.count == b.count && a.knots == b.knots;
}

static std::shared_ptr<const Topology> topology(int key)
{
    auto value = std::make_shared<Topology>();
    value->count = key;
    value->knots.push_back(static_cast<double>(key));
    return value;
}

static std::shared_ptr<const Preparation> preparation(
    const std::shared_ptr<const Topology>& owner, double weight, double rest)
{
    auto value = std::make_shared<Preparation>();
    value->topology = owner;
    value->weights.push_back(weight);
    value->rest.push_back(rest);
    return value;
}

static bool matches(const Preparation& candidate, const std::shared_ptr<const Topology>& owner,
                    double weight, double rest)
{
    return candidate.topology == owner && candidate.weights.size() == 1 && candidate.rest.size() == 1 &&
           std::memcmp(&candidate.weights[0], &weight, sizeof(double)) == 0 &&
           std::memcmp(&candidate.rest[0], &rest, sizeof(double)) == 0;
}

static void testExactKeys()
{
    Cache cache;
    auto plusZero = std::make_shared<Topology>();
    plusZero->count = 1;
    plusZero->knots.push_back(0.0);
    auto cached = cache.publishTopology(plusZero, 100, sameTopology);
    Topology minusZero = {1, std::vector<double>(1, -0.0)};
    assert(cache.findTopology(minusZero, sameTopology) == cached);

    Topology nanKey = {1, std::vector<double>(1, std::numeric_limits<double>::quiet_NaN())};
    assert(!cache.findTopology(nanKey, sameTopology));
    auto nanValue = std::make_shared<Topology>(nanKey);
    cache.publishTopology(nanValue, 100, sameTopology);
    assert(!cache.findTopology(nanKey, sameTopology));
}

static void testReuseAndEviction()
{
    Cache cache;
    auto a = cache.publishTopology(topology(1), 100, sameTopology);
    auto b = cache.publishTopology(topology(2), 100, sameTopology);
    auto first = preparation(a, 1.0, 10.0);
    auto second = preparation(a, 2.0, 10.0);
    cache.publishPreparation(a, first, 40, [&](const Preparation& p) { return matches(p, a, 1.0, 10.0); });
    cache.publishPreparation(a, second, 40, [&](const Preparation& p) { return matches(p, a, 2.0, 10.0); });
    assert(cache.preparationCount() == 2);
    assert(cache.findTopology(*a, sameTopology) == a);
    assert(cache.findPreparation(a, [&](const Preparation& p) { return matches(p, a, 1.0, 10.0); }) == first);
    assert(!cache.findPreparation(a, [&](const Preparation& p) { return matches(p, a, 1.0, 11.0); }));
    assert(!cache.findPreparation(b, [&](const Preparation& p) { return matches(p, b, 1.0, 10.0); }));
    auto third = preparation(a, 3.0, 10.0);
    cache.publishPreparation(a, third, 40, [&](const Preparation& p) { return matches(p, a, 3.0, 10.0); });
    assert(cache.preparationCount() == 2);
    assert(!cache.findPreparation(a, [&](const Preparation& p) { return matches(p, a, 2.0, 10.0); }));
    for (int key = 3; key <= 5; ++key)
        cache.publishTopology(topology(key), 100, sameTopology);
    assert(cache.topologyCount() == 4);
    assert(!cache.findTopology(*b, sameTopology));
    assert(cache.findTopology(*a, sameTopology) == a);
    // A borrowed immutable preparation remains valid after its entry is evicted.
    for (int key = 6; key <= 9; ++key)
        cache.publishTopology(topology(key), 100, sameTopology);
    assert(!cache.findTopology(*a, sameTopology));
    assert(first->weights[0] == 1.0 && first->topology == a);
}

static void testPayloadLimits()
{
    Cache cache;
    auto a = cache.publishTopology(topology(1), Cache::maxPayloadBytes - 100, sameTopology);
    auto b = cache.publishTopology(topology(2), 200, sameTopology);
    assert(cache.topologyCount() == 1);
    assert(cache.retainedPayloadBytes() == 200);
    assert(!cache.findTopology(*a, sameTopology));
    auto huge = cache.publishTopology(topology(3), Cache::maxPayloadBytes + 1, sameTopology);
    assert(huge && cache.topologyCount() == 1);
    assert(cache.retainedPayloadBytes() <= Cache::maxPayloadBytes);
    auto oversizedPrep = preparation(b, 1.0, 1.0);
    cache.publishPreparation(b, oversizedPrep, Cache::maxPayloadBytes,
                             [&](const Preparation& p) { return matches(p, b, 1.0, 1.0); });
    assert(cache.preparationCount() == 0);
    assert(cache.retainedPayloadBytes() == 200);
}

static void testConcurrentPublication()
{
    Cache cache;
    const int workers = 8;
    std::atomic<int> ready(0);
    std::atomic<bool> start(false);
    std::shared_ptr<const Topology> returned[workers];
    std::thread threads[workers];
    for (int i = 0; i < workers; ++i)
        threads[i] = std::thread([&, i] {
            auto built = topology(7);
            ++ready;
            while (!start.load()) {}
            returned[i] = cache.publishTopology(built, 100, sameTopology);
        });
    while (ready.load() != workers) {}
    start = true;
    for (auto& thread : threads)
        thread.join();
    assert(cache.topologyCount() == 1);
    for (const auto& value : returned)
        assert(value == returned[0]);

    auto owner = returned[0];
    std::shared_ptr<const Preparation> prepared[workers];
    ready = 0;
    start = false;
    for (int i = 0; i < workers; ++i)
        threads[i] = std::thread([&, i] {
            auto built = preparation(owner, 1.0, 2.0);
            ++ready;
            while (!start.load()) {}
            prepared[i] = cache.publishPreparation(owner, built, 30,
                [&](const Preparation& p) { return matches(p, owner, 1.0, 2.0); });
        });
    while (ready.load() != workers) {}
    start = true;
    for (auto& thread : threads)
        thread.join();
    assert(cache.preparationCount() == 1);
    for (const auto& value : prepared)
        assert(value == prepared[0]);
}

static void testNoPublicationOnFailure()
{
    Cache cache;
    auto accepted = cache.publishTopology(topology(1), 100, sameTopology);
    auto failed = topology(2);
    assert(failed && cache.topologyCount() == 1);
    assert(cache.findTopology(*accepted, sameTopology) == accepted);
    auto unretained = cache.publishTopology(topology(3), Cache::maxPayloadBytes + 1, sameTopology);
    auto prepared = preparation(unretained, 1.0, 2.0);
    cache.publishPreparation(unretained, prepared, 30,
        [&](const Preparation& p) { return matches(p, unretained, 1.0, 2.0); });
    assert(cache.preparationCount() == 0);
}

int main()
{
    testExactKeys();
    testReuseAndEviction();
    testPayloadLimits();
    testConcurrentPublication();
    testNoPublicationOnFailure();
}
