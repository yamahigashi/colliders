"""Exercise the workspace and pool code extracted from the current deformer source.

Run after Maya timing ends: python skirt_surface_workspace_test.py [--source path/to/skirtCollideDeformer.cpp]
The default source is ../sources/skirtCollideDeformer.cpp relative to this file.
"""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def extract(text, start, end):
    begin = text.index(start)
    return text[begin:text.index(end, begin)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path,
        default=Path(__file__).resolve().parents[1] / "sources/skirtCollideDeformer.cpp",
    )
    args = parser.parse_args()
    source = args.source.resolve()
    text = source.read_text()
    types = extract(text, "struct SurfaceContact\n", "\nnamespace\n{\n\nvoid warnInvalidBasis")
    methods = extract(text, "SkirtCollideDeformer::~SkirtCollideDeformer()", "MStatus SkirtCollideDeformer::initialize()")

    harness = r'''
#include <algorithm>
#include <emmintrin.h>
#define YDD_SURFACE_SSE2 1
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>
#include "skirtLegModel.h"
struct MPoint { double x = 0.0, y = 0.0, z = 0.0, w = 1.0; };
// This test concerns ownership and reset, not TBB allocator alignment.
namespace tbb { template<class T> using cache_aligned_allocator = std::allocator<T>; }
'''
    harness += types
    harness += r'''
class SkirtCollideDeformer {
public:
    SkirtCollideDeformer() : workspacePool(std::make_shared<SkirtCollideSurfaceWorkspacePool>()) {}
    ~SkirtCollideDeformer();
    static std::unique_ptr<SkirtCollideSurfaceWorkspace>
    acquireWorkspace(const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>&);
    static void releaseWorkspace(const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>&,
                                 std::unique_ptr<SkirtCollideSurfaceWorkspace>) noexcept;
    std::shared_ptr<SkirtCollideSurfaceWorkspacePool> workspacePool;
};
'''
    harness += methods
    harness += r'''
int main() {
    SkirtCollideDeformer node;
    auto pool = node.workspacePool;
    auto first = node.acquireWorkspace(pool);
    auto* address = first.get();
    first->reset(96, 400, 12, 64, 20, 200, 80, 3);
    first->restPoints.resize(96);
    first->restPoints[80] = MPoint{7.0, 8.0, 9.0, 2.0};
    first->restValid.resize(96, 1);
    first->members.resize(96, 1);
    first->membership.resize(96, 80);
    first->paint.resize(96, 0.2);
    first->rawPaint.resize(96, std::make_pair(80u, 0.2));
    const auto restCapacity = first->restPoints.capacity();
    const auto paintCapacity = first->paint.capacity();
    assert(first->capsuleColumns[0].p0x.empty());
    first->capsuleColumns[0].p0x.reserve(12);
    auto capsuleCapacity = first->capsuleColumns[0].p0x.capacity();
    first->samples.denominators[0] = 9;
    first->rows[0] = 9;
    first->beta[0] = 9;
    first->groups[0].paint = 0;
    first->groups[0].value = LegVec3(7, 8, 9);
    first->transports[0][0].valid[0] = true;
    first->chunkScratch[2].value.polynomial[0].roots.reserve(100);
    const auto rootsCapacity = first->chunkScratch[2].value.polynomial[0].roots.capacity();
    node.releaseWorkspace(pool, std::move(first));

    auto smaller = node.acquireWorkspace(pool);
    assert(smaller.get() == address);
    assert(smaller->restPoints.capacity() == restCapacity);
    assert(smaller->paint.capacity() == paintCapacity);
    smaller->restPoints.clear();
    smaller->restValid.clear();
    smaller->restPoints.resize(4);
    smaller->restValid.resize(4, 0);
    smaller->members.resize(4);
    std::fill(smaller->members.begin(), smaller->members.end(), 0);
    smaller->membership.clear();
    smaller->paint.resize(4);
    std::fill(smaller->paint.begin(), smaller->paint.end(), 1.0);
    smaller->rawPaint.clear();
    const MPoint defaultRest;
    assert(smaller->restPoints[2].x == defaultRest.x &&
           smaller->restPoints[2].y == defaultRest.y &&
           smaller->restPoints[2].z == defaultRest.z &&
           smaller->restPoints[2].w == defaultRest.w);
    assert(smaller->restValid[2] == 0 && smaller->members[2] == 0);
    assert(smaller->membership.empty() && smaller->rawPaint.empty());
    assert(smaller->paint[2] == 1.0);
    // A sparse rest input must not expose values left by a previous evaluation.
    smaller->restPoints.clear();
    smaller->restValid.clear();
    smaller->restPoints.resize(96);
    smaller->restValid.resize(96, 0);
    assert(smaller->restPoints[80].x == defaultRest.x &&
           smaller->restPoints[80].y == defaultRest.y &&
           smaller->restPoints[80].z == defaultRest.z &&
           smaller->restPoints[80].w == defaultRest.w);
    assert(smaller->restValid[80] == 0);
    assert(smaller->capsuleColumns[0].p0x.capacity() == capsuleCapacity);
    smaller->reset(1, 2, 1, 2, 1, 1, 1, 1);
    assert(smaller->samples.denominators[0] == 0);
    assert(smaller->capsuleColumns[0].p0x.capacity() == capsuleCapacity);
    assert(smaller->rows[0] == 0 && smaller->beta[0] == 0);
    assert(smaller->groups[0].paint == 1);
    assert(smaller->groups[0].value.x == 0 && smaller->groups[0].value.y == 0);
    assert(!smaller->transports[0][0].valid[0]);
    assert(smaller->chunkScratch.size() == 3);
    assert(smaller->chunkScratch[2].value.polynomial[0].roots.capacity() == rootsCapacity);
    smaller->reset(96, 400, 12, 64, 20, 200, 80, 3);
    assert(smaller->samples.denominators[0] == 0);
    assert(smaller->chunkScratch[2].value.polynomial[0].roots.capacity() == rootsCapacity);
    node.releaseWorkspace(pool, std::move(smaller));

    auto one = node.acquireWorkspace(pool);
    auto two = node.acquireWorkspace(pool);
    assert(one.get() != two.get());
    node.releaseWorkspace(pool, std::move(one));
    node.releaseWorkspace(pool, std::move(two));
    assert(pool->idle.size() == 2);
    auto three = node.acquireWorkspace(pool);
    auto four = node.acquireWorkspace(pool);
    auto five = node.acquireWorkspace(pool);
    assert(three.get() != four.get() && four.get() != five.get() && three.get() != five.get());
    node.releaseWorkspace(pool, std::move(three));
    node.releaseWorkspace(pool, std::move(four));
    node.releaseWorkspace(pool, std::move(five));
    assert(pool->idle.size() == 2);

    // The input vectors count toward the pool's 32 MiB retained-byte limit.
    SkirtCollideDeformer oversizedNode;
    auto oversizedPool = oversizedNode.workspacePool;
    auto oversized = oversizedNode.acquireWorkspace(oversizedPool);
    oversized->restPoints.reserve(32 * 1024 * 1024 / sizeof(MPoint) + 1);
    oversizedNode.releaseWorkspace(oversizedPool, std::move(oversized));
    assert(oversizedPool->idle.empty());

    try {
        auto returnWorkspace = [pool](SkirtCollideSurfaceWorkspace* p) {
            SkirtCollideDeformer::releaseWorkspace(pool, std::unique_ptr<SkirtCollideSurfaceWorkspace>(p));
        };
        std::unique_ptr<SkirtCollideSurfaceWorkspace, decltype(returnWorkspace)> lease(
            SkirtCollideDeformer::acquireWorkspace(pool).release(), returnWorkspace);
        throw 7;
    } catch (int) {}
    assert(pool->idle.size() == 2);

    std::vector<std::thread> workers;
    for (int i = 0; i < 8; ++i)
        workers.emplace_back([pool] {
            for (int j = 0; j < 100; ++j) {
                auto item = SkirtCollideDeformer::acquireWorkspace(pool);
                item->reset(32, 64, 4, 16, 4, 8, 8, 1);
                assert(item->samples.denominators[0] == 0);
                item->samples.denominators[0] = 11;
                SkirtCollideDeformer::releaseWorkspace(pool, std::move(item));
            }
        });
    for (auto& worker : workers) worker.join();
    assert(pool->idle.size() <= 2);

    auto shortLivedNode = std::unique_ptr<SkirtCollideDeformer>(new SkirtCollideDeformer);
    auto retainedPool = shortLivedNode->workspacePool;
    auto outstanding = SkirtCollideDeformer::acquireWorkspace(retainedPool);
    shortLivedNode.reset();
    SkirtCollideDeformer::releaseWorkspace(retainedPool, std::move(outstanding));
    assert(retainedPool->idle.size() == 1);
}
'''
    with tempfile.TemporaryDirectory(prefix="colliders-workspace-") as directory:
        cpp = Path(directory) / "workspace_test.cpp"
        binary = Path(directory) / "workspace_test"
        cpp.write_text(harness)
        subprocess.run(
            [os.environ.get("CXX", "g++"), "-std=c++17", "-pthread", "-I", str(source.parent),
             str(cpp), "-o", str(binary)],
            check=True,
        )
        subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    main()
