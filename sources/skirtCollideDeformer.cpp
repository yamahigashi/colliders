#include <maya/MArrayDataHandle.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MPointArray.h>
#include <maya/MFnData.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnGenericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MItGeometry.h>
#include <maya/MItMeshVertex.h>
#include <maya/MObjectHandle.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>

#include <tbb/blocked_range.h>
#include <tbb/cache_aligned_allocator.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#if defined(_M_X64) || defined(__SSE2__)
#include <emmintrin.h>
#define YDD_SURFACE_SSE2 1
#endif
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <type_traits>
#include <vector>
#include <utility>

#include "colliderInputValidation.h"
#include "pluginIdentity.h"
#include "skirtCollideDeformer.h"
#include "skirtLegBuild.h"
#include "utils.hpp"

struct SkirtCollideSurfaceTopology
{
    struct Sample
    {
        double u, v, quadratureWeight;
        int tensorBegin, tensorEnd, rowBegin, rowEnd, betaBegin, betaEnd;
    };
    struct Tensor
    {
        int cv, row;
        double basis;
    };
    struct Row
    {
        int cv, beta;
    };
    struct Group
    {
        int mass = 0, edgeBegin = 0, edgeEnd = 0, attachmentBegin = 0, attachmentEnd = 0;
    };
    struct EdgeSource
    {
        int a, b, edge;
    };
    struct Attachment
    {
        int sample, beta;
    };
    int nu, nv, du, dv, formU, formV;
    bool closedU;
    std::vector<double> ku, kv;
    std::vector<int> canonical, groupOf, betaGroups, edgeGroups;
    std::vector<Sample> samples;
    std::vector<Tensor> tensors;
    std::vector<Row> rows;
    std::vector<Group> groups;
    std::vector<EdgeSource> edgeSources;
    std::vector<Attachment> attachments;
};

struct SkirtCollideInputSnapshot
{
    std::vector<double> values;
    std::uint64_t hash = 14695981039346656037ULL;
    bool finite = true;

    void append(double value)
    {
        values.push_back(value);
        finite = finite && std::isfinite(value);
        // Equal floating-point values, including signed zero, must hash equally.
        hash ^= static_cast<std::uint64_t>(std::hash<double>{}(value));
        hash *= 1099511628211ULL;
    }

    void append(const MPoint& point)
    {
        append(point.x);
        append(point.y);
        append(point.z);
        append(point.w);
    }

    void append(const MMatrix& matrix)
    {
        for (unsigned int row = 0; row < 4; ++row)
            for (unsigned int column = 0; column < 4; ++column)
                append(matrix[row][column]);
    }

    bool matches(const SkirtCollideInputSnapshot& other) const
    {
        if (values.size() != other.values.size() || hash != other.hash)
            return false;
        for (std::size_t i = 0; i < values.size(); ++i)
            if (!(values[i] == other.values[i]))
                return false;
        return true;
    }
};

struct SkirtCollideEvaluation
{
    SkirtCollideInputSnapshot input;
    MPointArray output;
};

struct SkirtCollideSurfacePreparation
{
    std::shared_ptr<const SkirtCollideSurfaceTopology> topology;
    std::vector<double> weights, rowBasis, beta, denominators;
    std::vector<MPoint> restPoints, sampleRest;
    std::vector<double> restKey;
    std::vector<double> radii, edges;

    bool matches(const std::shared_ptr<const SkirtCollideSurfaceTopology>& candidate,
                 const MPointArray& cvs, const std::vector<MPoint>& rest,
                 const std::vector<double>& candidateRestKey) const
    {
        if (topology != candidate || restKey.size() != candidateRestKey.size() ||
            std::memcmp(restKey.data(), candidateRestKey.data(), restKey.size() * sizeof(double)) != 0 ||
            weights.size() != cvs.length() || restPoints.size() != rest.size())
            return false;
        for (unsigned int i = 0; i < cvs.length(); ++i)
        {
            const double weight = cvs[i].w;
            if (std::memcmp(&weights[i], &weight, sizeof(double)) != 0)
                return false;
        }
        for (std::size_t i = 0; i < rest.size(); ++i)
            if (std::memcmp(&restPoints[i].x, &rest[i].x, sizeof(double)) != 0 ||
                std::memcmp(&restPoints[i].y, &rest[i].y, sizeof(double)) != 0 ||
                std::memcmp(&restPoints[i].z, &rest[i].z, sizeof(double)) != 0 ||
                std::memcmp(&restPoints[i].w, &rest[i].w, sizeof(double)) != 0)
                return false;
        return true;
    }
};

MTypeId SkirtCollideDeformer::typeId(PluginIdentity::kSkirtCollideTypeId);

MObject SkirtCollideDeformer::attr_bellMatrix;
MObject SkirtCollideDeformer::attr_leftHipMatrix;
MObject SkirtCollideDeformer::attr_leftKneeMatrix;
MObject SkirtCollideDeformer::attr_leftHeelMatrix;
MObject SkirtCollideDeformer::attr_rightHipMatrix;
MObject SkirtCollideDeformer::attr_rightKneeMatrix;
MObject SkirtCollideDeformer::attr_rightHeelMatrix;

MObject SkirtCollideDeformer::attr_restGeometry;
MObject SkirtCollideDeformer::attr_restBellMatrix;
MObject SkirtCollideDeformer::attr_restLeftHipMatrix;
MObject SkirtCollideDeformer::attr_restLeftKneeMatrix;
MObject SkirtCollideDeformer::attr_restLeftHeelMatrix;
MObject SkirtCollideDeformer::attr_restRightHipMatrix;
MObject SkirtCollideDeformer::attr_restRightKneeMatrix;
MObject SkirtCollideDeformer::attr_restRightHeelMatrix;

MObject SkirtCollideDeformer::attr_skirtType;
MObject SkirtCollideDeformer::attr_ringScale;
MObject SkirtCollideDeformer::attr_thighRadiusX;
MObject SkirtCollideDeformer::attr_thighRadiusZ;
MObject SkirtCollideDeformer::attr_kneeRadiusX;
MObject SkirtCollideDeformer::attr_kneeRadiusZ;
MObject SkirtCollideDeformer::attr_calfRadiusX;
MObject SkirtCollideDeformer::attr_calfRadiusZ;
MObject SkirtCollideDeformer::attr_ankleRadiusX;
MObject SkirtCollideDeformer::attr_ankleRadiusZ;
MObject SkirtCollideDeformer::attr_thighPosition;
MObject SkirtCollideDeformer::attr_calfPosition;

MObject SkirtCollideDeformer::attr_leftRingAxis;
MObject SkirtCollideDeformer::attr_rightRingAxis;
MObject SkirtCollideDeformer::attr_falloff;

namespace
{
MObject attr_closedU;

double determinant(const double matrix[3][3])
{
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

MVector correction(const std::vector<LegContact> &input, const MVector &x)
{
    // Reused across the many calls per evaluation; deform may run on several threads.
    static thread_local std::vector<LegContact> constraints;
    constraints.clear();
    bool allSatisfied = true;
    for (const auto& c : input)
    {
        if (c.weight == 0.0)
            continue;
        const double h = c.distance - c.normal * legVector(x);
        allSatisfied = allSatisfied && h <= 0.0;
        constraints.push_back({c.normal, 1.25 * h, 4.0 * c.weight});
    }
    if (allSatisfied)
        return MVector(0.0, 0.0, 0.0);
    MVector best(0.0, 0.0, 0.0), bestAny(0.0, 0.0, 0.0);
    double bestObjective = std::numeric_limits<double>::infinity();
    double bestAnyObjective = std::numeric_limits<double>::infinity();
    const unsigned int count = static_cast<unsigned int>(constraints.size());
    for (unsigned int mask = 0; mask < (1u << count); mask++)
    {
        double matrix[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        double rhs[3] = {};
        for (unsigned int i = 0; i < count; i++)
        {
            if ((mask & (1u << i)) == 0)
                continue;
            const LegContact &constraint = constraints[i];
            for (int row = 0; row < 3; row++)
            {
                rhs[row] += constraint.weight * constraint.distance * constraint.normal[row];
                for (int column = 0; column < 3; column++)
                    matrix[row][column] += constraint.weight * constraint.normal[row] * constraint.normal[column];
            }
        }
        MVector delta(0.0, 0.0, 0.0);
        if (mask != 0)
        {
            const double det = determinant(matrix);
            for (int column = 0; column < 3; column++)
            {
                double replaced[3][3];
                for (int row = 0; row < 3; row++)
                    for (int j = 0; j < 3; j++)
                        replaced[row][j] = j == column ? rhs[row] : matrix[row][j];
                delta[column] = determinant(replaced) / det;
            }
        }
        bool consistent = true;
        double objective = delta * delta;
        for (unsigned int i = 0; i < count; i++)
        {
            const double projected = constraints[i].normal * legVector(delta);
            const double distance = constraints[i].distance;
            if ((mask & (1u << i)) != 0 ? projected > distance + 1.25e-9 : projected < distance - 1.25e-9)
                consistent = false;
            const double residual = std::max(0.0, distance - projected);
            objective += constraints[i].weight * residual * residual;
        }
        if (objective < bestAnyObjective)
        {
            bestAny = delta;
            bestAnyObjective = objective;
        }
        if (consistent && objective < bestObjective)
        {
            best = delta;
            bestObjective = objective;
        }
    }
    return std::isfinite(bestObjective) ? best : bestAny;
}

struct PointState
{
    int index;
    MPoint point, restPoint;
    double q = 0.0, radius = 0.0;
    bool finite = false;
    MVector delta = MVector(0.0, 0.0, 0.0), u = MVector(0.0, 0.0, 0.0);
    std::vector<LegContact> constraints;
    std::map<int, double> edges;
};

void warnUnsupportedGeometry(const MObject& node)
{
    // Handles distinguish node lifetimes without storing state in geometry data.
    static std::mutex mutex;
    static std::vector<MObjectHandle> warned;
    const std::lock_guard<std::mutex> lock(mutex);
    warned.erase(std::remove_if(warned.begin(), warned.end(),
                                [](const MObjectHandle& h) { return !h.isValid() || !h.isAlive(); }),
                 warned.end());
    for (const auto& h : warned)
        if (h.object() == node)
            return;
    warned.emplace_back(node);
    MGlobal::displayWarning(MFnDependencyNode(node).name() +
                            ": collision coupling supports only mesh and nurbsSurface geometry.");
}

MStatus buildAdjacency(const MDataHandle& geometry, const MObject& node, std::map<int, PointState>& points)
{
    const auto edge = [&](int i, int j) {
        auto a = points.find(i), b = points.find(j);
        if (a != points.end() && b != points.end() && a->second.finite && b->second.finite)
        {
            a->second.edges.emplace(j, 0.0);
            b->second.edges.emplace(i, 0.0);
        }
    };
    MStatus status;
    if (geometry.type() == MFnData::kMesh)
    {
        MObject mesh = geometry.asMesh();
        MItMeshVertex vertex(mesh, &status);
        if (!status)
        {
            warnUnsupportedGeometry(node);
            return MS::kSuccess;
        }
        for (; !vertex.isDone(); vertex.next())
        {
            MIntArray neighbours;
            status = vertex.getConnectedVertices(neighbours);
            if (!status)
                continue;
            for (unsigned int j = 0; j < neighbours.length(); ++j)
                edge(vertex.index(), neighbours[j]);
        }
    }
    else
        warnUnsupportedGeometry(node);
    return MS::kSuccess;
}

void couple(std::map<int, PointState> &points, const std::vector<LegSegment> &rest, double kappa)
{
    for (auto& entry : points)
        if (entry.second.finite)
            entry.second.radius = legLocalRadius(rest, legPoint(entry.second.restPoint));
    for (auto& entry : points)
    {
        PointState& p = entry.second;
        for (auto& edge : p.edges)
        {
            const PointState& other = points.at(edge.first);
            const double scale = kappa * (p.radius + other.radius) / 2.0;
            const MVector offset = p.restPoint - other.restPoint;
            edge.second = kappa > 0.0 && !rest.empty() ? scale * scale / (std::max)(offset * offset, 1e-12) : 0.0;
        }
    }
    for (int iteration = 0; iteration < 30; ++iteration)
    {
        for (auto& entry : points)
        {
            PointState& p = entry.second;
            if (p.q == 0.0)
                continue;
            MVector numerator(0.0, 0.0, 0.0);
            double total = 0.0;
            for (const auto& edge : p.edges)
            {
                total += edge.second;
                numerator += edge.second * points.at(edge.first).u;
            }
            const MVector x = numerator * (1.0 / (1.0 + total));
            p.delta = x + correction(p.constraints, x);
            p.u = p.q * p.delta;
        }
    }
    for (auto& entry : points)
    {
        PointState& p = entry.second;
        if (p.q > 0.0)
        {
            p.delta += correction(p.constraints, p.delta);
            p.u = p.q * p.delta;
        }
    }
}


} // namespace

struct SurfaceContact
{
    LegVec3 normal;
    double distance = 0.0, weight = 0.0, cache = 0.0;
};

struct SurfaceSamples
{
    void reset(std::size_t count)
    {
        points.resize(count);
        rests.resize(count);
        denominators.resize(count);
        std::fill(denominators.begin(), denominators.end(), 0.0);
        contacts.resize(count * 2);
    }

    std::size_t size() const { return points.size(); }
    SurfaceContact& contact(std::size_t sample, int side) { return contacts[sample * 2 + side]; }
    const SurfaceContact& contact(std::size_t sample, int side) const { return contacts[sample * 2 + side]; }

    std::vector<LegVec3> points, rests;
    std::vector<double> denominators;
    std::vector<SurfaceContact, tbb::cache_aligned_allocator<SurfaceContact>> contacts;
};

struct SurfaceGroup
{
    double paint = 1.0, diagonal = 0.0;
    LegVec3 value;
};

struct SurfacePoint
{
    LegVec3 point, restPoint;
    double q = 0.0, radius = 0.0, weight = 0.0;
    bool finite = false;
};

struct SurfaceTerms
{
    void resize(std::size_t capacity)
    {
        normalX.resize(capacity);
        normalY.resize(capacity);
        normalZ.resize(capacity);
        target.resize(capacity);
        weight.resize(capacity);
        basis.resize(capacity);
        contactIndex.resize(capacity);
        active.resize(capacity);
    }

    std::vector<double> normalX, normalY, normalZ, target, weight, basis;
    std::vector<std::size_t> contactIndex;
    std::vector<char> active;

    double dot(std::size_t i, const LegVec3& value) const
    {
        return normalX[i] * value.x + normalY[i] * value.y + normalZ[i] * value.z;
    }

    void dots(const LegVec3& value, std::size_t begin, std::size_t count, double* result) const
    {
        std::size_t j = 0;
#ifdef YDD_SURFACE_SSE2
        const __m128d x = _mm_set1_pd(value.x), y = _mm_set1_pd(value.y), z = _mm_set1_pd(value.z);
        for (; j + 1 < count; j += 2)
        {
            const std::size_t i = begin + j;
            const __m128d product = _mm_add_pd(
                _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(normalX.data() + i), x),
                           _mm_mul_pd(_mm_loadu_pd(normalY.data() + i), y)),
                _mm_mul_pd(_mm_loadu_pd(normalZ.data() + i), z));
            _mm_storeu_pd(result + j, product);
        }
#endif
        for (; j < count; ++j)
            result[j] = dot(begin + j, value);
    }

    void residuals(const LegVec3& value, std::size_t begin, std::size_t count,
                   double* result, bool positivePart) const
    {
        std::size_t j = 0;
#ifdef YDD_SURFACE_SSE2
        const __m128d x = _mm_set1_pd(value.x), y = _mm_set1_pd(value.y), z = _mm_set1_pd(value.z);
        const __m128d zero = _mm_setzero_pd();
        for (; j + 1 < count; j += 2)
        {
            const std::size_t i = begin + j;
            const __m128d product = _mm_add_pd(
                _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(normalX.data() + i), x),
                           _mm_mul_pd(_mm_loadu_pd(normalY.data() + i), y)),
                _mm_mul_pd(_mm_loadu_pd(normalZ.data() + i), z));
            const __m128d residual = _mm_sub_pd(_mm_loadu_pd(target.data() + i),
                                                 _mm_mul_pd(_mm_loadu_pd(basis.data() + i), product));
            _mm_storeu_pd(result + j, positivePart ? _mm_and_pd(_mm_cmpgt_pd(residual, zero), residual) : residual);
        }
#endif
        for (; j < count; ++j)
        {
            const std::size_t i = begin + j;
            const double residual = target[i] - basis[i] * dot(i, value);
            result[j] = positivePart ? (std::max)(0.0, residual) : residual;
        }
    }
};

constexpr std::size_t kSurfaceChunkSize = 32;

struct alignas(64) SurfaceChunkScratch
{
    LegScratch value;
    std::array<std::size_t, kSurfaceChunkSize> deferredSamples;
    std::size_t deferredCount = 0;
    bool invalid = false;
};

struct SkirtCollideSurfaceWorkspace
{
    std::vector<MPoint> restPoints;
    std::vector<char> restValid, members;
    std::vector<int> membership;
    std::vector<double> paint;
    std::vector<std::pair<unsigned int, double>> rawPaint;
    LegCapsuleColumns capsuleColumns[2];
    std::vector<SurfacePoint> points;
    std::vector<SurfaceGroup> groups;
    SurfaceSamples samples;
    SurfaceTerms terms;
    std::vector<double> edges, rows, beta, termScratch;
    std::vector<std::size_t> groupCounts;
    std::vector<std::array<LegSampleTransport, 2>> transports;
    std::vector<SurfaceChunkScratch, tbb::cache_aligned_allocator<SurfaceChunkScratch>> chunkScratch;
    std::vector<std::array<double, 2>> capsuleGaps;
    std::vector<std::array<char, 2>> deferred;

    void reset(std::size_t sampleCount, std::size_t slotCount, std::size_t groupCount,
               std::size_t maxGroupSlots, std::size_t edgeCount, std::size_t rowCount,
               std::size_t betaCount, std::size_t chunkCount)
    {
        samples.reset(sampleCount);
        terms.resize(slotCount);
        edges.resize(edgeCount);
        rows.resize(rowCount);
        std::fill(rows.begin(), rows.end(), 0.0);
        beta.resize(betaCount);
        std::fill(beta.begin(), beta.end(), 0.0);
        termScratch.resize(maxGroupSlots);
        groupCounts.resize(groupCount);
        groups.resize(groupCount);
        std::fill(groups.begin(), groups.end(), SurfaceGroup{});
        transports.resize(sampleCount);
        std::fill(transports.begin(), transports.end(), std::array<LegSampleTransport, 2>{});
        if (chunkScratch.size() < chunkCount)
            chunkScratch.resize(chunkCount);
        capsuleGaps.resize(sampleCount);
        deferred.resize(sampleCount);
    }
};

struct SkirtCollideSurfaceWorkspacePool
{
    std::mutex mutex;
    std::vector<std::unique_ptr<SkirtCollideSurfaceWorkspace>> idle;
};

namespace
{

void warnInvalidBasis(const MObject& node)
{
    static std::mutex mutex;
    static std::vector<MObjectHandle> warned;
    const std::lock_guard<std::mutex> lock(mutex);
    warned.erase(std::remove_if(warned.begin(), warned.end(),
                                [](const MObjectHandle& h) { return !h.isValid() || !h.isAlive(); }),
                 warned.end());
    for (const auto& h : warned)
        if (h.object() == node)
            return;
    warned.emplace_back(node);
    MGlobal::displayWarning(MFnDependencyNode(node).name() + ": invalid NURBS sample basis; geometry unchanged.");
}

std::vector<double> fullKnots(const MDoubleArray& knots)
{
    std::vector<double> result;
    if (knots.length() == 0)
        return result;
    result.reserve(knots.length() + 2);
    result.push_back(knots[0]);
    for (unsigned int i = 0; i < knots.length(); ++i)
        result.push_back(knots[i]);
    result.push_back(knots[knots.length() - 1]);
    return result;
}

std::vector<double> spanSites(const std::vector<double>& knots, int degree, int count, int divisions)
{
    std::vector<double> result;
    for (int i = degree; i < count; ++i)
        if (knots[i + 1] > knots[i])
            for (int k = 1; k <= divisions; ++k)
                result.push_back(knots[i] + (knots[i + 1] - knots[i]) * (k - 0.5) / divisions);
    return result;
}

struct QuadratureSite
{
    double value, weight;
};

std::vector<QuadratureSite> spanGaussSites(const std::vector<double>& knots, int degree, int count)
{
    std::vector<QuadratureSite> result;
    const double offset = 0.5 * std::sqrt(3.0 / 5.0);
    for (int i = degree; i < count; ++i)
        if (knots[i + 1] > knots[i])
        {
            const double start = knots[i], width = knots[i + 1] - start;
            result.push_back({start + width * (0.5 - offset), 10.0 / 9.0});
            result.push_back({start + width * 0.5, 16.0 / 9.0});
            result.push_back({start + width * (0.5 + offset), 10.0 / 9.0});
        }
    return result;
}

int splineBasis(const std::vector<double>& knots, int degree, int count, double t, std::vector<double>& values)
{
    int span = count - 1;
    if (t < knots[count])
        for (int i = degree; i < count; ++i)
            if (knots[i] <= t && t < knots[i + 1])
            {
                span = i;
                break;
            }
    std::fill(values.begin(), values.end(), 0.0);
    values[0] = 1.0;
    for (int j = 1; j <= degree; ++j)
    {
        double saved = 0.0;
        for (int r = 0; r < j; ++r)
        {
            const double right = knots[span + r + 1] - t;
            const double left = t - knots[span + 1 - j + r];
            const double temp = right + left != 0.0 ? values[r] / (right + left) : 0.0;
            values[r] = saved + right * temp;
            saved = left * temp;
        }
        values[j] = saved;
    }
    return span - degree;
}

bool buildSurfaceTopology(SkirtCollideSurfaceTopology& topology)
{
    auto& t = topology;
    const int nu = t.nu, nv = t.nv, du = t.du, dv = t.dv, count = nu * nv;
    const bool pu = t.formU == MFnNurbsSurface::kPeriodic, pv = t.formV == MFnNurbsSurface::kPeriodic;
    const int uu = pu ? nu - du : nu, vv = pv ? nv - dv : nv;
    if (du < 1 || dv < 1 || nu <= du || nv <= dv || t.ku.size() != static_cast<std::size_t>(nu + du + 1) ||
        t.kv.size() != static_cast<std::size_t>(nv + dv + 1))
        return false;
    t.canonical.resize(count);
    t.groupOf.resize(count);
    std::vector<int> parent(count), groupIndex(count, -1);
    for (int i = 0; i < count; ++i)
    {
        t.canonical[i] = (i / nv % uu) * nv + i % nv % vv;
        parent[i] = i;
    }
    const auto root = [&](int i) {
        while (parent[i] != i)
        {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    if (t.closedU && !pu)
        for (int v = 0; v < vv; ++v)
        {
            const int a = root(v), b = root((nu - 1) * nv + v);
            parent[(std::max)(a, b)] = (std::min)(a, b);
        }
    for (int i = 0; i < count; ++i)
        if (t.canonical[i] == i)
        {
            const int r = root(i);
            if (groupIndex[r] == -1)
            {
                groupIndex[r] = static_cast<int>(t.groups.size());
                t.groups.emplace_back();
            }
            const int g = groupIndex[r];
            t.groupOf[i] = g;
            ++t.groups[g].mass;
        }
    for (int i = 0; i < count; ++i)
        t.groupOf[i] = t.groupOf[t.canonical[i]];

    std::vector<std::pair<int, int>> edges;
    edges.reserve(count * 4);
    t.edgeSources.reserve(count * 4);
    for (int u = 0; u < uu; ++u)
        for (int v = 0; v < vv; ++v)
        {
            const int i = u * nv + v, g = t.groupOf[i];
            const int offsets[4][2] = {{u - 1, v}, {u + 1, v}, {u, v - 1}, {u, v + 1}};
            int neighbours[4], size = 0;
            for (const auto& offset : offsets)
            {
                int a = offset[0], b = offset[1];
                if (pu)
                    a = (a + uu) % uu;
                if (pv)
                    b = (b + vv) % vv;
                if (a >= 0 && a < uu && b >= 0 && b < vv)
                    neighbours[size++] = a * nv + b;
            }
            std::sort(neighbours, neighbours + size);
            size = static_cast<int>(std::unique(neighbours, neighbours + size) - neighbours);
            for (int n = 0; n < size; ++n)
            {
                const int j = neighbours[n], h = t.groupOf[j];
                if (g != h)
                {
                    edges.emplace_back(g, h);
                    t.edgeSources.push_back({i, j, 0});
                }
            }
        }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    for (const auto& edge : edges)
    {
        ++t.groups[edge.first].edgeEnd;
        t.edgeGroups.push_back(edge.second);
    }
    int edgeOffset = 0;
    for (auto& group : t.groups)
    {
        group.edgeBegin = edgeOffset;
        group.edgeEnd += edgeOffset;
        edgeOffset = group.edgeEnd;
    }
    for (auto& source : t.edgeSources)
        source.edge = static_cast<int>(
            std::lower_bound(edges.begin(), edges.end(), std::make_pair(t.groupOf[source.a], t.groupOf[source.b])) -
            edges.begin());

    auto us = spanSites(t.ku, du, nu, 2);
    const auto vs = spanGaussSites(t.kv, dv, nv);
    if (!pu)
    {
        us.push_back(t.ku[du]);
        us.push_back(t.ku[nu]);
    }
    const int support = (du + 1) * (dv + 1);
    t.samples.reserve(us.size() * vs.size());
    t.tensors.reserve(us.size() * vs.size() * support);
    t.rows.reserve(us.size() * vs.size() * support);
    t.betaGroups.reserve(us.size() * vs.size() * support);
    // Reuse degree-bounded scratch arrays; no sample owns a heap allocation.
    std::vector<double> bu(du + 1), bv(dv + 1);
    std::vector<int> rowIndices(support), betaIndices(support);
    for (double u : us)
    {
        const int firstU = splineBasis(t.ku, du, nu, u, bu);
        for (const auto& vSite : vs)
        {
            const int firstV = splineBasis(t.kv, dv, nv, vSite.value, bv);
            SkirtCollideSurfaceTopology::Sample sample;
            sample.u = u;
            sample.v = vSite.value;
            sample.quadratureWeight = vSite.weight;
            sample.tensorBegin = static_cast<int>(t.tensors.size());
            sample.rowBegin = static_cast<int>(t.rows.size());
            sample.betaBegin = static_cast<int>(t.betaGroups.size());
            int size = 0;
            for (int a = 0; a <= du; ++a)
                if (bu[a] != 0.0)
                    for (int b = 0; b <= dv; ++b)
                        if (bv[b] != 0.0)
                        {
                            const int cv = (firstU + a) * nv + firstV + b;
                            t.tensors.push_back({cv, 0, bu[a] * bv[b]});
                            rowIndices[size++] = t.canonical[cv];
                        }
            std::sort(rowIndices.begin(), rowIndices.begin() + size);
            size = static_cast<int>(std::unique(rowIndices.begin(), rowIndices.begin() + size) - rowIndices.begin());
            for (int i = 0; i < size; ++i)
                betaIndices[i] = t.groupOf[rowIndices[i]];
            std::sort(betaIndices.begin(), betaIndices.begin() + size);
            const int betaSize =
                static_cast<int>(std::unique(betaIndices.begin(), betaIndices.begin() + size) - betaIndices.begin());
            for (int i = 0; i < betaSize; ++i)
            {
                t.betaGroups.push_back(betaIndices[i]);
                ++t.groups[betaIndices[i]].attachmentEnd;
            }
            for (int i = 0; i < size; ++i)
            {
                const int beta = static_cast<int>(
                    std::lower_bound(betaIndices.begin(), betaIndices.begin() + betaSize, t.groupOf[rowIndices[i]]) -
                    betaIndices.begin());
                t.rows.push_back({rowIndices[i], sample.betaBegin + beta});
            }
            sample.tensorEnd = static_cast<int>(t.tensors.size());
            sample.rowEnd = static_cast<int>(t.rows.size());
            sample.betaEnd = static_cast<int>(t.betaGroups.size());
            for (int k = sample.tensorBegin; k < sample.tensorEnd; ++k)
                t.tensors[k].row =
                    sample.rowBegin + static_cast<int>(std::lower_bound(rowIndices.begin(), rowIndices.begin() + size,
                                                                        t.canonical[t.tensors[k].cv]) -
                                                       rowIndices.begin());
            t.samples.push_back(sample);
        }
    }
    int attachmentOffset = 0;
    for (auto& group : t.groups)
    {
        group.attachmentBegin = attachmentOffset;
        group.attachmentEnd += attachmentOffset;
        attachmentOffset = group.attachmentEnd;
    }
    t.attachments.resize(attachmentOffset);
    std::vector<int> next(t.groups.size());
    for (std::size_t g = 0; g < t.groups.size(); ++g)
        next[g] = t.groups[g].attachmentBegin;
    for (std::size_t s = 0; s < t.samples.size(); ++s)
        for (int b = t.samples[s].betaBegin; b < t.samples[s].betaEnd; ++b)
            t.attachments[next[t.betaGroups[b]]++] = {static_cast<int>(s), b};
    return true;
}

template<class T>
std::size_t vectorPayloadBytes(const std::vector<T>& values)
{
    return values.capacity() * sizeof(T);
}

std::size_t topologyPayloadBytes(const SkirtCollideSurfaceTopology& t)
{
    return sizeof(t) + vectorPayloadBytes(t.ku) + vectorPayloadBytes(t.kv) +
           vectorPayloadBytes(t.canonical) + vectorPayloadBytes(t.groupOf) +
           vectorPayloadBytes(t.betaGroups) + vectorPayloadBytes(t.edgeGroups) +
           vectorPayloadBytes(t.samples) + vectorPayloadBytes(t.tensors) + vectorPayloadBytes(t.rows) +
           vectorPayloadBytes(t.groups) + vectorPayloadBytes(t.edgeSources) +
           vectorPayloadBytes(t.attachments);
}

std::size_t preparationPayloadBytes(const SkirtCollideSurfacePreparation& p)
{
    return sizeof(p) + vectorPayloadBytes(p.weights) + vectorPayloadBytes(p.rowBasis) +
           vectorPayloadBytes(p.beta) + vectorPayloadBytes(p.denominators) +
           vectorPayloadBytes(p.restPoints) + vectorPayloadBytes(p.sampleRest) +
           vectorPayloadBytes(p.restKey) + vectorPayloadBytes(p.radii) + vectorPayloadBytes(p.edges);
}

bool topologyMatches(const SkirtCollideSurfaceTopology& a, const SkirtCollideSurfaceTopology& b)
{
    return a.nu == b.nu && a.nv == b.nv && a.du == b.du && a.dv == b.dv &&
           a.formU == b.formU && a.formV == b.formV && a.closedU == b.closedU &&
           a.ku == b.ku && a.kv == b.kv;
}

std::shared_ptr<const SkirtCollideSurfaceTopology> surfaceTopologyFor(
    MFnNurbsSurface& surface, bool closedU,
    SkirtSurfaceCache<SkirtCollideSurfaceTopology, SkirtCollideSurfacePreparation>& cache)
{
    MDoubleArray mayaU, mayaV;
    if (!surface.getKnotsInU(mayaU) || !surface.getKnotsInV(mayaV))
        return {};
    SkirtCollideSurfaceTopology key;
    key.nu = surface.numCVsInU();
    key.nv = surface.numCVsInV();
    key.du = surface.degreeU();
    key.dv = surface.degreeV();
    key.formU = surface.formInU();
    key.formV = surface.formInV();
    key.closedU = closedU;
    key.ku = fullKnots(mayaU);
    key.kv = fullKnots(mayaV);
    // Only degree, knots, form, CV counts and closedU invalidate topology. Rational
    // weights and all evaluation state stay outside the immutable, per-node cache.
    if (auto cached = cache.findTopology(key, topologyMatches))
        return cached;
    auto built = std::make_shared<SkirtCollideSurfaceTopology>(std::move(key));
    if (!buildSurfaceTopology(*built))
        return {};
    return cache.publishTopology(built, topologyPayloadBytes(*built), topologyMatches);
}

LegVec3 solveSurfaceGroup(const SurfaceGroup& group, const SkirtCollideSurfaceTopology::Group& layout,
                          const SkirtCollideSurfaceTopology& topology, const std::vector<double>& edges,
                          const std::vector<SurfaceGroup>& groups, SurfaceTerms& terms, std::size_t begin, std::size_t count, double* residualScratch)
{
    LegVec3 rhs0;
    for (int e = layout.edgeBegin; e < layout.edgeEnd; ++e)
        rhs0 += edges[e] * groups[topology.edgeGroups[e]].value;
    const auto objective = [&](const LegVec3& z) {
        double value = layout.mass * (z * z);
        for (int e = layout.edgeBegin; e < layout.edgeEnd; ++e)
        {
            const LegVec3 offset = z - groups[topology.edgeGroups[e]].value;
            value += edges[e] * (offset * offset);
        }
        terms.residuals(z, begin, count, residualScratch, true);
        for (std::size_t j = 0; j < count; ++j)
        {
            const std::size_t i = begin + j;
            value += terms.weight[i] * residualScratch[j] * residualScratch[j];
        }
        return value;
    };
    const auto consistent = [&](const LegVec3& z) {
        for (std::size_t j = 0; j < count; ++j)
        {
            const std::size_t i = begin + j;
            if (static_cast<bool>(terms.active[i]) != (terms.target[i] - terms.basis[i] * terms.dot(i, z) > 0.0))
                return false;
        }
        return true;
    };
    LegVec3 z = group.value;
    for (int iteration = 0; iteration < 50; ++iteration)
    {
        const double diagonal = group.diagonal;
        double matrix[3][3] = {{diagonal, 0.0, 0.0}, {0.0, diagonal, 0.0}, {0.0, 0.0, diagonal}};
        LegVec3 rhs = rhs0;
        terms.residuals(z, begin, count, residualScratch, false);
        for (std::size_t j = 0; j < count; ++j)
        {
            const std::size_t i = begin + j;
            terms.active[i] = residualScratch[j] > 0.0;
            if (!terms.active[i])
                continue;
            const double normal[3] = {terms.normalX[i], terms.normalY[i], terms.normalZ[i]};
            for (int r = 0; r < 3; ++r)
            {
                rhs[r] += terms.weight[i] * terms.basis[i] * terms.target[i] * normal[r];
                for (int c = 0; c < 3; ++c)
                    matrix[r][c] += terms.weight[i] * terms.basis[i] * terms.basis[i] * normal[r] * normal[c];
            }
        }
        const double det = determinant(matrix);
        LegVec3 candidate;
        for (int k = 0; k < 3; ++k)
        {
            double replaced[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    replaced[r][c] = c == k ? rhs[r] : matrix[r][c];
            candidate[k] = determinant(replaced) / det;
        }
        const double oldValue = objective(z);
        if (objective(candidate) >= oldValue)
        {
            const LegVec3 step = candidate - z;
            candidate = z;
            double scale = 0.5;
            for (int k = 1; k <= 30; ++k, scale *= 0.5)
            {
                const LegVec3 trial = z + scale * step;
                if (objective(trial) < oldValue)
                {
                    candidate = trial;
                    break;
                }
                if (k == 1 && step.x == 0.0 && step.y == 0.0 && step.z == 0.0)
                    break;
            }
        }
        const bool done = consistent(candidate) || (candidate - z).length() <= 1e-12 * (1.0 + z.length());
        z = candidate;
        if (done)
            break;
    }
    return z;
}

bool coupleSurface(const SkirtCollideSurfaceTopology &topology,
                   std::vector<SurfacePoint> &points, const Leg legs[2], const std::vector<LegSegment> &restCylinders,
                   const MPoint &waist, const MVector &axis, double kappa, std::vector<SurfaceGroup> &groups,
                   const SkirtCollideSurfacePreparation* prepared,
                   SkirtCollideSurfacePreparation* built, SkirtCollideSurfaceWorkspace& workspace)
{
    const auto& t = topology;
    const LegVec3 waistPosition = legPoint(waist);
    const LegVec3 axisDirection = legVector(axis);
    for (const auto& point : points)
        if (!std::isfinite(point.weight) || point.weight <= 0.0)
            return false;
    std::size_t maxGroupSlots = 0;
    for (const auto& layout : t.groups)
        maxGroupSlots = (std::max)(maxGroupSlots,
                                  static_cast<std::size_t>(layout.attachmentEnd - layout.attachmentBegin) * 2);
    constexpr std::size_t chunkSize = kSurfaceChunkSize;
    const std::size_t chunkCount = (t.samples.size() + chunkSize - 1) / chunkSize;
    workspace.reset(t.samples.size(), t.attachments.size() * 2, t.groups.size(), maxGroupSlots,
                    t.edgeGroups.size(), prepared ? 0 : t.rows.size(), t.betaGroups.size(), chunkCount);
    auto& samples = workspace.samples;
    auto& terms = workspace.terms;
    auto& edges = workspace.edges;
    auto& rows = workspace.rows;
    auto& beta = workspace.beta;
    auto& groupCounts = workspace.groupCounts;
    auto& termScratch = workspace.termScratch;
    auto& transports = workspace.transports;
    auto& chunkScratch = workspace.chunkScratch;
    auto& capsuleGaps = workspace.capsuleGaps;
    auto& deferred = workspace.deferred;
    // reset() has initialized each group's paint, diagonal and displacement.
    if (built)
        built->radii.resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        if (t.canonical[i] != static_cast<int>(i))
            continue;
        auto& point = points[i];
        auto& group = groups[t.groupOf[i]];
        group.paint = (std::min)(group.paint, point.q);
        if (point.finite)
        {
            point.radius = prepared ? prepared->radii[i] : legLocalRadius(restCylinders, point.restPoint);
            if (built)
                built->radii[i] = point.radius;
        }
    }
    if (prepared)
        edges = prepared->edges;
    else
        std::fill(edges.begin(), edges.end(), 0.0);
    if (!prepared)
        for (const auto& source : t.edgeSources)
        {
            const auto& point = points[source.a];
            const auto& other = points[source.b];
            if (!point.finite || !other.finite)
                continue;
            const double scale = kappa * (point.radius + other.radius) / 2.0;
            const LegVec3 offset = point.restPoint - other.restPoint;
            edges[source.edge] +=
                kappa > 0.0 && !restCylinders.empty() ? scale * scale / (std::max)(offset * offset, 1e-12) : 0.0;
        }
    if (built)
        built->edges = edges;
    for (std::size_t g = 0; g < groups.size(); ++g)
    {
        groups[g].diagonal = t.groups[g].mass;
        for (int e = t.groups[g].edgeBegin; e < t.groups[g].edgeEnd; ++e)
            groups[g].diagonal += edges[e];
    }
    if (prepared)
        beta = prepared->beta;
    if (built)
    {
        built->rowBasis.resize(t.rows.size());
        built->denominators.resize(t.samples.size());
        built->sampleRest.resize(t.samples.size());
    }
    const auto prepareSample = [&](std::size_t s) {
        const auto& layout = t.samples[s];
        double total = 0.0;
        if (!prepared)
        {
            for (int k = layout.tensorBegin; k < layout.tensorEnd; ++k)
            {
                const auto& coefficient = t.tensors[k];
                const double w = coefficient.basis * points[coefficient.cv].weight;
                rows[coefficient.row] += w;
                total += w;
            }
            if (!std::isfinite(total) || total <= 0.0)
                return false;
        }
        double sum = 0.0;
        LegVec3 p, r;
        for (int k = layout.rowBegin; k < layout.rowEnd; ++k)
        {
            const double b = prepared ? prepared->rowBasis[k] : rows[k] / total;
            if (!prepared)
            {
                if (!std::isfinite(b) || b < 0.0)
                    return false;
                sum += b;
                samples.denominators[s] += b * b;
                if (built)
                    built->rowBasis[k] = b;
            }
            const auto& point = points[t.rows[k].cv];
            p += b * point.point;
            if (!prepared)
            {
                r += b * point.restPoint;
                beta[t.rows[k].beta] += b;
            }
        }
        if (!prepared && std::abs(sum - 1.0) > 1e-9)
            return false;
        samples.points[s] = p;
        if (prepared)
        {
            samples.denominators[s] = prepared->denominators[s];
            samples.rests[s] = legPoint(prepared->sampleRest[s]);
        }
        else
        {
            samples.rests[s] = r;
            if (built)
                built->denominators[s] = samples.denominators[s];
        }
        return true;
    };
    bool mixedPaint = false;
    for (const auto& group : groups)
        if (group.paint > 0.0 && group.paint < 1.0)
            mixedPaint = true;
    const LegContactGeometry contactGeometry[2] = {legPrepareContactGeometry(legs[0], kappa),
                                                    legPrepareContactGeometry(legs[1], kappa)};
    workspace.capsuleColumns[0].prepare(contactGeometry[0]);
    workspace.capsuleColumns[1].prepare(contactGeometry[1]);
    bool anyDeferred = false;
    for (auto& gaps : capsuleGaps)
    {
        gaps[0] = gaps[1] = 0.0;
    }
    for (auto& flags : deferred)
    {
        flags[0] = flags[1] = 0;
    }
    const auto processSample = [&](std::size_t s, SurfaceChunkScratch& chunkData)
    {
        if (!prepareSample(s))
        {
            chunkData.invalid = true;
            return;
        }
        LegScratch& scratch = chunkData.value;
        const LegVec3 position = samples.points[s];
        const LegVec3 offset = samples.rests[s] - waistPosition;
        LegVec3 radial = offset - axisDirection * (offset * axisDirection);
        const bool valid = finiteVector(position) && finiteVector(samples.rests[s]) && radial.length() >= 1e-8;
        if (valid)
            radial.normalize();
        if (valid)
            for (int side = 0; side < 2; ++side)
                legSampleTransport(legs[side], samples.rests[s], radial, transports[s][side]);
        for (int side = 0; side < 2; ++side)
        {
            SurfaceContact& contact = samples.contact(s, side);
            contact = SurfaceContact();
            if (!valid)
                continue;
            bool inside = false;
            for (const auto& segment : legs[side].segments)
                if (legLocalPhi(segment.current, position) <= 0.0)
                {
                    inside = true;
                    break;
                }
            LegContact constraint;
            if (inside)
            {
                if (legConstraint(legs[side], position, transports[s][side], kappa, constraint,
                                  scratch, contactGeometry[side]) &&
                    constraint.weight != 0.0)
                {
                    contact.normal = constraint.normal;
                    contact.cache = constraint.normal * LegVec3();
                    contact.distance = constraint.distance + contact.cache;
                    contact.weight = 4.0 * t.samples[s].quadratureWeight * constraint.weight / samples.denominators[s];
                }
                continue;
            }
            const double gap = mixedPaint ? 0.0 : legCapsuleDistance(position, workspace.capsuleColumns[side]);
            capsuleGaps[s][side] = gap;
            const double scaled = gap * (1.0 - 1e-9);
            if (mixedPaint || !std::isfinite(gap) || !(1.25 * scaled > 0.0))
            {
                if (legConstraint(legs[side], position, transports[s][side], kappa, constraint,
                                  scratch, contactGeometry[side]) &&
                    constraint.weight != 0.0)
                {
                    contact.normal = constraint.normal;
                    contact.cache = constraint.normal * LegVec3();
                    contact.distance = constraint.distance + contact.cache;
                    contact.weight = 4.0 * t.samples[s].quadratureWeight * constraint.weight / samples.denominators[s];
                }
                continue;
            }
            deferred[s][side] = 1;
        }
        if (deferred[s][0] || deferred[s][1])
            chunkData.deferredSamples[chunkData.deferredCount++] = s;
    };
    const auto processChunk = [&](std::size_t chunk) {
        SurfaceChunkScratch& chunkData = chunkScratch[chunk];
        chunkData.deferredCount = 0;
        chunkData.invalid = false;
        const std::size_t end = (std::min)((chunk + 1) * chunkSize, samples.size());
        for (std::size_t s = chunk * chunkSize; s < end; ++s)
        {
            processSample(s, chunkData);
            if (chunkData.invalid)
                break;
        }
    };
    if (samples.size() < 64)
        for (std::size_t chunk = 0; chunk < chunkCount; ++chunk)
            processChunk(chunk);
    else
    {
        tbb::task_group_context context(tbb::task_group_context::isolated,
                                        tbb::task_group_context::default_traits |
                                            tbb::task_group_context::fp_settings);
        context.capture_fp_settings();
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, chunkCount, 1),
                          [&](const tbb::blocked_range<std::size_t>& range) {
                              for (std::size_t chunk = range.begin(); chunk < range.end(); ++chunk)
                                  processChunk(chunk);
                          }, context);
    }
    for (std::size_t chunk = 0; chunk < chunkCount; ++chunk)
        if (chunkScratch[chunk].invalid)
            return false;
    if (built)
    {
        if (!prepared)
            for (std::size_t s = 0; s < samples.size(); ++s)
            {
                const LegVec3& rest = samples.rests[s];
                built->sampleRest[s] = MPoint(rest.x, rest.y, rest.z, 1.0);
            }
        built->beta = beta;
    }
    for (std::size_t chunk = 0; chunk < chunkCount; ++chunk)
        anyDeferred = anyDeferred || chunkScratch[chunk].deferredCount != 0;
    double bound = 0.0;
    if (!anyDeferred)
        bound = std::numeric_limits<double>::quiet_NaN();
    else
    {
    // Initial group values are zero; a nonfinite edge used to make this bound NaN.
    for (double edge : edges)
        if (!std::isfinite(edge))
        {
            bound = std::numeric_limits<double>::quiet_NaN();
            break;
        }
    for (std::size_t s = 0; s < samples.size(); ++s)
        for (int side = 0; side < 2; ++side)
        {
            const SurfaceContact& contact = samples.contact(s, side);
            if (contact.weight == 0.0)
                continue;
            const double hinge = (std::max)(0.0, 1.25 * contact.distance - contact.cache);
            bound += contact.weight * hinge * hinge;
        }
    }
    const auto processDeferredChunk = [&](std::size_t chunk) {
        SurfaceChunkScratch& chunkData = chunkScratch[chunk];
        LegScratch& scratch = chunkData.value;
        for (std::size_t item = 0; item < chunkData.deferredCount; ++item)
        {
            const std::size_t s = chunkData.deferredSamples[item];
            double support = 0.0;
            bool supportFinite = std::isfinite(bound);
            for (int b = t.samples[s].betaBegin; b < t.samples[s].betaEnd; ++b)
            {
                const int g = t.betaGroups[b];
                const double mass = static_cast<double>(t.groups[g].mass);
                if (!(mass > 0.0) || !std::isfinite(beta[b]))
                {
                    supportFinite = false;
                    break;
                }
                support += beta[b] * beta[b] / mass;
            }
            const double radius = supportFinite && std::isfinite(support) && support >= 0.0
                                      ? std::sqrt(bound * support)
                                      : std::numeric_limits<double>::quiet_NaN();
            for (int side = 0; side < 2; ++side)
            {
                if (!deferred[s][side])
                    continue;
                deferred[s][side] = 0;
                const LegVec3 position = samples.points[s];
                const double scaled = capsuleGaps[s][side] * (1.0 - 1e-9);
                if (mixedPaint || !std::isfinite(scaled) || !std::isfinite(radius) ||
                    !(1.25 * scaled > radius))
                {
                    LegContact constraint;
                    if (legConstraint(legs[side], position, transports[s][side], kappa, constraint,
                                      scratch, contactGeometry[side]) && constraint.weight != 0.0)
                    {
                        SurfaceContact& contact = samples.contact(s, side);
                        contact.normal = constraint.normal;
                        contact.cache = constraint.normal * LegVec3();
                        contact.distance = constraint.distance + contact.cache;
                        contact.weight = 4.0 * t.samples[s].quadratureWeight * constraint.weight / samples.denominators[s];
                    }
                }
            }
        }
    };
    if (anyDeferred)
    {
        if (samples.size() < 64)
            for (std::size_t chunk = 0; chunk < chunkCount; ++chunk)
                processDeferredChunk(chunk);
        else
        {
            tbb::task_group_context context(tbb::task_group_context::isolated,
                                            tbb::task_group_context::default_traits |
                                                tbb::task_group_context::fp_settings);
            context.capture_fp_settings();
            tbb::parallel_for(tbb::blocked_range<std::size_t>(0, chunkCount, 1),
                              [&](const tbb::blocked_range<std::size_t>& range) {
                                  for (std::size_t chunk = range.begin(); chunk < range.end(); ++chunk)
                                      processDeferredChunk(chunk);
                              }, context);
        }
    }
    const auto fillGroup = [&](std::size_t g) {
        const auto& layout = t.groups[g];
        const std::size_t begin = static_cast<std::size_t>(layout.attachmentBegin) * 2;
        std::size_t count = 0;
        if (groups[g].paint != 0.0)
            for (int a = layout.attachmentBegin; a < layout.attachmentEnd; ++a)
            {
                const auto& attachment = t.attachments[a];
                for (int side = 0; side < 2; ++side)
                {
                    const std::size_t contactIndex = static_cast<std::size_t>(attachment.sample) * 2 + side;
                    const SurfaceContact& contact = samples.contacts[contactIndex];
                    if (contact.weight == 0.0)
                        continue;
                    const std::size_t i = begin + count++;
                    terms.normalX[i] = contact.normal.x;
                    terms.normalY[i] = contact.normal.y;
                    terms.normalZ[i] = contact.normal.z;
                    terms.basis[i] = beta[attachment.beta];
                    terms.weight[i] = contact.weight;
                    terms.contactIndex[i] = contactIndex;
                }
            }
        groupCounts[g] = count;
    };
    for (std::size_t g = 0; g < groups.size(); ++g)
        fillGroup(g);
    for (int sweep = 0; sweep < 3; ++sweep)
        for (std::size_t g = 0; g < groups.size(); ++g)
        {
            auto& group = groups[g];
            const auto& layout = t.groups[g];
            if (group.paint == 0.0)
                continue;
            const std::size_t begin = static_cast<std::size_t>(layout.attachmentBegin) * 2;
            const std::size_t count = groupCounts[g];
            terms.dots(group.value, begin, count, termScratch.data());
            for (std::size_t j = 0; j < count; ++j)
            {
                const std::size_t i = begin + j;
                const SurfaceContact& contact = samples.contacts[terms.contactIndex[i]];
                terms.target[i] =
                    1.25 * contact.distance - (contact.cache - terms.basis[i] * termScratch[j]);
            }
            const LegVec3 value =
                group.paint * solveSurfaceGroup(group, layout, t, edges, groups,
                                                terms, begin, count, termScratch.data());
            const LegVec3 change = value - group.value;
            group.value = value;
            terms.dots(change, begin, count, termScratch.data());
            for (std::size_t j = 0; j < count; ++j)
            {
                const std::size_t i = begin + j;
                SurfaceContact& contact = samples.contacts[terms.contactIndex[i]];
                contact.cache += terms.basis[i] * termScratch[j];
            }
        }
    return true;
}

} // namespace

SkirtCollideDeformer::SkirtCollideDeformer()
    : MPxDeformerNode(), restGeometryWarningIssued(false),
      workspacePool(std::make_shared<SkirtCollideSurfaceWorkspacePool>())
{
}

SkirtCollideDeformer::~SkirtCollideDeformer() = default;

std::unique_ptr<SkirtCollideSurfaceWorkspace> SkirtCollideDeformer::acquireWorkspace(
    const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>& pool)
{
    {
        std::lock_guard<std::mutex> lock(pool->mutex);
        if (!pool->idle.empty())
        {
            auto workspace = std::move(pool->idle.back());
            pool->idle.pop_back();
            return workspace;
        }
    }
    return std::unique_ptr<SkirtCollideSurfaceWorkspace>(new SkirtCollideSurfaceWorkspace);
}

void SkirtCollideDeformer::releaseWorkspace(
    const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>& pool,
    std::unique_ptr<SkirtCollideSurfaceWorkspace> workspace) noexcept
{
    constexpr std::size_t maxIdle = 2;
    constexpr std::size_t maxRetainedBytes = 32 * 1024 * 1024;
    const auto bytes = [](const auto& values) {
        return values.capacity() * sizeof(typename std::decay<decltype(values)>::type::value_type);
    };
    const auto scratchBytes = [&](const LegScratch& scratch) {
        std::size_t total = bytes(scratch.raySites) + bytes(scratch.intervals);
        for (const auto& level : scratch.polynomial)
            total += bytes(level.roots) + bytes(level.sites);
        return total;
    };
    std::size_t retained = bytes(workspace->restPoints) + bytes(workspace->restValid) +
        bytes(workspace->members) + bytes(workspace->membership) +
        bytes(workspace->paint) + bytes(workspace->rawPaint) +
        bytes(workspace->points) + bytes(workspace->groups) +
        bytes(workspace->samples.points) + bytes(workspace->samples.rests) +
        bytes(workspace->samples.denominators) + bytes(workspace->samples.contacts) +
        bytes(workspace->terms.normalX) + bytes(workspace->terms.normalY) +
        bytes(workspace->terms.normalZ) + bytes(workspace->terms.target) +
        bytes(workspace->terms.weight) + bytes(workspace->terms.basis) +
        bytes(workspace->terms.contactIndex) + bytes(workspace->terms.active) +
        bytes(workspace->edges) + bytes(workspace->rows) + bytes(workspace->beta) +
        bytes(workspace->termScratch) +
        bytes(workspace->groupCounts) + bytes(workspace->transports) +
        bytes(workspace->chunkScratch) + bytes(workspace->capsuleGaps) + bytes(workspace->deferred);
    for (const auto& chunk : workspace->chunkScratch)
        retained += scratchBytes(chunk.value);
    for (const auto& columns : workspace->capsuleColumns)
        retained += bytes(columns.p0x) + bytes(columns.p0y) + bytes(columns.p0z) +
                    bytes(columns.spanx) + bytes(columns.spany) + bytes(columns.spanz) +
                    bytes(columns.length2) + bytes(columns.radius);
    if (retained > maxRetainedBytes)
        return;
    try
    {
        std::lock_guard<std::mutex> lock(pool->mutex);
        if (pool->idle.size() < maxIdle)
            pool->idle.push_back(std::move(workspace));
    }
    catch (...)
    {
        // A failed cache insertion must not affect evaluation or exception unwinding.
    }
}

MStatus SkirtCollideDeformer::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnGenericAttribute gAttr;
    MStatus stat;

    attr_closedU = nAttr.create("closedU", "clu", MFnNumericData::kBoolean, false, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setStorable(true);
    nAttr.setReadable(true);
    nAttr.setWritable(true);
    nAttr.setKeyable(false);
    stat = addAttribute(attr_closedU);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_bellMatrix);

    attr_leftHipMatrix = mAttr.create("leftHipMatrix", "leftHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHipMatrix);

    attr_leftKneeMatrix = mAttr.create("leftKneeMatrix", "leftKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftKneeMatrix);

    attr_leftHeelMatrix = mAttr.create("leftHeelMatrix", "leftHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHeelMatrix);

    attr_rightHipMatrix = mAttr.create("rightHipMatrix", "rightHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHipMatrix);

    attr_rightKneeMatrix = mAttr.create("rightKneeMatrix", "rightKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightKneeMatrix);

    attr_rightHeelMatrix = mAttr.create("rightHeelMatrix", "rightHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHeelMatrix);

    attr_restGeometry = gAttr.create("restGeometry", "restGeometry", &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = gAttr.addDataAccept(MFnData::kMesh);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = gAttr.addDataAccept(MFnData::kNurbsSurface);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    gAttr.setHidden(true);
    gAttr.setStorable(false);
    gAttr.setReadable(false);
    addAttribute(attr_restGeometry);

    attr_restBellMatrix = mAttr.create("restBellMatrix", "restBellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restBellMatrix);

    attr_restLeftHipMatrix = mAttr.create("restLeftHipMatrix", "restLeftHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restLeftHipMatrix);

    attr_restLeftKneeMatrix = mAttr.create("restLeftKneeMatrix", "restLeftKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restLeftKneeMatrix);

    attr_restLeftHeelMatrix = mAttr.create("restLeftHeelMatrix", "restLeftHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restLeftHeelMatrix);

    attr_restRightHipMatrix = mAttr.create("restRightHipMatrix", "restRightHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restRightHipMatrix);

    attr_restRightKneeMatrix = mAttr.create("restRightKneeMatrix", "restRightKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restRightKneeMatrix);

    attr_restRightHeelMatrix = mAttr.create("restRightHeelMatrix", "restRightHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_restRightHeelMatrix);

    attr_skirtType = eAttr.create("skirtType", "skirtType", 1);
    eAttr.addField("Short", 0);
    eAttr.addField("Long", 1);
    eAttr.setKeyable(true);
    addAttribute(attr_skirtType);

    attr_ringScale = nAttr.create("ringScale", "ringScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.5, 1.0, 0.5);
    nAttr.setKeyable(true);
    addAttribute(attr_ringScale);

    attr_leftRingAxis = eAttr.create("leftRingAxis", "leftRingAxis", 0);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_leftRingAxis);

    attr_rightRingAxis = eAttr.create("rightRingAxis", "rightRingAxis", 3);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_rightRingAxis);

    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0.2f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    attr_thighRadiusX = nAttr.create("thighRadiusX", "thighRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_thighRadiusX);

    attr_thighRadiusZ = nAttr.create("thighRadiusZ", "thighRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_thighRadiusZ);

    attr_kneeRadiusX = nAttr.create("kneeRadiusX", "kneeRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_kneeRadiusX);

    attr_kneeRadiusZ = nAttr.create("kneeRadiusZ", "kneeRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_kneeRadiusZ);

    attr_calfRadiusX = nAttr.create("calfRadiusX", "calfRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_calfRadiusX);

    attr_calfRadiusZ = nAttr.create("calfRadiusZ", "calfRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_calfRadiusZ);

    attr_ankleRadiusX = nAttr.create("ankleRadiusX", "ankleRadiusX", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_ankleRadiusX);

    attr_ankleRadiusZ = nAttr.create("ankleRadiusZ", "ankleRadiusZ", MFnNumericData::kDouble, 1.0);
    nAttr.setMin(0.001);
    nAttr.setKeyable(true);
    addAttribute(attr_ankleRadiusZ);

    attr_thighPosition = nAttr.create("thighPosition", "thighPosition", MFnNumericData::kDouble, 0.5);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    addAttribute(attr_thighPosition);

    attr_calfPosition = nAttr.create("calfPosition", "calfPosition", MFnNumericData::kDouble, 0.5);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    addAttribute(attr_calfPosition);

    const MObject affects[] = {attr_thighRadiusX,        attr_thighRadiusZ,
                               attr_kneeRadiusX,         attr_kneeRadiusZ,
                               attr_calfRadiusX,         attr_calfRadiusZ,
                               attr_ankleRadiusX,        attr_ankleRadiusZ,
                               attr_thighPosition,       attr_calfPosition,
                               attr_bellMatrix,          attr_leftHipMatrix,
                               attr_leftKneeMatrix,      attr_leftHeelMatrix,
                               attr_rightHipMatrix,      attr_rightKneeMatrix,
                               attr_rightHeelMatrix,     attr_skirtType,
                               attr_ringScale,           attr_leftRingAxis,
                               attr_rightRingAxis,       attr_falloff,
                               attr_restGeometry,        attr_restBellMatrix,
                               attr_restLeftHipMatrix,   attr_restLeftKneeMatrix,
                               attr_restLeftHeelMatrix,  attr_restRightHipMatrix,
                               attr_restRightKneeMatrix, attr_restRightHeelMatrix,
                               attr_closedU};
    for (const MObject& attr : affects)
    {
        stat = attributeAffects(attr, outputGeom);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    return MS::kSuccess;
}

MStatus SkirtCollideDeformer::deform(MDataBlock& dataBlock, MItGeometry& iter, const MMatrix&, unsigned int multiIndex)
{
    MStatus stat;
    auto pool = workspacePool;
    auto returnWorkspace = [pool](SkirtCollideSurfaceWorkspace* value) {
        releaseWorkspace(pool, std::unique_ptr<SkirtCollideSurfaceWorkspace>(value));
    };
    std::unique_ptr<SkirtCollideSurfaceWorkspace, decltype(returnWorkspace)> workspace(
        acquireWorkspace(pool).release(), returnWorkspace);

    // Rest points are keyed by geometry index: on periodic surfaces the deformer iterator
    // skips the duplicated CVs, so its count differs from the rest data's iteration.
    auto& restPoints = workspace->restPoints;
    auto& restValid = workspace->restValid;
    restPoints.clear();
    restValid.clear();
    MDataHandle restHandle = dataBlock.inputValue(attr_restGeometry, &stat);
    if (stat && !restHandle.data().isNull())
    {
        if (restHandle.type() == MFnData::kNurbsSurface)
        {
            MFnNurbsSurface restSurface(restHandle.asNurbsSurface(), &stat);
            MPointArray cvs;
            if (stat && restSurface.getCVs(cvs, MSpace::kObject))
            {
                restPoints.resize(cvs.length());
                restValid.resize(cvs.length(), 1);
                for (unsigned int i = 0; i < cvs.length(); ++i)
                    restPoints[i] = cvs[i];
            }
        }
        else
        {
            MItGeometry restIter(restHandle, true, &stat);
            if (stat)
            {
                for (; !restIter.isDone(); restIter.next())
                {
                    const int restIndex = restIter.index();
                    if (restIndex < 0)
                        continue;
                    if (static_cast<std::size_t>(restIndex) >= restPoints.size())
                    {
                        restPoints.resize(restIndex + 1);
                        restValid.resize(restIndex + 1, 0);
                    }
                    restPoints[restIndex] = restIter.position();
                    restValid[restIndex] = 1;
                }
            }
        }
    }
    bool restMissing = restPoints.empty();
    if (!restMissing)
    {
        for (; !iter.isDone(); iter.next())
        {
            const int index = iter.index();
            if (index < 0 || static_cast<std::size_t>(index) >= restPoints.size() || !restValid[index])
            {
                restMissing = true;
                break;
            }
        }
        iter.reset();
    }
    if (restMissing)
    {
        if (!restGeometryWarningIssued.exchange(true))
        {
            MGlobal::displayWarning(MFnDependencyNode(thisMObject()).name() +
                                    ": restGeometry is missing or its points do not match the input geometry.");
        }
        return MS::kSuccess;
    }

    // The current bell matrix is pulled for its status only; its value is not used by the
    // computation and therefore is not part of the input snapshot below.
    dataBlock.inputValue(attr_bellMatrix, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short skirtType = dataBlock.inputValue(attr_skirtType, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MVector ringScale = dataBlock.inputValue(attr_ringScale, &stat).asVector();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short leftRingAxis = dataBlock.inputValue(attr_leftRingAxis, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short rightRingAxis = dataBlock.inputValue(attr_rightRingAxis, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double falloff = dataBlock.inputValue(attr_falloff, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double envelopeValue = dataBlock.inputValue(envelope, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    if (envelopeValue == 0.0 || !std::isfinite(envelopeValue) || !std::isfinite(falloff))
        return MS::kSuccess;

    if (!ColliderInput::validSkirtType(skirtType))
    {
        MGlobal::displayError("SkirtCollideDeformer: skirtType must be 0 or 1.");
        return MS::kInvalidParameter;
    }
    if (!ColliderInput::validAxis(leftRingAxis) || !ColliderInput::validAxis(rightRingAxis))
    {
        MGlobal::displayError("SkirtCollideDeformer: leftRingAxis and rightRingAxis must be between 0 and 5.");
        return MS::kInvalidParameter;
    }

    const MObject currentAttributes[2][3] = {{attr_leftHipMatrix, attr_leftKneeMatrix, attr_leftHeelMatrix},
                                             {attr_rightHipMatrix, attr_rightKneeMatrix, attr_rightHeelMatrix}};
    const MObject restAttributes[2][3] = {
        {attr_restLeftHipMatrix, attr_restLeftKneeMatrix, attr_restLeftHeelMatrix},
        {attr_restRightHipMatrix, attr_restRightKneeMatrix, attr_restRightHeelMatrix}};
    const short axes[] = {leftRingAxis, rightRingAxis};
    const double thighPosition = dataBlock.inputValue(attr_thighPosition, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double calfPosition = dataBlock.inputValue(attr_calfPosition, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MObject radiusAttributes[] = {attr_thighRadiusX, attr_thighRadiusZ, attr_kneeRadiusX,  attr_kneeRadiusZ,
                                        attr_calfRadiusX,  attr_calfRadiusZ,  attr_ankleRadiusX, attr_ankleRadiusZ};
    double radii[8];
    for (int i = 0; i < 8; ++i)
    {
        radii[i] = dataBlock.inputValue(radiusAttributes[i], &stat).asDouble();
        CHECK_MSTATUS_AND_RETURN_IT(stat);
    }
    MMatrix current[2][3], rest[2][3];
    for (int side = 0; side < 2; side++)
        for (int joint = 0; joint < 3; joint++)
        {
            current[side][joint] = dataBlock.inputValue(currentAttributes[side][joint], &stat).asMatrix();
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            rest[side][joint] = dataBlock.inputValue(restAttributes[side][joint], &stat).asMatrix();
            CHECK_MSTATUS_AND_RETURN_IT(stat);
        }
    Leg legs[2];
    std::vector<LegSegment> restCylinders;
    const MMatrix restBellMatrix = dataBlock.inputValue(attr_restBellMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MPoint restWaist = taxis(restBellMatrix);
    const MVector restAxis = getAxis(restBellMatrix, 1).normal();
    if (!isFinitePoint(restWaist) || !finiteVector(legVector(restAxis)) || restAxis.length() < 1e-8)
        return MS::kSuccess;

    MArrayDataHandle inputs = dataBlock.inputArrayValue(input, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = inputs.jumpToElement(multiIndex);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    MDataHandle element = inputs.inputValue(&stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MDataHandle geometry = element.child(inputGeom);
    if (geometry.type() == MFnData::kNurbsSurface)
    {
        MFnNurbsSurface surface(geometry.asNurbsSurface(), &stat);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        const bool closedU = dataBlock.inputValue(attr_closedU, &stat).asBool();
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        const auto topology = surfaceTopologyFor(surface, closedU, surfaceCache);
        MPointArray cvs;
        stat = surface.getCVs(cvs, MSpace::kObject);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        if (!topology || cvs.length() != topology->canonical.size())
        {
            warnInvalidBasis(thisMObject());
            return MS::kSuccess;
        }
        const auto& t = *topology;
        const std::size_t count = t.canonical.size();
        auto& members = workspace->members;
        members.resize(count);
        std::fill(members.begin(), members.end(), 0);
        auto& membership = workspace->membership;
        membership.clear();
        membership.reserve(count);
        bool orderedFullMembership = true;
        for (; !iter.isDone(); iter.next())
        {
            const int index = iter.index();
            if (index < 0 || static_cast<std::size_t>(index) >= count)
            {
                warnInvalidBasis(thisMObject());
                return MS::kSuccess;
            }
            orderedFullMembership = orderedFullMembership &&
                                    static_cast<std::size_t>(index) == membership.size();
            membership.push_back(index);
            members[t.canonical[index]] = 1;
        }
        // Surface CV indices include periodic overlap slots; the deformer weight
        // cache can use a compact index domain. Read logical paint indices directly.
        auto& paint = workspace->paint;
        paint.resize(count);
        std::fill(paint.begin(), paint.end(), 1.0);
        auto& rawPaint = workspace->rawPaint;
        rawPaint.clear();
        MArrayDataHandle weightLists = dataBlock.inputArrayValue(weightList, &stat);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        if (weightLists.jumpToElement(multiIndex) == MS::kSuccess)
        {
            MDataHandle weightListHandle = weightLists.inputValue(&stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            MArrayDataHandle surfaceWeights(weightListHandle.child(weights), &stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            const unsigned int size = surfaceWeights.elementCount(&stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            rawPaint.reserve(size);
            for (unsigned int entry = 0; entry < size; ++entry)
            {
                stat = surfaceWeights.jumpToArrayElement(entry);
                CHECK_MSTATUS_AND_RETURN_IT(stat);
                const unsigned int index = surfaceWeights.elementIndex(&stat);
                CHECK_MSTATUS_AND_RETURN_IT(stat);
                const MDataHandle value = surfaceWeights.inputValue(&stat);
                CHECK_MSTATUS_AND_RETURN_IT(stat);
                const double weight = value.asFloat();
                rawPaint.emplace_back(index, weight);
                if (index < count && t.canonical[index] == static_cast<int>(index))
                    paint[index] = weight;
            }
        }
        for (std::size_t i = 0; i < count; ++i)
        {
            if (i >= restPoints.size() || !restValid[i])
            {
                if (!restGeometryWarningIssued.exchange(true))
                {
                    MGlobal::displayWarning(MFnDependencyNode(thisMObject()).name() +
                                            ": restGeometry does not cover the surface sample support.");
                }
                return MS::kSuccess;
            }
        }

        SkirtCollideInputSnapshot snapshot;
        snapshot.values.reserve(4 * cvs.length() + 5 * restPoints.size() + t.ku.size() + t.kv.size() +
                                2 * rawPaint.size() + membership.size() + 240);
        // Lengths delimit variable-size fields in the scalar comparison sequence.
        snapshot.append(cvs.length());
        bool unitWeights = true;
        for (unsigned int i = 0; i < cvs.length(); ++i)
        {
            unitWeights = unitWeights && cvs[i].w == 1.0;
            snapshot.append(cvs[i]);
        }
        snapshot.append(t.nu);
        snapshot.append(t.nv);
        snapshot.append(t.du);
        snapshot.append(t.dv);
        snapshot.append(t.formU);
        snapshot.append(t.formV);
        snapshot.append(t.closedU);
        snapshot.append(t.ku.size());
        for (double knot : t.ku)
            snapshot.append(knot);
        snapshot.append(t.kv.size());
        for (double knot : t.kv)
            snapshot.append(knot);
        snapshot.append(restPoints.size());
        for (std::size_t i = 0; i < restPoints.size(); ++i)
        {
            snapshot.append(restPoints[i]);
            snapshot.append(restValid[i]);
        }
        for (int side = 0; side < 2; ++side)
            for (int joint = 0; joint < 3; ++joint)
            {
                snapshot.append(current[side][joint]);
                snapshot.append(rest[side][joint]);
            }
        snapshot.append(restBellMatrix);
        for (double radius : radii)
            snapshot.append(radius);
        snapshot.append(thighPosition);
        snapshot.append(calfPosition);
        snapshot.append(ringScale.x);
        snapshot.append(ringScale.y);
        snapshot.append(ringScale.z);
        snapshot.append(leftRingAxis);
        snapshot.append(rightRingAxis);
        snapshot.append(skirtType);
        snapshot.append(falloff);
        snapshot.append(envelopeValue);
        snapshot.append(rawPaint.size());
        for (const auto& entry : rawPaint)
        {
            snapshot.append(entry.first);
            snapshot.append(entry.second);
        }
        snapshot.append(membership.size());
        for (int index : membership)
            snapshot.append(index);
        const bool iteratorOutput = t.formU != MFnNurbsSurface::kPeriodic &&
                                    t.formV != MFnNurbsSurface::kPeriodic &&
                                    unitWeights && orderedFullMembership && membership.size() == count;

        std::shared_ptr<const SkirtCollideEvaluation> cached;
        if (snapshot.finite)
        {
            const std::lock_guard<std::mutex> lock(evaluationMutex);
            const auto found = evaluations.find(multiIndex);
            if (found != evaluations.end())
                cached = found->second;
        }
        const bool hit = cached && snapshot.matches(cached->input);
        if (!hit)
        {
            auto& points = workspace->points;
            points.resize(count);
            std::fill(points.begin(), points.end(), SurfacePoint{});
            const LegVec3 restWaistPosition = legPoint(restWaist);
            const LegVec3 restAxisDirection = legVector(restAxis);
            for (std::size_t i = 0; i < count; ++i)
            {
                auto& p = points[i];
                // getCVs returns Cartesian xyz and a separate rational weight in w.
                // The solver uses xyz; the basis reads the stored w and output restores it.
                p.point = legPoint(cvs[static_cast<unsigned int>(i)]);
                p.weight = cvs[static_cast<unsigned int>(i)].w;
                p.restPoint = legPoint(restPoints[i]);
                // Maya omits periodic duplicates from the deformer iterator. Their
                // paint and membership both follow the canonical CV.
                const double q = paint[t.canonical[i]];
                p.q = members[t.canonical[i]] && std::isfinite(q) ? clampValue(q, 0.0, 1.0) : 0.0;
                p.finite = finiteVector(p.point) && finiteVector(p.restPoint);
                const LegVec3 offset = p.restPoint - restWaistPosition;
                const LegVec3 radial = offset - restAxisDirection * (offset * restAxisDirection);
                if (!p.finite || radial.length() < 1e-8)
                    p.q = 0.0;
            }
            buildLegs(current, rest, axes, ringScale, radii, thighPosition, calfPosition, skirtType == 1, legs,
                      restCylinders);
            std::vector<double> restKey;
            restKey.reserve(2 * 3 * 16 + 2 + 3 + 8 + 4);
            for (int side = 0; side < 2; ++side)
                for (int joint = 0; joint < 3; ++joint)
                    for (unsigned int row = 0; row < 4; ++row)
                        for (unsigned int column = 0; column < 4; ++column)
                            restKey.push_back(rest[side][joint][row][column]);
            restKey.push_back(leftRingAxis);
            restKey.push_back(rightRingAxis);
            restKey.push_back(ringScale.x);
            restKey.push_back(ringScale.y);
            restKey.push_back(ringScale.z);
            for (double radius : radii)
                restKey.push_back(radius);
            restKey.push_back(thighPosition);
            restKey.push_back(calfPosition);
            restKey.push_back(skirtType);
            restKey.push_back(falloff);
            std::shared_ptr<const SkirtCollideSurfacePreparation> preparation;
            if (snapshot.finite)
                preparation = surfaceCache.findPreparation(topology, [&](const SkirtCollideSurfacePreparation& candidate) {
                    return candidate.matches(topology, cvs, restPoints, restKey);
                });
            std::shared_ptr<SkirtCollideSurfacePreparation> built;
            if (!preparation && snapshot.finite)
                built = std::make_shared<SkirtCollideSurfacePreparation>();
            auto& groups = workspace->groups;
            if (!coupleSurface(t, points, legs, restCylinders, restWaist, restAxis, falloff, groups,
                               preparation.get(), built.get(), *workspace))
            {
                warnInvalidBasis(thisMObject());
                return MS::kSuccess;
            }
            if (built)
            {
                built->topology = topology;
                built->weights.reserve(cvs.length());
                for (unsigned int i = 0; i < cvs.length(); ++i)
                    built->weights.push_back(cvs[i].w);
                built->restPoints = restPoints;
                built->restKey = std::move(restKey);
                surfaceCache.publishPreparation(topology, built, preparationPayloadBytes(*built),
                    [&](const SkirtCollideSurfacePreparation& candidate) {
                        return candidate.matches(topology, cvs, restPoints, built->restKey);
                    });
            }
            for (std::size_t i = 0; i < count; ++i)
            {
                const auto& group = groups[t.groupOf[i]];
                if (group.paint == 0.0)
                    continue;
                const LegVec3 resultPosition = points[i].point + envelopeValue * group.value;
                MPoint result(resultPosition.x, resultPosition.y, resultPosition.z,
                              cvs[static_cast<unsigned int>(i)].w);
                if (isFinitePoint(result))
                    cvs[static_cast<unsigned int>(i)] = result;
            }
        }
        const MPointArray& outputCVs = hit ? cached->output : cvs;
        if (iteratorOutput)
            stat = iter.setAllPositions(outputCVs, MSpace::kObject);
        else
        {
            // The iterator omits periodic duplicates and cannot preserve rational w.
            MArrayDataHandle outputs = dataBlock.outputArrayValue(outputGeom, &stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            stat = outputs.jumpToElement(multiIndex);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            MDataHandle outputGeometry = outputs.outputValue(&stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            MFnNurbsSurface outputSurface(outputGeometry.asNurbsSurface(), &stat);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            stat = outputSurface.setCVs(outputCVs, MSpace::kObject);
        }
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        if (!hit && snapshot.finite)
        {
            auto result = std::make_shared<SkirtCollideEvaluation>();
            result->input = std::move(snapshot);
            // The surface is already written; a failed copy only means this result is not reusable.
            if (result->output.copy(cvs) == MS::kSuccess)
            {
                std::shared_ptr<const SkirtCollideEvaluation> completed = std::move(result);
                const std::lock_guard<std::mutex> lock(evaluationMutex);
                evaluations[multiIndex].swap(completed);
            }
        }
        return MS::kSuccess;
    }
    buildLegs(current, rest, axes, ringScale, radii, thighPosition, calfPosition, skirtType == 1, legs, restCylinders);
    const LegContactGeometry contactGeometry[2] = {legPrepareContactGeometry(legs[0], falloff),
                                                    legPrepareContactGeometry(legs[1], falloff)};
    std::map<int, PointState> points;
    LegScratch scratch;
    for (; !iter.isDone(); iter.next())
    {
        PointState p;
        p.index = iter.index();
        p.point = iter.position();
        p.restPoint = restPoints[p.index];
        const double paint = weightValue(dataBlock, multiIndex, p.index);
        p.q = std::isfinite(paint) ? clampValue(paint, 0.0, 1.0) : 0.0;
        p.finite = isFinitePoint(p.point) && isFinitePoint(p.restPoint);
        if (!p.finite)
            p.q = 0.0;
        if (p.finite)
        {
            const MVector offset = p.restPoint - restWaist;
            MVector waistRadial = offset - restAxis * (offset * restAxis);
            if (waistRadial.length() < 1e-8)
                p.q = 0.0;
            if (p.q > 0.0)
            {
                waistRadial.normalize();
                for (int side = 0; side < 2; ++side)
                {
                    LegContact constraint;
                    if (legConstraint(legs[side], legPoint(p.point), legPoint(p.restPoint), legVector(waistRadial),
                                      falloff, constraint, scratch, contactGeometry[side]))
                        p.constraints.push_back(constraint);
                }
            }
        }
        points.emplace(p.index, p);
    }
    stat = buildAdjacency(geometry, thisMObject(), points);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    couple(points, restCylinders, falloff);
    iter.reset();
    for (; !iter.isDone(); iter.next())
    {
        const PointState& p = points.at(iter.index());
        if (p.q == 0.0)
            continue;
        MPoint result = p.point + envelopeValue * p.u;
        if (isFinitePoint(result))
        {
            stat = iter.setPosition(result);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
        }
    }
    return MS::kSuccess;
}
