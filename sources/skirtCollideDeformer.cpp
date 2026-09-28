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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
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
        double u, v;
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
    std::size_t termCapacity = 0;
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


struct SurfaceContact
{
    MVector normal = MVector(0.0, 0.0, 0.0);
    double distance = 0.0, weight = 0.0, cache = 0.0;
};

struct SurfaceSample
{
    MPoint point = MPoint(0.0, 0.0, 0.0), rest = MPoint(0.0, 0.0, 0.0);
    double denominator = 0.0;
    SurfaceContact contacts[2];
};

struct SurfaceGroup
{
    double paint = 1.0, diagonal = 0.0;
    MVector value = MVector(0.0, 0.0, 0.0);
};

struct SurfaceTerm
{
    MVector normal;
    double target, weight, basis;
    bool active;
};

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
    const auto vs = spanSites(t.kv, dv, nv, 4);
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
        for (double v : vs)
        {
            const int firstV = splineBasis(t.kv, dv, nv, v, bv);
            SkirtCollideSurfaceTopology::Sample sample;
            sample.u = u;
            sample.v = v;
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
        t.termCapacity = (std::max)(t.termCapacity, static_cast<std::size_t>(group.attachmentEnd) * 2);
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

std::shared_ptr<const SkirtCollideSurfaceTopology> surfaceTopologyFor(
    MFnNurbsSurface& surface, bool closedU, std::mutex& mutex,
    std::shared_ptr<const SkirtCollideSurfaceTopology>& cache)
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
    const std::lock_guard<std::mutex> lock(mutex);
    if (cache && cache->nu == key.nu && cache->nv == key.nv && cache->du == key.du && cache->dv == key.dv &&
        cache->formU == key.formU && cache->formV == key.formV && cache->closedU == key.closedU &&
        cache->ku == key.ku && cache->kv == key.kv)
        return cache;
    auto built = std::make_shared<SkirtCollideSurfaceTopology>(std::move(key));
    if (!buildSurfaceTopology(*built))
        return {};
    cache = built;
    return cache;
}

MVector solveSurfaceGroup(const SurfaceGroup& group, const SkirtCollideSurfaceTopology::Group& layout,
                          const SkirtCollideSurfaceTopology& topology, const std::vector<double>& edges,
                          const std::vector<SurfaceGroup>& groups, SurfaceTerm* terms, std::size_t count)
{
    MVector rhs0(0.0, 0.0, 0.0);
    for (int e = layout.edgeBegin; e < layout.edgeEnd; ++e)
        rhs0 += edges[e] * groups[topology.edgeGroups[e]].value;
    const auto objective = [&](const MVector& z) {
        double value = layout.mass * (z * z);
        for (int e = layout.edgeBegin; e < layout.edgeEnd; ++e)
        {
            const MVector offset = z - groups[topology.edgeGroups[e]].value;
            value += edges[e] * (offset * offset);
        }
        for (std::size_t i = 0; i < count; ++i)
        {
            const SurfaceTerm& term = terms[i];
            const double hinge = (std::max)(0.0, term.target - term.basis * (term.normal * z));
            value += term.weight * hinge * hinge;
        }
        return value;
    };
    const auto consistent = [&](const MVector& z) {
        for (std::size_t i = 0; i < count; ++i)
            if (terms[i].active != (terms[i].target - terms[i].basis * (terms[i].normal * z) > 0.0))
                return false;
        return true;
    };
    MVector z = group.value;
    for (int iteration = 0; iteration < 50; ++iteration)
    {
        const double diagonal = group.diagonal;
        double matrix[3][3] = {{diagonal, 0.0, 0.0}, {0.0, diagonal, 0.0}, {0.0, 0.0, diagonal}};
        MVector rhs = rhs0;
        for (std::size_t i = 0; i < count; ++i)
        {
            SurfaceTerm& term = terms[i];
            term.active = term.target - term.basis * (term.normal * z) > 0.0;
            if (!term.active)
                continue;
            for (int r = 0; r < 3; ++r)
            {
                rhs[r] += term.weight * term.basis * term.target * term.normal[r];
                for (int c = 0; c < 3; ++c)
                    matrix[r][c] += term.weight * term.basis * term.basis * term.normal[r] * term.normal[c];
            }
        }
        const double det = determinant(matrix);
        MVector candidate;
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
            const MVector step = candidate - z;
            candidate = z;
            double scale = 0.5;
            for (int k = 1; k <= 30; ++k, scale *= 0.5)
            {
                const MVector trial = z + scale * step;
                if (objective(trial) < oldValue)
                {
                    candidate = trial;
                    break;
                }
            }
        }
        const bool done = consistent(candidate) || (candidate - z).length() <= 1e-12 * (1.0 + z.length());
        z = candidate;
        if (done)
            break;
    }
    return z;
}

struct SurfacePoint
{
    MPoint point, restPoint;
    double q = 0.0, radius = 0.0;
    bool finite = false;
};

bool coupleSurface(const SkirtCollideSurfaceTopology &topology, const MPointArray &cvs,
                   std::vector<SurfacePoint> &points, const Leg legs[2], const std::vector<LegSegment> &restCylinders,
                   const MPoint &waist, const MVector &axis, double kappa, std::vector<SurfaceGroup> &groups,
                   int sweeps, int regenerations)
{
    const auto& t = topology;
    for (unsigned int i = 0; i < cvs.length(); ++i)
        if (!std::isfinite(cvs[i].w) || cvs[i].w <= 0.0)
            return false;
    groups.resize(t.groups.size());
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        if (t.canonical[i] != static_cast<int>(i))
            continue;
        auto& point = points[i];
        auto& group = groups[t.groupOf[i]];
        group.paint = (std::min)(group.paint, point.q);
        if (point.finite)
            point.radius = legLocalRadius(restCylinders, legPoint(point.restPoint));
    }
    std::vector<double> edges(t.edgeGroups.size(), 0.0);
    for (const auto& source : t.edgeSources)
    {
        const auto& point = points[source.a];
        const auto& other = points[source.b];
        if (!point.finite || !other.finite)
            continue;
        const double scale = kappa * (point.radius + other.radius) / 2.0;
        const MVector offset = point.restPoint - other.restPoint;
        edges[source.edge] +=
            kappa > 0.0 && !restCylinders.empty() ? scale * scale / (std::max)(offset * offset, 1e-12) : 0.0;
    }
    for (std::size_t g = 0; g < groups.size(); ++g)
    {
        groups[g].diagonal = t.groups[g].mass;
        for (int e = t.groups[g].edgeBegin; e < t.groups[g].edgeEnd; ++e)
            groups[g].diagonal += edges[e];
    }
    std::vector<SurfaceSample> samples(t.samples.size());
    std::vector<double> rows(t.rows.size(), 0.0), beta(t.betaGroups.size(), 0.0);
    for (std::size_t s = 0; s < samples.size(); ++s)
    {
        auto& sample = samples[s];
        const auto& layout = t.samples[s];
        double total = 0.0;
        for (int k = layout.tensorBegin; k < layout.tensorEnd; ++k)
        {
            const auto& coefficient = t.tensors[k];
            const double w = coefficient.basis * cvs[coefficient.cv].w;
            rows[coefficient.row] += w;
            total += w;
        }
        if (!std::isfinite(total) || total <= 0.0)
            return false;
        double sum = 0.0;
        MVector p(0.0, 0.0, 0.0), r(0.0, 0.0, 0.0);
        for (int k = layout.rowBegin; k < layout.rowEnd; ++k)
        {
            const double b = rows[k] / total;
            if (!std::isfinite(b) || b < 0.0)
                return false;
            sum += b;
            sample.denominator += b * b;
            const auto& point = points[t.rows[k].cv];
            p += b * MVector(point.point.x, point.point.y, point.point.z);
            r += b * MVector(point.restPoint.x, point.restPoint.y, point.restPoint.z);
            beta[t.rows[k].beta] += b;
        }
        if (std::abs(sum - 1.0) > 1e-9)
            return false;
        sample.point = MPoint(p);
        sample.rest = MPoint(r);
    }
    std::vector<SurfaceTerm> terms(t.termCapacity);
    std::vector<std::array<LegSampleTransport, 2>> transports(samples.size());
    LegScratch scratch;
    bool mixedPaint = false;
    for (const auto& group : groups)
        if (group.paint > 0.0 && group.paint < 1.0)
            mixedPaint = true;
    std::vector<MVector> displacements(samples.size());
    std::vector<double> startNorms(samples.size(), 0.0);
    std::vector<char> sampleValid(samples.size(), 0);
    std::vector<std::array<double, 2>> capsuleGaps(samples.size());
    std::vector<std::array<char, 2>> deferred(samples.size());
    for (int block = 0; block <= regenerations; ++block)
    {
        for (auto& gaps : capsuleGaps)
        {
            gaps[0] = gaps[1] = 0.0;
        }
        for (auto& flags : deferred)
        {
            flags[0] = flags[1] = 0;
        }
        for (std::size_t s = 0; s < samples.size(); ++s)
        {
            auto& sample = samples[s];
            MVector displacement(0.0, 0.0, 0.0);
            for (int b = t.samples[s].betaBegin; b < t.samples[s].betaEnd; ++b)
                displacement += beta[b] * groups[t.betaGroups[b]].value;
            displacements[s] = displacement;
            const MPoint position = sample.point + displacement;
            const MVector offset = sample.rest - waist;
            MVector radial = offset - axis * (offset * axis);
            const bool valid = isFinitePoint(position) && isFinitePoint(sample.rest) && radial.length() >= 1e-8;
            sampleValid[s] = valid ? 1 : 0;
            startNorms[s] = valid ? displacement.length() : 0.0;
            if (valid)
                radial.normalize();
            if (block == 0 && valid)
                for (int side = 0; side < 2; ++side)
                    legSampleTransport(legs[side], legPoint(sample.rest), legVector(radial), transports[s][side]);
            for (int side = 0; side < 2; ++side)
            {
                SurfaceContact& contact = sample.contacts[side];
                contact = SurfaceContact();
                if (!valid)
                    continue;
                bool inside = false;
                for (const auto& segment : legs[side].segments)
                    if (legLocal(segment.current, legPoint(position)).phi <= 0.0)
                    {
                        inside = true;
                        break;
                    }
                LegContact constraint;
                if (inside)
                {
                    if (legConstraint(legs[side], legPoint(position), transports[s][side], kappa, constraint,
                                      scratch) &&
                        constraint.weight != 0.0)
                    {
                        contact.normal = mayaLegVector(constraint.normal);
                        contact.cache = constraint.normal * legVector(displacement);
                        contact.distance = constraint.distance + contact.cache;
                        contact.weight = 4.0 * constraint.weight / sample.denominator;
                    }
                    continue;
                }
                const double gap = mixedPaint ? 0.0 : legCapsuleDistance(legs[side], legPoint(position));
                capsuleGaps[s][side] = gap;
                const double scaled = gap * (1.0 - 1e-9);
                const double start = startNorms[s];
                if (mixedPaint || !std::isfinite(gap) || !std::isfinite(start) || !(1.25 * scaled > 1.25 * start))
                {
                    if (legConstraint(legs[side], legPoint(position), transports[s][side], kappa, constraint,
                                      scratch) &&
                        constraint.weight != 0.0)
                    {
                        contact.normal = mayaLegVector(constraint.normal);
                        contact.cache = constraint.normal * legVector(displacement);
                        contact.distance = constraint.distance + contact.cache;
                        contact.weight = 4.0 * constraint.weight / sample.denominator;
                    }
                    continue;
                }
                deferred[s][side] = 1;
            }
        }
        double bound = 0.0;
        for (std::size_t g = 0; g < groups.size(); ++g)
            bound += static_cast<double>(t.groups[g].mass) * (groups[g].value * groups[g].value);
        for (std::size_t g = 0; g < groups.size(); ++g)
            for (int e = t.groups[g].edgeBegin; e < t.groups[g].edgeEnd; ++e)
            {
                const MVector offset = groups[g].value - groups[topology.edgeGroups[e]].value;
                bound += 0.5 * edges[e] * (offset * offset);
            }
        for (std::size_t s = 0; s < samples.size(); ++s)
            for (int side = 0; side < 2; ++side)
            {
                const SurfaceContact& contact = samples[s].contacts[side];
                if (contact.weight == 0.0)
                    continue;
                const double hinge = (std::max)(0.0, 1.25 * contact.distance - contact.cache);
                bound += contact.weight * hinge * hinge;
            }
        for (std::size_t s = 0; s < samples.size(); ++s)
        {
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
                auto& sample = samples[s];
                const MPoint position = sample.point + displacements[s];
                const double scaled = capsuleGaps[s][side] * (1.0 - 1e-9);
                const double start = startNorms[s];
                if (mixedPaint || !std::isfinite(scaled) || !std::isfinite(start) || !std::isfinite(radius) ||
                    !(1.25 * scaled > 1.25 * start + radius))
                {
                    LegContact constraint;
                    if (sampleValid[s] &&
                        legConstraint(legs[side], legPoint(position), transports[s][side], kappa, constraint,
                                      scratch) &&
                        constraint.weight != 0.0)
                    {
                        SurfaceContact& contact = sample.contacts[side];
                        contact.normal = mayaLegVector(constraint.normal);
                        contact.cache = constraint.normal * legVector(displacements[s]);
                        contact.distance = constraint.distance + contact.cache;
                        contact.weight = 4.0 * constraint.weight / sample.denominator;
                    }
                }
            }
        }
        for (int sweep = block * sweeps / (regenerations + 1); sweep < (block + 1) * sweeps / (regenerations + 1); ++sweep)
            for (std::size_t g = 0; g < groups.size(); ++g)
            {
                auto& group = groups[g];
                const auto& layout = t.groups[g];
                if (group.paint == 0.0)
                    continue;
                std::size_t count = 0;
                for (int a = layout.attachmentBegin; a < layout.attachmentEnd; ++a)
                {
                    const auto& attachment = t.attachments[a];
                    for (const auto& contact : samples[attachment.sample].contacts)
                        if (contact.weight != 0.0)
                        {
                            SurfaceTerm& term = terms[count++];
                            term.normal = contact.normal;
                            term.basis = beta[attachment.beta];
                            term.weight = contact.weight;
                            term.target =
                                1.25 * contact.distance - (contact.cache - term.basis * (contact.normal * group.value));
                        }
                }
                const MVector value =
                    group.paint * solveSurfaceGroup(group, layout, t, edges, groups, terms.data(), count);
                const MVector change = value - group.value;
                group.value = value;
                for (int a = layout.attachmentBegin; a < layout.attachmentEnd; ++a)
                {
                    const auto& attachment = t.attachments[a];
                    for (auto& contact : samples[attachment.sample].contacts)
                        contact.cache += beta[attachment.beta] * (contact.normal * change);
                }
            }
    }
    return true;
}

} // namespace

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

    // Rest points are keyed by geometry index: on periodic surfaces the deformer iterator
    // skips the duplicated CVs, so its count differs from the rest data's iteration.
    std::vector<MPoint> restPoints;
    std::vector<char> restValid;
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
        const auto topology = surfaceTopologyFor(surface, closedU, surfaceTopologyMutex, surfaceTopology);
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
        std::vector<char> members(count, 0);
        std::vector<int> membership;
        membership.reserve(count);
        for (; !iter.isDone(); iter.next())
        {
            const int index = iter.index();
            if (index < 0 || static_cast<std::size_t>(index) >= count)
            {
                warnInvalidBasis(thisMObject());
                return MS::kSuccess;
            }
            membership.push_back(index);
            members[t.canonical[index]] = 1;
        }
        // Surface CV indices include periodic overlap slots; the deformer weight
        // cache can use a compact index domain. Read logical paint indices directly.
        std::vector<double> paint(count, 1.0);
        std::vector<std::pair<unsigned int, double>> rawPaint;
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
        for (unsigned int i = 0; i < cvs.length(); ++i)
            snapshot.append(cvs[i]);
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
            std::vector<SurfacePoint> points(count);
            for (std::size_t i = 0; i < count; ++i)
            {
                auto& p = points[i];
                p.point = cvs[static_cast<unsigned int>(i)];
                p.restPoint = restPoints[i];
                // getCVs returns Cartesian xyz and a separate rational weight in w.
                p.point.w = p.restPoint.w = 1.0;
                // Maya omits periodic duplicates from the deformer iterator. Their
                // paint and membership both follow the canonical CV.
                const double q = paint[t.canonical[i]];
                p.q = members[t.canonical[i]] && std::isfinite(q) ? clampValue(q, 0.0, 1.0) : 0.0;
                p.finite = isFinitePoint(p.point) && isFinitePoint(p.restPoint);
                const MVector offset = p.restPoint - restWaist;
                const MVector radial = offset - restAxis * (offset * restAxis);
                if (!p.finite || radial.length() < 1e-8)
                    p.q = 0.0;
            }
            buildLegs(current, rest, axes, ringScale, radii, thighPosition, calfPosition, skirtType == 1, legs,
                      restCylinders);
            std::vector<SurfaceGroup> groups;
            if (!coupleSurface(t, cvs, points, legs, restCylinders, restWaist, restAxis, falloff, groups, 10, 1))
            {
                warnInvalidBasis(thisMObject());
                return MS::kSuccess;
            }
            for (std::size_t i = 0; i < count; ++i)
            {
                const auto& group = groups[t.groupOf[i]];
                if (group.paint == 0.0)
                    continue;
                MPoint result = points[i].point + envelopeValue * group.value;
                result.w = cvs[static_cast<unsigned int>(i)].w;
                if (isFinitePoint(result))
                    cvs[static_cast<unsigned int>(i)] = result;
            }
        }
        // Write the output surface directly: iterator writes discard rational w
        // and omit periodic duplicates. Input geometry remains read-only.
        MArrayDataHandle outputs = dataBlock.outputArrayValue(outputGeom, &stat);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        stat = outputs.jumpToElement(multiIndex);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        MDataHandle outputGeometry = outputs.outputValue(&stat);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        MFnNurbsSurface outputSurface(outputGeometry.asNurbsSurface(), &stat);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
        const MPointArray& outputCVs = hit ? cached->output : cvs;
        stat = outputSurface.setCVs(outputCVs, MSpace::kObject);
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
                for (const auto& leg : legs)
                {
                    LegContact constraint;
                    if (legConstraint(leg, legPoint(p.point), legPoint(p.restPoint), legVector(waistRadial), falloff,
                                      constraint, scratch))
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
