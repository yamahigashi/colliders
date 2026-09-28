#pragma once

#include <maya/MMatrix.h>
#include <maya/MObject.h>
#include <maya/MPxDeformerNode.h>
#include <maya/MStatus.h>
#include <maya/MTypeId.h>

#include <atomic>
#include <map>
#include <memory>
#include <mutex>

#include "skirtSurfaceCache.h"

struct SkirtCollideSurfaceTopology;
struct SkirtCollideEvaluation;
struct SkirtCollideSurfacePreparation;
struct SkirtCollideSurfaceWorkspace;
struct SkirtCollideSurfaceWorkspacePool;

class SkirtCollideDeformer : public MPxDeformerNode
{
public:
    static MTypeId typeId;

    static MObject attr_bellMatrix;
    static MObject attr_leftHipMatrix;
    static MObject attr_leftKneeMatrix;
    static MObject attr_leftHeelMatrix;
    static MObject attr_rightHipMatrix;
    static MObject attr_rightKneeMatrix;
    static MObject attr_rightHeelMatrix;

    static MObject attr_restGeometry;
    static MObject attr_restBellMatrix;
    static MObject attr_restLeftHipMatrix;
    static MObject attr_restLeftKneeMatrix;
    static MObject attr_restLeftHeelMatrix;
    static MObject attr_restRightHipMatrix;
    static MObject attr_restRightKneeMatrix;
    static MObject attr_restRightHeelMatrix;

    static MObject attr_skirtType;
    static MObject attr_ringScale;
    static MObject attr_thighRadiusX;
    static MObject attr_thighRadiusZ;
    static MObject attr_kneeRadiusX;
    static MObject attr_kneeRadiusZ;
    static MObject attr_calfRadiusX;
    static MObject attr_calfRadiusZ;
    static MObject attr_ankleRadiusX;
    static MObject attr_ankleRadiusZ;
    static MObject attr_thighPosition;
    static MObject attr_calfPosition;

    static MObject attr_leftRingAxis;
    static MObject attr_rightRingAxis;
    static MObject attr_falloff;

    SkirtCollideDeformer();
    virtual ~SkirtCollideDeformer() override;

    static void* creator() { return new SkirtCollideDeformer(); }
    static MStatus initialize();

    virtual MStatus deform(MDataBlock& dataBlock, MItGeometry& iter, const MMatrix&, unsigned int multiIndex) override;

private:
    std::atomic<bool> restGeometryWarningIssued;
    SkirtSurfaceCache<SkirtCollideSurfaceTopology, SkirtCollideSurfacePreparation> surfaceCache;
    std::mutex evaluationMutex;
    std::map<unsigned int, std::shared_ptr<const SkirtCollideEvaluation>> evaluations;
    std::shared_ptr<SkirtCollideSurfaceWorkspacePool> workspacePool;

    static std::unique_ptr<SkirtCollideSurfaceWorkspace>
    acquireWorkspace(const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>& pool);
    static void releaseWorkspace(const std::shared_ptr<SkirtCollideSurfaceWorkspacePool>& pool,
                                 std::unique_ptr<SkirtCollideSurfaceWorkspace> workspace) noexcept;
};
