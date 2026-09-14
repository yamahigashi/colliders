#include <set>
#include <map>

#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MIntArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnDependencyNode.h>

#include <maya/MTransformationMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MEulerRotation.h>

#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MMeshIntersector.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MFnMatrixData.h>

#include <maya/MArrayDataBuilder.h>
#include <maya/MArrayDataHandle.h>

#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>

#include <tbb/parallel_for.h>

#include "bellCollider.h"
#include "colliderInputValidation.h"
#include "pluginIdentity.h"
#include "utils.hpp"

using namespace std;

MTypeId BellCollider::typeId(PluginIdentity::kBellTypeId);

MObject BellCollider::attr_bellMatrix;
MObject BellCollider::attr_ringMatrix;
MObject BellCollider::attr_bellSubdivision;
MObject BellCollider::attr_ringSubdivision;
MObject BellCollider::attr_bellBottomRadius;
MObject BellCollider::attr_falloff;
MObject BellCollider::attr_collision;
MObject BellCollider::attr_drawColor;
MObject BellCollider::attr_drawOpacity;
MObject BellCollider::attr_outputCurve;
MObject BellCollider::attr_outputBellMesh;

MString BellCollider::drawDbClassification = MString("drawdb/geometry/") + PluginIdentity::kBellNodeName;
MString BellCollider::drawRegistrantId = PluginIdentity::kDrawRegistrant;

MStatus BellCollider::compute(const MPlug &plug, MDataBlock &dataBlock)
{
    if (plug != attr_outputCurve && plug != attr_outputBellMesh)
        return MS::kUnknownParameter;

    // Extract inputs
    BellColliderInputs inputs;
    inputs.bellSubdivision = dataBlock.inputValue(attr_bellSubdivision).asInt();
    if (!ColliderInput::validSubdivision(inputs.bellSubdivision))
    {
        MGlobal::displayError("BellCollider: bellSubdivision must be between 3 and 4096.");
        return MS::kInvalidParameter;
    }
    inputs.bellMatrix = dataBlock.inputValue(attr_bellMatrix).asMatrix();
    
    MStatus stat;
    auto ringMatrixHandle = dataBlock.inputArrayValue(attr_ringMatrix, &stat);
    if (!stat) return stat;
    const unsigned int ringCount = ringMatrixHandle.elementCount(&stat);
    if (!stat) return stat;
    inputs.rings.reserve(ringCount);
    for (unsigned int i = 0; i < ringCount; ++i)
    {
        stat = ringMatrixHandle.jumpToArrayElement(i);
        if (!stat) return stat;
        const MDataHandle ringValue = ringMatrixHandle.inputValue(&stat);
        if (!stat) return stat;
        inputs.rings.emplace_back(ringValue.asMatrix());
    }

    const float bottomRadius = dataBlock.inputValue(attr_bellBottomRadius).asFloat();
    inputs.falloff = dataBlock.inputValue(attr_falloff).asFloat();
    inputs.collision = dataBlock.inputValue(attr_collision).asFloat();

    MPointArray baseBellPoints = BellColliderSolver::makeBellPoints(inputs.bellMatrix, 1, inputs.bellSubdivision, 1, bottomRadius, 1);
    BellColliderSolver::roundMeshPoints(baseBellPoints);

    BellColliderOutputs outputs;
    stat = BellColliderSolver::solve(inputs, baseBellPoints, outputs);
    if (stat != MS::kSuccess)
        return stat;

    dataBlock.outputValue(attr_outputCurve).setMObject(BellColliderSolver::makeBellCurve(outputs.points, inputs.bellSubdivision));
    dataBlock.outputValue(attr_outputBellMesh).setMObject(BellColliderSolver::makeBellMesh(outputs.points, inputs.bellSubdivision));

    dataBlock.setClean(attr_outputCurve);
    dataBlock.setClean(attr_outputBellMesh);

    return MS::kSuccess;
}

MStatus BellCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;
    
    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_bellMatrix);

    attr_ringMatrix = mAttr.create("ringMatrix", "ringMatrix");
    mAttr.setArray(true);
    mAttr.setHidden(true);
    addAttribute(attr_ringMatrix);

    attr_bellSubdivision = nAttr.create("bellSubdivision", "bellSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_bellSubdivision);

    attr_ringSubdivision = nAttr.create("ringSubdivision", "ringSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_ringSubdivision);

    attr_bellBottomRadius = nAttr.create("bellBottomRadius", "bellBottomRadius", MFnNumericData::kFloat, 0.8);
    nAttr.setMin(0);
    nAttr.setKeyable(true);
    addAttribute(attr_bellBottomRadius);

    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0);
    nAttr.setMin(-1);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    attr_collision = nAttr.create("collision", "collision", MFnNumericData::kFloat, 0);
    nAttr.setMin(0);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_collision);

    attr_drawColor = nAttr.create("drawColor", "drawColor", MFnNumericData::k3Double);
    nAttr.setDefault(0.0, 0.01, 0.11);
    nAttr.setMin(0, 0, 0);
    nAttr.setMax(1, 1, 1);
    nAttr.setKeyable(true);
    addAttribute(attr_drawColor);

    attr_drawOpacity = nAttr.create("drawOpacity", "drawOpacity", MFnNumericData::kFloat, 0.3);
    nAttr.setMin(0);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_drawOpacity);

    attr_outputCurve = tAttr.create("outputCurve", "outputCurve", MFnData::kNurbsCurve);
    tAttr.setHidden(true);
    addAttribute(attr_outputCurve);

    attr_outputBellMesh = tAttr.create("outputBellMesh", "outputBellMesh", MFnData::kMesh);
    tAttr.setHidden(true);
    addAttribute(attr_outputBellMesh);

    attributeAffects(attr_bellMatrix, attr_outputCurve);
    attributeAffects(attr_ringMatrix, attr_outputCurve);
    attributeAffects(attr_bellSubdivision, attr_outputCurve);
    attributeAffects(attr_bellBottomRadius, attr_outputCurve);
    attributeAffects(attr_falloff, attr_outputCurve);
    attributeAffects(attr_collision, attr_outputCurve);

    attributeAffects(attr_bellMatrix, attr_outputBellMesh);
    attributeAffects(attr_ringMatrix, attr_outputBellMesh);
    attributeAffects(attr_bellSubdivision, attr_outputBellMesh);
    attributeAffects(attr_bellBottomRadius, attr_outputBellMesh);
    attributeAffects(attr_falloff, attr_outputBellMesh);
    attributeAffects(attr_collision, attr_outputBellMesh);

    return MS::kSuccess;
}

MUserData* BellColliderDrawOverride::prepareForDraw(
    const MDagPath& objPath, 
    const MDagPath& cameraPath, 
    const MHWRender::MFrameContext& frameContext, 
    MUserData* oldData)
{
    MStatus stat;
    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto* data = dynamic_cast<BellColliderDrawData*>(oldData);
    if (!data)
        data = new BellColliderDrawData();

    MPlug bellMatrixPlug(obj, BellCollider::attr_bellMatrix);
    MPlug ringMatrixPlug(obj, BellCollider::attr_ringMatrix);
    MPlug ringSubdivisionPlug(obj, BellCollider::attr_ringSubdivision);
    MPlug drawColorPlug(obj, BellCollider::attr_drawColor);
    MPlug drawOpacityPlug(obj, BellCollider::attr_drawOpacity);
    MPlug outputBellMeshPlug(obj, BellCollider::attr_outputBellMesh);

    MMatrix bellMatrix;
    MObject matrixObj;
    if (bellMatrixPlug.getValue(matrixObj) == MS::kSuccess)
        bellMatrix = MFnMatrixData(matrixObj).matrix();

    const MMatrix bellMatrixInverse = bellMatrix.inverse();

    int ringSubdivision = 16;
    ringSubdivisionPlug.getValue(ringSubdivision);

    double rVal = 0.0, gVal = 0.01, bVal = 0.11;
    if (drawColorPlug.numChildren() == 3)
    {
        drawColorPlug.child(0).getValue(rVal);
        drawColorPlug.child(1).getValue(gVal);
        drawColorPlug.child(2).getValue(bVal);
    }

    float drawOpacity = 0.3f;
    drawOpacityPlug.getValue(drawOpacity);

    MObject bellMesh;
    outputBellMeshPlug.getValue(bellMesh);

    data->drawData.color = MColor(rVal, gVal, bVal, drawOpacity);
    data->drawData.bellMesh.update(bellMesh);
    std::vector<MMatrix> ringMatrices;
    data->drawData.collisionPointBellList.clear();
    data->drawData.collisionPointRingList.clear();
    data->drawData.directionLines.clear();

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();
    const Plane bellPlane(bell_translate, bellNormal);

    unsigned int numRings = ringMatrixPlug.numElements();
    for (unsigned int i = 0; i < numRings; i++)
    {
        MPlug ringMatrixElementPlug = ringMatrixPlug.elementByPhysicalIndex(i);
        MObject ringMatrixObj;
        if (ringMatrixElementPlug.getValue(ringMatrixObj) != MS::kSuccess)
            continue;
        MMatrix ringMatrix = MFnMatrixData(ringMatrixObj).matrix();
        const PreparedBellRing ring(ringMatrix);
        data->drawData.directionLines.append(ring.translation);
        data->drawData.directionLines.append(ring.translation + bellPlane.projectVector(ring.direction));

        MPoint collisionPointBell, collisionPointRing, linePoint;
        if (BellColliderSolver::collisionPoints(bellMatrix, bellMatrixInverse, bellPlane, ring, collisionPointBell, collisionPointRing, linePoint))
        {
            data->drawData.collisionPointRingList.push_back(linePoint);
            data->drawData.collisionPointRingList.push_back(collisionPointRing);
            data->drawData.collisionPointBellList.push_back(collisionPointBell);
        }

        ringMatrices.push_back(ringMatrix);
    }

    data->drawData.rings.update(ringMatrices, ringSubdivision);
    return data;
}

void BellColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath, 
    MHWRender::MUIDrawManager& drawManager, 
    const MHWRender::MFrameContext& frameContext, 
    const MUserData* data)
{
    auto* bellColliderData = dynamic_cast<const BellColliderDrawData*>(data);

    if (bellColliderData)
    {
        drawManager.beginDrawable();

        const auto& drawData = bellColliderData->drawData;
        if (drawData.bellMesh.points.length())
        {
            drawData.bellMesh.geometry.draw(drawManager, drawData.color, MColor(0, 0, 0));
            drawData.rings.geometry.draw(drawManager, drawData.color * 0.5, MColor(0, 0, 0));

            drawManager.setColor(MColor(0, 0, 0));
            drawManager.setPointSize(5);

            for (const auto& p : drawData.collisionPointBellList)
                drawManager.point(p);

            for (const auto& p : drawData.collisionPointRingList)
                drawManager.point(p);

            if (drawData.directionLines.length())
                drawManager.lineList(drawData.directionLines, false);
        }

        drawManager.endDrawable();
    }
}
