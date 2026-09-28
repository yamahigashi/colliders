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

#include <maya/MArrayDataBuilder.h>
#include <maya/MArrayDataHandle.h>

#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>

#include <tbb/parallel_for.h>

#include <maya/MFnMatrixData.h>

#include "planeCollider.h"
#include "colliderInputValidation.h"
#include "pluginIdentity.h"
#include "utils.hpp"

using namespace std;

MTypeId PlaneCollider::typeId(PluginIdentity::kPlaneTypeId);

MObject PlaneCollider::attr_planeMatrix;
MObject PlaneCollider::attr_normalAxis;
MObject PlaneCollider::attr_inputPosition;
MObject PlaneCollider::attr_drawColor;
MObject PlaneCollider::attr_drawOpacity;
MObject PlaneCollider::attr_outputPosition;

MString PlaneCollider::drawDbClassification = MString("drawdb/geometry/") + PluginIdentity::kPlaneNodeName;
MString PlaneCollider::drawRegistrantId = PluginIdentity::kDrawRegistrant;

MStatus PlaneCollider::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    if (plug != attr_outputPosition && !(plug.isChild() && plug.parent() == attr_outputPosition))
        return MS::kUnknownParameter;

    const MMatrix planeMatrix = dataBlock.inputValue(attr_planeMatrix).asMatrix();
    const short normalAxis = dataBlock.inputValue(attr_normalAxis).asShort();
    if (!ColliderInput::validAxis(normalAxis))
    {
        MGlobal::displayError("PlaneCollider: normalAxis must be between 0 and 5.");
        return MS::kInvalidParameter;
    }

    const MVector inputPosition = dataBlock.inputValue(attr_inputPosition).asVector();

    Plane plane(taxis(planeMatrix), getAxis(planeMatrix, normalAxis));

    MPoint outputPosition = inputPosition;
    if (plane.distance(inputPosition) < 0)
        outputPosition = plane.projectPoint(inputPosition);

    MDataHandle outputHandle = dataBlock.outputValue(attr_outputPosition);
    outputHandle.setMVector(outputPosition);
    outputHandle.setClean();
    dataBlock.setClean(plug);

    return MS::kSuccess;
}

MStatus PlaneCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;

    attr_planeMatrix = mAttr.create("planeMatrix", "planeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_planeMatrix);

    attr_normalAxis = eAttr.create("normalAxis", "normalAxis", 1); // Y by default
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setChannelBox(true);
    addAttribute(attr_normalAxis);

    attr_inputPosition = nAttr.create("inputPosition", "inputPosition", MFnNumericData::k3Double);
    nAttr.setKeyable(true);
    addAttribute(attr_inputPosition);

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

    attr_outputPosition = nAttr.create("outputPosition", "outputPosition", MFnNumericData::k3Double);
    addAttribute(attr_outputPosition);

    attributeAffects(attr_planeMatrix, attr_outputPosition);
    attributeAffects(attr_normalAxis, attr_outputPosition);
    attributeAffects(attr_inputPosition, attr_outputPosition);
    attributeAffects(attr_drawColor, attr_outputPosition);
    attributeAffects(attr_drawOpacity, attr_outputPosition);

    return MS::kSuccess;
}

MUserData* PlaneColliderDrawOverride::prepareForDraw(
    const MDagPath& objPath,
    const MDagPath& cameraPath,
    const MHWRender::MFrameContext& frameContext,
    MUserData* oldData)
{
    MStatus stat;

    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto* data = dynamic_cast<PlaneColliderDrawData*>(oldData);

    if (!data)
        data = new PlaneColliderDrawData();

    MPlug planeMatrixPlug(obj, PlaneCollider::attr_planeMatrix);
    MPlug normalAxisPlug(obj, PlaneCollider::attr_normalAxis);
    MPlug drawColorPlug(obj, PlaneCollider::attr_drawColor);
    MPlug drawOpacityPlug(obj, PlaneCollider::attr_drawOpacity);

    MMatrix planeMatrix;
    MObject matrixObj;
    if (planeMatrixPlug.getValue(matrixObj) == MS::kSuccess)
    {
        MFnMatrixData matrixData(matrixObj);
        planeMatrix = matrixData.matrix();
    }

    short normalAxis = 1;
    normalAxisPlug.getValue(normalAxis);
    if (!ColliderInput::validAxis(normalAxis))
    {
        data->drawData.size = 0.0;
        return data;
    }

    double rVal = 0.0, gVal = 0.01, bVal = 0.11;
    if (drawColorPlug.numChildren() == 3)
    {
        drawColorPlug.child(0).getValue(rVal);
        drawColorPlug.child(1).getValue(gVal);
        drawColorPlug.child(2).getValue(bVal);
    }

    float drawOpacity = 0.3f;
    drawOpacityPlug.getValue(drawOpacity);

    const MVector normal = getAxis(planeMatrix, normalAxis);
    Plane plane(taxis(planeMatrix), normal);

    data->drawData.color = MColor(rVal, gVal, bVal, drawOpacity);
    data->drawData.planeCenter = plane.orig;
    data->drawData.planeNormal = plane.normal;
    data->drawData.size = normal.length();

    return data;
}

void PlaneColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath,
    MHWRender::MUIDrawManager& drawManager,
    const MHWRender::MFrameContext& frameContext,
    const MUserData* data)
{
    auto* PlaneColliderData = dynamic_cast<const PlaneColliderDrawData*>(data);

    if (PlaneColliderData && PlaneColliderData->drawData.size > 0.0)
    {
        drawManager.beginDrawable();
        drawManager.setColor(PlaneColliderData->drawData.color);
        drawManager.circle(PlaneColliderData->drawData.planeCenter, PlaneColliderData->drawData.planeNormal, PlaneColliderData->drawData.size, true);
        drawManager.setLineWidth(3);
        drawManager.setColor(PlaneColliderData->drawData.color * 0.25);
        drawManager.line(PlaneColliderData->drawData.planeCenter, PlaneColliderData->drawData.planeCenter + PlaneColliderData->drawData.planeNormal * PlaneColliderData->drawData.size);
        drawManager.endDrawable();
    }
}
