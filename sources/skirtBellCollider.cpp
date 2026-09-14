#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MPlug.h>
#include <maya/MPointArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MVector.h>
#include <maya/MPoint.h>
#include <maya/MGlobal.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MRampAttribute.h>
#include <maya/MUIDrawManager.h>
#include <maya/MFnMatrixData.h>

#include <vector>
#include <cmath>

#include "skirtBellCollider.h"
#include "colliderInputValidation.h"
#include "skirtLegProfile.h"
#include "pluginIdentity.h"
#include "bellColliderSolver.h"
#include "utils.hpp"
#include "skirtRingFrames.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MTypeId SkirtBellCollider::typeId(PluginIdentity::kSkirtBellTypeId);
MString SkirtBellCollider::typeName(PluginIdentity::kSkirtBellNodeName);

MString SkirtBellCollider::drawDbClassification = MString("drawdb/geometry/") + PluginIdentity::kSkirtBellNodeName;
MString SkirtBellCollider::drawRegistrantId = PluginIdentity::kDrawRegistrant;

MObject SkirtBellCollider::attr_bellMatrix;
MObject SkirtBellCollider::attr_leftHipMatrix;
MObject SkirtBellCollider::attr_leftKneeMatrix;
MObject SkirtBellCollider::attr_leftHeelMatrix;
MObject SkirtBellCollider::attr_rightHipMatrix;
MObject SkirtBellCollider::attr_rightKneeMatrix;
MObject SkirtBellCollider::attr_rightHeelMatrix;

MObject SkirtBellCollider::attr_skirtType;
MObject SkirtBellCollider::attr_height;
MObject SkirtBellCollider::attr_ringScale;
MObject SkirtBellCollider::attr_thighRadiusX;
MObject SkirtBellCollider::attr_thighRadiusZ;
MObject SkirtBellCollider::attr_kneeRadiusX;
MObject SkirtBellCollider::attr_kneeRadiusZ;
MObject SkirtBellCollider::attr_calfRadiusX;
MObject SkirtBellCollider::attr_calfRadiusZ;
MObject SkirtBellCollider::attr_ankleRadiusX;
MObject SkirtBellCollider::attr_ankleRadiusZ;
MObject SkirtBellCollider::attr_thighPosition;
MObject SkirtBellCollider::attr_calfPosition;

MObject SkirtBellCollider::attr_bellScale;
MObject SkirtBellCollider::attr_bellSubdivision;
MObject SkirtBellCollider::attr_ringSubdivision;
MObject SkirtBellCollider::attr_falloff;
MObject SkirtBellCollider::attr_collision;
MObject SkirtBellCollider::attr_tightness;
MObject SkirtBellCollider::attr_smoothness;
MObject SkirtBellCollider::attr_follow;
MObject SkirtBellCollider::attr_bellScaleRamp;
MObject SkirtBellCollider::attr_leftRingAxis;
MObject SkirtBellCollider::attr_rightRingAxis;
MObject SkirtBellCollider::attr_bellAxis;

MObject SkirtBellCollider::attr_outputSurface;

static MPointArray getCurvePoints(const MPointArray& points, int bellSubdivision, bool use_bottom)
{
    MPointArray cvs;
    int start = use_bottom ? 1 : bellSubdivision + 1;
    int count = bellSubdivision;
    cvs.setLength(count + 3);
    for (int i = 0; i < count; i++)
    {
        cvs.set(points[start + i], i);
    }
    // Overlap the first 3 CVs for cubic periodic continuity
    cvs.set(cvs[0], count);
    cvs.set(cvs[1], count + 1);
    cvs.set(cvs[2], count + 2);
    return cvs;
}

void SkirtBellCollider::postConstructor()
{
    MObject thisObj = thisMObject();
    MRampAttribute rampAttr(thisObj, attr_bellScaleRamp);

    MFloatArray positions;
    MFloatArray values;
    MIntArray interpolations;

    positions.append(0.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    positions.append(1.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    rampAttr.addEntries(positions, values, interpolations);
}

MStatus SkirtBellCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;
    MStatus stat;

    // Joint matrices
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

    // Skirt Type
    attr_skirtType = eAttr.create("skirtType", "skirtType", 1);
    eAttr.addField("Short", 0);
    eAttr.addField("Long", 1);
    eAttr.setKeyable(true);
    addAttribute(attr_skirtType);


    // Height
    attr_height = nAttr.create("height", "height", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.01f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_height);

    // Ring Scale
    attr_ringScale = nAttr.create("ringScale", "ringScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.5, 1.0, 0.5);
    nAttr.setKeyable(true);
    addAttribute(attr_ringScale);

    // Bell Scale
    attr_bellScale = nAttr.create("bellScale", "bellScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.8, 1.0, 0.8);
    nAttr.setKeyable(true);
    addAttribute(attr_bellScale);

    // Bell Subdivision
    attr_bellSubdivision = nAttr.create("bellSubdivision", "bellSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_bellSubdivision);

    // Ring Subdivision
    attr_ringSubdivision = nAttr.create("ringSubdivision", "ringSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(ColliderInput::kMinSubdivision);
    nAttr.setMax(ColliderInput::kMaxSubdivision);
    nAttr.setKeyable(true);
    addAttribute(attr_ringSubdivision);

    // Falloff
    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(-1.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    // Collision
    attr_collision = nAttr.create("collision", "collision", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_collision);

    // Tightness
    attr_tightness = nAttr.create("tightness", "tightness", MFnNumericData::kFloat, 0.5f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_tightness);

    attr_smoothness = nAttr.create("smoothness", "smoothness", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_smoothness);

    attr_follow = nAttr.create("follow", "follow", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_follow);

    // Bell Scale Ramp (Curve)
    attr_bellScaleRamp = MRampAttribute::createCurveRamp("bellScaleRamp", "bellScaleRamp");
    addAttribute(attr_bellScaleRamp);

    // Left Leg Axis
    attr_leftRingAxis = eAttr.create("leftRingAxis", "leftRingAxis", 0);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_leftRingAxis);

    // Right Leg Axis
    attr_rightRingAxis = eAttr.create("rightRingAxis", "rightRingAxis", 3);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_rightRingAxis);

    // Skirt Axis
    attr_bellAxis = eAttr.create("bellAxis", "bellAxis", 1);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_bellAxis);

    // Output NURBS Surface (Visible to users)
    attr_outputSurface = tAttr.create("outputSurface", "outputSurface", MFnData::kNurbsSurface);
    tAttr.setWritable(false);
    tAttr.setStorable(false);
    addAttribute(attr_outputSurface);

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

    // Set up attribute affects relationships
    const MObject affects[] = {
        attr_thighRadiusX, attr_thighRadiusZ, attr_kneeRadiusX, attr_kneeRadiusZ,
        attr_calfRadiusX, attr_calfRadiusZ, attr_ankleRadiusX, attr_ankleRadiusZ,
        attr_thighPosition, attr_calfPosition,
        attr_bellMatrix, attr_leftHipMatrix, attr_leftKneeMatrix, attr_leftHeelMatrix,
        attr_rightHipMatrix, attr_rightKneeMatrix, attr_rightHeelMatrix, attr_skirtType,
        attr_height, attr_ringScale, attr_bellScale, attr_bellSubdivision,
        attr_falloff, attr_collision, attr_tightness, attr_smoothness, attr_follow,
        attr_bellScaleRamp, attr_leftRingAxis, attr_rightRingAxis, attr_bellAxis
    };
    for (const MObject& attr : affects) {
        attributeAffects(attr, attr_outputSurface);
    }

    return MS::kSuccess;
}

MStatus SkirtBellCollider::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    if (plug != attr_outputSurface)
        return MS::kUnknownParameter;

    MStatus stat;
    const double kFollowDamping = 0.5;

    // Get input values
    const MMatrix inputBellMatrix = dataBlock.inputValue(attr_bellMatrix).asMatrix();
    const MMatrix leftHipMatrix = dataBlock.inputValue(attr_leftHipMatrix).asMatrix();
    const MMatrix leftKneeMatrix = dataBlock.inputValue(attr_leftKneeMatrix).asMatrix();
    const MMatrix leftHeelMatrix = dataBlock.inputValue(attr_leftHeelMatrix).asMatrix();
    const MMatrix rightHipMatrix = dataBlock.inputValue(attr_rightHipMatrix).asMatrix();
    const MMatrix rightKneeMatrix = dataBlock.inputValue(attr_rightKneeMatrix).asMatrix();
    const MMatrix rightHeelMatrix = dataBlock.inputValue(attr_rightHeelMatrix).asMatrix();

    const short skirtType = dataBlock.inputValue(attr_skirtType).asShort();
    const float height = dataBlock.inputValue(attr_height).asFloat();
    const MVector ringScale = dataBlock.inputValue(attr_ringScale).asVector();
    const MVector bellScale = dataBlock.inputValue(attr_bellScale).asVector();
    const int bellSubdivision = dataBlock.inputValue(attr_bellSubdivision).asInt();
    const float falloff = dataBlock.inputValue(attr_falloff).asFloat();
    const float collision = dataBlock.inputValue(attr_collision).asFloat();
    const float tightness = dataBlock.inputValue(attr_tightness).asFloat();
    const float smoothness = dataBlock.inputValue(attr_smoothness).asFloat();
    const float follow = dataBlock.inputValue(attr_follow).asFloat();
    const short leftRingAxis = dataBlock.inputValue(attr_leftRingAxis).asShort();
    const short rightRingAxis = dataBlock.inputValue(attr_rightRingAxis).asShort();
    const short bellAxis = dataBlock.inputValue(attr_bellAxis).asShort();

    if (!ColliderInput::validSubdivision(bellSubdivision))
    {
        MGlobal::displayError("SkirtBellCollider: bellSubdivision must be between 3 and 4096.");
        return MS::kInvalidParameter;
    }
    if (!ColliderInput::validSkirtType(skirtType))
    {
        MGlobal::displayError("SkirtBellCollider: skirtType must be 0 or 1.");
        return MS::kInvalidParameter;
    }
    if (!ColliderInput::validAxis(leftRingAxis) || !ColliderInput::validAxis(rightRingAxis) ||
        !ColliderInput::validAxis(bellAxis))
    {
        MGlobal::displayError("SkirtBellCollider: leftRingAxis, rightRingAxis and bellAxis must be between 0 and 5.");
        return MS::kInvalidParameter;
    }

    // Ramp attribute
    MRampAttribute rampAttr(thisMObject(), attr_bellScaleRamp);

    // Midpoints
    const MPoint LH = taxis(leftHipMatrix);
    const MPoint RH = taxis(rightHipMatrix);

    const MPoint H = (LH + RH) * 0.5;

    const MPoint W = taxis(inputBellMatrix);

    const SkirtRingFrames ringFrames(leftHipMatrix, leftKneeMatrix, leftHeelMatrix,
                                    rightHipMatrix, rightKneeMatrix, rightHeelMatrix,
                                    ringScale, leftRingAxis, rightRingAxis, skirtType == 1);
    const double L_thigh = ringFrames.thighLength;
    const double L_calf = ringFrames.calfLength;

    const SkirtLegProfile profile(
        dataBlock.inputValue(attr_thighRadiusX).asDouble(),
        dataBlock.inputValue(attr_thighRadiusZ).asDouble(),
        dataBlock.inputValue(attr_kneeRadiusX).asDouble(),
        dataBlock.inputValue(attr_kneeRadiusZ).asDouble(),
        dataBlock.inputValue(attr_calfRadiusX).asDouble(),
        dataBlock.inputValue(attr_calfRadiusZ).asDouble(),
        dataBlock.inputValue(attr_ankleRadiusX).asDouble(),
        dataBlock.inputValue(attr_ankleRadiusZ).asDouble(),
        dataBlock.inputValue(attr_thighPosition).asDouble(),
        dataBlock.inputValue(attr_calfPosition).asDouble(),
        L_thigh, L_calf);

    const double d_hip = (H - W).length();
    const double d_knee = d_hip + L_thigh;
    const double d_heel = d_knee + L_calf;
    const double d_mid = (d_hip + d_knee) * 0.5;

    double h_param = (double)height;
    if (h_param < 0.0) h_param = 0.0;
    if (h_param > 1.0) h_param = 1.0;

    double h_val = 0.0;
    if (skirtType == 1) // Long skirt: bottom interpolates from Knee to Heel
    {
        h_val = d_knee + (d_heel - d_knee) * h_param;
    }
    else // Short skirt: bottom interpolates from Hip to Knee
    {
        h_val = d_hip + (d_knee - d_hip) * h_param;
    }

    const double raw_defaultHeight = (skirtType == 1) ? d_heel : d_knee;
    const double defaultHeight = raw_defaultHeight < 1e-5 ? 1.0 : raw_defaultHeight;
    const double s = h_val / defaultHeight;

    // Get direction vector based on the selected bell matrix axis
    const MVector raw_dir_vector = getAxis(inputBellMatrix, bellAxis);
    const MVector dir_vector = raw_dir_vector.length() < 1e-4 ? MVector(0, -1, 0) : raw_dir_vector.normal();

    // Start point of the skirt aligned rigidly along the chosen waist axis
    const MPoint P_start = W;

    const int N = (skirtType == 0) ? 2 : 3;

    // Setup level distances along the skirt axis
    vector<double> levelDistances;
    levelDistances.push_back(0.0);
    if (skirtType == 1) // Long skirt: height controls the heel level only
    {
        levelDistances.push_back(d_mid);
        levelDistances.push_back(d_knee);
        levelDistances.push_back(h_val);
    }
    else // Short skirt: height scales all levels
    {
        levelDistances.push_back(d_mid * s);
        levelDistances.push_back(d_knee * s);
    }

    vector<vector<PreparedBellRing>> bellRings(N + 1);
    vector<vector<PreparedBellRing>> ringsWithoutKnee(N + 1);
    vector<vector<PreparedBellRing>> ringsWithKnee(N + 1);
    for (int r = 1; r <= N; ++r)
    {
        const double s_level = profile.levelParameter(levelDistances[r], d_hip);
        auto prepare = [&](const MMatrix& base, SkirtLegProfile::Ring kind) {
            const auto radius = profile.forRing(s_level, kind, ringScale.y);
            MMatrix matrix = base;
            for (unsigned int c = 0; c < 4; ++c)
            {
                matrix[0][c] *= radius.x;
                matrix[2][c] *= radius.z;
            }
            return PreparedBellRing(matrix);
        };
        bellRings[r].push_back(prepare(ringFrames.leftKnee, SkirtLegProfile::Ring::Knee));
        bellRings[r].push_back(prepare(ringFrames.rightKnee, SkirtLegProfile::Ring::Knee));
        if (skirtType == 1)
        {
            ringsWithoutKnee[r].push_back(prepare(ringFrames.leftHeel, SkirtLegProfile::Ring::Heel));
            ringsWithoutKnee[r].push_back(prepare(ringFrames.rightHeel, SkirtLegProfile::Ring::Heel));
            bellRings[r].insert(bellRings[r].end(), ringsWithoutKnee[r].begin(), ringsWithoutKnee[r].end());
            ringsWithKnee[r] = ringsWithoutKnee[r];
            ringsWithKnee[r].push_back(prepare(ringFrames.leftExtended, SkirtLegProfile::Ring::Extended));
            ringsWithKnee[r].push_back(prepare(ringFrames.rightExtended, SkirtLegProfile::Ring::Extended));
        }
    }

    vector<MPointArray> rows(N + 1);
    vector<MVector> deltaDirect(N, MVector(0, 0, 0));
    MPointArray controlPoints;
    MDoubleArray uKnots;
    MDoubleArray vKnots;

    // Setup uKnots: uniform for kPeriodic form (degree 3)
    int numKnotsInU = bellSubdivision + 5;
    for (int j = 0; j < numKnotsInU; j++)
    {
        uKnots.append((double)(j - 2));
    }

    const double h_safe = h_val < 1e-5 ? 1e-5 : h_val;

    // Bell orientation depends only on the waist matrix and the skirt axis, so
    // every stage of this evaluation shares one frame (bulletproof against
    // zero/unconnected inputs).
    const MVector dir_y = dir_vector;

    const MVector raw_waist_X = xaxis(inputBellMatrix);
    const MVector waist_X = raw_waist_X.length() < 1e-4 ? MVector(1, 0, 0) : raw_waist_X.normal();

    const MVector raw_X = (waist_X - (waist_X * dir_y) * dir_y);
    const MVector X = [&]() {
        if (raw_X.length() < 1e-4) {
            const MVector raw_waist_Z = zaxis(inputBellMatrix);
            const MVector waist_Z = raw_waist_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_waist_Z.normal();
            const MVector crossed = dir_y ^ waist_Z;
            if (crossed.length() < 1e-4) {
                return MVector(1, 0, 0);
            }
            return crossed.normal();
        }
        return raw_X.normal();
    }();

    const MVector raw_Z = X ^ dir_y;
    const MVector Z = raw_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_Z.normal();

    // Every stage uses the same subdivision, so the circle samples are shared.
    const BellCircleTable circle(bellSubdivision);

    // Solve for each bell
    for (int i = 0; i < N; i++)
    {
        const double dist = levelDistances[i+1] - levelDistances[i];
        const double safeDist = dist < 1e-4 ? 1e-4 : dist;
        const MPoint P_bell = P_start + dir_vector * levelDistances[i];

        // Query Ramp values using actual distance ratios
        const float t_bottom = (float)(levelDistances[i] / h_safe);
        const float t_top = (float)(levelDistances[i+1] / h_safe);
        float raw_scale_bottom = 1.0f;
        float raw_scale_top = 1.0f;
        rampAttr.getValueAtPosition(t_bottom, raw_scale_bottom);
        rampAttr.getValueAtPosition(t_top, raw_scale_top);

        const float scale_bottom = raw_scale_bottom;
        const float scale_top = raw_scale_top < 1e-5f ? 1e-5f : raw_scale_top;

        // Scale axes using bellScale attributes
        const MVector scaled_X = X * (bellScale.x * scale_top);
        const MVector Y = dir_y * (safeDist * bellScale.y);
        const MVector scaled_Z = Z * (bellScale.z * scale_top);

        const double m[4][4] = {
            {scaled_X.x, scaled_X.y, scaled_X.z, 0.0},
            {Y.x, Y.y, Y.z, 0.0},
            {scaled_Z.x, scaled_Z.y, scaled_Z.z, 0.0},
            {P_bell.x, P_bell.y, P_bell.z, 1.0}
        };
        const MMatrix bellMatrix(m);

        MPointArray baseBellPoints = BellColliderSolver::makeBellPoints(bellMatrix, 1, circle, 1, scale_bottom / scale_top, 1);
        BellColliderSolver::roundMeshPoints(baseBellPoints);

        auto solveForRings = [&](const vector<PreparedBellRing>& rings, MPointArray& outBottom, MPointArray& outTop, MVector& outMeanDisplacement) -> MStatus {
            BellColliderInputs inputs;
            inputs.bellMatrix = bellMatrix;
            inputs.rings = rings;
            inputs.bellSubdivision = bellSubdivision;
            inputs.falloff = falloff;
            inputs.collision = collision;
            inputs.smoothness = smoothness;
            inputs.followGain = follow * (levelDistances[i + 1] / h_safe) * kFollowDamping;

            BellColliderOutputs outputs;
            const MStatus solveStat = BellColliderSolver::solve(inputs, baseBellPoints, outputs);
            if (solveStat != MS::kSuccess)
            {
                MGlobal::displayError("SkirtBellCollider: Solver failed at bell index " + MString(to_string(i).c_str()) + " with: " + solveStat.errorString());
                return solveStat;
            }

            MPointArray& meshPoints = outputs.points;
            BellColliderSolver::roundMeshPoints(meshPoints);

            if (i == 0) outBottom = getCurvePoints(meshPoints, bellSubdivision, true);
            outTop = getCurvePoints(meshPoints, bellSubdivision, false);
            outMeanDisplacement = outputs.meanDisplacement;

            return MS::kSuccess;
        };

        if (skirtType == 1 && i == 2)
        {
            MPointArray topWith, topWithout;
            MVector directWith(0, 0, 0);
            MVector directWithout(0, 0, 0);

            if (tightness < 1.0f) {
                MPointArray dummy;
                MStatus s = solveForRings(ringsWithKnee[i + 1], dummy, topWith, directWith);
                if (s != MS::kSuccess) return s;
            }
            if (tightness > 0.0f) {
                MPointArray dummy;
                MStatus s = solveForRings(ringsWithoutKnee[i + 1], dummy, topWithout, directWithout);
                if (s != MS::kSuccess) return s;
            }

            if (tightness <= 0.0f) {
                rows[i + 1] = topWith;
                deltaDirect[i] = directWith;
            } else if (tightness >= 1.0f) {
                rows[i + 1] = topWithout;
                deltaDirect[i] = directWithout;
            } else {
                MPointArray blendedTop;
                blendedTop.setLength(topWith.length());
                for (unsigned int j = 0; j < topWith.length(); j++) {
                    blendedTop.set(topWith[j] * (1.0 - tightness) + topWithout[j] * tightness, j);
                }
                rows[i + 1] = blendedTop;
                deltaDirect[i] = directWith * (1.0 - tightness) + directWithout * tightness;
            }
        }
        else
        {
            MPointArray bottom, top;
            MVector directDisplacement;
            MStatus s = solveForRings(bellRings[i + 1], bottom, top, directDisplacement);
            if (s != MS::kSuccess) return s;
            if (i == 0) rows[0] = bottom;
            rows[i + 1] = top;
            deltaDirect[i] = directDisplacement;
        }
    }

    // Layer 3 runs only when it adds displacement; re-relaxing untouched rows would deepen the residual to (1-c)^2.
    if (follow > 0.0f)
    {
        MVector cumulativeDirect(0, 0, 0);
        for (int r = 2; r <= N; r++)
        {
            cumulativeDirect += deltaDirect[r - 2];
            const double rowTightness = skirtType == 1 && r == N ? tightness : 0.0;
            const MVector propagation = cumulativeDirect * (follow * (1.0 - rowTightness));

            for (int j = 0; j < bellSubdivision; j++)
                rows[r][j] += propagation;

            if (skirtType == 1 && r == N)
            {
                BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithoutKnee[r][0], collision, 0, bellSubdivision);
                BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithoutKnee[r][1], collision, 0, bellSubdivision);

                if (tightness < 1.0f)
                {
                    const double kneeCollision = collision * (1.0 - tightness);
                    BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithKnee[r][2], kneeCollision, 0, bellSubdivision);
                    BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithKnee[r][3], kneeCollision, 0, bellSubdivision);
                }
            }
            else
            {
                BellColliderSolver::relaxTowardRingBoundary(rows[r], bellRings[r][0], collision, 0, bellSubdivision);
                BellColliderSolver::relaxTowardRingBoundary(rows[r], bellRings[r][1], collision, 0, bellSubdivision);

                if (skirtType == 1)
                {
                    BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithoutKnee[r][0], collision, 0, bellSubdivision);
                    BellColliderSolver::relaxTowardRingBoundary(rows[r], ringsWithoutKnee[r][1], collision, 0, bellSubdivision);
                }
            }

            // Solver mesh top CVs start at bellSubdivision + 1; surface rows start at 0 and keep three periodic duplicates.
            for (int j = 0; j < 3; j++)
                rows[r].set(rows[r][j], bellSubdivision + j);
        }
    }

    // Populate vKnots matching actual distance levels
    for (int i = 0; i <= N; i++)
    {
        vKnots.append(levelDistances[i]);
    }

    // Populate controlPoints using the transposed indexing layout required by Maya:
    // index = (numCVsInV * uIndex) + vIndex
    int numU = bellSubdivision + 3;
    int numV = N + 1;
    controlPoints.setLength(numU * numV);
    for (int u = 0; u < numU; u++)
    {
        for (int v = 0; v < numV; v++)
        {
            controlPoints.set(rows[v][numU - 1 - u], u * numV + v);
        }
    }

    // Create NURBS surface
    MFnNurbsSurfaceData surfaceDataFn;
    MObject surfaceData = surfaceDataFn.create(&stat);
    if (stat != MS::kSuccess)
    {
        MGlobal::displayError("SkirtBellCollider: Failed to create NURBS surface data object: " + stat.errorString());
        return stat;
    }

    MFnNurbsSurface surfaceFn;
    surfaceFn.create(
        controlPoints,
        uKnots,
        vKnots,
        3, // uDegree (cubic periodic)
        1, // vDegree (linear open)
        MFnNurbsSurface::kPeriodic, // uForm
        MFnNurbsSurface::kOpen, // vForm
        false, // createRational
        surfaceData,
        &stat
    );
    if (stat != MS::kSuccess)
    {
        MGlobal::displayError("SkirtBellCollider: MFnNurbsSurface::create failed: " + stat.errorString());
        return stat;
    }

    // Set to output plug
    MDataHandle outputHandle = dataBlock.outputValue(attr_outputSurface, &stat);
    if (stat != MS::kSuccess)
        return stat;

    outputHandle.setMObject(surfaceData);
    dataBlock.setClean(plug);

    return MS::kSuccess;
}

MUserData* SkirtBellColliderDrawOverride::prepareForDraw(
    const MDagPath& objPath,
    const MDagPath& cameraPath,
    const MHWRender::MFrameContext& frameContext,
    MUserData* oldData)
{
    MStatus stat;
    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto* data = dynamic_cast<SkirtBellColliderDrawData*>(oldData);
    if (!data)
        data = new SkirtBellColliderDrawData();

    // Extract attributes
    auto getMatrix = [&obj](const MObject& attr, MMatrix& outMat) {
        MPlug plug(obj, attr);
        MObject tempObj;
        if (plug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
            outMat = MFnMatrixData(tempObj).matrix();
    };

    MMatrix leftHipMatrix, leftKneeMatrix, leftHeelMatrix;
    MMatrix rightHipMatrix, rightKneeMatrix, rightHeelMatrix;

    getMatrix(SkirtBellCollider::attr_leftHipMatrix, leftHipMatrix);
    getMatrix(SkirtBellCollider::attr_leftKneeMatrix, leftKneeMatrix);
    getMatrix(SkirtBellCollider::attr_leftHeelMatrix, leftHeelMatrix);
    getMatrix(SkirtBellCollider::attr_rightHipMatrix, rightHipMatrix);
    getMatrix(SkirtBellCollider::attr_rightKneeMatrix, rightKneeMatrix);
    getMatrix(SkirtBellCollider::attr_rightHeelMatrix, rightHeelMatrix);

    short skirtType = 1;
    MPlug(obj, SkirtBellCollider::attr_skirtType).getValue(skirtType);

    MVector ringScale(0.5, 1.0, 0.5);
    MPlug ringScalePlug(obj, SkirtBellCollider::attr_ringScale);
    if (ringScalePlug.numChildren() == 3)
    {
        double sx = 0.5, sy = 1.0, sz = 0.5;
        ringScalePlug.child(0).getValue(sx);
        ringScalePlug.child(1).getValue(sy);
        ringScalePlug.child(2).getValue(sz);
        ringScale = MVector(sx, sy, sz);
    }

    int ringSubdivision = 16;
    MPlug(obj, SkirtBellCollider::attr_ringSubdivision).getValue(ringSubdivision);

    short leftRingAxis = 0, rightRingAxis = 0;
    MPlug(obj, SkirtBellCollider::attr_leftRingAxis).getValue(leftRingAxis);
    MPlug(obj, SkirtBellCollider::attr_rightRingAxis).getValue(rightRingAxis);

    if (!ColliderInput::validSkirtType(skirtType) ||
        !ColliderInput::validAxis(leftRingAxis) || !ColliderInput::validAxis(rightRingAxis))
    {
        data->drawData.rings.update({}, 0);
        data->drawData.curves.update(MPointArray(), 0, 0);
        return data;
    }

    const SkirtRingFrames ringFrames(leftHipMatrix, leftKneeMatrix, leftHeelMatrix,
                                    rightHipMatrix, rightKneeMatrix, rightHeelMatrix,
                                    ringScale, leftRingAxis, rightRingAxis, skirtType == 1);
    auto getDouble = [&obj](const MObject& attr, double value) {
        MPlug(obj, attr).getValue(value);
        return value;
    };
    const SkirtLegProfile profile(
        getDouble(SkirtBellCollider::attr_thighRadiusX, 1.0),
        getDouble(SkirtBellCollider::attr_thighRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_kneeRadiusX, 1.0),
        getDouble(SkirtBellCollider::attr_kneeRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_calfRadiusX, 1.0),
        getDouble(SkirtBellCollider::attr_calfRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_ankleRadiusX, 1.0),
        getDouble(SkirtBellCollider::attr_ankleRadiusZ, 1.0),
        getDouble(SkirtBellCollider::attr_thighPosition, 0.5),
        getDouble(SkirtBellCollider::attr_calfPosition, 0.5),
        ringFrames.thighLength, ringFrames.calfLength);
    std::vector<std::array<double, 2>> farMultipliers;
    const auto matrices = ringFrames.visibleMatrices(profile, farMultipliers);
    data->drawData.rings.update(matrices, ringSubdivision, farMultipliers);

    MObject surfaceData;
    MPointArray cvs;
    unsigned int numU = 0, numV = 0;
    if (MPlug(obj, SkirtBellCollider::attr_outputSurface).getValue(surfaceData) == MS::kSuccess && !surfaceData.isNull()) {
        MFnNurbsSurface surfaceFn(surfaceData, &stat);
        if (stat == MS::kSuccess && surfaceFn.getCVs(cvs, MSpace::kObject) == MS::kSuccess) {
            numU = surfaceFn.numCVsInU();
            numV = surfaceFn.numCVsInV();
        }
    }
    data->drawData.curves.update(cvs, numU, numV);

    // Default transparent cyan color for drawing collider rings
    data->drawData.color = MColor(0.0f, 0.6f, 1.0f, 0.25f);

    return data;
}

void SkirtBellColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath,
    MHWRender::MUIDrawManager& drawManager,
    const MHWRender::MFrameContext& frameContext,
    const MUserData* data)
{
    auto* skirtDrawData = dynamic_cast<const SkirtBellColliderDrawData*>(data);
    if (!skirtDrawData)
        return;

    const auto& drawData = skirtDrawData->drawData;

    drawManager.beginDrawable();

    drawData.rings.geometry.draw(drawManager, drawData.color, MColor(0.0f, 0.1f, 0.2f, 1.0f));
    drawManager.setColor(MColor(1.0f, 0.75f, 0.0f, 1.0f));
    drawManager.setLineWidth(2.0f);
    if (drawData.curves.lines.length())
        drawManager.lineList(drawData.curves.lines, false);

    drawManager.endDrawable();
}
