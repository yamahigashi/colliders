#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>

#include <cstddef>
#include <cmath>
#include <vector>

#include "skirtCollideDeformer.h"
#include "pluginIdentity.h"
#include "utils.hpp"

MTypeId SkirtCollideDeformer::typeId(PluginIdentity::kSkirtCollideTypeId);

MObject SkirtCollideDeformer::attr_bellMatrix;
MObject SkirtCollideDeformer::attr_leftHipMatrix;
MObject SkirtCollideDeformer::attr_leftKneeMatrix;
MObject SkirtCollideDeformer::attr_leftHeelMatrix;
MObject SkirtCollideDeformer::attr_rightHipMatrix;
MObject SkirtCollideDeformer::attr_rightKneeMatrix;
MObject SkirtCollideDeformer::attr_rightHeelMatrix;

MObject SkirtCollideDeformer::attr_skirtType;
MObject SkirtCollideDeformer::attr_ringScale;
MObject SkirtCollideDeformer::attr_leftRingAxis;
MObject SkirtCollideDeformer::attr_rightRingAxis;
MObject SkirtCollideDeformer::attr_collision;
MObject SkirtCollideDeformer::attr_falloff;
MObject SkirtCollideDeformer::attr_endFade;

struct RingVolume
{
    MMatrix matrix;
    MMatrix inverse;
};

static bool isFiniteMatrix(const MMatrix& matrix)
{
    for (unsigned int row = 0; row < 4; row++)
    {
        for (unsigned int column = 0; column < 4; column++)
        {
            if (!std::isfinite(matrix[row][column]))
                return false;
        }
    }
    return true;
}

static bool isFinitePoint(const MPoint& point)
{
    return std::isfinite(point.x) && std::isfinite(point.y)
        && std::isfinite(point.z) && std::isfinite(point.w);
}

static double smooth01(double t)
{
    if (t <= 0.0)
        return 0.0;
    if (t >= 1.0)
        return 1.0;
    return t * t * (3.0 - 2.0 * t);
}

static void appendRingVolume(const MMatrix& ringMatrix, std::vector<RingVolume>& ringVolumes)
{
    const double determinant = ringMatrix.det4x4();
    if (!std::isfinite(determinant) || std::fabs(determinant) < 1e-8)
        return;

    const MMatrix ringMatrixInverse = ringMatrix.inverse();
    if (!isFiniteMatrix(ringMatrix) || !isFiniteMatrix(ringMatrixInverse))
        return;

    RingVolume ringVolume;
    ringVolume.matrix = ringMatrix;
    ringVolume.inverse = ringMatrixInverse;
    ringVolumes.push_back(ringVolume);
}

MStatus SkirtCollideDeformer::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MStatus stat;

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

    attr_collision = nAttr.create("collision", "collision", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(2.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_collision);

    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0.2f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    attr_endFade = nAttr.create("endFade", "endFade", MFnNumericData::kFloat, 0.1f);
    nAttr.setMin(0.0f);
    nAttr.setMax(0.5f);
    nAttr.setKeyable(true);
    addAttribute(attr_endFade);

    const MObject affects[] = {
        attr_bellMatrix,
        attr_leftHipMatrix, attr_leftKneeMatrix, attr_leftHeelMatrix,
        attr_rightHipMatrix, attr_rightKneeMatrix, attr_rightHeelMatrix,
        attr_skirtType, attr_ringScale, attr_leftRingAxis, attr_rightRingAxis,
        attr_collision, attr_falloff, attr_endFade
    };
    for (const MObject& attr : affects)
    {
        stat = attributeAffects(attr, outputGeom);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    return MS::kSuccess;
}

MStatus SkirtCollideDeformer::deform(MDataBlock& dataBlock, MItGeometry& iter,
    const MMatrix& localToWorldMatrix, unsigned int multiIndex)
{
    MStatus stat;

    const MMatrix bellMatrix = dataBlock.inputValue(attr_bellMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix leftHipMatrix = dataBlock.inputValue(attr_leftHipMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix leftKneeMatrix = dataBlock.inputValue(attr_leftKneeMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix leftHeelMatrix = dataBlock.inputValue(attr_leftHeelMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix rightHipMatrix = dataBlock.inputValue(attr_rightHipMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix rightKneeMatrix = dataBlock.inputValue(attr_rightKneeMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix rightHeelMatrix = dataBlock.inputValue(attr_rightHeelMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    const short skirtType = dataBlock.inputValue(attr_skirtType, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MVector ringScale = dataBlock.inputValue(attr_ringScale, &stat).asVector();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short leftRingAxis = dataBlock.inputValue(attr_leftRingAxis, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short rightRingAxis = dataBlock.inputValue(attr_rightRingAxis, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const float collision = dataBlock.inputValue(attr_collision, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const float falloff = dataBlock.inputValue(attr_falloff, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const float endFade = dataBlock.inputValue(attr_endFade, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const float envelopeValue = dataBlock.inputValue(envelope, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    if (collision == 0.0f || envelopeValue == 0.0f)
        return MS::kSuccess;

    const MPoint LH = taxis(leftHipMatrix);
    const MPoint LK = taxis(leftKneeMatrix);
    const MPoint LHe = taxis(leftHeelMatrix);
    const MPoint RH = taxis(rightHipMatrix);
    const MPoint RK = taxis(rightKneeMatrix);
    const MPoint RHe = taxis(rightHeelMatrix);

    const double L_thigh = ((LK - LH).length() + (RK - RH).length()) * 0.5;
    const double L_calf = ((LHe - LK).length() + (RHe - RK).length()) * 0.5;

    const MVector thighScale(ringScale.x, L_thigh * ringScale.y, ringScale.z);
    const MMatrix leftHipToKnee = createRingMatrix(leftHipMatrix, thighScale, leftRingAxis, &LK);
    const MMatrix rightHipToKnee = createRingMatrix(rightHipMatrix, thighScale, rightRingAxis, &RK);

    std::vector<RingVolume> ringVolumes;
    ringVolumes.reserve(6);
    appendRingVolume(leftHipToKnee, ringVolumes);
    appendRingVolume(rightHipToKnee, ringVolumes);

    if (skirtType == 1)
    {
        const double legLen = L_thigh + L_calf;
        const MVector heelScale(ringScale.x, legLen * ringScale.y, ringScale.z);
        const MMatrix leftHipToHeel = createRingMatrix(leftHipMatrix, heelScale, leftRingAxis, &LHe);
        const MMatrix rightHipToHeel = createRingMatrix(rightHipMatrix, heelScale, rightRingAxis, &RHe);
        const MMatrix leftHipToKneeExtended = createRingMatrix(leftHipMatrix, heelScale, leftRingAxis, &LK);
        const MMatrix rightHipToKneeExtended = createRingMatrix(rightHipMatrix, heelScale, rightRingAxis, &RK);

        appendRingVolume(leftHipToKneeExtended, ringVolumes);
        appendRingVolume(rightHipToKneeExtended, ringVolumes);
        appendRingVolume(leftHipToHeel, ringVolumes);
        appendRingVolume(rightHipToHeel, ringVolumes);
    }

    const MPoint bellPosition = taxis(bellMatrix);
    const MVector rawBellAxis = getAxis(bellMatrix, 1);
    const MVector bellAxisDirection = rawBellAxis.length() > 1e-8 ? rawBellAxis.normal() : MVector(0.0, -1.0, 0.0);

    const MMatrix worldToLocalMatrix = localToWorldMatrix.inverse();
    for (; !iter.isDone(); iter.next())
    {
        const float pointWeight = weightValue(dataBlock, multiIndex, iter.index());
        const double strength = (double)collision * (double)envelopeValue * (double)pointWeight;
        if (!std::isfinite(strength) || strength == 0.0)
            continue;

        MPoint pointWorld = iter.position() * localToWorldMatrix;
        bool changed = false;

        // Push direction reference: the skirt cone's radial direction from the bell
        // axis (bellMatrix Y). Using the point's own leg-radial direction flips the
        // push to the far side once a point tunnels past the leg axis.
        const MVector bellRadialWorld = [&]() {
            const MVector offset = pointWorld - bellPosition;
            const MVector radial = offset - bellAxisDirection * (offset * bellAxisDirection);
            return radial.length() > 1e-8 ? radial.normal() : MVector(0.0, 0.0, 0.0);
        }();

        for (std::size_t i = 0; i < ringVolumes.size(); i++)
        {
            const RingVolume& ringVolume = ringVolumes[i];
            MPoint pointRing = pointWorld * ringVolume.inverse;
            // The ring frame spans hip (y=0) to target (y=1); the mirrored y<0 lobe
            // is above the hip and must never collide (it grabbed the waist rows).
            // Both ends fade smoothly (crossing a hard y gate popped the correction
            // on and off); the far end fades OUTSIDE (0,1] so the hem keeps full
            // coverage at y <= 1.
            const double fade = (double)endFade;
            const double y = pointRing.y;
            if (y <= 0.0)
                continue;
            double axialWeight = 1.0;
            if (fade > 1e-8)
            {
                if (y >= 1.0 + fade)
                    continue;
                if (y < fade)
                    axialWeight = smooth01(y / fade);
                else if (y > 1.0)
                    axialWeight = 1.0 - smooth01((y - 1.0) / fade);
            }
            else if (y > 1.0)
                continue;

            // Cylinder volume (constant unit radius along the segment): the previous
            // spherical taper vanished at y=1, leaving the hem uncorrected.
            const double r_needed = 1.0;
            const double band = (double)falloff;

            const double r_xz = std::sqrt(pointRing.x * pointRing.x + pointRing.z * pointRing.z);
            if (r_xz >= r_needed + band)
                continue;
            // Outside the cylinder the bulge fades to zero at the radial gate, so
            // entering/leaving the band region never steps.
            double radialWeight = 1.0;
            if (r_xz > r_needed)
            {
                if (band <= 1e-8)
                    continue;
                radialWeight = 1.0 - smooth01((r_xz - r_needed) / band);
            }

            // Direction in ring local XZ, from the bell-radial reference; fall back to
            // the point's own radial direction when the reference degenerates.
            MVector pushDirection = bellRadialWorld * ringVolume.inverse;
            pushDirection.y = 0.0;
            if (pushDirection.length() > 1e-8)
                pushDirection.normalize();
            else if (r_xz > 1e-8)
                pushDirection = MVector(pointRing.x / r_xz, 0.0, pointRing.z / r_xz);
            else
                pushDirection = MVector(1.0, 0.0, 0.0);

            // Decompose the XZ offset into the push direction and its perpendicular.
            // The exit lies on the CHORD at s_boundary = sqrt(R^2 - q^2), not at the
            // full radius: pushing every point to s = R shoved grazing contacts
            // (|q| near R) a full radius sideways and hurled near-outside points
            // placed perpendicular to the push direction - the boundary rattling.
            const double s_current = pointRing.x * pushDirection.x + pointRing.z * pushDirection.z;
            const double q_x = pointRing.x - pushDirection.x * s_current;
            const double q_z = pointRing.z - pushDirection.z * s_current;
            const double q_sq = q_x * q_x + q_z * q_z;
            const double s_boundary = q_sq >= r_needed * r_needed
                ? 0.0
                : std::sqrt(r_needed * r_needed - q_sq);

            // Signed extent along the push direction; tunneled points are negative and
            // land back on the near-side chord exit. C1 soft clamp as before.
            double s_target;
            if (band > 1e-8)
            {
                if (s_current >= s_boundary + band)
                    continue;
                if (s_current <= s_boundary - band)
                    s_target = s_boundary;
                else
                {
                    const double t = s_current - s_boundary + band;
                    s_target = s_boundary + (t * t) / (4.0 * band);
                }
            }
            else
            {
                if (s_current >= s_boundary)
                    continue;
                s_target = s_boundary;
            }

            const double shift = (s_target - s_current) * strength * axialWeight * radialWeight;
            pointRing.x += pushDirection.x * shift;
            pointRing.z += pushDirection.z * shift;

            const MPoint pushedWorld = pointRing * ringVolume.matrix;
            if (!isFinitePoint(pushedWorld))
                continue;

            pointWorld = pushedWorld;
            changed = true;
        }

        if (changed)
        {
            const MPoint pointObject = pointWorld * worldToLocalMatrix;
            if (!isFinitePoint(pointObject))
                continue;

            stat = iter.setPosition(pointObject);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
        }
    }

    return MS::kSuccess;
}
