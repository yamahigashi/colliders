#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnData.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnGenericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MQuaternion.h>
#include <maya/MVector.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "colliderInputValidation.h"
#include "pluginIdentity.h"
#include "skirtCollideDeformer.h"
#include "skirtLegProfile.h"
#include "utils.hpp"

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
struct Station
{
    double z, radiusX, radiusZ;
};

struct Cylinder
{
    MPoint origin;
    MVector axis, x, z;
    double length;
    std::vector<Station> stations;
};

struct CylinderPair
{
    Cylinder current, rest;
};

struct Constraint
{
    MVector normal;
    double distance;
};

bool isFinitePoint(const MPoint& point)
{
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && std::isfinite(point.w);
}

bool validMatrix(const MMatrix& matrix)
{
    for (unsigned int row = 0; row < 4; row++)
        for (unsigned int column = 0; column < 4; column++)
            if (!std::isfinite(matrix[row][column]))
                return false;
    const double determinant = matrix.det4x4();
    return std::isfinite(determinant) && std::fabs(determinant) >= 1e-8;
}

std::vector<double> stationParameters(const SkirtLegProfile& profile)
{
    const double candidates[] = {0.0, profile.knee * profile.thighPosition, profile.knee,
        profile.knee + (1.0 - profile.knee) * profile.calfPosition, 1.0};
    const int precedence[] = {4, 2, 1, 3, 0};
    std::vector<double> parameters;
    for (int index : precedence)
    {
        bool coincident = false;
        for (double s : parameters)
            coincident = coincident || std::fabs(candidates[index] - s) <= 1e-9;
        if (!coincident)
            parameters.push_back(candidates[index]);
    }
    std::sort(parameters.begin(), parameters.end());
    return parameters;
}

bool makeCylinder(const MMatrix& joint, const MMatrix& target, short ringAxis, const MVector& scale,
    const SkirtLegProfile& profile, bool calf, Cylinder& cylinder)
{
    if (!validMatrix(joint) || !validMatrix(target))
        return false;
    cylinder.origin = taxis(joint);
    const MPoint end = taxis(target);
    const MVector direction = end - cylinder.origin;
    cylinder.length = direction.length() * scale.y;
    if (!std::isfinite(cylinder.length) || cylinder.length < 1e-6)
        return false;
    cylinder.axis = direction.normal();
    const MMatrix frame = createRingMatrix(joint, MVector(1.0, 1.0, 1.0), ringAxis, &end);
    if (!validMatrix(frame))
        return false;
    cylinder.x = getAxis(frame, 0).normal();
    cylinder.z = getAxis(frame, 2).normal();
    const double start = calf ? profile.knee : 0.0;
    const double finish = calf ? 1.0 : profile.knee;
    for (double s : stationParameters(profile))
    {
        if (s < start || s > finish)
            continue;
        const SkirtLegProfile::Radius radius = profile.evaluate(s, cylinder.length);
        const Station station = {(s - start) * profile.legLength * scale.y, scale.x * radius.x, scale.z * radius.z};
        if (!std::isfinite(station.z) || !std::isfinite(station.radiusX) || !std::isfinite(station.radiusZ))
            return false;
        cylinder.stations.push_back(station);
    }
    return cylinder.stations.size() >= 2;
}

double radialSupport(const Cylinder& cylinder, const Station& station, const MVector& normal)
{
    const double x = station.radiusX * (normal * cylinder.x);
    const double z = station.radiusZ * (normal * cylinder.z);
    return std::sqrt(x * x + z * z);
}

double support(const Cylinder& cylinder, const MVector& normal)
{
    double result = -std::numeric_limits<double>::infinity();
    for (const Station& station : cylinder.stations)
        result = (std::max)(result, station.z * (normal * cylinder.axis) + radialSupport(cylinder, station, normal));
    return normal * MVector(cylinder.origin) + result;
}

double activationWeight(const Cylinder& cylinder, const MPoint& point, double falloff)
{
    const MVector offset = point - cylinder.origin;
    const double z = offset * cylinder.axis;
    const MVector radial = offset - z * cylinder.axis;
    const double radius = radial.length();
    const MVector direction = radius < 1e-8 ? cylinder.x : radial / radius;
    Station section = cylinder.stations.back();
    if (z <= cylinder.stations.front().z)
        section = cylinder.stations.front();
    else
    {
        for (std::size_t j = 1; j < cylinder.stations.size(); j++)
        {
            const Station& first = cylinder.stations[j - 1];
            const Station& second = cylinder.stations[j];
            if (z <= second.z)
            {
                const double t = (z - first.z) / (second.z - first.z);
                section.radiusX = first.radiusX + t * (second.radiusX - first.radiusX);
                section.radiusZ = first.radiusZ + t * (second.radiusZ - first.radiusZ);
                break;
            }
        }
    }
    const double surfaceRadius = radialSupport(cylinder, section, direction);
    const double distance = (std::max)(radius - surfaceRadius, (std::max)(-z, z - cylinder.length));
    if (distance <= 0.0)
        return 1.0;
    const double band = falloff * surfaceRadius;
    if (band < 1e-8 || distance >= band)
        return 0.0;
    const double t = distance / band;
    return 1.0 - t * t * (3.0 - 2.0 * t);
}

bool pointConstraint(const CylinderPair& pair, const MPoint& restPoint, const MVector& waistRadial, const MPoint& point,
    double falloff, Constraint& constraint)
{
    const Cylinder& rest = pair.rest;
    const Cylinder& current = pair.current;
    const double weight = activationWeight(current, point, falloff);
    if (weight == 0.0)
        return false;
    MVector outward = waistRadial - rest.axis * (waistRadial * rest.axis);
    if (outward.length() < 1e-6)
    {
        const MVector offset = restPoint - rest.origin;
        outward = offset - rest.axis * (offset * rest.axis);
        if (outward.length() < 1e-8)
            return false;
    }
    outward.normalize();
    double radius = 0.0;
    double maxRadius = 0.0;
    for (const Station& station : rest.stations)
    {
        radius = (std::max)(radius, radialSupport(rest, station, outward));
        maxRadius = (std::max)(maxRadius, (std::max)(station.radiusX, station.radiusZ));
    }
    const double clearance = (std::max)(outward * (restPoint - rest.origin) - radius, 1e-3 * maxRadius);
    const MVector transported =
        rest.axis * current.axis < -1.0 + 1e-6 ? outward : outward.rotateBy(MQuaternion(rest.axis, current.axis));
    const double z = current.axis * (point - current.origin);
    if (z <= 0.0)
        constraint.normal = -current.axis;
    else if (z >= current.length)
        constraint.normal = current.axis;
    else
        constraint.normal =
            (z * (current.length - z) * transported + clearance * (2.0 * z - current.length) * current.axis).normal();
    constraint.distance = weight * (support(current, constraint.normal) - constraint.normal * MVector(point));
    return std::isfinite(constraint.distance);
}

double determinant(const double matrix[3][3], int size)
{
    if (size == 1)
        return matrix[0][0];
    if (size == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

MVector minimumDisplacement(const std::vector<Constraint>& constraints)
{
    MVector best(0.0, 0.0, 0.0);
    double bestSquared = std::numeric_limits<double>::infinity();
    const unsigned int count = static_cast<unsigned int>(constraints.size());
    for (unsigned int mask = 0; mask < (1u << count); mask++)
    {
        int selected[3];
        int size = 0;
        for (unsigned int i = 0; i < count; i++)
        {
            if ((mask & (1u << i)) == 0)
                continue;
            if (size < 3)
                selected[size] = static_cast<int>(i);
            size++;
        }
        if (size > 3)
            continue;
        MVector delta(0.0, 0.0, 0.0);
        bool valid = true;
        if (size > 0)
        {
            double gram[3][3] = {};
            for (int row = 0; row < size; row++)
                for (int column = 0; column < size; column++)
                    gram[row][column] = constraints[selected[row]].normal * constraints[selected[column]].normal;
            const double det = determinant(gram, size);
            if (!std::isfinite(det) || std::fabs(det) < 1e-12)
                continue;
            for (int column = 0; column < size; column++)
            {
                double replaced[3][3];
                for (int row = 0; row < size; row++)
                    for (int j = 0; j < size; j++)
                        replaced[row][j] = j == column ? constraints[selected[row]].distance : gram[row][j];
                const double lambda = determinant(replaced, size) / det;
                if (!std::isfinite(lambda) || lambda < 0.0)
                {
                    valid = false;
                    break;
                }
                delta += lambda * constraints[selected[column]].normal;
            }
        }
        if (!valid)
            continue;
        for (const Constraint& constraint : constraints)
            if (constraint.normal * delta < constraint.distance - 1e-9)
                valid = false;
        const double squared = delta * delta;
        if (valid && std::isfinite(squared) && squared < bestSquared)
        {
            best = delta;
            bestSquared = squared;
        }
    }
    if (std::isfinite(bestSquared))
        return best;
    const Constraint& largest = *std::max_element(constraints.begin(), constraints.end(),
        [](const Constraint& a, const Constraint& b) { return a.distance < b.distance; });
    return largest.distance * largest.normal;
}
}

MStatus SkirtCollideDeformer::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnGenericAttribute gAttr;
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

    const MObject affects[] = {attr_thighRadiusX, attr_thighRadiusZ, attr_kneeRadiusX, attr_kneeRadiusZ,
        attr_calfRadiusX, attr_calfRadiusZ, attr_ankleRadiusX, attr_ankleRadiusZ, attr_thighPosition, attr_calfPosition,
        attr_bellMatrix, attr_leftHipMatrix, attr_leftKneeMatrix, attr_leftHeelMatrix, attr_rightHipMatrix,
        attr_rightKneeMatrix, attr_rightHeelMatrix, attr_skirtType, attr_ringScale, attr_leftRingAxis,
        attr_rightRingAxis, attr_falloff, attr_restGeometry, attr_restBellMatrix,
        attr_restLeftHipMatrix, attr_restLeftKneeMatrix, attr_restLeftHeelMatrix, attr_restRightHipMatrix,
        attr_restRightKneeMatrix, attr_restRightHeelMatrix};
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
        if (!restGeometryWarningIssued)
        {
            MGlobal::displayWarning(MFnDependencyNode(thisMObject()).name() +
                ": restGeometry is missing or its points do not match the input geometry.");
            restGeometryWarningIssued = true;
        }
        return MS::kSuccess;
    }

    const MMatrix bellMatrix = dataBlock.inputValue(attr_bellMatrix, &stat).asMatrix();
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
    const float envelopeValue = dataBlock.inputValue(envelope, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    if (envelopeValue == 0.0f)
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
    const MObject restAttributes[2][3] = {{attr_restLeftHipMatrix, attr_restLeftKneeMatrix, attr_restLeftHeelMatrix},
        {attr_restRightHipMatrix, attr_restRightKneeMatrix, attr_restRightHeelMatrix}};
    const short axes[] = {leftRingAxis, rightRingAxis};
    const double thighPosition = dataBlock.inputValue(attr_thighPosition).asDouble();
    const double calfPosition = dataBlock.inputValue(attr_calfPosition).asDouble();
    const double radii[] = {dataBlock.inputValue(attr_thighRadiusX).asDouble(),
        dataBlock.inputValue(attr_thighRadiusZ).asDouble(), dataBlock.inputValue(attr_kneeRadiusX).asDouble(),
        dataBlock.inputValue(attr_kneeRadiusZ).asDouble(), dataBlock.inputValue(attr_calfRadiusX).asDouble(),
        dataBlock.inputValue(attr_calfRadiusZ).asDouble(), dataBlock.inputValue(attr_ankleRadiusX).asDouble(),
        dataBlock.inputValue(attr_ankleRadiusZ).asDouble()};
    std::vector<CylinderPair> cylinders;
    for (int side = 0; side < 2; side++)
    {
        MMatrix current[3], rest[3];
        for (int joint = 0; joint < 3; joint++)
        {
            current[joint] = dataBlock.inputValue(currentAttributes[side][joint], &stat).asMatrix();
            CHECK_MSTATUS_AND_RETURN_IT(stat);
            rest[joint] = dataBlock.inputValue(restAttributes[side][joint], &stat).asMatrix();
            CHECK_MSTATUS_AND_RETURN_IT(stat);
        }
        bool finitePositions = true;
        for (int joint = 0; joint < 3; joint++)
            finitePositions =
                finitePositions && isFinitePoint(taxis(current[joint])) && isFinitePoint(taxis(rest[joint]));
        if (!finitePositions)
            continue;
        const auto makeProfile = [&](const MMatrix matrices[3])
        {
            return SkirtLegProfile(radii[0], radii[1], radii[2], radii[3], radii[4], radii[5], radii[6], radii[7],
                thighPosition, calfPosition, (taxis(matrices[1]) - taxis(matrices[0])).length(),
                (taxis(matrices[2]) - taxis(matrices[1])).length());
        };
        const SkirtLegProfile currentProfile = makeProfile(current);
        const SkirtLegProfile restProfile = makeProfile(rest);
        for (int segment = 0; segment < (skirtType == 1 ? 2 : 1); segment++)
        {
            CylinderPair pair;
            if (makeCylinder(current[segment], current[segment + 1], axes[side], ringScale, currentProfile,
                    segment == 1, pair.current) &&
                makeCylinder(
                    rest[segment], rest[segment + 1], axes[side], ringScale, restProfile, segment == 1, pair.rest))
                cylinders.push_back(pair);
        }
    }
    const MMatrix restBellMatrix = dataBlock.inputValue(attr_restBellMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MPoint restWaist = taxis(restBellMatrix);
    const MVector restAxis = getAxis(restBellMatrix, 1).normal();
    if (!isFinitePoint(restWaist) || restAxis.length() < 1e-8)
        return MS::kSuccess;
    const MVector rawUp = -getAxis(bellMatrix, 1);
    const MVector up = rawUp.length() >= 1e-8 ? rawUp.normal() : MVector(0.0, 1.0, 0.0);

    for (; !iter.isDone(); iter.next())
    {
        const double strength =
            static_cast<double>(envelopeValue) * weightValue(dataBlock, multiIndex, iter.index());
        if (!std::isfinite(strength) || strength == 0.0)
            continue;
        const int index = iter.index();
        if (index < 0 || static_cast<std::size_t>(index) >= restPoints.size() || !restValid[index])
            continue;
        const MPoint& restPoint = restPoints[index];
        const MPoint point = iter.position();
        if (!isFinitePoint(restPoint) || !isFinitePoint(point))
            continue;
        const MVector offset = restPoint - restWaist;
        MVector waistRadial = offset - restAxis * (offset * restAxis);
        if (waistRadial.length() < 1e-8)
            continue;
        waistRadial.normalize();
        std::vector<Constraint> constraints;
        double upward = 0.0;
        for (const CylinderPair& cylinder : cylinders)
        {
            Constraint constraint;
            if (pointConstraint(cylinder, restPoint, waistRadial, point, falloff, constraint))
            {
                constraints.push_back(constraint);
                upward =
                    (std::max)(upward, (std::max)(0.0, constraint.distance) * (std::max)(0.0, up * constraint.normal));
            }
        }
        if (upward > 0.0)
            constraints.push_back({up, upward});
        const MPoint result = point + strength * minimumDisplacement(constraints);
        if (isFinitePoint(result))
        {
            stat = iter.setPosition(result);
            CHECK_MSTATUS_AND_RETURN_IT(stat);
        }
    }
    return MS::kSuccess;
}
