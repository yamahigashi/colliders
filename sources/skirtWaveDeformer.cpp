#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFloatArray.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MRampAttribute.h>
#include <maya/MVector.h>

#include <cmath>
#include <cstdint>
#include <vector>

#include "skirtWaveDeformer.h"
#include "pluginIdentity.h"
#include "utils.hpp"

MTypeId SkirtWaveDeformer::typeId(PluginIdentity::kSkirtWaveTypeId);

MObject SkirtWaveDeformer::attr_bellMatrix;
MObject SkirtWaveDeformer::attr_evaluationToWorldRotation;
MObject SkirtWaveDeformer::attr_amplitude;
MObject SkirtWaveDeformer::attr_amplitudeRamp;
MObject SkirtWaveDeformer::attr_idleAmplitude;
MObject SkirtWaveDeformer::attr_idleAmplitudeV;
MObject SkirtWaveDeformer::attr_idleAmplitudeU;
MObject SkirtWaveDeformer::attr_idleDirectionality;
MObject SkirtWaveDeformer::attr_idleDirectionX;
MObject SkirtWaveDeformer::attr_idleDirectionZ;
MObject SkirtWaveDeformer::attr_idleDirectionSpace;
MObject SkirtWaveDeformer::attr_phaseSpread;
MObject SkirtWaveDeformer::attr_impulseAmount;
MObject SkirtWaveDeformer::attr_wavePhaseV;
MObject SkirtWaveDeformer::attr_wavePhaseU;
MObject SkirtWaveDeformer::attr_idleComplexity;
MObject SkirtWaveDeformer::attr_waveCountV;
MObject SkirtWaveDeformer::attr_waveCountU;
MObject SkirtWaveDeformer::attr_impulseX;
MObject SkirtWaveDeformer::attr_impulseZ;
MObject SkirtWaveDeformer::attr_impulsePosition;
MObject SkirtWaveDeformer::attr_impulseWidth;
MObject SkirtWaveDeformer::attr_impulseSpace;
MObject SkirtWaveDeformer::attr_directionality;
MObject SkirtWaveDeformer::attr_noiseAmplitude;
MObject SkirtWaveDeformer::attr_noisePhase;
MObject SkirtWaveDeformer::attr_noiseFrequencyV;
MObject SkirtWaveDeformer::attr_noiseFrequencyU;
MObject SkirtWaveDeformer::attr_skew;
MObject SkirtWaveDeformer::attr_sharpness;

namespace
{
const double kMatrixTolerance = 1e-8;
// evaluationToWorldRotation must be a finite pure rotation: orthonormal 3x3
// block with determinant +1, no translation, no perspective.
const double kRotationOrthonormalTolerance = 1e-6;
const double kRotationAffineTolerance = 1e-8;
const double kRadialTolerance = 1e-8;
const double kImpulseTolerance = 1e-5;
const double kMinimumHemHeight = 1e-5;
const double kPi = 3.14159265358979323846;

struct WavePoint
{
    MPoint evaluationPoint;
    MPoint bellLocalPoint;
    float weight;
    bool finite;
};

bool isFiniteMatrix(const MMatrix& matrix)
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

bool isPureRotation(const MMatrix& matrix)
{
    if (!isFiniteMatrix(matrix))
        return false;
    for (unsigned int index = 0; index < 3; index++)
    {
        if (std::fabs(matrix[index][3]) > kRotationAffineTolerance)
            return false;
        if (std::fabs(matrix[3][index]) > kRotationAffineTolerance)
            return false;
    }
    if (std::fabs(matrix[3][3] - 1.0) > kRotationAffineTolerance)
        return false;
    for (unsigned int row = 0; row < 3; row++)
    {
        for (unsigned int column = 0; column < 3; column++)
        {
            double dot = 0.0;
            for (unsigned int k = 0; k < 3; k++)
                dot += matrix[row][k] * matrix[column][k];
            const double expected = row == column ? 1.0 : 0.0;
            if (std::fabs(dot - expected) > kRotationOrthonormalTolerance)
                return false;
        }
    }
    const double determinant =
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return std::fabs(determinant - 1.0) <= kRotationOrthonormalTolerance;
}

// Brings World direction inputs into the evaluation space. The inverse lives
// in this frame only; deform() keeps the two direction vectors as before.
void rotateWorldDirections(const MMatrix& evaluationToWorldRotation,
    bool rotateImpulse, MVector& impulseVector, bool rotateIdle, MVector& idleDirection)
{
    const MMatrix worldToEvaluationRotation = evaluationToWorldRotation.inverse();
    if (rotateImpulse)
        impulseVector = impulseVector * worldToEvaluationRotation;
    if (rotateIdle)
        idleDirection = idleDirection * worldToEvaluationRotation;
}

bool isFinitePoint(const MPoint& point)
{
    return std::isfinite(point.x) && std::isfinite(point.y)
        && std::isfinite(point.z) && std::isfinite(point.w);
}

bool isFiniteVector(const MVector& vector)
{
    return std::isfinite(vector.x) && std::isfinite(vector.y)
        && std::isfinite(vector.z);
}

double shapedSine(double phase, double skew, double sharpness)
{
    const double distortedPhase = phase + skew * std::sin(phase);
    const double sine = std::sin(distortedPhase);
    if (sine > 0.0)
        return std::pow(sine, 1.0 + 3.0 * sharpness);
    if (sine < 0.0)
        return -std::pow(-sine, 1.0 + 3.0 * sharpness);
    return 0.0;
}

double impulseKernel(double v, double position, double width)
{
    const double offset = v - position;
    if (offset <= -width || offset >= width)
        return 0.0;
    const double lobe = std::cos(kPi * offset / (2.0 * width));
    return lobe * lobe;
}

double latticeValue(std::uint32_t x, std::uint32_t y, std::uint32_t z)
{
    std::uint32_t h = x * 73856093u
        ^ y * 19349663u
        ^ z * 83492791u;
    h = (h ^ (h >> 13)) * 1274126177u;
    h = h ^ (h >> 16);
    return ((double)(h & 0xffffffu) / (double)0xffffffu) * 2.0 - 1.0;
}

std::uint32_t latticeCoordinate(double cell)
{
    constexpr double kUint32Period = 4294967296.0;
    double wrapped = std::fmod(cell, kUint32Period);
    if (wrapped < 0.0)
        wrapped += kUint32Period;
    return static_cast<std::uint32_t>(wrapped);
}

double fadeCurve(double t)
{
    return t * t * (3.0 - 2.0 * t);
}

// Deterministic hash-lattice value noise in [-1, 1]; theta is fed in as a
// cos/sin circle embedding by the caller, so the U seam never shows.
double valueNoise(double x, double y, double z)
{
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        return 0.0;

    const double fx = std::floor(x);
    const double fy = std::floor(y);
    const double fz = std::floor(z);
    const std::uint32_t x0 = latticeCoordinate(fx);
    const std::uint32_t y0 = latticeCoordinate(fy);
    const std::uint32_t z0 = latticeCoordinate(fz);
    const double tx = fadeCurve(x - fx);
    const double ty = fadeCurve(y - fy);
    const double tz = fadeCurve(z - fz);

    double result = 0.0;
    for (int dz = 0; dz <= 1; dz++)
    {
        for (int dy = 0; dy <= 1; dy++)
        {
            for (int dx = 0; dx <= 1; dx++)
            {
                const double wx = dx ? tx : 1.0 - tx;
                const double wy = dy ? ty : 1.0 - ty;
                const double wz = dz ? tz : 1.0 - tz;
                result += wx * wy * wz
                    * latticeValue(x0 + dx, y0 + dy, z0 + dz);
            }
        }
    }
    return result;
}
}

SkirtWaveDeformer::SkirtWaveDeformer()
    : MPxDeformerNode(), bellMatrixWarningIssued(false), rotationWarningIssued(false)
{
}

void SkirtWaveDeformer::postConstructor()
{
    MStatus stat;
    MRampAttribute ramp(thisMObject(), attr_amplitudeRamp, &stat);
    if (!stat)
    {
        MGlobal::displayError(MString(PluginIdentity::kSkirtWaveNodeName) + ": failed to initialize amplitudeRamp.");
        return;
    }

    MFloatArray positions;
    MFloatArray values;
    MIntArray interpolations;

    stat = positions.append(0.0f);
    if (!stat)
        return;
    stat = values.append(0.0f);
    if (!stat)
        return;
    stat = interpolations.append(MRampAttribute::kLinear);
    if (!stat)
        return;

    stat = positions.append(1.0f);
    if (!stat)
        return;
    stat = values.append(1.0f);
    if (!stat)
        return;
    stat = interpolations.append(MRampAttribute::kLinear);
    if (!stat)
        return;

    ramp.addEntries(positions, values, interpolations, &stat);
    if (!stat)
        MGlobal::displayError(MString(PluginIdentity::kSkirtWaveNodeName) + ": failed to add amplitudeRamp defaults.");
}

MStatus SkirtWaveDeformer::initialize()
{
    MFnEnumAttribute eAttr;
    MFnMatrixAttribute mAttr;
    MFnNumericAttribute nAttr;
    MStatus stat;

    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix", MFnMatrixAttribute::kDouble, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    mAttr.setHidden(true);
    stat = addAttribute(attr_bellMatrix);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    // Rotation from the evaluation space of the geometry and bellMatrix to the
    // world the World direction modes refer to. Identity keeps existing scenes
    // bitwise unchanged; rig construction connects it, animators never see it.
    attr_evaluationToWorldRotation = mAttr.create("evaluationToWorldRotation", "etwr",
        MFnMatrixAttribute::kDouble, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    mAttr.setHidden(true);
    mAttr.setKeyable(false);
    mAttr.setStorable(true);
    mAttr.setConnectable(true);
    stat = addAttribute(attr_evaluationToWorldRotation);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_amplitude = nAttr.create("amplitude", "amplitude", MFnNumericData::kDouble, 1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setMax(3.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_amplitude);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_amplitudeRamp = MRampAttribute::createCurveRamp("amplitudeRamp", "amplitudeRamp", &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    stat = addAttribute(attr_amplitudeRamp);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleAmplitude = nAttr.create("idleAmplitude", "idleAmplitude", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleAmplitude);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleAmplitudeV = nAttr.create("idleAmplitudeV", "idleAmplitudeV", MFnNumericData::kDouble, 1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(2.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleAmplitudeV);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleAmplitudeU = nAttr.create("idleAmplitudeU", "idleAmplitudeU", MFnNumericData::kDouble, 1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(2.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleAmplitudeU);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleDirectionality = nAttr.create("idleDirectionality", "idleDirectionality", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleDirectionality);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleDirectionX = nAttr.create("idleDirectionX", "idleDirectionX", MFnNumericData::kDouble, -1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-1.0);
    nAttr.setSoftMax(1.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleDirectionX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleDirectionZ = nAttr.create("idleDirectionZ", "idleDirectionZ", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-1.0);
    nAttr.setSoftMax(1.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_idleDirectionZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_phaseSpread = nAttr.create("phaseSpread", "phaseSpread", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(0.5);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_phaseSpread);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulseAmount = nAttr.create("impulseAmount", "impulseAmount", MFnNumericData::kDouble, 1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(2.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_impulseAmount);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_idleDirectionSpace = eAttr.create("idleDirectionSpace", "idleDirectionSpace", 0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    eAttr.addField("World", 0);
    eAttr.addField("Bell Local", 1);
    eAttr.setKeyable(false);
    eAttr.setChannelBox(true);
    stat = addAttribute(attr_idleDirectionSpace);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_wavePhaseV = nAttr.create("wavePhaseV", "wavePhaseV", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-5.0);
    nAttr.setSoftMax(5.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_wavePhaseV);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_wavePhaseU = nAttr.create("wavePhaseU", "wavePhaseU", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-5.0);
    nAttr.setSoftMax(5.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_wavePhaseU);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    // Material/stylization attributes are rigger-tuned (ADR-0004): visible in
    // the channel box but excluded from the keyable set so Key All and anim
    // layers do not capture them.
    attr_idleComplexity = nAttr.create("idleComplexity", "idleComplexity", MFnNumericData::kDouble, 0.35, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_idleComplexity);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_waveCountV = nAttr.create("waveCountV", "waveCountV", MFnNumericData::kDouble, 1.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMax(1.5);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_waveCountV);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_waveCountU = nAttr.create("waveCountU", "waveCountU", MFnNumericData::kInt, 2, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_waveCountU);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulseX = nAttr.create("impulseX", "impulseX", MFnNumericData::kDouble, -0.5, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-2.0);
    nAttr.setSoftMax(2.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_impulseX);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulseZ = nAttr.create("impulseZ", "impulseZ", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-2.0);
    nAttr.setSoftMax(2.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_impulseZ);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulsePosition = nAttr.create("impulsePosition", "impulsePosition", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(0.0);
    nAttr.setSoftMax(1.5);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_impulsePosition);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulseWidth = nAttr.create("impulseWidth", "impulseWidth", MFnNumericData::kDouble, 0.25, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.01);
    nAttr.setMax(1.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_impulseWidth);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_impulseSpace = eAttr.create("impulseSpace", "impulseSpace", 0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    eAttr.addField("World", 0);
    eAttr.addField("Bell Local", 1);
    eAttr.setKeyable(false);
    eAttr.setChannelBox(true);
    stat = addAttribute(attr_impulseSpace);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_noiseAmplitude = nAttr.create("noiseAmplitude", "noiseAmplitude", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_noiseAmplitude);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_noisePhase = nAttr.create("noisePhase", "noisePhase", MFnNumericData::kDouble, 0.0, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setSoftMin(-5.0);
    nAttr.setSoftMax(5.0);
    nAttr.setKeyable(true);
    stat = addAttribute(attr_noisePhase);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_noiseFrequencyV = nAttr.create("noiseFrequencyV", "noiseFrequencyV", MFnNumericData::kDouble, 1.5, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(5.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_noiseFrequencyV);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_noiseFrequencyU = nAttr.create("noiseFrequencyU", "noiseFrequencyU", MFnNumericData::kDouble, 1.5, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(5.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_noiseFrequencyU);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_directionality = nAttr.create("directionality", "directionality", MFnNumericData::kDouble, 0.8, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_directionality);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_skew = nAttr.create("skew", "skew", MFnNumericData::kDouble, 0.3, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(-1.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_skew);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    attr_sharpness = nAttr.create("sharpness", "sharpness", MFnNumericData::kDouble, 0.3, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(false);
    nAttr.setChannelBox(true);
    stat = addAttribute(attr_sharpness);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    const MObject affects[] = {
        attr_bellMatrix,
        attr_evaluationToWorldRotation,
        attr_amplitude,
        attr_amplitudeRamp,
        attr_idleAmplitude,
        attr_idleAmplitudeV,
        attr_idleAmplitudeU,
        attr_idleDirectionality,
        attr_idleDirectionX,
        attr_idleDirectionZ,
        attr_idleDirectionSpace,
        attr_phaseSpread,
        attr_impulseAmount,
        attr_wavePhaseV,
        attr_wavePhaseU,
        attr_idleComplexity,
        attr_waveCountV,
        attr_waveCountU,
        attr_impulseX,
        attr_impulseZ,
        attr_impulsePosition,
        attr_impulseWidth,
        attr_impulseSpace,
        attr_directionality,
        attr_noiseAmplitude,
        attr_noisePhase,
        attr_noiseFrequencyV,
        attr_noiseFrequencyU,
        attr_skew,
        attr_sharpness
    };
    for (const MObject& attribute : affects)
    {
        stat = attributeAffects(attribute, outputGeom);
        CHECK_MSTATUS_AND_RETURN_IT(stat);
    }

    return MS::kSuccess;
}

MStatus SkirtWaveDeformer::deform(MDataBlock& dataBlock, MItGeometry& iter,
    const MMatrix&, unsigned int multiIndex)
{
    MStatus stat;
    const MMatrix bellMatrix = dataBlock.inputValue(attr_bellMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const MMatrix evaluationToWorldRotation =
        dataBlock.inputValue(attr_evaluationToWorldRotation, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double amplitude = dataBlock.inputValue(attr_amplitude, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleAmplitude = dataBlock.inputValue(attr_idleAmplitude, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleAmplitudeV = dataBlock.inputValue(attr_idleAmplitudeV, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleAmplitudeU = dataBlock.inputValue(attr_idleAmplitudeU, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleDirectionality = dataBlock.inputValue(attr_idleDirectionality, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleDirectionX = dataBlock.inputValue(attr_idleDirectionX, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleDirectionZ = dataBlock.inputValue(attr_idleDirectionZ, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short idleDirectionSpace = dataBlock.inputValue(attr_idleDirectionSpace, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double phaseSpread = dataBlock.inputValue(attr_phaseSpread, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double impulseAmount = dataBlock.inputValue(attr_impulseAmount, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double wavePhaseV = dataBlock.inputValue(attr_wavePhaseV, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double wavePhaseU = dataBlock.inputValue(attr_wavePhaseU, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleComplexity = dataBlock.inputValue(attr_idleComplexity, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double waveCountV = dataBlock.inputValue(attr_waveCountV, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const int waveCountU = dataBlock.inputValue(attr_waveCountU, &stat).asInt();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double impulseX = dataBlock.inputValue(attr_impulseX, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double impulseZ = dataBlock.inputValue(attr_impulseZ, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double impulsePosition = dataBlock.inputValue(attr_impulsePosition, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double impulseWidth = dataBlock.inputValue(attr_impulseWidth, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const short impulseSpace = dataBlock.inputValue(attr_impulseSpace, &stat).asShort();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double directionality = dataBlock.inputValue(attr_directionality, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double noiseAmplitude = dataBlock.inputValue(attr_noiseAmplitude, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double noisePhase = dataBlock.inputValue(attr_noisePhase, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double noiseFrequencyV = dataBlock.inputValue(attr_noiseFrequencyV, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double noiseFrequencyU = dataBlock.inputValue(attr_noiseFrequencyU, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double skew = dataBlock.inputValue(attr_skew, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double sharpness = dataBlock.inputValue(attr_sharpness, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const float envelopeValue = dataBlock.inputValue(envelope, &stat).asFloat();
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    const bool finiteBellMatrix = isFiniteMatrix(bellMatrix);
    const double determinant = finiteBellMatrix ? bellMatrix.det4x4() : 0.0;
    const MVector rawBellAxis = finiteBellMatrix
        ? maxis(bellMatrix, 1) : MVector(0.0, 0.0, 0.0);
    const double globalScale = rawBellAxis.length();
    const bool invalidBellMatrix = !finiteBellMatrix || !std::isfinite(determinant)
        || std::fabs(determinant) < kMatrixTolerance
        || !std::isfinite(globalScale) || globalScale < kMatrixTolerance;
    const MMatrix bellInverse = invalidBellMatrix ? MMatrix() : bellMatrix.inverse();
    if (invalidBellMatrix || !isFiniteMatrix(bellInverse))
    {
        if (!bellMatrixWarningIssued)
        {
            MFnDependencyNode nodeFn(thisMObject());
            MGlobal::displayWarning(MString(PluginIdentity::kSkirtWaveNodeName) + " " + nodeFn.name()
                + ": degenerate bellMatrix; passing geometry through unchanged.");
            bellMatrixWarningIssued = true;
        }
        return MS::kSuccess;
    }

    if (amplitude == 0.0 || envelopeValue == 0.0f)
        return MS::kSuccess;

    const MPoint bellEvaluationPosition = taxis(bellMatrix);
    const MVector bellAxis = rawBellAxis / globalScale;

    // World (default): X/Z are global axes, projected off the cone axis so the
    // keyed direction matches viewport intuition. Bell Local: the waist frame.
    // When the geometry is evaluated in a space whose orientation differs from
    // the world, World directions are first brought into that space with the
    // inverse of evaluationToWorldRotation. An exact identity skips the
    // multiplication so existing scenes keep their operation order.
    const bool impulseUseWorld = (impulseSpace == 0);
    const bool idleUseWorld = (idleDirectionSpace == 0);
    const bool rotateDirections = (impulseUseWorld || idleUseWorld)
        && !(evaluationToWorldRotation == MMatrix::identity);
    if (rotateDirections && !isPureRotation(evaluationToWorldRotation))
    {
        if (!rotationWarningIssued)
        {
            MFnDependencyNode nodeFn(thisMObject());
            MGlobal::displayWarning(MString(PluginIdentity::kSkirtWaveNodeName) + " " + nodeFn.name()
                + ": evaluationToWorldRotation is not a finite pure rotation; passing geometry through unchanged.");
            rotationWarningIssued = true;
        }
        return MS::kSuccess;
    }

    MVector impulseVector(0.0, 0.0, 0.0);
    if (impulseAmount != 0.0)
        impulseVector = MVector(impulseX, 0.0, impulseZ);
    MVector idleDirection(idleDirectionX, 0.0, idleDirectionZ);
    const double idleDirectionLength = std::hypot(idleDirectionX, idleDirectionZ);
    if (idleDirectionLength < kImpulseTolerance)
        idleDirection = MVector(0.0, 0.0, 0.0);
    else
        idleDirection /= idleDirectionLength;
    if (rotateDirections)
        rotateWorldDirections(evaluationToWorldRotation, impulseUseWorld, impulseVector,
            idleUseWorld, idleDirection);

    if (impulseUseWorld)
        impulseVector -= bellAxis * (impulseVector * bellAxis);
    const double impulseLength = impulseVector.length();

    if (idleAmplitude == 0.0 && impulseLength < kImpulseTolerance
        && noiseAmplitude < kImpulseTolerance)
        return MS::kSuccess;

    if (idleUseWorld && idleDirectionLength >= kImpulseTolerance)
        idleDirection -= bellAxis * (idleDirection * bellAxis);

    std::vector<WavePoint> points;
    double maximumHeight = -kMinimumHemHeight;
    for (; !iter.isDone(); iter.next())
    {
        WavePoint point;
        point.evaluationPoint = iter.position();
        point.bellLocalPoint = point.evaluationPoint * bellInverse;
        point.weight = weightValue(dataBlock, multiIndex, iter.index());
        point.finite = isFinitePoint(point.evaluationPoint) && isFinitePoint(point.bellLocalPoint);

        if (point.finite && point.bellLocalPoint.y > maximumHeight)
            maximumHeight = point.bellLocalPoint.y;
        points.push_back(point);
    }
    const double hemHeight = maximumHeight < kMinimumHemHeight
        ? kMinimumHemHeight : maximumHeight;

    MStatus rampStatus;
    MRampAttribute amplitudeRamp(thisMObject(), attr_amplitudeRamp, &rampStatus);
    CHECK_MSTATUS_AND_RETURN_IT(rampStatus);

    const double goldenRatio = (1.0 + std::sqrt(5.0)) * 0.5;

    MPointArray outputPoints;
    stat = outputPoints.setLength(static_cast<unsigned int>(points.size()));
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    for (unsigned int index = 0; index < points.size(); index++)
    {
        const WavePoint& point = points[index];
        MPoint outputPoint = point.evaluationPoint;
        if (!point.finite || point.weight == 0.0f)
        {
            outputPoints.set(outputPoint, index);
            continue;
        }

        const double v = clamp(point.bellLocalPoint.y / hemHeight, 0.0, 1.0);
        const double theta = std::atan2(point.bellLocalPoint.z, point.bellLocalPoint.x);

        MVector radialLocal(point.bellLocalPoint.x, 0.0, point.bellLocalPoint.z);
        if (radialLocal.length() < kRadialTolerance)
        {
            outputPoints.set(outputPoint, index);
            continue;
        }
        radialLocal.normalize();

        const MVector offset = point.evaluationPoint - bellEvaluationPosition;
        MVector radialEvaluation = offset - bellAxis * (offset * bellAxis);
        if (radialEvaluation.length() < kRadialTolerance)
        {
            outputPoints.set(outputPoint, index);
            continue;
        }
        radialEvaluation.normalize();

        double idleWave = 0.0;
        if (idleAmplitude != 0.0)
        {
            const double directionDot = idleUseWorld
                ? radialEvaluation * idleDirection : radialLocal * idleDirection;
            const double verticalGain = idleAmplitudeV
                * ((1.0 - idleDirectionality) + idleDirectionality * directionDot);
            double localPhaseV = wavePhaseV;
            if (phaseSpread != 0.0 && verticalGain != 0.0)
                localPhaseV += phaseSpread * valueNoise(
                    std::cos(theta) * noiseFrequencyU + 7.31,
                    std::sin(theta) * noiseFrequencyU + 3.17,
                    noisePhase);
            const auto combineWaves = [&](double phaseV, double phaseU)
            {
                const double vertical = verticalGain == 0.0
                    ? 0.0 : verticalGain * shapedSine(phaseV, skew, sharpness);
                if (waveCountU == 0)
                    return vertical;
                const double around = idleAmplitudeU == 0.0
                    ? 0.0 : idleAmplitudeU * shapedSine(phaseU, skew, sharpness);
                return 0.5 * (vertical + around);
            };
            double basePrimary = 0.0;
            if (idleComplexity != 1.0)
            {
                const double phiV = 2.0 * kPi * (localPhaseV - waveCountV * v);
                const double phiU = waveCountU * (theta - 2.0 * kPi * wavePhaseU);
                basePrimary = combineWaves(phiV, phiU);
            }
            double baseSecondary = 0.0;
            if (idleComplexity != 0.0)
            {
                const double phiV2 = 2.0 * kPi
                    * (goldenRatio * localPhaseV - waveCountV * v)
                    + 0.5 * kPi;
                const double phiU2 = waveCountU
                    * (theta - 2.0 * kPi * goldenRatio * wavePhaseU)
                    + 0.5 * kPi;
                baseSecondary = combineWaves(phiV2, phiU2);
            }
            idleWave = idleAmplitude
                * ((1.0 - idleComplexity) * basePrimary
                    + idleComplexity * baseSecondary);
        }

        const double impulseDot = impulseUseWorld
            ? (radialEvaluation * impulseVector)
            : (radialLocal * impulseVector);
        const double impulseWave = impulseAmount == 0.0 || impulseLength < kImpulseTolerance
            ? 0.0
            : impulseAmount * impulseKernel(v, impulsePosition, impulseWidth)
                * ((1.0 - directionality) * impulseLength
                    + directionality * impulseDot);

        const double noiseWave = noiseAmplitude < kImpulseTolerance
            ? 0.0
            : noiseAmplitude * valueNoise(
                std::cos(theta) * noiseFrequencyU + 7.31,
                std::sin(theta) * noiseFrequencyU + 3.17,
                v * noiseFrequencyV + noisePhase);

        // Skip ramp sampling and preserve the input when the wave contributes nothing.
        const double waveSum = idleWave + impulseWave + noiseWave;
        if (waveSum == 0.0)
        {
            outputPoints.set(outputPoint, index);
            continue;
        }

        float rampValue = 0.0f;
        amplitudeRamp.getValueAtPosition(static_cast<float>(v), rampValue, &rampStatus);
        CHECK_MSTATUS_AND_RETURN_IT(rampStatus);

        const double displacementScale = globalScale * (double)envelopeValue
            * (double)point.weight * amplitude * (double)rampValue * waveSum;
        const MVector displacement = radialEvaluation * displacementScale;
        if (displacementScale == 0.0 || !std::isfinite(displacementScale)
            || !isFiniteVector(displacement))
        {
            outputPoints.set(outputPoint, index);
            continue;
        }

        const MPoint displacedEvaluation = point.evaluationPoint + displacement;
        if (isFinitePoint(displacedEvaluation))
            outputPoint = displacedEvaluation;
        outputPoints.set(outputPoint, index);
    }

    stat = iter.setAllPositions(outputPoints, MSpace::kObject);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    return MS::kSuccess;
}
