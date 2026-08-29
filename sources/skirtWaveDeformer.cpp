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
#include <vector>

#include "skirtWaveDeformer.h"
#include "utils.hpp"

MTypeId SkirtWaveDeformer::typeId(1274438);

MObject SkirtWaveDeformer::attr_bellMatrix;
MObject SkirtWaveDeformer::attr_amplitude;
MObject SkirtWaveDeformer::attr_amplitudeRamp;
MObject SkirtWaveDeformer::attr_idleAmplitude;
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
const double kRadialTolerance = 1e-8;
const double kImpulseTolerance = 1e-5;
const double kMinimumHemHeight = 1e-5;
const double kPi = 3.14159265358979323846;

struct WavePoint
{
    MPoint objectPoint;
    MPoint worldPoint;
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

double latticeValue(int x, int y, int z)
{
    unsigned int h = (unsigned int)x * 73856093u
        ^ (unsigned int)y * 19349663u
        ^ (unsigned int)z * 83492791u;
    h = (h ^ (h >> 13)) * 1274126177u;
    h = h ^ (h >> 16);
    return ((double)(h & 0xffffffu) / (double)0xffffffu) * 2.0 - 1.0;
}

double fadeCurve(double t)
{
    return t * t * (3.0 - 2.0 * t);
}

// Deterministic hash-lattice value noise in [-1, 1]; theta is fed in as a
// cos/sin circle embedding by the caller, so the U seam never shows.
double valueNoise(double x, double y, double z)
{
    const double fx = std::floor(x);
    const double fy = std::floor(y);
    const double fz = std::floor(z);
    const int x0 = (int)fx;
    const int y0 = (int)fy;
    const int z0 = (int)fz;
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
    : MPxDeformerNode(), bellMatrixWarningIssued(false)
{
}

void SkirtWaveDeformer::postConstructor()
{
    MStatus stat;
    MRampAttribute ramp(thisMObject(), attr_amplitudeRamp, &stat);
    if (!stat)
    {
        MGlobal::displayError("skirtWaveDeformer: failed to initialize amplitudeRamp.");
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
        MGlobal::displayError("skirtWaveDeformer: failed to add amplitudeRamp defaults.");
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
        attr_amplitude,
        attr_amplitudeRamp,
        attr_idleAmplitude,
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
    const MMatrix& localToWorldMatrix, unsigned int multiIndex)
{
    MStatus stat;
    const MMatrix bellMatrix = dataBlock.inputValue(attr_bellMatrix, &stat).asMatrix();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double amplitude = dataBlock.inputValue(attr_amplitude, &stat).asDouble();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    const double idleAmplitude = dataBlock.inputValue(attr_idleAmplitude, &stat).asDouble();
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
            MGlobal::displayWarning(MString("skirtWaveDeformer ") + nodeFn.name()
                + ": degenerate bellMatrix; passing geometry through unchanged.");
            bellMatrixWarningIssued = true;
        }
        return MS::kSuccess;
    }

    const MPoint bellPosition = taxis(bellMatrix);
    const MVector bellAxis = rawBellAxis / globalScale;

    std::vector<WavePoint> points;
    double maximumHeight = -kMinimumHemHeight;
    for (; !iter.isDone(); iter.next())
    {
        WavePoint point;
        point.objectPoint = iter.position();
        point.worldPoint = point.objectPoint * localToWorldMatrix;
        point.bellLocalPoint = point.worldPoint * bellInverse;
        point.weight = weightValue(dataBlock, multiIndex, iter.index());
        point.finite = isFinitePoint(point.objectPoint)
            && isFinitePoint(point.worldPoint) && isFinitePoint(point.bellLocalPoint);

        if (point.finite && point.bellLocalPoint.y > maximumHeight)
            maximumHeight = point.bellLocalPoint.y;
        points.push_back(point);
    }
    const double hemHeight = maximumHeight < kMinimumHemHeight
        ? kMinimumHemHeight : maximumHeight;

    MStatus rampStatus;
    MRampAttribute amplitudeRamp(thisMObject(), attr_amplitudeRamp, &rampStatus);
    CHECK_MSTATUS_AND_RETURN_IT(rampStatus);

    const MMatrix worldToLocalMatrix = localToWorldMatrix.inverse();
    const double goldenRatio = (1.0 + std::sqrt(5.0)) * 0.5;

    // World (default): X/Z are global axes, projected off the cone axis so the
    // keyed direction matches viewport intuition. Bell Local: the waist frame.
    const bool impulseUseWorld = (impulseSpace == 0);
    MVector impulseVector(impulseX, 0.0, impulseZ);
    if (impulseUseWorld)
        impulseVector -= bellAxis * (impulseVector * bellAxis);
    const double impulseLength = impulseVector.length();

    MPointArray outputPoints;
    stat = outputPoints.setLength(static_cast<unsigned int>(points.size()));
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    for (unsigned int index = 0; index < points.size(); index++)
    {
        const WavePoint& point = points[index];
        MPoint outputPoint = point.objectPoint;
        if (!point.finite)
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

        const MVector offset = point.worldPoint - bellPosition;
        MVector radialWorld = offset - bellAxis * (offset * bellAxis);
        if (radialWorld.length() < kRadialTolerance)
        {
            outputPoints.set(outputPoint, index);
            continue;
        }
        radialWorld.normalize();

        const double phiV = 2.0 * kPi * (wavePhaseV - waveCountV * v);
        const double phiU = waveCountU * (theta - 2.0 * kPi * wavePhaseU);
        const double phiV2 = 2.0 * kPi
            * (goldenRatio * wavePhaseV - waveCountV * v)
            + 0.5 * kPi;
        const double phiU2 = waveCountU
            * (theta - 2.0 * kPi * goldenRatio * wavePhaseU)
            + 0.5 * kPi;

        const double basePrimary = waveCountU > 0
            ? 0.5 * (shapedSine(phiV, skew, sharpness)
                + shapedSine(phiU, skew, sharpness))
            : shapedSine(phiV, skew, sharpness);
        const double baseSecondary = waveCountU > 0
            ? 0.5 * (shapedSine(phiV2, skew, sharpness)
                + shapedSine(phiU2, skew, sharpness))
            : shapedSine(phiV2, skew, sharpness);
        const double idleWave = idleAmplitude
            * ((1.0 - idleComplexity) * basePrimary
                + idleComplexity * baseSecondary);

        const double impulseDot = impulseUseWorld
            ? (radialWorld * impulseVector)
            : (radialLocal * impulseVector);
        const double impulseWave = impulseLength < kImpulseTolerance
            ? 0.0
            : impulseKernel(v, impulsePosition, impulseWidth)
                * ((1.0 - directionality) * impulseLength
                    + directionality * impulseDot);

        const double noiseWave = noiseAmplitude < kImpulseTolerance
            ? 0.0
            : noiseAmplitude * valueNoise(
                std::cos(theta) * noiseFrequencyU + 7.31,
                std::sin(theta) * noiseFrequencyU + 3.17,
                v * noiseFrequencyV + noisePhase);

        // Exact no-op at rest: skip the ramp sample and the lossy
        // world/object round trip whenever the wave contributes nothing.
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
        const MVector displacement = radialWorld * displacementScale;
        if (displacementScale == 0.0 || !std::isfinite(displacementScale)
            || !isFiniteVector(displacement))
        {
            outputPoints.set(outputPoint, index);
            continue;
        }

        const MPoint displacedWorld = point.worldPoint + displacement;
        const MPoint displacedObject = displacedWorld * worldToLocalMatrix;
        if (isFinitePoint(displacedWorld) && isFinitePoint(displacedObject))
            outputPoint = displacedObject;
        outputPoints.set(outputPoint, index);
    }

    stat = iter.setAllPositions(outputPoints, MSpace::kObject);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    return MS::kSuccess;
}
