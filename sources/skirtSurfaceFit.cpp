#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MPointArray.h>

#include <cmath>
#include <vector>

#include "skirtSurfaceFit.h"
#include "nurbsRefitKernel.h"
#include "pluginIdentity.h"

MTypeId SkirtSurfaceFit::typeId(PluginIdentity::kSkirtSurfaceFitTypeId);
MObject SkirtSurfaceFit::attr_inputSurface;
MObject SkirtSurfaceFit::attr_spansV;
MObject SkirtSurfaceFit::attr_outputSurface;

namespace
{
MStatus invalidInput(const char* reason)
{
    MGlobal::displayError(MString("yddSkirtSurfaceFit: ") + reason);
    return MS::kInvalidParameter;
}
}

MStatus SkirtSurfaceFit::initialize()
{
    MStatus stat;
    MFnTypedAttribute inputAttr;
    attr_inputSurface = inputAttr.create("inputSurface", "is", MFnData::kNurbsSurface,
        MObject::kNullObj, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    inputAttr.setWritable(true);
    inputAttr.setStorable(false);
    inputAttr.setConnectable(true);
    inputAttr.setDisconnectBehavior(MFnAttribute::kReset);
    stat = addAttribute(attr_inputSurface);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    MFnNumericAttribute nAttr;
    attr_spansV = nAttr.create("spansV", "sv", MFnNumericData::kInt, 4, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    nAttr.setKeyable(false);
    nAttr.setStorable(true);
    stat = addAttribute(attr_spansV);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    MFnTypedAttribute outputAttr;
    attr_outputSurface = outputAttr.create("outputSurface", "os", MFnData::kNurbsSurface,
        MObject::kNullObj, &stat);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    outputAttr.setWritable(false);
    outputAttr.setStorable(false);
    stat = addAttribute(attr_outputSurface);
    CHECK_MSTATUS_AND_RETURN_IT(stat);

    stat = attributeAffects(attr_inputSurface, attr_outputSurface);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    return attributeAffects(attr_spansV, attr_outputSurface);
}

MStatus SkirtSurfaceFit::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    if (plug != attr_outputSurface)
        return MS::kUnknownParameter;

    MStatus stat;
    const int spansV = dataBlock.inputValue(attr_spansV, &stat).asInt();
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    if (spansV < 1 || spansV > 256)
        return invalidInput("spansV must be between 1 and 256.");

    const MObject input = dataBlock.inputValue(attr_inputSurface, &stat).asNurbsSurface();
    if (stat != MS::kSuccess || input.isNull())
        return invalidInput("inputSurface must be a readable NURBS surface.");
    MFnNurbsSurface inputFn(input, &stat);
    if (stat != MS::kSuccess)
        return invalidInput("inputSurface must be a readable NURBS surface.");
    if (inputFn.formInV() != MFnNurbsSurface::kOpen)
        return invalidInput("inputSurface V form must be open.");

    MPointArray cvs;
    stat = inputFn.getCVs(cvs, MSpace::kObject);
    if (stat != MS::kSuccess)
        return invalidInput("inputSurface CVs must be readable.");
    for (unsigned int i = 0; i < cvs.length(); ++i)
        if (!(std::fabs(cvs[i].w - 1.0) <= 1e-9))
            return invalidInput("inputSurface CV weights must be 1 within 1e-9.");

    const int nu = inputFn.numCVsInU();
    const int nv = inputFn.numCVsInV();
    const int pu = inputFn.degreeU();
    const int q = inputFn.degreeV();
    if (pu < 1 || q < 1 || nu < pu + 1 || nv < q + 1)
        return invalidInput("inputSurface CV counts must be at least degree + 1 in U and V.");

    MDoubleArray uKnots, vKnots;
    stat = inputFn.getKnotsInU(uKnots);
    if (stat != MS::kSuccess)
        return invalidInput("inputSurface U knots must be readable.");
    stat = inputFn.getKnotsInV(vKnots);
    if (stat != MS::kSuccess)
        return invalidInput("inputSurface V knots must be readable.");
    if (uKnots.length() != static_cast<unsigned int>(nu + pu - 1) ||
        vKnots.length() != static_cast<unsigned int>(nv + q - 1))
        return invalidInput("inputSurface knot counts must equal CV count + degree - 1.");

    for (unsigned int i = 0; i < cvs.length(); ++i)
        if (!std::isfinite(cvs[i].x) || !std::isfinite(cvs[i].y) ||
            !std::isfinite(cvs[i].z) || !std::isfinite(cvs[i].w))
            return invalidInput("inputSurface CVs and knots must be finite.");
    for (unsigned int i = 0; i < uKnots.length(); ++i)
        if (!std::isfinite(uKnots[i]))
            return invalidInput("inputSurface CVs and knots must be finite.");
    for (unsigned int i = 0; i < vKnots.length(); ++i)
        if (!std::isfinite(vKnots[i]))
            return invalidInput("inputSurface CVs and knots must be finite.");

    for (unsigned int i = 1; i < uKnots.length(); ++i)
        if (uKnots[i] < uKnots[i - 1])
            return invalidInput("inputSurface U knots must be nondecreasing.");
    for (unsigned int i = 1; i < vKnots.length(); ++i)
        if (vKnots[i] < vKnots[i - 1])
            return invalidInput("inputSurface V knots must be nondecreasing.");

    std::vector<double> inputVKnots(vKnots.length());
    for (unsigned int i = 0; i < vKnots.length(); ++i)
        inputVKnots[i] = vKnots[i];
    const int m = spansV + 3;
    std::vector<double> samples;
    if (!NurbsRefit::sampleMatrix(nv, q, inputVKnots, m, samples))
        return invalidInput("inputSurface V parameter domain length must be greater than 1e-12.");

    std::vector<double> points(nu * nv * 3);
    for (int i = 0; i < nu * nv; ++i)
    {
        points[i * 3] = cvs[i].x;
        points[i * 3 + 1] = cvs[i].y;
        points[i * 3 + 2] = cvs[i].z;
    }
    std::vector<double> outputPoints;
    NurbsRefit::refitColumns(nu, nv, points, m, samples, outputPoints);
    for (double value : outputPoints)
        if (!std::isfinite(value))
            return invalidInput("outputSurface CVs must be finite.");

    MPointArray outputCVs;
    outputCVs.setLength(nu * m);
    for (int i = 0; i < nu * m; ++i)
        outputCVs.set(MPoint(outputPoints[i * 3], outputPoints[i * 3 + 1],
            outputPoints[i * 3 + 2], 1.0), i);
    const std::vector<double> outputVKnots = NurbsRefit::clampedUniformKnotsMaya(spansV, 3);
    MDoubleArray outputKnots;
    for (double value : outputVKnots)
        outputKnots.append(value);

    MFnNurbsSurfaceData surfaceDataFn;
    MObject surfaceData = surfaceDataFn.create(&stat);
    if (stat != MS::kSuccess)
        return stat;
    MFnNurbsSurface surfaceFn;
    surfaceFn.create(outputCVs, uKnots, outputKnots, pu, 3,
        inputFn.formInU(), MFnNurbsSurface::kOpen, false, surfaceData, &stat);
    if (stat != MS::kSuccess)
        return stat;

    MDataHandle outputHandle = dataBlock.outputValue(attr_outputSurface, &stat);
    if (stat != MS::kSuccess)
        return stat;
    outputHandle.setMObject(surfaceData);
    return dataBlock.setClean(plug);
}
