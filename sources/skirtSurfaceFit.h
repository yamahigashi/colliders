#pragma once

#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <maya/MTypeId.h>

class SkirtSurfaceFit : public MPxNode
{
public:
    static MTypeId typeId;
    static MObject attr_inputSurface;
    static MObject attr_spansV;
    static MObject attr_outputSurface;

    static void* creator() { return new SkirtSurfaceFit(); }
    static MStatus initialize();
    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock) override;
};
