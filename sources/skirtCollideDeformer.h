#pragma once

#include <maya/MPxDeformerNode.h>
#include <maya/MMatrix.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <maya/MTypeId.h>

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

    static MObject attr_skirtType;
    static MObject attr_ringScale;
    static MObject attr_leftRingAxis;
    static MObject attr_rightRingAxis;
    static MObject attr_collision;
    static MObject attr_falloff;
    static MObject attr_endFade;

    SkirtCollideDeformer() : MPxDeformerNode() {}
    virtual ~SkirtCollideDeformer() override {}

    static void* creator() { return new SkirtCollideDeformer(); }
    static MStatus initialize();

    virtual MStatus deform(MDataBlock& dataBlock, MItGeometry& iter,
        const MMatrix& localToWorldMatrix, unsigned int multiIndex) override;
};
