#pragma once

#include <maya/MPxDeformerNode.h>
#include <maya/MMatrix.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <maya/MTypeId.h>

class SkirtWaveDeformer : public MPxDeformerNode
{
public:
    static MTypeId typeId;

    static MObject attr_bellMatrix;
    static MObject attr_amplitude;
    static MObject attr_amplitudeRamp;
    static MObject attr_idleAmplitude;
    static MObject attr_wavePhaseV;
    static MObject attr_wavePhaseU;
    static MObject attr_idleComplexity;
    static MObject attr_waveCountV;
    static MObject attr_waveCountU;
    static MObject attr_impulseX;
    static MObject attr_impulseZ;
    static MObject attr_impulsePosition;
    static MObject attr_impulseWidth;
    static MObject attr_impulseSpace;
    static MObject attr_directionality;
    static MObject attr_noiseAmplitude;
    static MObject attr_noisePhase;
    static MObject attr_noiseFrequencyV;
    static MObject attr_noiseFrequencyU;
    static MObject attr_skew;
    static MObject attr_sharpness;

    SkirtWaveDeformer();
    virtual ~SkirtWaveDeformer() override {}

    static void* creator() { return new SkirtWaveDeformer(); }
    static MStatus initialize();

    virtual void postConstructor() override;
    virtual MStatus deform(MDataBlock& dataBlock, MItGeometry& iter,
        const MMatrix& localToWorldMatrix, unsigned int multiIndex) override;

private:
    bool bellMatrixWarningIssued;
};
