#include <maya/MFnPlugin.h>

#include "pluginIdentity.h"
#include "bellCollider.h"
#include "planeCollider.h"
#include "skirtBellCollider.h"
#include "skirtCollideDeformer.h"
#include "skirtWaveDeformer.h"
#include "skirtSurfaceFit.h"

MStatus initializePlugin(MObject plugin)
{
	MStatus stat;

	MFnPlugin pluginFn(plugin, PluginIdentity::kVendor, PluginIdentity::kVersion, "Any");
	stat = pluginFn.registerNode(PluginIdentity::kBellNodeName, BellCollider::typeId, BellCollider::creator, BellCollider::initialize, MPxNode::kLocatorNode, &BellCollider::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::registerDrawOverrideCreator(BellCollider::drawDbClassification, BellCollider::drawRegistrantId, BellColliderDrawOverride::creator);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode(PluginIdentity::kPlaneNodeName, PlaneCollider::typeId, PlaneCollider::creator, PlaneCollider::initialize, MPxNode::kLocatorNode, &PlaneCollider::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::registerDrawOverrideCreator(PlaneCollider::drawDbClassification, PlaneCollider::drawRegistrantId, PlaneColliderDrawOverride::creator);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode(PluginIdentity::kSkirtBellNodeName, SkirtBellCollider::typeId, SkirtBellCollider::creator, SkirtBellCollider::initialize, MPxNode::kLocatorNode, &SkirtBellCollider::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::registerDrawOverrideCreator(SkirtBellCollider::drawDbClassification, SkirtBellCollider::drawRegistrantId, SkirtBellColliderDrawOverride::creator);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode(PluginIdentity::kSkirtCollideNodeName, SkirtCollideDeformer::typeId, SkirtCollideDeformer::creator, SkirtCollideDeformer::initialize, MPxNode::kDeformerNode);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode(PluginIdentity::kSkirtWaveNodeName, SkirtWaveDeformer::typeId, SkirtWaveDeformer::creator, SkirtWaveDeformer::initialize, MPxNode::kDeformerNode);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode(PluginIdentity::kSkirtSurfaceFitNodeName, SkirtSurfaceFit::typeId, SkirtSurfaceFit::creator, SkirtSurfaceFit::initialize, MPxNode::kDependNode);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return MS::kSuccess;
}

MStatus uninitializePlugin(MObject plugin)
{
	MStatus stat;

	MFnPlugin pluginFn(plugin);

	stat = pluginFn.deregisterNode(SkirtSurfaceFit::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(SkirtWaveDeformer::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(SkirtCollideDeformer::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(SkirtBellCollider::drawDbClassification, SkirtBellCollider::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(SkirtBellCollider::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(PlaneCollider::drawDbClassification, PlaneCollider::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(PlaneCollider::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(BellCollider::drawDbClassification, BellCollider::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(BellCollider::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return MS::kSuccess;
}
