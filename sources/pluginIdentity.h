#pragma once

// Fixed local, internal-use IDs; see docs/adr/0003-ydd-plugin-identity.md.
namespace PluginIdentity
{
constexpr char kVendor[] = "yamahigashi";
constexpr char kVersion[] = "4.0.0";
constexpr char kDrawRegistrant[] = "yddCollidersPlugin";
constexpr char kBellNodeName[] = "yddBellCollider";
constexpr unsigned int kBellTypeId = 0x0007DD00;
constexpr char kPlaneNodeName[] = "yddPlaneCollider";
constexpr unsigned int kPlaneTypeId = 0x0007DD01;
constexpr char kSkirtBellNodeName[] = "yddSkirtBellCollider";
constexpr unsigned int kSkirtBellTypeId = 0x0007DD02;
constexpr char kSkirtCollideNodeName[] = "yddSkirtCollideDeformer";
constexpr unsigned int kSkirtCollideTypeId = 0x0007DD03;
constexpr char kSkirtWaveNodeName[] = "yddSkirtWaveDeformer";
constexpr unsigned int kSkirtWaveTypeId = 0x0007DD04;
}
