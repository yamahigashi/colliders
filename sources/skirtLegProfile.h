#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

struct SkirtLegProfile {
  enum class Ring { Knee, Heel, Extended };
  struct Radius {
    double x, z;
  };
  struct Node {
    double s;
    Radius radius;
  };

  double thighPosition, calfPosition;
  double thighLength, legLength, knee;
  std::array<Radius, 5> stations;

  SkirtLegProfile(double thighX, double thighZ, double kneeX, double kneeZ,
                  double calfX, double calfZ, double ankleX, double ankleZ,
                  double thighPos, double calfPos, double thigh, double calf)
      : thighPosition(thighPos), calfPosition(calfPos), thighLength(thigh),
        legLength(thigh + calf),
        knee(legLength < 1e-6 ? 0.0 : thigh / legLength),
        stations{{{1.0, 1.0},
                  {thighX, thighZ},
                  {kneeX, kneeZ},
                  {calfX, calfZ},
                  {ankleX, ankleZ}}} {
    if (legLength < 1e-6)
      return;
    const std::array<Node, 5> candidates{
        {{0.0, stations[0]},
         {knee * thighPosition, stations[1]},
         {knee, stations[2]},
         {knee + (1.0 - knee) * calfPosition, stations[3]},
         {1.0, stations[4]}}};
    // Resolve tolerance chains by priority before sorting the surviving nodes.
    const int precedence[] = {4, 2, 1, 3, 0};
    for (int index : precedence) {
      const Node &candidate = candidates[index];
      bool coincident = false;
      for (std::size_t i = 0; i < count; ++i)
        coincident = coincident || std::abs(candidate.s - nodes[i].s) <= 1e-9;
      if (!coincident)
        nodes[count++] = candidate;
    }
    std::sort(nodes.begin(), nodes.begin() + count,
              [](const Node &a, const Node &b) { return a.s < b.s; });
  }

  double effectiveLength(Ring ring, double scaleY) const {
    return (ring == Ring::Knee ? thighLength : legLength) * scaleY;
  }

  double parameter(double y, Ring ring, double scaleY) const {
    if (legLength < 1e-6 || effectiveLength(ring, scaleY) < 1e-6)
      return 0.0;
    if (ring == Ring::Knee)
      return y * scaleY * thighLength / legLength;
    return ring == Ring::Extended ? (std::min)(y * scaleY, knee) : y * scaleY;
  }

  double levelParameter(double distance, double hipDistance) const {
    return legLength < 1e-6 ? 0.0 : (distance - hipDistance) / legLength;
  }

  Radius evaluate(double s, double length) const {
    if (legLength < 1e-6 || length < 1e-6)
      return {1.0, 1.0};
    if (s < 0.0)
      return {1.0, 1.0};
    if (s > 1.0)
      return stations[4];
    if (s <= nodes[0].s)
      return nodes[0].radius;
    for (std::size_t i = 1; i < count; ++i) {
      if (s <= nodes[i].s) {
        const Node &a = nodes[i - 1], &b = nodes[i];
        const double t = (s - a.s) / (b.s - a.s);
        return {a.radius.x + (b.radius.x - a.radius.x) * t,
                a.radius.z + (b.radius.z - a.radius.z) * t};
      }
    }
    return nodes[count - 1].radius;
  }

  double fX(double s, double length) const { return evaluate(s, length).x; }
  double fZ(double s, double length) const { return evaluate(s, length).z; }

  Radius forRing(double s, Ring ring, double scaleY) const {
    return evaluate(ring == Ring::Extended ? (std::min)(s, knee) : s,
                    effectiveLength(ring, scaleY));
  }

private:
  std::array<Node, 5> nodes{};
  std::size_t count = 0;
};
