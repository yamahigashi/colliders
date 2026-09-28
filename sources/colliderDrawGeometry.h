#pragma once

#include <array>
#include <cmath>
#include <cstring>
#include <maya/MFloatPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MUIDrawManager.h>
#include "colliderInputValidation.h"
#include <vector>

namespace ColliderDraw {
inline bool samePoints(const MPointArray &a, const MPointArray &b) {
  if (a.length() != b.length())
    return false;
  for (unsigned int i = 0; i < a.length(); ++i)
    if (a[i] != b[i])
      return false;
  return true;
}
inline bool sameMatrices(const std::vector<MMatrix> &a,
                         const std::vector<MMatrix> &b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    for (unsigned int r = 0; r < 4; ++r)
      for (unsigned int c = 0; c < 4; ++c)
        if (a[i][r][c] != b[i][r][c])
          return false;
  return true;
}
struct Geometry {
  MFloatPointArray triangles;
  MPointArray lines;
  void draw(MHWRender::MUIDrawManager &manager, const MColor &fill,
            const MColor &wire) const {
    manager.setColor(fill);
    if (triangles.length())
      manager.mesh(MHWRender::MUIDrawManager::kTriangles, triangles);
    manager.setColor(wire);
    if (lines.length())
      manager.lineList(lines, false);
  }
};
inline void appendLines(MPointArray &lines, const MPointArray &points,
                        unsigned int sides) {
  for (unsigned int i = 0; i < sides; ++i) {
    const unsigned int a = i + 1, b = (i + 1) % sides + 1;
    lines.append(points[a]);
    lines.append(points[b]);
    lines.append(points[a + sides]);
    lines.append(points[b + sides]);
    lines.append(points[a]);
    lines.append(points[a + sides]);
  }
}
struct Rings {
  Geometry geometry;
  std::vector<MMatrix> matrices;
  std::vector<std::array<double, 2>> farMultipliers;
  MPointArray unitPoints;
  MIntArray indices;
  int subdivision = 0;
  bool clearInvalidCache() {
    const bool changed = subdivision != 0 || !matrices.empty() ||
                         !farMultipliers.empty() || unitPoints.length() != 0 ||
                         indices.length() != 0 || geometry.triangles.length() != 0 ||
                         geometry.lines.length() != 0;
    geometry.triangles.clear();
    geometry.lines.clear();
    matrices.clear();
    farMultipliers.clear();
    unitPoints.clear();
    indices.clear();
    subdivision = 0;
    return changed;
  }
  bool update(const std::vector<MMatrix> &next, int sides,
              const std::vector<std::array<double, 2>> &nextFar = {}) {
    if (!ColliderInput::validSubdivision(sides) ||
        (!nextFar.empty() && nextFar.size() != next.size()))
      return clearInvalidCache();
    if (subdivision == sides && sameMatrices(matrices, next) &&
        farMultipliers == nextFar)
      return false;
    if (subdivision != sides) {
      subdivision = sides;
      unitPoints.clear();
      indices.clear();
      unitPoints.append(MPoint(0, 0, 0));
      for (int row = 0; row < 2; ++row)
        for (int i = 0; i < sides; ++i) {
          const double angle = i * (6.28318530717958647692 / sides);
          unitPoints.append(MPoint(std::cos(angle), row, std::sin(angle)));
        }
      for (int i = 0; i < sides; ++i) {
        const int a = i + 1, b = (i + 1) % sides + 1;
        const int face[] = {0, a, b, a, a + sides, b + sides, a, b + sides, b};
        for (int index : face)
          indices.append(index);
      }
    }
    matrices = next;
    farMultipliers = nextFar;
    geometry.triangles.clear();
    geometry.lines.clear();
    for (size_t m = 0; m < matrices.size(); ++m) {
      const auto &matrix = matrices[m];
      const std::array<double, 2> endRadius =
          farMultipliers.empty() ? std::array<double, 2>{{1.0, 1.0}}
                                 : farMultipliers.at(m);
      MPointArray points;
      for (unsigned int i = 0; i < unitPoints.length(); ++i) {
        MPoint unit = unitPoints[i];
        if (i > static_cast<unsigned int>(sides)) {
          unit.x *= endRadius[0];
          unit.z *= endRadius[1];
        }
        const MPoint p = unit * matrix;
        points.append(MPoint(static_cast<float>(p.x), static_cast<float>(p.y),
                             static_cast<float>(p.z)));
      }
      for (unsigned int i = 0; i < indices.length(); ++i)
        geometry.triangles.append(MFloatPoint(points[indices[i]]));
      appendLines(geometry.lines, points, sides);
    }
    return true;
  }
};
struct BellMesh {
  Geometry geometry;
  MPointArray points;
  bool update(const MObject &mesh) {
    MPointArray next;
    MStatus status;
    MFnMesh fn(mesh, &status);
    if (status)
      fn.getPoints(next);
    if (samePoints(points, next))
      return false;
    points = next;
    geometry.triangles.clear();
    geometry.lines.clear();
    if (!points.length())
      return true;
    MIntArray counts, indices;
    fn.getTriangles(counts, indices);
    for (unsigned int i = 0; i < indices.length(); ++i)
      geometry.triangles.append(MFloatPoint(points[indices[i]]));
    appendLines(geometry.lines, points, (points.length() - 1) / 2);
    return true;
  }
};
struct Curves {
  struct Element {
    unsigned int panelId = 0;
    MPointArray points;
    unsigned int numU = 0, numV = 0;
    MPointArray lines;
  };
  std::vector<Element> elements;
  MPointArray lines;

  static bool samePointBits(const MPointArray &a, const MPointArray &b) {
    if (a.length() != b.length())
      return false;
    for (unsigned int i = 0; i < a.length(); ++i)
      if (std::memcmp(&a[i].x, &b[i].x, sizeof(double)) != 0 ||
          std::memcmp(&a[i].y, &b[i].y, sizeof(double)) != 0 ||
          std::memcmp(&a[i].z, &b[i].z, sizeof(double)) != 0 ||
          std::memcmp(&a[i].w, &b[i].w, sizeof(double)) != 0)
        return false;
    return true;
  }
  bool clearInvalidCache() {
    const bool changed = !elements.empty();
    elements.clear();
    lines.clear();
    return changed;
  }
  bool update(const std::vector<Element> &next) {
    if (next.empty())
      return clearInvalidCache();
    for (const auto &element : next) {
      if (element.numU == 0 || element.numV == 0 ||
          element.points.length() / element.numV != element.numU ||
          element.points.length() % element.numV != 0)
        return clearInvalidCache();
      for (unsigned int i = 0; i < element.points.length(); ++i) {
        const auto &point = element.points[i];
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z) || !std::isfinite(point.w))
          return clearInvalidCache();
      }
    }
    const size_t previousSize = elements.size();
    bool changed = previousSize != next.size();
    elements.resize(next.size());
    for (size_t i = 0; i < next.size(); ++i) {
      auto &element = elements[i];
      const auto &input = next[i];
      if (i < previousSize && element.panelId == input.panelId &&
          element.numU == input.numU && element.numV == input.numV &&
          samePointBits(element.points, input.points))
        continue;
      changed = true;
      element.panelId = input.panelId;
      element.points = input.points;
      element.numU = input.numU;
      element.numV = input.numV;
      element.lines.clear();
      for (unsigned int v = 0; v < element.numV; ++v)
        for (unsigned int u = 1; u < element.numU; ++u) {
          element.lines.append(element.points[(u - 1) * element.numV + v]);
          element.lines.append(element.points[u * element.numV + v]);
        }
    }
    if (changed) {
      lines.clear();
      for (const auto &element : elements)
        for (unsigned int i = 0; i < element.lines.length(); ++i)
          lines.append(element.lines[i]);
    }
    return changed;
  }
  bool update(const MPointArray &next, unsigned int uCount,
              unsigned int vCount) {
    return update(std::vector<Element>{{0, next, uCount, vCount, {}}});
  }
};
} // namespace ColliderDraw
