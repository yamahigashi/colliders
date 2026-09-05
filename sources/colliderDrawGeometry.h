#pragma once

#include <cmath>
#include <maya/MFloatPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MUIDrawManager.h>
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
  MPointArray unitPoints;
  MIntArray indices;
  int subdivision = 0;
  bool update(const std::vector<MMatrix> &next, int sides) {
    if (sides < 3)
      sides = 3;
    if (subdivision == sides && sameMatrices(matrices, next))
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
    geometry.triangles.clear();
    geometry.lines.clear();
    for (const auto &matrix : matrices) {
      MPointArray points;
      for (unsigned int i = 0; i < unitPoints.length(); ++i) {
        const MPoint p = unitPoints[i] * matrix;
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
  MPointArray points, lines;
  unsigned int numU = 0, numV = 0;
  bool update(const MPointArray &next, unsigned int uCount,
              unsigned int vCount) {
    if (numU == uCount && numV == vCount && samePoints(points, next))
      return false;
    points = next;
    numU = uCount;
    numV = vCount;
    lines.clear();
    if (points.length() != numU * numV)
      return true;
    for (unsigned int v = 0; v < numV; ++v)
      for (unsigned int u = 1; u < numU; ++u) {
        lines.append(points[(u - 1) * numV + v]);
        lines.append(points[u * numV + v]);
      }
    return true;
  }
};
} // namespace ColliderDraw
