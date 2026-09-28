#pragma once

namespace ColliderInput {
constexpr int kMinSubdivision = 3;
// Bound geometry allocation and evaluation work even for connected inputs.
constexpr int kMaxSubdivision = 4096;

inline bool validSubdivision(int value) {
    return value >= kMinSubdivision && value <= kMaxSubdivision;
}

inline bool validAxis(short value) {
    return value >= 0 && value <= 5;
}

inline bool validSkirtType(short value) {
    return value == 0 || value == 1;
}
} // namespace ColliderInput
