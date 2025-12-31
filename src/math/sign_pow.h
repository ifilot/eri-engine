#pragma once

namespace eri::math {

inline double sign_pow(int n) noexcept {
    return (n & 1) ? -1.0 : 1.0;
}

} // namespace eri::math
