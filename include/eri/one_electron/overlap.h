#pragma once

#include <eri/cgf.h>

namespace eri::one_electron {

enum class OverlapMethod {
    Huzinaga,
    Hellsing
};

double overlap(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               OverlapMethod method = OverlapMethod::Huzinaga);

} // namespace eri::one_electron