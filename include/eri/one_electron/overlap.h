#pragma once

#include <eri/enums.h>
#include <eri/basis/cgf.h>

namespace eri::one_electron {

double overlap(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               eri::enums::OverlapMethod method = eri::enums::OverlapMethod::Huzinaga);

} // namespace eri::one_electron