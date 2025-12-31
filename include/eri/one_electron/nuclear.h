#pragma once

#include <eri/enums.h>
#include <eri/basis/cgf.h>

namespace eri::one_electron {

double nuclear(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               const std::array<double, 3>& C,
               eri::enums::NuclearMethod method = eri::enums::NuclearMethod::Huzinaga);

} // namespace eri::one_electron