#pragma once

#include <eri/enums.h>
#include <eri/basis/cgf.h>

namespace eri::one_electron {

double kinetic(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               eri::enums::KineticMethod method = eri::enums::KineticMethod::Huzinaga);

} // namespace eri::one_electron
