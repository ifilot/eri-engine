#pragma once

#include <eri/enums.h>
#include <eri/basis/cgf.h>

namespace eri::two_electron {

double eri(const eri::basis::CGF& a,
           const eri::basis::CGF& b,
           const eri::basis::CGF& c,
           const eri::basis::CGF& d,
           eri::enums::ERIMethod method = eri::enums::ERIMethod::Huzinaga);

} // namespace eri::one_electron