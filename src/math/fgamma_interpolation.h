#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include <eri/math/fgamma_interpolation.h>
#include "math/fgamma.h"

#include <cstddef>

namespace eri::math {

// Interpolated Boys function (falls back if outside table range)
double Fgamma_interp(int nu, double T);

// Compute F[0..nu_max] using interpolation + stable recurrence
void Fgamma_block_interp(int nu_max, double T, double* F);

} // namespace eri::math