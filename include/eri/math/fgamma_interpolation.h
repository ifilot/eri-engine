#pragma once

namespace eri::math {

// Initialize interpolation tables up to nu_table_max
// MUST be called once before using Fgamma_interp / Fgamma_block_interp
void init_Fgamma_interp_table(int nu_table_max);

} // namespace eri::math
