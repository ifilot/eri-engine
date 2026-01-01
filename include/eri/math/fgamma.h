#pragma once

namespace eri::math {

double Fgamma_acc(int a, double x);

double Fgamma_nr(int a, double x);

// Initialize interpolation tables up to nu_table_max
// MUST be called once before using Fgamma_interp / Fgamma_block_interp
void init_Fgamma_interp_table(int nu_table_max);

double Fgamma_interp(int nu, double T);

void Fgamma_block_interp(int nu_max, double T, double* F);

} // namespace eri::math