#include "helmholtz_discrete_equation_y.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteEquationY, symY, LS_STAG_R_PD, _sol, const int, _N, const int, _NM)
	LS_STAG_GR_KERNEL_2D(HelmholtzDiscreteEquationY, correctionY, LS_STAG_R_PD, _right, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _Number, LS_STAG_CR_PCT, _type,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _theta_N, LS_STAG_CR_PC, _curV_C, LS_STAG_CR_PI, _CellOnAirfoill, const int, _hs,
		const int, _hf, const int, _N, const int, _NM)
}