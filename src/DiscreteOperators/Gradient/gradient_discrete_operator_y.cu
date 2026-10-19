#include "gradient_discrete_operator_y.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(GradientDiscreteOperatorY, calc, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell, LS_STAG_R_PC, _theta, 
		LS_STAG_R_PC, _delta, LS_STAG_CR_PCT, _type, const int, _N, const int, _NM)
	LS_STAG_G_S_KERNEL(GradientDiscreteOperatorY, portraitBuildGy, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, baseN,
		const int, NM, const int, m_K)
}