#include "gradient_discrete_operator_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(GradientDiscreteOperatorX, calc, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, LS_STAG_R_PC, _delta, LS_STAG_R_PC, _theta, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_S_KERNEL(GradientDiscreteOperatorX, portraitBuildGx, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, m_M,
		const int, m_N, const int, m_K)
}