#include "divergence_discrete_operator_y.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(DivergenceDiscreteOperatorY, calcDy, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, const coordinate*, _delta, const coordinate*, _theta, const int, _N)
	LS_STAG_G_S_KERNEL(DivergenceDiscreteOperatorY, numBuildDy, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK)
	LS_STAG_G_S_KERNEL(DivergenceDiscreteOperatorY, portraitBuildDy, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type,
		const int, NM, const int, m_N, const int, m_K)
}