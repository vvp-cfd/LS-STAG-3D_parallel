#include "divergence_discrete_operator_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(DivergenceDiscreteOperatorX, calcDx, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, const coordinate*, _delta, const coordinate*, _theta)
	LS_STAG_G_S_KERNEL(DivergenceDiscreteOperatorX, numBuildDx, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK)
	LS_STAG_G_S_KERNEL(DivergenceDiscreteOperatorX, portraitBuildDx, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type,
		const int, m_M, const int, m_N, const int, m_K)
}