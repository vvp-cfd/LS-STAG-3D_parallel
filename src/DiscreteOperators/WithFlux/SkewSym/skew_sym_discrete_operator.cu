#include "skew_sym_discrete_operator.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(SkewSymDiscreteOperator, numBuildS, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NM)
	LS_STAG_G_S_KERNEL(SkewSymDiscreteOperator, portraitBuildS, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type,
		const int, NMK, const int, m_N, const int, m_NM)
}