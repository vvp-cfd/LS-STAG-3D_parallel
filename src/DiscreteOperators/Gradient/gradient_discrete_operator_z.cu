#include "gradient_discrete_operator_z.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(GradientDiscreteOperatorZ, calc, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _delta, LS_STAG_CR_PCT, _type)
	LS_STAG_G_S_KERNEL(GradientDiscreteOperatorZ, portraitBuildGz, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, baseNM, const int, NMK)
}