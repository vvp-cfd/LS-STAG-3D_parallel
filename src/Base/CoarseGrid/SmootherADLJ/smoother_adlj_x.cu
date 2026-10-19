#include "smoother_adlj_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(SmootherADLJx, tdma, LS_STAG_R_PD, _solution, LS_STAG_CR_PD, _d, LS_STAG_CR_PD, _c,
		LS_STAG_CR_PD, _e, LS_STAG_CR_PD, _right_member, const int, _N, const int, _NM)
}