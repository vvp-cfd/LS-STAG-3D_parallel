#include "smoother_adlj_y.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(SmootherADLJy, tdma, LS_STAG_R_PD, _solution, LS_STAG_R_PD, _aux_sol, LS_STAG_CR_PD, _d,
		LS_STAG_CR_PD, _c, LS_STAG_CR_PD, _e, LS_STAG_CR_PD, _right_member, LS_STAG_CR_PI, _ho, const int, _M2, const int, _NM)
	LS_STAG_G_KERNEL_2D(SmootherADLJy, initOrder, LS_STAG_R_PI, _ho, LS_STAG_R_PI, _vo, const int, _N, const int, _M, const int, _NM)
}