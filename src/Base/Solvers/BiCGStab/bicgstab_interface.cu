#include "bicgstab_interface.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(BiCGStabInterface, solUpd, LS_STAG_CR_PD, _s, LS_STAG_R_PD, _y, LS_STAG_R_PD, _z,
		LS_STAG_R_PD, _sol, LS_STAG_R_PD, _t, LS_STAG_R_PD, _r, const double, alpha, const double, w)
	LS_STAG_G_KERNEL_R(double, BiCGStabInterface, pUpd, LS_STAG_R_PD, _p, LS_STAG_CR_PD, _nu_v, LS_STAG_CR_PD, _r,
		const double, _coef, const double, _beta)
}