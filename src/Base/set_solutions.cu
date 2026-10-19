#include "set_solutions.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(SetSolutions, corr, LS_STAG_R_PD, _p, LS_STAG_R_PD, _predP, LS_STAG_CR_PPI, _OnAirfoil,
		const double, _boundP, const double, _coefCorrectionPredP)
	LS_STAG_GR_KERNEL_2D(SetSolutions, corr, LS_STAG_R_PD, _u, LS_STAG_R_PD, _prevU, LS_STAG_CR_PCT, _xTypes,
		const double, value, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL_2D(SetSolutions, corr, LS_STAG_R_PD, _u, LS_STAG_R_PD, _prevU, LS_STAG_CR_PCT, _xTypes, LS_STAG_CR_PD, values, const int, _N, const int, _M, const int, _NM)
	LS_STAG_GR_KERNEL_2D(SetSolutions, corrOutFlow, LS_STAG_R_PD, _u, LS_STAG_R_PD, _v, LS_STAG_R_PD, _w, const double, vBase, const int, _N, const int, _M, const int, _NM)
	LS_STAG_GR_KERNEL_2D(SetSolutions, corrOutFlowTB, LS_STAG_R_PD, _w, const int, _N, const int, _NMK2)
	LS_STAG_G_KERNEL_2D(SetSolutions, copyOutflow, LS_STAG_R_PD, _u, LS_STAG_CR_PD, _prevU, const int, _N, const int, _NM)
}