#include "loads_calculator.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(LoadsCalculator, calc, LS_STAG_CR_PCT, _type, LS_STAG_CR_PPI, _OnAirfoil, LS_STAG_R_PD, _ForceX,
		LS_STAG_R_PD, _ForceY, LS_STAG_R_PD, _ForceZ, LS_STAG_CR_PC, _Vib, LS_STAG_R_PC, _VibYX, LS_STAG_R_PC, _VibZY,
		LS_STAG_CR_PC, _delta, LS_STAG_R_PC, _inv_delta_ib_ne, LS_STAG_R_PC, _inv_delta_ib_tn, 
		LS_STAG_CR_PC, _delta_ib_max, LS_STAG_CR_PC, _theta, LS_STAG_CR_PD, _Ux, LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _Uz,
		LS_STAG_CR_PD, _P, LS_STAG_CR_PD, _invV, LS_STAG_CR_PD, _nu, double, nu, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_KERNEL_R(double, LoadsCalculator, calcMoment, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceX, LS_STAG_CR_PD, _ForceY, LS_STAG_CR_PC, _curC)
	LS_STAG_G_KERNEL_R(double, LoadsCalculator, calcMomentX, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceY, LS_STAG_CR_PD, _ForceZ, LS_STAG_CR_PC, _curC)
	LS_STAG_G_KERNEL_R(double, LoadsCalculator, calcMomentY, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceX, LS_STAG_CR_PD, _ForceZ, LS_STAG_CR_PC, _curC)
}