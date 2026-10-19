#include "laplace_discrete_operator.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(LaplaceDiscreteOperator, calcLap, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _delta, LS_STAG_R_PD, _inv_Mx, const double*, _inv_My,
		const double*, _inv_Mz, LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM)
	LS_STAG_GR_KERNEL_2D(LaplaceDiscreteOperator, ibFlux, LS_STAG_R_PD, _flux, const coordinate*, _theta,
		LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL_2D(LaplaceDiscreteOperator, initFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PC, _delta, 
		LS_STAG_CR_PD, _U_W, const double, _U_E, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_KERNEL_2D(LaplaceDiscreteOperator, setBaseFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PC, _delta, 
		LS_STAG_CR_PD, _Ux, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_KERNEL(LaplaceDiscreteOperator, addFluxNSp, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, V_S, const double, V_N,const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(LaplaceDiscreteOperator, addFluxTBp, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_PC, m_theta, LS_STAG_R_PC, m_delta, 
		 const double, W_B, const double, W_T, const int, _N, const int, _NMK2)
	LS_STAG_G_S_KERNEL(LaplaceDiscreteOperator, fixPointFlow, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, LS_STAG_CR_PI, Portrait, const int, numA)
	LS_STAG_G_S_KERNEL(LaplaceDiscreteOperator, fixPointFlowBase, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, const int, NM)
	LS_STAG_G_S_KERNEL(LaplaceDiscreteOperator, fixPointDef, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, const int, N2)
}