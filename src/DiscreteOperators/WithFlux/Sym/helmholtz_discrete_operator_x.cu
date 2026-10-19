#include "helmholtz_discrete_operator_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorX, calcHelmX, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, const coordinate*, _alpha, const coordinate*, _delta, const coordinate*, _inv_delta_ib, const coordinate*, _delta_ib_en,
		const coordinate*, _delta_ib_ws, LS_STAG_CR_PD, _inv_V, LS_STAG_CR_PD, _Mx, LS_STAG_CR_PCT, _type,
		const int, _N, const int, _M, const int, _NM, const double, _nu, const double, _coefIntegr, const bool, _isFlow)
	LS_STAG_GR_KERNEL_2D(HelmholtzDiscreteOperatorX, ibFluxHelmX, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType,
		const coordinate*, _theta, const coordinate*, _alpha, LS_STAG_CR_PC, _delta, const coordinate*, _inv_delta_ib, 
		const coordinate*, _delta_ib_en, const coordinate*, _delta_ib_ws, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PC, _VibZY, 
		 LS_STAG_CR_PC, _VibYX, LS_STAG_CR_PD, _inv_V, const int, _N, const int, _M, const int, _NM, const double, _nu)
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorX, setPlusOutFlow, LS_STAG_R_PD, _PlusOutFlow, LS_STAG_CR_PC, _delta,
		const double, _nu, const int, _M, const int, _N3)
	LS_STAG_G_KERNEL(HelmholtzDiscreteOperatorX, addFluxNS, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, U_S, const double, U_N, const int, _N, const int, _M, const int, _NM, const double, _nu)
	LS_STAG_G_KERNEL(HelmholtzDiscreteOperatorX, addFluxTB, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta, 
		const double, U_B, const double, U_T, const int, _N, const int, _M, const int, _K2, const int, _NMK2, const double, _nu)
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorX, inFlowFluxHelmX, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		LS_STAG_CR_PD, _U_W, const double, _U_E, const bool, isArtery, const bool, isCavity, const int, _N, const int, _M, 
		const int, _NM, const int, _vs, const int, _vf, const int, _zs, const int, _zf, const double, _nu)
}