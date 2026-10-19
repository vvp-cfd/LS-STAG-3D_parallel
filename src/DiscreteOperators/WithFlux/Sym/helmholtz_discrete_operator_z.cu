#include "helmholtz_discrete_operator_z.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorZ, calcHelmZ, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, const coordinate*, _alpha, const coordinate*, _delta, 
		const coordinate*, _inv_delta_ib_tn, const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, 
		LS_STAG_CR_PD, _inv_V, LS_STAG_CR_PD, _Mz, LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, 
		const double, _nu, const double, _coefIntegr, const bool, _isFlow)
	LS_STAG_GR_KERNEL_2D(HelmholtzDiscreteOperatorZ, ibFluxHelmZ, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType,
		const coordinate*, _theta, const coordinate*, _alpha, LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _inv_delta_ib_tn, 
		const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PC, _VibZY, 
		 LS_STAG_CR_PC, _VibYX, LS_STAG_CR_PD, _inv_V, const int, _N, const int, _NM, const double, _nu)
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorZ, setPlusOutFlowZ, LS_STAG_R_PD, _PlusOutFlow, LS_STAG_CR_PC, _delta,
		const double, _nu, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_KERNEL(HelmholtzDiscreteOperatorZ, addFluxNShz, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, V_S, const double, V_N, const int, _N, const int, _NM, const double, _nu)
	LS_STAG_G_KERNEL(HelmholtzDiscreteOperatorZ, addFluxTBhz, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, V_B, const double, V_T, const int, _N, const int, _NMK2, const double, _nu)
	LS_STAG_G_KERNEL_2D(HelmholtzDiscreteOperatorZ, inFlowFluxhz, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta,
		const double, _V_W, const double, _V_E, const bool, isCavity, const int, _N, const int, _NM, const double, _nu)
}