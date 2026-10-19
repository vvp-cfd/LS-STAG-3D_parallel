#include "convective_discrete_operator_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorX, baseFluxX, LS_STAG_R_PD, _flux, LS_STAG_R_PC, _delta,
		LS_STAG_CR_PD, _Ux, LS_STAG_CR_PD, _U_W, const double, _U_E, const bool, isArtery, const bool, isCavity, 
		const int, _N, const int, _M, const int, _NM, const int, _vs, const int, _vf, const int, _zs, const int, _zf)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorX, calcConvX, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, LS_STAG_R_PC, _delta, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz,
		LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N3)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorX, ibFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType,
		LS_STAG_CR_PCT, _baseType, coordinate*, _theta, LS_STAG_R_PC, _delta, LS_STAG_CR_PC, _Vib,
		LS_STAG_CR_PC, _VibENT, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PD, _Uib,
		LS_STAG_CR_PD, _Ug, const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N3)
	LS_STAG_G_KERNEL(ConvectiveDiscreteOperatorX, addFluxNS, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, UV_S, const double, UV_N, const int, _N, const int, _M, const int, _NM)
	LS_STAG_G_KERNEL(ConvectiveDiscreteOperatorX, addFluxTB, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta, 
		const double, UW_B, const double, UW_T, const int, _N, const int, _NMK2)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorX, addFluxTBOut, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_theta, LS_STAG_CR_PC, m_delta, LS_STAG_CR_PD, _Ux, LS_STAG_CR_PD, _Uz,
		const int, _N, const int, _NMx, const int, _K, const int, _Nbase, const int, _NMbase, const int, _Kw, const int, _M2)
}