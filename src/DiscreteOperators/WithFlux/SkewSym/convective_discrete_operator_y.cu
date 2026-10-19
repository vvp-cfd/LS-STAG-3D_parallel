#include "convective_discrete_operator_y.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorY, baseFluxY, LS_STAG_R_PD, _flux, LS_STAG_R_PC, _delta,
		const double*, _Ux, LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _U_W, const double, _V_W, const double, _UV_E,
		const bool, isArtery, const bool, isCavity, const int, _N, const int, _M3, const int, _NM, const int, _vs, 
		const int, _vf, const int, _zs, const int, _zf)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorY, calcConvY, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, coordinate*, _delta, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz,
		LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N2)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorY, ibFluxConvY, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType,
		LS_STAG_CR_PCT, _baseType, coordinate*, _theta, coordinate*, _delta, LS_STAG_CR_PC, _Vib,
		LS_STAG_CR_PC, _VibENT, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PD, _Uib,
		LS_STAG_CR_PD, _Ug, const int, _N, const int, _M3, const int, _NM, const int, baseNM)
	LS_STAG_G_KERNEL(ConvectiveDiscreteOperatorY, addFluxNSy, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		LS_STAG_CR_PD, m_Uy, const double, V_S, const double, V_N, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(ConvectiveDiscreteOperatorY, addFluxTBy, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta, 
		const double, VW_B, const double, VW_T, const int, _N, const int, _NMK2)
	LS_STAG_G_KERNEL_2D(ConvectiveDiscreteOperatorY, addFluxTBOutY, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_theta, LS_STAG_CR_PC, m_delta, LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _Uz,
		const int, _N, const int, _NMy, const int, _K, const int, _Nbase, const int, _NMbase, const int, _Kw)
}