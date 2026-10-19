#include "ls_stag_mesh_mov.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvuMov, const CellTypes*, _baseType, LS_STAG_R_PCT, _x_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshM, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvvMov, const CellTypes*, _baseType, LS_STAG_R_PCT, _y_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvw, const CellTypes*, _baseType, LS_STAG_R_PCT, _z_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvtxyMov, const CellTypes*, _y_Type, LS_STAG_R_PCT, _xy_Type, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvtxzMov, const CellTypes*, _z_Type, LS_STAG_R_PCT, _xz_Type, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, cvtyzMov, const CellTypes*, _z_Type, LS_STAG_R_PCT, _yz_Type, const int, _baseMeshN,
		const int, _baseMeshNM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, buildInvDeltaIb, LS_STAG_R_PC, invDib, LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, 
			LS_STAG_R_PC, _delta, const int, _N, const int, _NM, const int, _hsMov, const int, _hfMov)
	LS_STAG_GR_KERNEL_2D(LSSTAGMeshMov, buildInvBaseDeltaIb, LS_STAG_R_PC, invDibNE, LS_STAG_R_PC, invDibTN, 
		LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM,
	    const int, _hsMov, const int, _hfMov)
	LS_STAG_G_KERNEL(LSSTAGMeshMov, linVel, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, coordinate*, _v_ib, LS_STAG_CR_PC, _curV_C,
		coordinate*, _v_ib_yx, coordinate*, _v_ib_zy, coordinate*, _v_ib_ENT)
	LS_STAG_G_KERNEL(LSSTAGMeshMov, angVel, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, coordinate*, _v_ib, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _alpha,
		LS_STAG_CR_PC, P_C, coordinate*, _v_ib_yx, coordinate*, _v_ib_zy, coordinate*, _v_ib_ENT, const double, Omega, const int, _N)
}