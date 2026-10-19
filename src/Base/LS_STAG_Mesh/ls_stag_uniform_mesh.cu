#include "ls_stag_uniform_mesh.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, corr, LS_STAG_R_PD, _vec, LS_STAG_CR_PCT, _type)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, cvu, const CellTypes*, _baseType, LS_STAG_R_PCT, _x_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshM, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, cvv, const CellTypes*, _baseType, LS_STAG_R_PCT, _y_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, cvw, const CellTypes*, _baseType, LS_STAG_R_PCT, _z_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, cvtxy, const CellTypes*, _y_Type, LS_STAG_R_PCT, _xy_Type, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, cvtxz, const CellTypes*, _z_Type, LS_STAG_R_PCT, _xz_Type, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, cvtyz, const CellTypes*, _z_Type, LS_STAG_R_PCT, _yz_Type, const int, _baseMeshN,
		const int, _baseMeshNM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, initIB, LS_STAG_R_PC, _curV_C, LS_STAG_R_PC, _cur_C)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, calcDeltaIB, LS_STAG_CR_PCT, _baseType, LS_STAG_R_PC, _theta, LS_STAG_R_PC, _delta,
		LS_STAG_R_PC, _delta_ib_ws, LS_STAG_R_PC, _delta_ib_en, LS_STAG_R_PC, _delta_ib_nb, LS_STAG_R_PC, _delta_ib_st, LS_STAG_R_PC, _delta_ib_max, LS_STAG_CR_P(coord<SideTypes>), _SideType, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, buildInvDeltaIb, LS_STAG_R_PC, invDib, LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, 
			LS_STAG_R_PC, _delta, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, buildInvBaseDeltaIb, LS_STAG_R_PC, invDibNE, LS_STAG_R_PC, invDibTN, 
		LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, buildDeltaNodes, LS_STAG_R_PC, _nodes, LS_STAG_R_PC, _delta,
		const double, deltaX, const double, deltaY, const double, deltaZ, const double, a, const double, b, 
		const double, c, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, mx, LS_STAG_R_PD, _x_VOF, LS_STAG_CR_PD, _VOF, const int, _baseMeshN, 
		const int, _baseMeshM, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, my, LS_STAG_R_PD, _y_VOF, LS_STAG_CR_PD, _VOF, const int, _baseMeshN, 
		const int, _baseMeshNM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, mz, LS_STAG_R_PD, _z_VOF, LS_STAG_CR_PD, _VOF, const int, _NM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, v, LS_STAG_R_PD, _VOF, LS_STAG_CR_PCT, _baseType, LS_STAG_R_PC, _theta,
		LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, vXY, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, vXZ, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, vYZ, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, mXY, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, mXZ, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM)
	LS_STAG_G_KERNEL_2D(LSSTAGUniformMesh, mYZ, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, _baseMeshNM)
	LS_STAG_G_KERNEL(LSSTAGUniformMesh, projectionVib, const double, fullV, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PC, _delta,
		LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_yx, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT, const int, _N)
	LS_STAG_G_S_KERNEL(LSSTAGUniformMesh, calcTypeOfCellBaseBox, LS_STAG_R_PCT, _baseType, const int, xN, const int, yM, const int, zK)
	LS_STAG_G_S_KERNEL(LSSTAGUniformMesh, calcTypeOfSideBaseBox, LS_STAG_R_P(coord<SideTypes>), _SideType, const int, xN, const int, yM, const int, zK)
	LS_STAG_G_S_KERNEL(LSSTAGUniformMesh, calcThetaBaseBox, LS_STAG_R_PC, _theta, const int, xN, const int, yM, const int, zK)
	LS_STAG_G_S_KERNEL(LSSTAGUniformMesh, calcAlphaBaseBox, LS_STAG_R_PC, _alpha, const int, xN, const int, yM, const int, zK)
}