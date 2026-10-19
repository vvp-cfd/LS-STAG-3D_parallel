#include "flow_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowTaskDescription, calcAlphaMovBlock, LS_STAG_R_PC, _alpha, LS_STAG_R_PD, _phi, const int, _hsMov, 
		const int, _hfMov, const int, _N, const int, _NM, const int, _zK)
	LS_STAG_GR_KERNEL_2D(FlowTaskDescription, calcThetaMovBlock, LS_STAG_R_PC, _theta, LS_STAG_R_PC, _alpha, 
		LS_STAG_CR_P(coord<SideTypes>), _baseType, const int, _hsMov, const int, _hfMov, const int, _N, const int, _NM, const int, _zK)
	LS_STAG_GR_KERNEL_2D(FlowTaskDescription, rebuildMovBlock, const int, hsMov, const int, hfMov, LS_STAG_R_PCT, _baseType,
		LS_STAG_R_P(coord<SideTypes>), _baseSideType, LS_STAG_R_PD, _phi, LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_R_PI, _cutcellCounts,
		LS_STAG_R_P(tripletInt*), _airfoilNums, const int, _N, const int, _M, const int, _NM, const int, _zK)
	LS_STAG_G_S_KERNEL(FlowTaskDescription, calcAlpha, LS_STAG_R_PC, _alpha, const int, hs_loc, const int, hf_loc,
		const int, vs_loc, const int, vf_loc, LS_STAG_R_PD, _phi, const int, xN, const int, yM, const int, zK, 
		const int, N, const int, NM)
	LS_STAG_G_S_KERNEL(FlowTaskDescription, calcTheta, LS_STAG_R_PC, _theta,  LS_STAG_R_PC, _alpha, 
		LS_STAG_CR_P(coord<SideTypes>), _baseType, const int, hs_loc, const int, hf_loc,
		const int, vs_loc, const int, vf_loc, const int, xN, const int, yM, const int, zK, 
		const int, N, const int, NM)
	LS_STAG_G_S_KERNEL(FlowTaskDescription, calcTypeOfCells, LS_STAG_R_PD, _phi, LS_STAG_R_PCT, _baseType,
		LS_STAG_R_P(coord<SideTypes>), _baseSideType, const int, hs_loc, const int, hf_loc, const int, vs_loc, const int, vf_loc,
		const int, xN, const int, yM, const int, zK, const int, N, const int, M, const int, NM, LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_R_PI, _cutcellCounts, 
		LS_STAG_R_P(tripletInt*), _airfoilNums)
}