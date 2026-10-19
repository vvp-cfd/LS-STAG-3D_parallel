#include "flow_cube_task_description.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(FlowCubeTaskDescription, calcAlphaCube, LS_STAG_R_PC, _alpha, const int, hs, const int, hf,
		const int, vs, const int, vf, const int, xN, const int, yM, const int, zK, const int, zs, const int, zf)
	LS_STAG_G_S_KERNEL(FlowCubeTaskDescription, calcThetaCube, LS_STAG_R_PC, _theta, const int, hs, const int, hf,
		const int, vs, const int, vf, const int, xN, const int, yM, const int, zK, const int, zs, const int, zf)
	LS_STAG_G_S_KERNEL(FlowCubeTaskDescription, calcTypeOfCellsCube, LS_STAG_R_PCT, _baseType, LS_STAG_R_P(coord<SideTypes>), _SideType,
		const int, hs, const int, hf, const int, vs, const int, vf, const int, zK, const int, xN, const int, yM, const int, N,
		const int, zs, const int, zf, LS_STAG_R_PI, _cutcellCounts, LS_STAG_R_P(tripletInt*), _airfoilNums)
}
