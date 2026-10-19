#include "ls_stag_mesh.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(LSSTAGMesh, buildNodeX, LS_STAG_R_PC, _nodes, LS_STAG_R_PC, _delta, const double, curDx, 
		const double, curDy, const double, curDz, const double, curShiftX, const double, curShiftY, 
		const double, curShiftZ, const int, curCountX, const int, xStartNum, const int, yStartNum, const int, zStartNum,
		const int, _N, const int, _NM)
}