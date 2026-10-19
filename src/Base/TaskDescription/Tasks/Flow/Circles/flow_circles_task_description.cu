#include "flow_circles_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowCircleSTaskDescription, cellsOnAir, LS_STAG_R_PI, _cellOnAirfoil, const int, airNum, const int, _xs,
		const int, _xf, const int, _N, const int, _NM)
	LS_STAG_G_KERNEL(FlowCircleSTaskDescription, vib, LS_STAG_CR_PD, _radiusA, LS_STAG_CR_PD, _OmegaF,
		LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PC, _delta, LS_STAG_R_PC, _alpha,
		LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_yx, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT, const int, _N)
	LS_STAG_GR_KERNEL_2D(FlowCircleSTaskDescription, calcLevelSetCircles, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, LS_STAG_CR_PI, _cellOnAirfoil, const coordinate*, _cur_C,
		const int, _N, const int, _NM, const double*, _radiusA)
}