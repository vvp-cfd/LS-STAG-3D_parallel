#include "flow_circles_mov_task_description.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(FlowCircleSMovTaskDescription, build, LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT,
		LS_STAG_R_PC, _v_ib_yx, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PI, _OnAirfoil, LS_STAG_CR_PC, _curV_C)
	LS_STAG_GR_KERNEL_2D(FlowCircleSMovTaskDescription, calcLevelSetCirclesMov, const int, hs, const int, hf, LS_STAG_R_PD, _phi, 
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, LS_STAG_CR_PI, _cellOnAirfoil, const coordinate*, _cur_C, const int, _N, const int, _NM, 
		const double, rad)
}