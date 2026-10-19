#include "flow_sphere_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowSphereTaskDescription, calcLevelSetSphere, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM,
		const double, rad)
}
