#include "flow_finite_cylinder_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowFiniteCylinderTaskDescription, calcLevelSetFiniteCylinder, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM,
		const double, rad, const double, halfLength, const double, ex, const double, ey, const double, ez)
}
