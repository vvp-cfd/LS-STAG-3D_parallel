#include "flow_rotor_darrieus_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowRotorDarrieusTaskDescription, calcLevelSetRotorDarrieus, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, BaseRadius,
		const int, numBlade, const double, radius, const double, step2, const double, step, const double, ANGLE,
		const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor, const double, len, const double, inv_len, 
		const double, delta, const double, _LE, const int, cU_L2, LS_STAG_CR_PD, cU_coef, const int, cdU_L2, 
		LS_STAG_CR_PD, cdU_coef, const int, cD_L2, LS_STAG_CR_PD, cD_coef, const int, cdD_L2, LS_STAG_CR_PD, cdD_coef, 
		const double, _eps)
}