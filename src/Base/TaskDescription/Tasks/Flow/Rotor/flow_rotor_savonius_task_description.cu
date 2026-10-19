#include "flow_rotor_savonius_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowRotorSavoniusTaskDescription, calcLevelSetRotorSavonius, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, AxisRadius,
		const int, numBlade, const double, BladeRadius, const double, BladeWidth, const double, step2, const double, step,
		const double, ANGLE, const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor)
}