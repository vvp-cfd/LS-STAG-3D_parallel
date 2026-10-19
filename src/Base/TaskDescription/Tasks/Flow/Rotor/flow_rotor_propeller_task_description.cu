#include "flow_rotor_propeller_task_description.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(FlowRotorPropellerTaskDescription, calcLevelSetRotorPropeller, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, PropellerRadius,
		const double, PropellerWidth, const int, numBlade, const double, L_BASE, const double, step2, const double, step,
		const double, ANGLE, const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor)
}