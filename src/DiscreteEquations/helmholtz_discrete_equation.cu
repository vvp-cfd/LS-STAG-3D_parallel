#include "helmholtz_discrete_equation.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(HelmholtzDiscreteEquation, initRM, LS_STAG_R_PD, _right, LS_STAG_CR_PD, _lastSol, LS_STAG_CR_PD, _prev_sol,
		LS_STAG_CR_PD, _VOF, const double, _coef2)
}