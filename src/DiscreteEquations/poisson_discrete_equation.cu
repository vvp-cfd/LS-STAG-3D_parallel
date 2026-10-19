#include "poisson_discrete_equation.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(PoissonDiscreteEquation, fixPoint, LS_STAG_R_PD, m_right, const int, num, const bool, isFlow,
		const int, num2)
}