#include "helmholtz_discrete_operator_base.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(HelmholtzDiscreteOperatorBase, addPlusOutFlow, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _numOutFlow,
		LS_STAG_CR_PD, _PlusOutFlow)
	LS_STAG_G_S_KERNEL(HelmholtzDiscreteOperatorBase, numPlusOutFlow, LS_STAG_R_PI, numOutFlow, LS_STAG_CR_PI, Number,
		const int, m_N, const int, m_M, const int, m_K)
}