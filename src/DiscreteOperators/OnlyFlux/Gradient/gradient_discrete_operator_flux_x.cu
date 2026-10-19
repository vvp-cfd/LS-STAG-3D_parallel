#include "gradient_discrete_operator_flux_x.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(GradientDiscreteOperatorFluxX, computeFluxGx, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_inv_delta_ib,
		LS_STAG_CR_PC, m_theta, const int, m_N, const int, m_M, const int, NM)
}