#include "gradient_discrete_operator_flux_y.h"

namespace LS_STAG
{
	LS_STAG_GR_KERNEL_2D(GradientDiscreteOperatorFluxY, computeFluxGy, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_inv_delta_ib,
		LS_STAG_CR_PC, m_theta, const int, m_N, const int, NM)
}