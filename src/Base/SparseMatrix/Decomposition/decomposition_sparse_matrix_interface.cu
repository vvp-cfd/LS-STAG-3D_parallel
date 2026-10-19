#include "decomposition_sparse_matrix_interface.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(DecompositionSparseMatrixInterface, calcDirectGauss, LS_STAG_R_PD, Linv_b,
		LS_STAG_CR_PD, L_Cell, LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PD, _b,
		const int, NM)
	LS_STAG_G_S_KERNEL(DecompositionSparseMatrixInterface, calcInverseGauss, LS_STAG_R_PD, m_x,
		LS_STAG_CR_PD, Linv_b, LS_STAG_CR_PD, U_Cell, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, U_Portrait,
		const int, NM2)
}