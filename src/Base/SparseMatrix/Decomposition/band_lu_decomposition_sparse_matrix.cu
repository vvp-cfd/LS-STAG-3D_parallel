#include "band_lu_decomposition_sparse_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(BandLUDecompositionSparseMatrix, portraitLU, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, L_Number,
		LS_STAG_R_PI, U_Portrait, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PD, L_Cell, const int, L_m_N)
	LS_STAG_G_S_KERNEL(BandLUDecompositionSparseMatrix, numBuildBLU, LS_STAG_R_PI, L_Number, LS_STAG_R_PI, U_Number,
		const int, L_m_N, const int, width_band, const int, NM)
	LS_STAG_G_S_KERNEL(BandLUDecompositionSparseMatrix, calcLU, LS_STAG_R_PD, L_Cell, LS_STAG_R_PD, U_Cell,
		LS_STAG_CR_PD, LU_Cell, LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, LU_Number,
		LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PI, LU_Portrait, const int, width_band, const int, NM,
		const int, L_m_N)
}