#include "lu_decomposition_sparse_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(LUDecompositionSparseMatrix, portraitLU, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, L_Number,
		LS_STAG_R_PI, U_Portrait, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PD, L_Cell)
	LS_STAG_G_KERNEL(LUDecompositionSparseMatrix, initScalePerm, LS_STAG_CR_PI, LU_Number, LS_STAG_CR_PI, LU_Portrait, LS_STAG_CR_PD, LU_Cell, LS_STAG_R_PD, scale, LS_STAG_R_PI, perm)
	LS_STAG_G_S_KERNEL(LUDecompositionSparseMatrix, numBuildLU, LS_STAG_R_PI, L_Number, LS_STAG_R_PI, U_Number,
		const int, NM)
	LS_STAG_G_S_KERNEL(LUDecompositionSparseMatrix, calcLU, LS_STAG_CR_PD, scale, LS_STAG_R_PI, perm, LS_STAG_R_PD, col, LS_STAG_R_PD, L_Cell, LS_STAG_R_PD, U_Cell,
		LS_STAG_CR_PD, LU_Cell, LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, LU_Number,
		LS_STAG_CR_PI, LU_Portrait, const int, NM)
}