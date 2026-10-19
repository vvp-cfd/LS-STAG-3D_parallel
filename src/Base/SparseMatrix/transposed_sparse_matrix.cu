#include "transposed_sparse_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransposedSparseMatrix, transpose, LS_STAG_R_PI, CNT, LS_STAG_R_PI, Number,
		LS_STAG_R_PI, Portrait, LS_STAG_R_PD, Cell, LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Portrait,
		LS_STAG_CR_PD, m_Cell, const int, elemCount, const int, NM, const int, CNM)
}