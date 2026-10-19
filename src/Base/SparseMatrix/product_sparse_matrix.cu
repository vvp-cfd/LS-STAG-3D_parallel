#include "product_sparse_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(ProductSparseMatrix, init, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell)
	LS_STAG_G_S_KERNEL(ProductSparseMatrix, numSpMM, LS_STAG_R_PI, Number, LS_STAG_R_PI, mark,
		LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Number2, LS_STAG_CR_PI, m_Portrait, LS_STAG_CR_PI, m_Portrait2, const int, NMK)
	LS_STAG_G_KERNEL(ProductSparseMatrix, spMM, LS_STAG_R_PI, Number, LS_STAG_R_PI, Portrait, LS_STAG_R_PD, Cell,
		LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Number2, LS_STAG_CR_PI, m_Portrait, LS_STAG_CR_PI, m_Portrait2,
		LS_STAG_CR_PD, m_Cell, LS_STAG_CR_PD, m_Cell2)
}