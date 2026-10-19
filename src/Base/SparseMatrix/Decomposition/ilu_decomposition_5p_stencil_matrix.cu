#include "ilu_decomposition_5p_stencil_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(ILUDecomposition5pStencilMatrix, calcInverseGauss5p, LS_STAG_R_PD, m_x, LS_STAG_CR_PD, Linv_b, 
		LS_STAG_CR_PD, U_Cell, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, U_Portrait, const int, levelNum, 
		const int, _N, const int, _N2)
	LS_STAG_G_KERNEL(ILUDecomposition5pStencilMatrix, calcDirectGauss5p, LS_STAG_R_PD, Linv_b, LS_STAG_CR_PD, L_Cell, 
		LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PD, _b, const int, levelNum, const int, _N, 
		const int, _N2)
}