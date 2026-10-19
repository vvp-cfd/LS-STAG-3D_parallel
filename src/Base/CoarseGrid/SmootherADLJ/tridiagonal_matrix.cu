#include "tridiagonal_matrix.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TriDiagonalMatrix, triLU, LS_STAG_R_PD, d, LS_STAG_R_PD, c, LS_STAG_R_PD, e, const int, NM)
}