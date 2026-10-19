#include "transfer_operators_even_even_even.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenEvenEven, numP, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenEvenEven, numR, LS_STAG_R_PI, _Number, const int, NMR)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenEvenEven, portraitP, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR, const int, NMR)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenEvenEven, portraitR, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR,
		const int, NP, const int, NR, const int, KR, const int, NMP)
}