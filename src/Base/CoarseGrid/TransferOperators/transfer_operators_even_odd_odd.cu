#include "transfer_operators_even_odd_odd.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddOdd, numPeoo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddOdd, numReoo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddOdd, portraitPeoo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddOdd, portraitReoo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2, const int, K2,
		const int, NP, const int, NR, const int, NMP)
}