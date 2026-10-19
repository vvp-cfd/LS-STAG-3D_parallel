#include "transfer_operators_odd_even_odd.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenOdd, numPoeo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenOdd, numRoeo, LS_STAG_R_PI, _Number, const int, N2, const int, MR, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenOdd, portraitPoeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenOdd, portraitRoeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR,
		const int, NP, const int, K2, const int, nn)
}