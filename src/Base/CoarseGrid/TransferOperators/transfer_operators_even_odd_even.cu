#include "transfer_operators_even_odd_even.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddEven, numPeoe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddEven, numReoe, LS_STAG_R_PI, _Number, const int, NR, const int, M2, const int, KR)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddEven, portraitPeoe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR)
	LS_STAG_G_S_KERNEL(TransferOperatorsEvenOddEven, portraitReoe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2, const int, K2,
		const int, NP, const int, NR, const int, MR, const int, KR, const int, NMP)
}