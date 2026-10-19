#include "transfer_operators_odd_even_even.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenEven, numPoee, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenEven, numRoee, LS_STAG_R_PI, _Number, const int, nn, const int, MR, const int, KR)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenEven, portraitPoee, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR, const int, NMR)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddEvenEven, portraitRoee, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR, const int, KR,
		const int, NP, const int, NR, const int, NMP)
}