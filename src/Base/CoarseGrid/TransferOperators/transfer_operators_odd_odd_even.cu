#include "transfer_operators_odd_odd_even.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddEven, numPooe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddEven, portraitPooe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddEven, numRooe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddEven, portraitRooe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2,
		const int, NP, const int, NMP, const int, KR, const int, nn)
}