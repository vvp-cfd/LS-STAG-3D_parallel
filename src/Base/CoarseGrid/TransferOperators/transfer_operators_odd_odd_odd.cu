#include "transfer_operators_odd_odd_odd.h"

namespace LS_STAG
{
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddOdd, numPooo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2)
	LS_STAG_G_S_KERNEL(TransferOperatorsOddOddOdd, portraitPooo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2)
}