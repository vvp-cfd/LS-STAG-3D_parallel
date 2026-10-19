#include "smoother_adlj_z.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_2D(SmootherADLJz, initOrderZ, LS_STAG_R_PI, _ho, LS_STAG_R_PI, _vo, const int, _N, const int, _K, const int, _NM, const int, _NK)
}