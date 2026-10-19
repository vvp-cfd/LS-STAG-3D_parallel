#include "smoother_adlj.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(SmootherADLJ, build, const int*, _Number, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
		LS_STAG_CR_PI, _ADLJ_yvo, LS_STAG_CR_PI, _ADLJ_zvo, LS_STAG_R_PD, _y_c, LS_STAG_R_PD, _x_c, LS_STAG_R_PD, _z_c, 
		LS_STAG_R_PD, _y_d, LS_STAG_R_PD, _x_d, LS_STAG_R_PD, _z_d,
		LS_STAG_R_PD, _y_e, LS_STAG_R_PD, _x_e, LS_STAG_R_PD, _z_e, const int, _N, const int, _NM)
}