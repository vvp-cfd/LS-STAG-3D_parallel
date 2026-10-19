#include "polynomial_derivative.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(PolynomialDerivative, calc, LS_STAG_R_PD, _coef, LS_STAG_CR_PD, _antiderivative_coef, const int, _length)
}