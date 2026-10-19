#include "mesh_function.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, copy, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, copy, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector,
		const ElemType, factor)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, scale, LS_STAG_R_P(ElemType), _Vector, const ElemType, factor)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, add, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, add, const ElemType, value, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, scaleAdd, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector, 
	    const ElemType, scale)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, subtract, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, set, const ElemType, value, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, set, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, set, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector, const ElemType, factor)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, add, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector,
		const ElemType, factor)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, diff, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, diff, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector, const ElemType, scale)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, setSqr, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, setSqr, LS_STAG_CR_P(ElemType), _values,
		LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, setSqr, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, MeshFunction<ElemType>, setSqrScale, LS_STAG_R_P(ElemType), _Vector, const ElemType, scale);
	LS_STAG_G_KERNEL_RT(ElemType, MeshFunction<ElemType>, dot, LS_STAG_R_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_RT(ElemType, MeshFunction<ElemType>, sum, LS_STAG_CR_P(ElemType), _Vector)
	LS_STAG_G_KERNEL_RT(ElemType, MeshFunction<ElemType>, dot, LS_STAG_CR_P(ElemType), _Vector, LS_STAG_CR_P(ElemType), _Vector2)

	template MeshFunction<double>;
	template MeshFunction<int>;
}