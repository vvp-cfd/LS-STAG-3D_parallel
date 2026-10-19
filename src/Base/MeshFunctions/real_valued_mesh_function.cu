#include "real_valued_mesh_function.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait,
		LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProdDiff, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait,
		LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _Vector2)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait,
		LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, const double, initVal)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait,
		LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, initVal)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait,
		LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, initVal, LS_STAG_CR_PD, initVal2)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProdScaleAdd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, 
		LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, const double, _scale, const double, _add)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProdScaleP, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
			LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _componentwiseMultiplier, 
			LS_STAG_CR_PD, _addComponentwiseMultiplier)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, mvProdScalePM, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
			LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _componentwiseMultiplier, 
			LS_STAG_CR_PD, _addComponentwiseMultiplier, const double, value)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, compProd, LS_STAG_CR_PD, _componentwiseMultiplier,
		LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, copyCompProd, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values, 
		LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, copyCompProdPerm, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_R_PD, _Vector, LS_STAG_CR_PI, _perm)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, copyCompProdPermTo, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_R_PD, _Vector, LS_STAG_CR_PI, _perm)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, copyCompProdSum, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values, 
		LS_STAG_CR_PD, _componentwiseMultiplier2, LS_STAG_CR_PD, _values2, LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, compProdScale, LS_STAG_CR_PD, _componentwiseMultiplier, const double, scale, LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, compProdScaleP, LS_STAG_CR_PD, _componentwiseMultiplier, const double, scale, 
			LS_STAG_CR_PD, _addValues, const double, factor, LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, compProdPM, LS_STAG_CR_PD, _componentwiseMultiplier,
		LS_STAG_CR_PD, _addValues, LS_STAG_CR_PD, _subtractValues, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _subtractValues2)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, compProdPM, LS_STAG_CR_PD, _componentwiseMultiplier,
		LS_STAG_CR_PD, _addValues, LS_STAG_CR_PD, _subtractValues, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _subtractValues2,
		LS_STAG_CR_PD, _subtractValues3)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, addSqrt, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, values, LS_STAG_CR_PD, values2)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, setSumSqrt, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, values, 
		LS_STAG_CR_PD, values2, LS_STAG_CR_PD, values3)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, subtract, LS_STAG_CR_PD, values, LS_STAG_R_PD, _Vector,
		const double, scale)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, subtract, LS_STAG_CR_PD, values, LS_STAG_CR_PD, values2,
		LS_STAG_CR_PD, values3, LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL(RealValuedMeshFunction, subtractAdd, LS_STAG_CR_PD, addValues, LS_STAG_CR_PD, subtractValues,
		LS_STAG_R_PD, _Vector)
	LS_STAG_G_KERNEL_R(double, RealValuedMeshFunction, mvRes, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _Vector2,
		LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier)
}