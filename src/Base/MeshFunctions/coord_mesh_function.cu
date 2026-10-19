#include "coord_mesh_function.h"

namespace LS_STAG
{
	LS_STAG_G_KERNEL_T(ElemType, CoordMeshFunction<ElemType>, set, const ElemType, value, LS_STAG_R_P(coord<ElemType>), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, CoordMeshFunction<ElemType>, copy, LS_STAG_CR_P(coord<ElemType>), _values,
		LS_STAG_R_P(coord<ElemType>), _Vector)
	LS_STAG_G_KERNEL_T(ElemType, CoordMeshFunction<ElemType>, add, const ElemType, _x, const ElemType, _y,
		const ElemType, _z, LS_STAG_R_P(coord<ElemType>), _Vector)
		
	template CoordMeshFunction<double>;
	template CoordMeshFunction<int>;
}