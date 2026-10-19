/*
* \copyright (C) 2013 Puzikova V.V.
* \file      coord_mesh_function.h
* \author    Puzikova Valeria (vvp) 
* \date      06.10.2013 11:42:00
* \brief     coord_mesh_function.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_COORD_MESH_FUNCTION_H_
#define _LS_STAG_COORD_MESH_FUNCTION_H_

// ================================================================== INCLUDES
#include <stdio.h>
#include "array.h"
#include "../profiler.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса CoordMeshFunction и определениями типов coordinate и tripletInt

Данный файл содержит в себе реализацию класса CoordMeshFunction и определения типов coordinate и tripletInt
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с двумерными координатами
*/
	template<typename ElemType>
	struct coord
	{
		LS_STAG_INLINE void set(const ElemType value)// если не double/int --- переопределить
		{
			x = value;
			y = value;
			z = value;
		};
		LS_STAG_INLINE void set(const ElemType _x, const ElemType _y, const ElemType _z)// если не double/int --- переопределить
		{
			x = _x;
			y = _y;
			z = _z;
		};
		LS_STAG_INLINE void set(const coord* values)
		{
			x = values->x;
			y = values->y;
			z = values->z;
		};
		LS_STAG_INLINE void setDiffXY(const coord& a, const coord& b)
		{
			x = a.x - b.x;
			y = a.y - b.y;
			z = 0.0;
		}
		LS_STAG_INLINE void setDiff(const coord& a, const coord& b)
		{
			x = a.x - b.x;
			y = a.y - b.y;
			z = a.z - b.z;
		}
		LS_STAG_INLINE void diff(const ElemType _x, const ElemType _y, const ElemType _z, const ElemType scale)
		{
			x -= scale * _x;
			y -= scale * _y;
			z -= scale * _z;
		}
		LS_STAG_INLINE void add(const ElemType _x, const ElemType _y, const ElemType _z)// если не double/int --- переопределить
		{
			x += _x;
			y += _y;
			z += _z;
		};
		LS_STAG_INLINE ElemType getMax() 
		{ 
			const ElemType aux = (x > y) ? x : y;
			return (aux > z) ? aux : z;
		}
		LS_STAG_INLINE bool inRange(const ElemType _min, const ElemType _max) const
		{
			return (((x > _min) && (x < _max)) || ((y > _min) && (y < _max)) || ((z > _min) && (z < _max)));
		}
		LS_STAG_INLINE ElemType getProdMax(ElemType _x, ElemType _y, ElemType _z) 
		{ 
			_x *= x; _y *= y; _z *= z; 
			const ElemType aux = (_x > _y) ? _x : _y;
			return (aux > _z) ? aux : _z;
		}
		LS_STAG_INLINE ElemType getProdMax(const coord* _cm) { return getProdMax(_cm->x, _cm->y, _cm->z); }
		LS_STAG_INLINE double norm() { return sqrt(x * x + y * y + z * z); }
		ElemType x;
		ElemType y;
		ElemType z;
	};
	LS_STAG_DECLARATION_KERNEL_T(ElemType, set, const ElemType, value, LS_STAG_R_P(coord<ElemType>), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, copy, LS_STAG_CR_P(coord<ElemType>), _values, LS_STAG_R_P(coord<ElemType>), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, add, const ElemType, _x, const ElemType, _y, const ElemType, _я, LS_STAG_R_P(coord<ElemType>), _Vector);
	/*!
	\brief Вывод в файл данных в координатном формате
	*/
	template<typename ElemType>
	std::ostream& operator<<(std::ofstream& os, const coord<ElemType>& p) { return os << p.x << " " << p.y << " " << p.z; }
	/*!
	\brief Класс, реализующий работу с координатными сеточными функциями

	Данный класс реализует работу с сеточными функциями, элементы которых представляют собой тройку значений по аналогии с координатами:
	в нем происходит вся работа с памятью (выделение, удаление,
	изменение числа контрольных обьемов), копирование, присвоение значений и т.д.
	*/
	template<typename ElemType>// double or int
	struct CoordMeshFunction : Array<coord<ElemType>>
	{
		LS_STAG_KERNEL(set, const ElemType value, LS_STAG_R_P(coord<ElemType>) _Vector)
		{
			_Vector[i].set(value);
		}
		LS_STAG_KERNEL(copy, LS_STAG_CR_P(coord<ElemType>) _values, LS_STAG_R_P(coord<ElemType>) _Vector)
		{
			_Vector[i].set(_values + i);
		}
		LS_STAG_KERNEL(add, const ElemType _x, const ElemType _y, const ElemType _z, LS_STAG_R_P(coord<ElemType>) _Vector)
		{
			_Vector[i].add(_x, _y, _z);
		}
		void set(const ElemType value)// если не double/int --- переопределить
		{
			LS_STAG_PROFILE("setCoord");
			LS_STAG_PAR_FOR(CoordMeshFunction<ElemType>, set, size, value, Vector)
		};
		void add(const ElemType _x, const ElemType _y, const ElemType _z)// если не double/int --- переопределить
		{
			LS_STAG_PROFILE("addCoord");
			LS_STAG_PAR_FOR(CoordMeshFunction<ElemType>, add, size, _x, _y, _z, Vector)
		};
		void copy(const coord<ElemType>* values)// ? copy
		{
			LS_STAG_PROFILE("copyCoord");
			LS_STAG_PAR_FOR(CoordMeshFunction<ElemType>, copy, size, values, Vector)
		};
	};
	typedef coord<double> coordinate;
	typedef coord<int> tripletInt;
}

#endif //_LS_STAG_COORD_MESH_FUNCTION_H_
