/*
* \copyright (C) 2013 Puzikova V.V.
* \file      mesh_function.h
* \author    Puzikova Valeria (vvp) 
* \date      05.10.2013 18:31:16
* \brief     mesh_function.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_MESH_FUNCTION_H_
#define _LS_STAG_MESH_FUNCTION_H_

// ================================================================== INCLUDES
#include "math.h"
#include "array.h"
#include "../profiler.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса MeshFunction и типа extDouble

Данный файл содержит в себе реализацию класса MeshFunction и определение типа extDouble
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_T(ElemType, copy, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, copy, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector,
		const ElemType, factor);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, scale, LS_STAG_R_P(ElemType), _Vector, const ElemType, factor);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, add, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, scaleAdd, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector, 
	    const ElemType, scale);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, add, const ElemType, value, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, subtract, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, set, const ElemType, value, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, set, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, set, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector, const ElemType, factor);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, add, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector,
		const ElemType, factor);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, diff, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, diff, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector, const ElemType, scale);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, setSqr, LS_STAG_CR_P(ElemType), _values, LS_STAG_CR_P(ElemType), _values2,
		LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, setSqr, LS_STAG_CR_P(ElemType), _values, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, setSqr, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_T(ElemType, setSqrScale, LS_STAG_R_P(ElemType), _Vector, const ElemType, scale);
	LS_STAG_DECLARATION_KERNEL_RT(ElemType, dot, LS_STAG_R_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_RT(ElemType, sum, LS_STAG_CR_P(ElemType), _Vector);
	LS_STAG_DECLARATION_KERNEL_RT(ElemType, dot, LS_STAG_CR_P(ElemType), _Vector, LS_STAG_CR_P(ElemType), _Vector2);
/*!
\brief Класс, реализующий работу с сеточными функциями

Данный класс реализует работу с сеточными функциями: в нем происходит вся работа с памятью (выделение, удаление,
изменение числа контрольных обьемов), выполнение различных арифметических операций, в т.ч. с другими сеточными функциями,
расчет нормы и т.д.
*/
	template<typename ElemType>// double or int
	struct MeshFunction : Array<ElemType>
	{
		LS_STAG_KERNEL(copy, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] = _values[i];
		}
		LS_STAG_KERNEL(copy, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector, const ElemType factor)
		{
			_Vector[i] = factor * _values[i];
		}
		LS_STAG_KERNEL(scale, LS_STAG_R_P(ElemType) _Vector, const ElemType factor)
		{
			_Vector[i] *= factor;
		}
		LS_STAG_KERNEL(add, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] += _values[i];
		}
		LS_STAG_KERNEL(scaleAdd, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector, const ElemType scale)
		{
			ElemType aux = _Vector[i];
			_Vector[i] = scale * aux + _values[i];
		}
		LS_STAG_KERNEL(add, const ElemType value, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] += value;
		}
		LS_STAG_KERNEL(subtract, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] -= _values[i];
		}
		LS_STAG_KERNEL(set, const ElemType value, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] = value;
		}
		LS_STAG_KERNEL(set, LS_STAG_CR_P(ElemType) _values, LS_STAG_CR_P(ElemType) _values2, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] = _values[i] + _values2[i];
		}
		LS_STAG_KERNEL(set, LS_STAG_CR_P(ElemType) _values, LS_STAG_CR_P(ElemType) _values2, LS_STAG_R_P(ElemType) _Vector,
			const ElemType factor)
		{
			_Vector[i] = factor * (_values[i] + _values2[i]);
		}
		LS_STAG_KERNEL(add, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector, const ElemType factor)
		{
			_Vector[i] += (factor * _values[i]);
		}
		LS_STAG_KERNEL(diff, LS_STAG_CR_P(ElemType) _values, LS_STAG_CR_P(ElemType) _values2, LS_STAG_R_P(ElemType) _Vector)
		{
			_Vector[i] = _values[i] - _values2[i];
		}
		LS_STAG_KERNEL(diff, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _values2,
			LS_STAG_R_P(ElemType) _Vector, const ElemType scale)
		{
			_values2[i] *= scale;
			_Vector[i] = _values[i] - _values2[i];
		}
		LS_STAG_KERNEL(setSqr, LS_STAG_CR_P(ElemType) _values, LS_STAG_CR_P(ElemType) _values2, LS_STAG_R_P(ElemType) _Vector)
		{
			ElemType aux = _values[i] + _values2[i];
			_Vector[i] = aux * aux;
		}
		LS_STAG_KERNEL(setSqr, LS_STAG_CR_P(ElemType) _values, LS_STAG_R_P(ElemType) _Vector)
		{
			ElemType aux = _values[i];
			_Vector[i] = aux * aux;
		}
		LS_STAG_KERNEL(setSqr, LS_STAG_R_P(ElemType) _Vector)
		{
			ElemType aux = _Vector[i];
			_Vector[i] = aux * aux;
		}
		LS_STAG_KERNEL(setSqrScale, LS_STAG_R_P(ElemType) _Vector, const ElemType scale)
		{
			ElemType aux = _Vector[i];
			_Vector[i] = scale * aux * aux;
		}
		LS_STAG_KERNEL_R(ElemType, dot, LS_STAG_R_P(ElemType) _Vector)
		{
			ElemType curElem = _Vector[i];
			return curElem * curElem;
		}
		LS_STAG_KERNEL_R(ElemType, sum, LS_STAG_CR_P(ElemType) _Vector)
		{
			return  _Vector[i];
		}
		LS_STAG_KERNEL_R(ElemType, dot, LS_STAG_CR_P(ElemType) _Vector, LS_STAG_CR_P(ElemType) _Vector2)
		{
			return _Vector[i] * _Vector2[i];
		}
		void set(const ElemType value)// если не double/int --- переопределить
		{
			LS_STAG_PROFILE("vSet");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, set, size, value, Vector)
		};
		void set(const ElemType value, const double part)// если не double/int --- переопределить
		{
			LS_STAG_PROFILE("vSetPart");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, set, (int)floor(size * part), value, Vector)
		};
		void setSum(const ElemType* values, const ElemType* values2)
		{
			LS_STAG_PROFILE("setSum");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, set, size, values, values2, Vector)
		};
		void setSumWithSqr(const ElemType* values, const ElemType* values2)
		{
			LS_STAG_PROFILE("setSumSqr");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, setSqr, size, values, values2, Vector)
		};
		void setSumWithScale(const ElemType* values, const ElemType* values2, const ElemType factor)
		{
			LS_STAG_PROFILE("setSumScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, set, size, values, values2, Vector, factor)
		};
		void addWithScale(const ElemType* values, const ElemType factor)
		{
			LS_STAG_PROFILE("addScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, add, size, values, Vector, factor)
		};
		void addWithScale(const MeshFunction<ElemType>* values, const ElemType factor)
		{
			ElemType* src;
			values->getCopy(src);
			addWithScale(src, factor);
		};
		void setDiff(const ElemType* values, const ElemType* values2)
		{
			LS_STAG_PROFILE("setDiff");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, diff, size, values, values2, Vector)
		};
		void setDiffWithScale(const ElemType* values, ElemType* values2, const ElemType scale)
		{
			LS_STAG_PROFILE("setDiffScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, diff, size, values, values2, Vector, scale)
		};
		void setSquare()
		{
			LS_STAG_PROFILE("setSqr");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, setSqr, size, Vector)
		};
		void setSquareScale(const ElemType scale)
		{
			LS_STAG_PROFILE("setSqrScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, setSqrScale, size, Vector, scale)
		};
		void setSquare(const ElemType* values)
		{
			LS_STAG_PROFILE("setSqr");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, setSqr, size, values, Vector)
		};
		void setSquare(MeshFunction<ElemType>* values)
		{
			LS_STAG_PROFILE("setSqr");
			ElemType* src;
			values->getCopy(src);
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, setSqr, size, src, Vector)
		}; 
		virtual double getNorm()
		{
			LS_STAG_PROFILE("norm");
			LS_STAG_PAR_SUM(MeshFunction<ElemType>, dot, ElemType, size, Vector)
			return sqrt(1.0 * LS_STAG_PAR_ACC_VALUE);
		};
		virtual ElemType getSum()
		{
			LS_STAG_PROFILE("sum");
			LS_STAG_RETURN_PAR_SUM(MeshFunction<ElemType>, sum, ElemType, size, Vector)
		}; 
		virtual void copy(const ElemType* values)//? copy
		{
			LS_STAG_PROFILE("vCopy");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, copy, size, values, Vector)
		};
		void copy(const MeshFunction<ElemType>* values)//? copy
		{
			ElemType* src;
			values->getCopy(src);
			copy(src);
		};
		void copyWithScale(const ElemType* values, const ElemType factor)//? copy
		{
			LS_STAG_PROFILE("vCopyScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, copy, size, values, Vector, factor)
		};
		void copyWithScale(const MeshFunction<ElemType>* values, const ElemType factor)//? copy
		{
			ElemType* src;
			values->getCopy(src);
			copyWithScale(src, factor);
		};
		void copyTo(ElemType* values, double divisor)//? copy
		{
			LS_STAG_PROFILE("vCopyTo");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, copy, size, Vector, values, (ElemType)(1.0 / divisor))
		};
		virtual void operator*=(const ElemType factor)
		{
			LS_STAG_PROFILE("vScale");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, scale, size, Vector, factor)
		};
		virtual void operator+=(const ElemType* values)
		{
			LS_STAG_PROFILE("vvP");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, add, size, values, Vector)
		};
		void scaleAdd(const ElemType* values, const ElemType scale)
		{
			LS_STAG_PROFILE("vScaleAdd");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, scaleAdd, size, values, Vector, scale)
		}; 
		void scaleAdd(const MeshFunction<ElemType>* values, const ElemType scale)
		{
			LS_STAG_PROFILE("vScaleAdd");
			ElemType* src;
			values->getCopy(src);
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, scaleAdd, size, src, Vector, scale)
		}; 
		void operator+=(const MeshFunction<ElemType>* values)
		{
			ElemType* src;
			values->getCopy(src);
			operator+=(src);
		};
		void operator+=(const ElemType value)
		{
			LS_STAG_PROFILE("vP");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, add, size, value, Vector)
		};
		virtual void operator-=(const ElemType* values)
		{
			LS_STAG_PROFILE("vvM");
			LS_STAG_PAR_FOR(MeshFunction<ElemType>, subtract, size, values, Vector)
		};
		void operator-=(const MeshFunction<ElemType>* values)
		{
			ElemType* src;
			values->getCopy(src);
			operator-=(src);
		};
		virtual ElemType operator*(const ElemType* Vector2)
		{
			LS_STAG_PROFILE("dot");
			LS_STAG_RETURN_PAR_SUM(MeshFunction<ElemType>, dot, ElemType, size, Vector, Vector2)
		}
	};
	typedef MeshFunction<double> extDouble;
}

#endif //_LS_STAG_MESH_FUNCTION_H_
