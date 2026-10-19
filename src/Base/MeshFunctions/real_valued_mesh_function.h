/*
* \copyright (C) 2013 Puzikova V.V.
* \file      real_valued_mesh_function.h
* \author    Puzikova Valeria (vvp) 
* \date      07.10.2013 9:32:02
* \brief     real_valued_mesh_function.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_REAL_VALUED_MESH_FUNCTION_H_
#define _LS_STAG_REAL_VALUED_MESH_FUNCTION_H_

// ================================================================== INCLUDES
#include "mesh_function.h"
#include "../SparseMatrix/sparse_matrix.h"


// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса RealValuedMeshFunction

Данный файл содержит в себе реализацию класса RealValuedMeshFunction
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, 
		_multiplier, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(mvProdDiff, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
		LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _Vector2);
	LS_STAG_DECLARATION_KERNEL(mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, 
		_multiplier, LS_STAG_R_PD, _Vector, const double, initVal);
	LS_STAG_DECLARATION_KERNEL(mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, 
		_multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, initVal);
	LS_STAG_DECLARATION_KERNEL(mvProd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, 
		_multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, initVal, LS_STAG_CR_PD, initVal2);
	LS_STAG_DECLARATION_KERNEL(mvProdScaleAdd, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, 
		LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, const double, _scale, const double, _add);
	LS_STAG_DECLARATION_KERNEL(mvProdScaleP, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
			LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _componentwiseMultiplier, 
			LS_STAG_CR_PD, _addComponentwiseMultiplier);
	LS_STAG_DECLARATION_KERNEL(mvProdScalePM, LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
			LS_STAG_CR_PD, _multiplier, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _componentwiseMultiplier, 
			LS_STAG_CR_PD, _addComponentwiseMultiplier, const double, value);
	LS_STAG_DECLARATION_KERNEL(addSqrt, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, values, LS_STAG_CR_PD, values2);
	LS_STAG_DECLARATION_KERNEL(setSumSqrt, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, values, LS_STAG_CR_PD, values2, 
		LS_STAG_CR_PD, values3);
	LS_STAG_DECLARATION_KERNEL(compProd, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(copyCompProd, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(copyCompProdPerm, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_R_PD, _Vector, LS_STAG_CR_PI, _perm);
	LS_STAG_DECLARATION_KERNEL(copyCompProdPermTo, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_R_PD, _Vector, LS_STAG_CR_PI, _perm);
	LS_STAG_DECLARATION_KERNEL(copyCompProdSum, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _values,
		LS_STAG_CR_PD, _componentwiseMultiplier2, LS_STAG_CR_PD, _values2, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(compProdScale, LS_STAG_CR_PD, _componentwiseMultiplier, const double, scale, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(compProdScaleP, LS_STAG_CR_PD, _componentwiseMultiplier, const double, scale, 
			LS_STAG_CR_PD, _addValues, const double, factor, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL(compProdPM, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _addValues, LS_STAG_CR_PD, 
		_subtractValues, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _subtractValues2);
	LS_STAG_DECLARATION_KERNEL(compProdPM, LS_STAG_CR_PD, _componentwiseMultiplier, LS_STAG_CR_PD, _addValues, LS_STAG_CR_PD, 
		_subtractValues, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _subtractValues2, LS_STAG_CR_PD, _subtractValues3);
	LS_STAG_DECLARATION_KERNEL(subtract, LS_STAG_CR_PD, values, LS_STAG_R_PD, _Vector, const double, scale);
	LS_STAG_DECLARATION_KERNEL(subtract, LS_STAG_CR_PD, values, LS_STAG_CR_PD, values2, LS_STAG_CR_PD, values3, LS_STAG_R_PD, 
		_Vector);
	LS_STAG_DECLARATION_KERNEL(subtractAdd, LS_STAG_CR_PD, addValues, LS_STAG_CR_PD, subtractValues, LS_STAG_R_PD, _Vector);
	LS_STAG_DECLARATION_KERNEL_R(double, mvRes, LS_STAG_R_PD, _Vector, LS_STAG_CR_PD, _Vector2,
		LS_STAG_CR_PI, _Num, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell, LS_STAG_CR_PD, _multiplier);
/*!
\brief Класс, реализующий работу с сеточными функциями с элементами типа double

Данный производный класс реализует работу с сеточными действительнозначными функциями.
Помимо функционала, определенного в базовом классе MeshFunction, в RealValuedMeshFunction 
определены операции матрично-векторного умножения и арифметические операции над несколькими сеточными функциями.
*/
	struct RealValuedMeshFunction : extDouble
	{
		LS_STAG_KERNEL(mvProd, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell, 
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, const double initVal = 0.0)
		{
			double aux = initVal;
			for (int j = _Num[i], num = _Num[i + 1]; j < num; j++)
				aux += _Cell[j] * _multiplier[_Portrait[j]];
			_Vector[i] = aux;
		}
		LS_STAG_KERNEL(mvProdScaleAdd, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell, 
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, const double _scale, const double _add)
		{
			double aux = 0.0;
			for (int j = _Num[i], num = _Num[i + 1]; j < num; j++)
				aux += _Cell[j] * _multiplier[_Portrait[j]];
			_Vector[i] = _scale * aux + _add;
		}
		LS_STAG_KERNEL(mvProd, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell,
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, LS_STAG_CR_PD initVal)
		{
			RealValuedMeshFunction::_k_mvProd(i, _Num, _Portrait, _Cell, _multiplier, _Vector, initVal[i]);
		}
		LS_STAG_KERNEL(mvProd, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell,
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, LS_STAG_CR_PD initVal, LS_STAG_CR_PD initVal2)
		{
			RealValuedMeshFunction::_k_mvProd(i, _Num, _Portrait, _Cell, _multiplier, _Vector, initVal[i] + initVal2[i]);
		}
		LS_STAG_KERNEL(mvProdScaleP, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell,
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _componentwiseMultiplier, 
			LS_STAG_CR_PD _addComponentwiseMultiplier)
		{
			RealValuedMeshFunction::_k_mvProdScaleAdd(i, _Num, _Portrait, _Cell, _multiplier, _Vector, 
				_componentwiseMultiplier[i],  _addComponentwiseMultiplier[i] * _multiplier[i]);
		}
		LS_STAG_KERNEL(mvProdScalePM, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell,
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _componentwiseMultiplier, 
			LS_STAG_CR_PD _addComponentwiseMultiplier, const double value)
		{
			RealValuedMeshFunction::_k_mvProdScaleAdd(i, _Num, _Portrait, _Cell, _multiplier, _Vector, 
				_componentwiseMultiplier[i],  _addComponentwiseMultiplier[i] * (_multiplier[i] - value));
		}
		LS_STAG_KERNEL(compProd, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_R_PD _Vector)
		{
			_Vector[i] *= _componentwiseMultiplier[i];
		}
		LS_STAG_KERNEL(copyCompProd, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _values, LS_STAG_R_PD _Vector)
		{
			_Vector[i] = _componentwiseMultiplier[i] * _values[i];
		}
		LS_STAG_KERNEL(copyCompProdPerm, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _values, LS_STAG_R_PD _Vector, LS_STAG_CR_PI _perm)
		{
			const int iPerm = _perm[i];
			_Vector[i] = _componentwiseMultiplier[iPerm] * _values[iPerm];
		}
		LS_STAG_KERNEL(copyCompProdPermTo, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _values, LS_STAG_R_PD _Vector, LS_STAG_CR_PI _perm)
		{
			const int iPerm = _perm[i];
			_Vector[iPerm] = _componentwiseMultiplier[iPerm] * _values[i];
		}
		LS_STAG_KERNEL(copyCompProdSum, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _values, 
			LS_STAG_CR_PD _componentwiseMultiplier2, LS_STAG_CR_PD _values2, LS_STAG_R_PD _Vector)
		{
			_Vector[i] = _componentwiseMultiplier[i] * _values[i] + _componentwiseMultiplier2[i] * _values2[i];
		}
		LS_STAG_KERNEL(compProdPM, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _addValues,
			LS_STAG_CR_PD _subtractValues, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _subtractValues2)
		{
			_Vector[i] *= _componentwiseMultiplier[i];
			double aux = _addValues[i] - _subtractValues[i] - _subtractValues2[i];
			_Vector[i] += aux;
		}
		LS_STAG_KERNEL(compProdScaleP, LS_STAG_CR_PD _componentwiseMultiplier, const double scale, 
			LS_STAG_CR_PD _addValues, const double factor, LS_STAG_R_PD _Vector)
		{
			double aux = _componentwiseMultiplier[i] * scale;
			_Vector[i] *= aux;
			_Vector[i] += (_addValues[i] * factor);
		}
		LS_STAG_KERNEL(compProdScale, LS_STAG_CR_PD _componentwiseMultiplier, const double scale, LS_STAG_R_PD _Vector)
		{
			double aux = _componentwiseMultiplier[i] * scale;
			_Vector[i] *= aux;
		}
		LS_STAG_KERNEL(compProdPM, LS_STAG_CR_PD _componentwiseMultiplier, LS_STAG_CR_PD _addValues,
			LS_STAG_CR_PD _subtractValues, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _subtractValues2, LS_STAG_CR_PD _subtractValues3)
		{
			_Vector[i] *= _componentwiseMultiplier[i];
			double aux = _addValues[i] - _subtractValues[i] - _subtractValues2[i] - _subtractValues3[i];
			_Vector[i] += aux;
		}
		LS_STAG_KERNEL(subtract, LS_STAG_CR_PD values, LS_STAG_R_PD _Vector, const double scale)
		{
			double aux = values[i] * scale;
			_Vector[i] -= aux;
		}
		LS_STAG_KERNEL(subtract, LS_STAG_CR_PD values, LS_STAG_CR_PD values2, LS_STAG_CR_PD values3, LS_STAG_R_PD _Vector)
		{
			double aux = values[i] + values2[i] + values3[i];
			_Vector[i] -= aux;
		}
		LS_STAG_KERNEL(addSqrt, LS_STAG_R_PD _Vector, LS_STAG_CR_PD values, LS_STAG_CR_PD values2)
		{
			double aux = _Vector[i] + values[i] + values2[i];
			_Vector[i] = sqrt(aux);
		}
		LS_STAG_KERNEL(setSumSqrt, LS_STAG_R_PD _Vector, LS_STAG_CR_PD values, LS_STAG_CR_PD values2, LS_STAG_CR_PD values3)
		{
			double aux = values3[i] + values[i] + values2[i];
			_Vector[i] = sqrt(aux);
		}
		LS_STAG_KERNEL(subtractAdd, LS_STAG_CR_PD addValues, LS_STAG_CR_PD subtractValues, LS_STAG_R_PD _Vector)
		{
			double aux = addValues[i] - subtractValues[i];
			_Vector[i] += aux;
		}
		LS_STAG_KERNEL(mvProdDiff, LS_STAG_CR_PI _Num, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell,
			LS_STAG_CR_PD _multiplier, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _Vector2)
		{
			double aux = _Vector2[i];
			for (int j = _Num[i], num = _Num[i + 1]; j < num; j++)
				aux -= _Cell[j] * _multiplier[_Portrait[j]];
			_Vector[i] = aux;
		}
		LS_STAG_KERNEL_R(double, mvRes, LS_STAG_R_PD _Vector, LS_STAG_CR_PD _Vector2, LS_STAG_CR_PI _Num,
			LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell, LS_STAG_CR_PD _multiplier)
		{
			double aux = _Vector2[i];
			for (int j = _Num[i], num = _Num[i + 1]; j < num; j++)
				aux -= _Cell[j] * _multiplier[_Portrait[j]];
			_Vector[i] = aux;
			return aux * aux;
		}
		void setProduct()
		{
			LS_STAG_PROFILE("SpMV");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProd, size, m_Num, m_Portrait, m_Cell, multiplier, Vector)
			/*
			char transa = 'N';
			char matdescra[4] = {'G','L','N','C'};
			double val = 1.0;
			int size2 = m_Portrait[m_Num[size]-1];
			 mkl_dcsrmv(&transa, &size, &size2, &val, matdescra, m_Cell, m_Portrait, m_Num, m_Num + 1, multiplier, &val, Vector );
			 */
		 };
		#ifdef LS_STAG_USE_MKL
			void setProduct(const double* _multiplier, double* _Vector, const int _size)// для FGMRES из MKL
			{
				LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProd, _size, m_Num, m_Portrait, m_Cell, _multiplier, _Vector)
			 };
		#endif
		void setProductWithSum(double* Vector2)
		{
			LS_STAG_PROFILE("SpMV_P");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProd, size, m_Num, m_Portrait, m_Cell, multiplier, Vector, Vector2)
		};
		void setProductDiff(double* _multiplier, double* Vector2)
		{
			LS_STAG_PROFILE("SpMV_Diff");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProdDiff, size, m_Num, m_Portrait, m_Cell, _multiplier, Vector, Vector2)
		};
		void setProductWithSum(double* Vector2, double* Vector3)
		{
			LS_STAG_PROFILE("SpMV_PP");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProd, size, m_Num, m_Portrait, m_Cell, multiplier, Vector, Vector2, Vector3)
		};
		 double setProductWithResidual(double* _multiplier, double* Vector2)
		 {
			 LS_STAG_PROFILE("SpMV_Res");
			 LS_STAG_PAR_SUM(RealValuedMeshFunction, mvRes, double, size, Vector, Vector2, m_Num, m_Portrait, m_Cell, _multiplier)
			 return sqrt(LS_STAG_PAR_ACC_VALUE);
		 };
		 double setProductWithResidual(double* Vector2)
		 {
			 LS_STAG_PROFILE("SpMV_Res");
			 LS_STAG_PAR_SUM(RealValuedMeshFunction, mvRes, double, size, Vector, Vector2, m_Num, m_Portrait, m_Cell, multiplier)
				 return sqrt(LS_STAG_PAR_ACC_VALUE);
		 };
		void setProduct(SparseMatrix* Matrix, double* Vector2)
		{
			m_A = Matrix;
			multiplier = Vector2;
			setProduct();
		 };
		void setProduct(SparseMatrix* Matrix)
		{
			m_A = Matrix;
			setMatrComp();
			setProduct();
		 };
		void setProduct(double* Vector2)
		{
			multiplier = Vector2;
			setProduct();
		 };
		void setOperator(SparseMatrix* _A)
		{
			m_A = _A;
			setMatrComp();
		};
		void setMultiplier(extDouble* _multiplier)
		{
			_multiplier->getCopy(multiplier);
		};
		void setMultipliers(SparseMatrix* _A, double* _multiplier)
		{
			m_A = _A;
			multiplier = _multiplier;
			setMatrComp();
		};
		void setMultipliers(SparseMatrix* _A, extDouble* _multiplier)
		{
			m_A = _A;
			_multiplier->getCopy(multiplier);
			setMatrComp();
		};
		void setMultipliers(SparseMatrix* _A, double* _multiplier, double* _componentwiseMultiplier)
		{
			m_A = _A;
			multiplier = _multiplier;
			setMatrComp();
			m_componentwiseMultiplier = _componentwiseMultiplier;
		};
		void setComponentwiseMultiplier(double* _componentwiseMultiplier)
		{
			m_componentwiseMultiplier = _componentwiseMultiplier;
		};
		void setComponentwiseMultiplier(extDouble* _componentwiseMultiplier)
		{
			_componentwiseMultiplier->getCopy(m_componentwiseMultiplier);
		};
		void copyR(RealValuedMeshFunction* values)
		{
			double* src;
			values->getCopy(src);
			extDouble::copy(src);
		}
		void componentwiseProduct()
		{
			LS_STAG_PROFILE("CompProd");
			#ifdef LS_STAG_USE_MKL
				vdMul(size, Vector, m_componentwiseMultiplier, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, compProd, size, m_componentwiseMultiplier, Vector)
			#endif
		};
		void componentwiseProductScale(const double scale)
		{
			LS_STAG_PROFILE("CompProdScale");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdScale, size, m_componentwiseMultiplier, scale, Vector)
		};
		void componentwiseProductScale(double* _componentwiseMultiplier, const double scale)
		{
			LS_STAG_PROFILE("CompProdScale");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdScale, size, _componentwiseMultiplier, scale, Vector)
		};
		void componentwiseProductScaleAdd(const double scale, const double factor, const double* addValues)
		{
			LS_STAG_PROFILE("CompProdScaleP");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdScaleP, size, m_componentwiseMultiplier, scale, addValues, factor, Vector)
		};
		void componentwiseProductScaleAdd(double* _componentwiseMultiplier, const double scale, const double factor, const double* addValues)
		{
			LS_STAG_PROFILE("CompProdScaleP");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdScaleP, size, _componentwiseMultiplier, scale, addValues, factor, Vector)
		};
		void setProductScaleAdd(const double _scale, const double _add)
		{
			LS_STAG_PROFILE("SpMV_scaleAdd");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProdScaleAdd, size, m_Num, m_Portrait, m_Cell, multiplier, Vector, 
			_scale, _add)
		};
		void setProductScaleP(double* _addComponentwiseMultiplier)
		{
			LS_STAG_PROFILE("SpMV_scaleP");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProdScaleP, size, m_Num, m_Portrait, m_Cell, multiplier, Vector, 
			m_componentwiseMultiplier, _addComponentwiseMultiplier)
		};
		void setProductScalePM(double* _addComponentwiseMultiplier, const double value)
		{
			LS_STAG_PROFILE("SpMV_scalePM");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, mvProdScalePM, size, m_Num, m_Portrait, m_Cell, multiplier, Vector, 
			m_componentwiseMultiplier, _addComponentwiseMultiplier, value)
		};
		void componentwiseProductAddAndSubtract(const double* addValues, const double* subtractValues, const double* subtractValues2)
		{
			LS_STAG_PROFILE("CompProdPM");
			#ifdef LS_STAG_USE_MKL
				vdMul(size, Vector, m_componentwiseMultiplier, Vector);
				vdAdd(size, Vector, addValues, Vector);
				vdSub(size, Vector, subtractValues, Vector);
				vdSub(size, Vector, subtractValues2, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdPM, size, m_componentwiseMultiplier, addValues, subtractValues, Vector, subtractValues2)
			#endif
		};
		void componentwiseProductAddAndSubtract(const double* addValues, const double* subtractValues, const double* subtractValues2, const double* subtractValues3)
		{
			LS_STAG_PROFILE("CompProdPMM");
			#ifdef LS_STAG_USE_MKL
				vdMul(size, Vector, m_componentwiseMultiplier, Vector);
				vdAdd(size, Vector, addValues, Vector);
				vdSub(size, Vector, subtractValues, Vector);
				vdSub(size, Vector, subtractValues2, Vector);
				vdSub(size, Vector, subtractValues3, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, compProdPM, size, m_componentwiseMultiplier, addValues, subtractValues, Vector, subtractValues2, subtractValues3)
			#endif
		};
		void componentwiseProduct(double* _componentwiseMultiplier)
		{
			LS_STAG_PROFILE("CompProd");
			#ifdef LS_STAG_USE_MKL
				vdMul(size, Vector, _componentwiseMultiplier, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, compProd, size, _componentwiseMultiplier, Vector)
			#endif
		};
		void copyComponentwiseProduct(const double* _values)
		{
			LS_STAG_PROFILE("CopyCompProd");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProd, size, m_componentwiseMultiplier, _values, Vector)
		};
		void copyComponentwiseProductTo(double* _values)
		{
			LS_STAG_PROFILE("CopyCompProdTo");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProd, size, m_componentwiseMultiplier, Vector, _values)
		};
		void copyComponentwiseProductPerm(const double* _values, const int* _perm)
		{
			LS_STAG_PROFILE("CopyCompProdPerm");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProdPerm, size, m_componentwiseMultiplier, _values, Vector, _perm)
		};
		void copyComponentwiseProductPermTo(double* _values, const int* _perm)
		{
			LS_STAG_PROFILE("CopyCompProdPermTo");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProdPermTo, size, m_componentwiseMultiplier, Vector, _values, _perm)
		};
		void copyComponentwiseProductSum(const double* _cm, const double* _values, const double* _cm2, 
			const double* _values2)
		{
			LS_STAG_PROFILE("CopyCompProdSum");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProdSum, size, _cm, _values, _cm2, _values2, Vector)
		};
		void copyComponentwiseProduct(extDouble* _values)
		{
			LS_STAG_PROFILE("CopyCompProd");
			double* src;
			_values->getCopy(src);
				LS_STAG_PAR_FOR(RealValuedMeshFunction, copyCompProd, size, m_componentwiseMultiplier, src, Vector)
		};
		void componentwiseProductR(RealValuedMeshFunction* _componentwiseMultiplier)
		{
			double* src;
			_componentwiseMultiplier->getCopy(src);
			componentwiseProduct(src);
		}
		void subtract(const double* values, const double scale)
		{
			LS_STAG_PROFILE("substractScale");
			LS_STAG_PAR_FOR(RealValuedMeshFunction, subtract, size, values, Vector, scale)
		}
		void setSumSqrt(const double* values, const double* values2, const double* values3)
		{
			LS_STAG_PROFILE("setSumSqrt");
				LS_STAG_PAR_FOR(RealValuedMeshFunction, setSumSqrt, size, Vector, values, values2, values3)
		};
		void addSqrt(const double* values, const double* values2)
		{
			LS_STAG_PROFILE("addSqrt");
				LS_STAG_PAR_FOR(RealValuedMeshFunction, addSqrt, size, Vector, values, values2)
		};
		void subtract(const double* values, const double* values2, const double* values3)
		{
			LS_STAG_PROFILE("substract");
			#ifdef LS_STAG_USE_MKL
				vdSub(size, Vector, values, Vector);
				vdSub(size, Vector, values2, Vector);
				vdSub(size, Vector, values3, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, subtract, size, values, values2, values3, Vector)
			#endif
		};
		void addAndSubtract(const double* addValues, const double* subtractValues)
		{
			LS_STAG_PROFILE("subtractAdd");
			#ifdef LS_STAG_USE_MKL
				vdAdd(size, Vector, addValues, Vector);
				vdSub(size, Vector, subtractValues, Vector);
			#else
				LS_STAG_PAR_FOR(RealValuedMeshFunction, subtractAdd, size, addValues, subtractValues, Vector)
			#endif
		}
		#ifdef LS_STAG_USE_MKL
			void setSum(const double* values, const double* values2) { vdAdd(size, values, values2, Vector); };
			void setDiff(const double* values, const double* values2) { vdSub(size, values, values2, Vector); };
			void setDiffWithScale(const double* values, double* values2, const double scale)
			{
				cblas_dscal(size, scale, values2, 1);
				vdSub(size, values, values2, Vector);
			};
			void setSquare() { vdSqr(size, Vector, Vector); };
			double getNorm() { return cblas_dnrm2(size, Vector, 1); };
			double getSum() { return cblas_dasum(size, Vector, 1); };
			void copy(const double* values) { cblas_dcopy(size, values, 1, Vector, 1); };
			void operator*=(const double factor) { cblas_dscal(size, factor, Vector, 1); };
			void operator+=(const double* values) { vdAdd(size, Vector, values, Vector); };
			void operator-=(const double* values) { vdSub(size, Vector, values, Vector); };
			double operator*(const double* Vector2) { return cblas_ddot(size, Vector, 1, Vector2, 1); }
		#endif
			void setMatr(SparseMatrix *_A)
			{
				m_A = _A;
				setMatrComp();
			}
	  private:
		void setMatrComp()
		{
			m_Cell = m_A->Cell;
			m_Portrait = m_A->Portrait;
			m_Num = m_A->Number;
		};
		SparseMatrix *m_A;///< Указатель на разреженную матрицу, с которой для данной сеточной функции производится матрично-векторное умножение
		double *m_Cell;// test 08.12.14
		int *m_Portrait, *m_Num;// test 08.12.14
		double *multiplier, *m_componentwiseMultiplier;
	};
}

#endif //_LS_STAG_REAL_VALUED_MESH_FUNCTION_H_
