/*
* \copyright (C) 2013 Puzikova V.V.
* \file      bezier_approximation.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 20:54:59
* \brief     bezier_approximation.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BEZIER_APPROXIMATION_H_
#define _LS_STAG_BEZIER_APPROXIMATION_H_

// ================================================================== INCLUDES
#include "math.h"
#include "base_for_basis_bezier.h"
#include "../Base/MeshFunctions/coord_mesh_function.h"
#include "Polynomial/polynomial_distance_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BezierApproximation

Данный файл содержит в себе реализацию класса BezierApproximation
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий аппроксимацию границы кривой Безье

Данный производный класс реализует построение аппроксимации погруженной границы кривой Безье. 
Расчет расстояния от точек расчетной области до данной кривой реализован в базовом классе PolynomialDistanceFunction.

 В моделировании обтекания тел с постоянным сечением вдоль оси Oz (бесконечный цилиндр, выдавливание профиля) 
 функция уровня (расстояние до границы) не зависит от координаты z. Поэтому в 3D изменений в иерархии этих классов почти нет
(координату z не учитываем).
*/
	struct BezierApproximation : PolynomialDistanceFunction
	{
		~BezierApproximation()
		{
			if(length > 0) delete[] coord;
		};
	  private:
		  /*!
		  Строит аппроксимацию кпогруженной границы кривой Безье
		  */
		void build()
		{
			for(int j = 0; j < length; j++)
			{
				curve.pCoef[j] = 0.0;
				for(int k = 0; k < length; k++)
					curve.pCoef[j] += coord[k].y * abs(m_base->BaseBezier[k][num]) * m_base->BaseBezier[j][num+k];
			}
			curve.m_coef.sync(curve.pCoef);
			PolynomialDistanceFunction::build();
		};
	  public:
		void set(BaseForBasisBezier* _base, const double _LE, const double yc)
		{
			// к этому моменту массив координат coord и его длина _length д.б. считаны !
			m_base = _base;
			for (int i = 0; i < length; i++)
				coord[i].add(_LE, yc, 0.0);
			init();
			num = m_base->LIMIT_BEZIER - length + 1;
			build();
		};
		  /*!
		  Загружает массив координат точек погруженной границы из текстового файла
		  */
		void load(std::ifstream &in)
		{
			in >> length;
			coord = new coordinate[length];
			coordinate* curCoord = coord;
			for(int j = 0; j < length; j++, curCoord++)
			{
				in >> curCoord->x;
				in >> curCoord->y;
				in >> curCoord->z;
			}
		};
		void print(std::ofstream &out)
		{
			out << length << "\n";
			coordinate* curCoord = coord;
			for(int j = 0; j < length; j++, curCoord++)
			{
				out << *curCoord << "\n";
			}
		};
		coordinate* coord;         ///< Массив координат точек погруженной границы, по которым строится аппроксимация
	  private:
		BaseForBasisBezier* m_base;///< Указатель на базисную матрицу Безье
		int num;
	};
}

#endif //_LS_STAG_BEZIER_APPROXIMATION_H_
