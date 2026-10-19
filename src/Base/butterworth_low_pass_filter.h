/*
* \copyright (C) 2013 Puzikova V.V.
* \file      butterworth_low_pass_filter.h
* \author    Puzikova Valeria (vvp) 
* \date      26.10.2013 0:28:45
* \brief     butterworth_low_pass_filter.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BUTTERWORTH_LOW_PASS_FILTER_H_
#define _LS_STAG_BUTTERWORTH_LOW_PASS_FILTER_H_

// ================================================================== INCLUDES
#include "math.h"
#include "ls_stag_real_constants.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ButterworthLowPassFilter

Данный файл содержит в себе реализацию класса ButterworthLowPassFilter
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу ФНЧ Баттерворта

Данный класс реализует построение фильтра низких частот Баттерворта и 
обработку сигнала данным фильтром в режиме реального времени
*/
	struct ButterworthLowPassFilter
	{
		static const int Nc = 3;//< Порядок фильтра + 1
		ButterworthLowPassFilter():
		thresh(FILTER_THRESH),
		eps(EQ_DOUBLE),
		res_n(0.0),
		res_n_1(0.0),
		res_filt_n(0.0),
		res_filt_n_1(0.0)
		{// < 5 Hz
			b[0] = 0.0004128228951732185;
			b[1] = 0.000825645790346437;
			b[2] = 0.0004128228951732185;
			a[0] = -1.94171217995202;
			a[1] = 0.9433634715327129;
		}
		/*!
		Расчет коэффициентов фильтра низких частот Батерворта по шагу дискретизации, частоте среза, уровню подавления
		*/
		void init(const double t, const double fs, const double Rs)
		{
			double aux2 = pow(10.0, -0.025 * Rs);
			double aux = aux2 * aux2;
			double xi = pow(1 - aux * aux, 0.25) / (aux2 * tan(t * fs * PI));
			aux2 = xi * xi;
			double aux3 = 1.0 + aux2;
			double aux4 = xi * TWO_ROOT;
			aux = 1.0 / (aux3 + aux4);
			b[0] = aux;
			b[2] = aux;
			a[1] = (aux3 - aux4) * aux;
			aux *= 2.0;
			b[1] = aux;
			a[0] = aux * (1.0 - aux2);
		}
		/*!
		Производит фильтрацию сигнала
		\param newValue новое значение сигнала
		\return отфильтрованное значение сигнала
		*/
		double getValue(const double newValue)
		{
			res = newValue;
			if((fabs(res - res_n) > thresh) && (fabs(res_n) > eps)) res = res_n;
			res_filt = b[0] * res + b[1] * res_n + b[2] * res_n_1 - (a[0] * res_filt_n + a[1] * res_filt_n_1);
			res_filt_n_1 = res_filt_n;
			res_filt_n = res_filt;
			res_n_1 = res_n;
			res_n = res;
			return res_filt;
		};
		double b[Nc], a[Nc-1];///< Коэффициенты фильтра низких частот Батерворта 
		double thresh, eps, res_n, res_n_1, res_filt_n, res_filt_n_1, res, res_filt;
	};
}

#endif //_LS_STAG_BUTTERWORTH_LOW_PASS_FILTER_H_
