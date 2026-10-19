/*
* \copyright (C) 2014 Puzikova V.V.
* \file      timer.h
* \author    Puzikova Valeria (vvp) 
* \date      06.01.2014 16:56:22
* \brief     timer.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TIMER_H_
#define _LS_STAG_TIMER_H_

// ================================================================== INCLUDES
#include <ctime>
#include "tstring.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса Timer

Данный файл содержит в себе реализацию класса Timer
*/
// ===========================================================================


namespace LS_STAG
{
/*!
\brief Класс для работы с датами и временем

Данный класс реализует конвертирование времени из типов time_t и tm в описанный в данном классе строковый аналог tm,
а также реализует некоторые операции над веменными метками.
*/
	struct Timer
	{
		struct tmString
		{
			static const int shiftMonth = 1;
			static const int shiftYear = 1900;
			void copy(const tm * _timeinfo)
			{
				convertT(_timeinfo->tm_sec, tm_secT);
				convertT(_timeinfo->tm_min, tm_minT);
				convertT(_timeinfo->tm_hour, tm_hourT);
				convertT(_timeinfo->tm_mday, tm_mdayT);
				convertT(_timeinfo->tm_mon + shiftMonth, tm_monT);
				tm_yearT = strConv.toTString(_timeinfo->tm_year + shiftYear);
				labelT = tm_yearT + tm_monT + tm_mdayT + _T("-") + tm_hourT + tm_minT + tm_secT;
				dateT = tm_mdayT + _T(".") + tm_monT + _T(".") + tm_yearT + _T(", ") + tm_hourT + _T(":") + tm_minT + _T(":") + tm_secT;
				convert(_timeinfo->tm_sec, tm_sec);
				convert(_timeinfo->tm_min, tm_min);
				convert(_timeinfo->tm_hour, tm_hour);
				convert(_timeinfo->tm_mday, tm_mday);
				convert(_timeinfo->tm_mon + shiftMonth, tm_mon);
				tm_year = strConv.toString(_timeinfo->tm_year + shiftYear);
				label = tm_year + tm_mon + tm_mday + "-" + tm_hour + tm_min + tm_sec;
				date = tm_mday + "." + tm_mon + "." + tm_year + ", " + tm_hour + ":" + tm_min + ":" + tm_sec;
			};
			void copy(const time_t timeSec)
			{
				int mdiff = (int)(timeSec / 60.0);
				int hh = (int)(mdiff / 60.0);
				int mm = mdiff - hh * 60;
				int ss = (int)timeSec - mdiff * 60;
				convertT(ss, tm_secT);
				convertT(mm, tm_minT);
				convertT(hh, tm_hourT);
				tm_mdayT = _T("");
				tm_monT = _T("");
				tm_yearT = _T("");
				labelT = strConv.toTString(timeSec);
				dateT = tm_hourT + _T(":") + tm_minT + _T(":") + tm_secT;
				convert(ss, tm_sec);
				convert(mm, tm_min);
				convert(hh, tm_hour);
				tm_mday = "";
				tm_mon = "";
				tm_year = "";
				label = strConv.toString(timeSec);
				date = tm_hour + ":" + tm_min + ":" + tm_sec;
			};
		  private:
			void convertT(const int _value, std::tstring &value)
			{
				if(_value < 10)
					value = _T("0") + strConv.toTString(_value);
				else
					value = strConv.toTString(_value);
			};
			void convert(const int _value, std::string &value)
			{
				if(_value < 10)
					value = "0" + strConv.toString(_value);
				else
					value = strConv.toString(_value);
			};
		  public:
			std::tstring labelT, dateT;
			std::tstring tm_secT;     ///< Cекунды после минут [0,59]
			std::tstring tm_minT;     ///< Минуты после часов [0,59]
			std::tstring tm_hourT;    ///< Часы после полуночи [0,23]
			std::tstring tm_mdayT;    ///< День месяца [1,31]
			std::tstring tm_monT;     ///< Месяц года (январь = 0) [0,11]
			std::tstring tm_yearT;    ///< День года (1 января = 0) [0,365]
			std::string label, date;
			std::string tm_sec;     /* seconds after the minute - [0,59] */
			std::string tm_min;     /* minutes after the hour - [0,59] */
			std::string tm_hour;    /* hours since midnight - [0,23] */
			std::string tm_mday;    /* day of the month - [1,31] */
			std::string tm_mon;     /* months since January - [0,11] */
			std::string tm_year;    /* years since 1900 */
			//    остальные поля пока не нужны
		  private:
			StringConvertor strConv;
		};
		void setTime()
		{
			time( &rawtime );                              
			timeinfo = localtime( &rawtime );
			timeInfoString.copy(timeinfo);
		};
		/*!
		Сохраняет разницу между временем rawtime и time2 в разных форматах
		*/
		void operator-=(const time_t time2)
		{
			rawtime -= time2;
			timeInfoString.copy(rawtime);
		};
		tmString timeInfoString;///< Конвертер времени из time_t и tm в строковый аналог tm 
		time_t rawtime;///< Время в секундах начиная с 1 января 1970 г.
		tm * timeinfo;///< Указатель на структуру, содержащую компоненты календарного времени
	};
}

#endif //_LS_STAG_TIMER_H_
