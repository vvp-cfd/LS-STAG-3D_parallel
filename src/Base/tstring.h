/*
* \copyright (C) 2014 Puzikova V.V.
* \file      tstring.h
* \author    Puzikova Valeria (vvp) 
* \date      07.01.2014 1:14:13
* \brief     tstring.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TSTRING_H_
#define _LS_STAG_TSTRING_H_

// ================================================================== INCLUDES
#include <string>
#include <sstream>
#if !defined(__linux__)
#include <tchar.h>
#else
#define _T(x) x
#endif
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса StringConvertor и типа tstring

Данный файл содержит в себе реализацию класса StringConvertor и определения типов tstring, tostringstream
*/
// ===========================================================================

namespace std
{
	#ifdef UNICODE
		typedef wstring tstring;
		typedef wostringstream tostringstream;
	#else
		typedef string tstring;
		typedef ostringstream tostringstream;
	#endif
}
namespace LS_STAG
{
/*!
\brief Класс, реализующий конвертирование переменных в строки

Данный класс реализует конвертирование переменных любого типа в строки типов string и tstring
*/
	struct StringConvertor
	{
		template <typename T>
		std::tstring toTString(T val)
		{
			ossT.str(std::tstring());
			ossT << val;
			return ossT.str();
		};
		template <typename T>
		std::string toString(T val)
		{
			oss.str(std::string());
			oss << val;
			return oss.str();
		};
		std::string toString(const std::tstring& val)
		{// 30.09.2026: узкая конверсия tstring (для имён файлов ArrayInterface::print)
		#if defined UNICODE && _WIN32
			if (val.empty()) return std::string();
			std::string out;
			out.resize(val.size());
			WideCharToMultiByte(1251, 0, val.c_str(), (int)val.size(), &out[0], (int)out.size(), NULL, NULL);
			return out;
		#else
			return val;
		#endif
		};
	  private:
		std::tostringstream ossT;
		std::ostringstream oss;
	};
}

#endif //_LS_STAG_TSTRING_H_
