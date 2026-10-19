/*
* \copyright (C) 2025 Puzikova V.V.
* \file      array_interface.h
* \author    Puzikova Valeria (vvp)
* \date      04.02.2025 17:33:15
* \brief     array_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ARRAY_INTERFACE_H_
#define _LS_STAG_ARRAY_INTERFACE_H_

// ================================================================== INCLUDES
#include <string.h>
#include <fstream>
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием абстрактного класса ArrayInterface

Данный файл содержит в себе реализацию абстрактного класса ArrayInterface
*/
// ===========================================================================

namespace LS_STAG
{
	/*!
	\brief Интерфейс для массивов

	Данный абстрактный класс реализует интерфейс массивов  
	*/
	template<typename ElemType>// double/int/coord
	struct ArrayInterface
	{
		ArrayInterface() :
			size(0),
			capasity(0)
		{
			;
		}
		virtual ~ArrayInterface() { ; }
		virtual ElemType getElem(const int index) = 0;
		virtual void setElem(const int index, const ElemType value) = 0;
		virtual void setCopy(ElemType*& cpuCopy) = 0;
		virtual void sync(const ElemType* cpuCopy) = 0;
		virtual void print(const std::string& name) = 0;
		virtual void loadState(std::ifstream& in) = 0;
		virtual void saveState(std::ofstream& out) = 0;
	protected:
		virtual void alloc() = 0;
	public:
		void resize(const int _size)
		{
			size = _size;
			capasity = size * sizeof(ElemType);
		};
		void init(const int _size)
		{
			if (!size)
			{
				resize(_size);
				alloc();
			}
		}; 
		void getCopy(ElemType*& m_copy) const { m_copy = Vector; }
		void init(const int _size, ElemType*& m_copy)
		{
			init(_size);
			getCopy(m_copy);
		};
	protected:
		ElemType* Vector;///< Массив, в котором хранится значения элементов сеточной функции
		int size;///< Размерность пространства, которому принадлежит сеточная функция
		size_t capasity;
	};
}

#endif //_LS_STAG_ARRAY_INTERFACE_H_
