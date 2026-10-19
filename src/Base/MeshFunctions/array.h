/*
* \copyright (C) 2025 Puzikova V.V.
* \file      array.h
* \author    Puzikova Valeria (vvp)
* \date      04.02.2025 19:48:33
* \brief     array.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ARRAY_H_
#define _LS_STAG_ARRAY_H_

// ================================================================== INCLUDES
#include <stdio.h>
#include "../ls_stag_parallel.h"
#include "array_interface.h"
#include "../profiler.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием классов ArrayCPU, ArrayGPU (для LS_STAG_USE_CUDA)

Данный файл содержит в себе реализации классов ArrayCPU, ArrayGPU, определение Array
*/
// ===========================================================================

namespace LS_STAG
{
	/*!
	\brief Класс, реализующий базовую работу с массивами на CPU

	Данный класс реализует работу с массивами на CPU: в нем происходит вся работа с памятью (выделение, удаление), 
	загрузка из файла, сохранение в файл, печать в текстовый файл, и т.д. 
	*/
	template<typename ElemType>// double/int/coord
	struct ArrayCPU : ArrayInterface<ElemType>
	{
		~ArrayCPU()
		{
			if (size > 0)
			{
				delete[] Vector;
				size = 0;
			}
		}
		ElemType getElem(const int index) { return Vector[index]; }
		void setElem(const int index, const ElemType value) { Vector[index] = value; }
		void setCopy(ElemType*& cpuCopy) { getCopy(cpuCopy); } // на CPU память заново не выделяем
		void sync(const ElemType* cpuCopy) { ; } // на CPU копия совпадает с исходным вектором, значения те же
		void print(const std::string& name)
		{
			std::ofstream f;
			std::string fileName = name + ".dat";
			f.open(fileName.c_str());
			for (int i = 0; i < size; i++)
				f << Vector[i] << std::endl;
			f.close();
		}
		void loadState(std::ifstream& in) { in.read((char*)(Vector), capasity); }
		void saveState(std::ofstream& out) { out.write((char*)(Vector), capasity); }
	private:
		void alloc() { Vector = new ElemType[size]; }
	};
#if LS_STAG_USE_CUDA
	/*!
	\brief Класс, реализующий базовую работу с массивами на GPU

	Данный класс реализует работу с массивами на GPU: в нем происходит вся работа с памятью (выделение, удаление),
	загрузка из файла, сохранение в файл, печать в текстовый файл, и т.д.
	*/
	template<typename ElemType>// double/int/coord
	struct ArrayGPU : ArrayInterface<ElemType>
	{
		ArrayGPU() :
			copy(NULL)
		{
			cudaMallocHost(&pRes, sizeof(ElemType)); // Pinned-память
		}
		~ArrayGPU()
		{
			cudaFreeHost(pRes);
			if (size > 0)
			{
				cudaFree(Vector);
				size = 0;
				copy = NULL;
			}
		}
		ElemType getElem(const int index) 
		{ 
			LS_STAG_PROFILE("getElemFromGPU");
			cudaMemcpy(pRes, Vector + index, sizeof(ElemType), cudaMemcpyDeviceToHost); // из GPU копируем на CPU
			return *pRes;
		}
		void setElem(const int index, const ElemType value) 
		{
			LS_STAG_PROFILE("setElemToGPU");
			cudaMemcpy(Vector + index, &value, sizeof(ElemType), cudaMemcpyHostToDevice);
		}
		void setCopy(ElemType*& cpuCopy) 
		{
			LS_STAG_PROFILE("setCopyFromGPU");
			createCopy();
			cudaMemcpy(copy, Vector, capasity, cudaMemcpyDeviceToHost); // из GPU копируем на CPU
			cpuCopy = copy;
		} 
		void sync(const ElemType* cpuCopy) 
		{ // из CPU копируем на GPU
			LS_STAG_PROFILE("syncToGPU");
			cudaMemcpy(Vector, cpuCopy, capasity, cudaMemcpyHostToDevice);
		} 
		void print(const std::string& name)
		{
			setCopy(copy);
			m_copy.print(name);
		}
		void loadState(std::ifstream& in) 
		{ 
			createCopy();
			m_copy.loadState(in);
			sync(copy);
		}
		void saveState(std::ofstream& out) 
		{ 
			setCopy(copy);
			m_copy.saveState(out);
		}
	private:
		void alloc() { cudaMalloc(&Vector, capasity); }
		void createCopy() { if(!copy) m_copy.init(size, copy); }
		ArrayCPU<ElemType> m_copy;
		ElemType* copy;
		ElemType* pRes;
	};

	template <typename ElemType>
	using Array = ArrayGPU<ElemType>;
#else
	template <typename ElemType>
	using Array = ArrayCPU<ElemType>;
#endif
}

#endif //_LS_STAG_ARRAY_H_
