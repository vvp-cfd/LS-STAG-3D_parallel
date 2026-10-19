/*
* \copyright (C) 2019 Puzikova V.V.
* \file      domain_decompositor.h
* \author    Puzikova Valeria (vvp) 
* \date      06.04.2019 20:26:47
* \brief     domain_decompositor.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DOMAIN_DECOMPOSITOR_H_
#define _LS_STAG_DOMAIN_DECOMPOSITOR_H_

// ================================================================== INCLUDES
#include "shur_interface_task.h"
#include "shur_subdomain_2.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DomainDecompositor

Данный файл содержит в себе реализацию класса DomainDecompositor
*/
// пока храним и полную матрицу, и матрицы подобластей
// но затем построение всех этих матриц нужно перенести непосредственно
// в построение глобальной матрицы --- это нужно хорошо обдумать....
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий декомпозицию области

Данный класс реализует построение декомпозиции расчетной области и хранит интерфейсную задачу и задачи в подобластях
*/
	struct DomainDecompositor
	{
		DomainDecompositor():
		domainCount(0)
		{
			;
		}
		~DomainDecompositor()
		{
			if(domainCount)
			{
				for(int i = 0; i < domainCount; i++)
					delete subTasks[i];
				delete [] subTasks;
				domainCount = 0;
			}
		}
		void set(SparseMatrix* _A, double* x, double* b, const int _p)
		{// + build
			domainCount = _p;
			subTasks = new ShurSubdomain * [domainCount];
			subTasks[0] = new ShurSubdomain;
			int lastNum = domainCount - 1;
			for (int i = 1; i < lastNum; i++)
				subTasks[i] = (ShurSubdomain*)(new ShurSubdomain2);
			int k = _A->m_K / _p;
			if (lastNum)
			{
				subTasks[lastNum] = new ShurSubdomain;
				interTask.set(_A, x, b, k, _p);
				for (int i = 0; i < domainCount; i++)
					subTasks[i]->set(_A, &(interTask.Ass), x, b, interTask.right, interTask.sol, i, k);
			}
			else
				subTasks[0]->set(_A, NULL, x, b, NULL, NULL, 0, k);
		}
		void rebuild()
		{
			for(int i = 0; i < domainCount; i++)
				subTasks[i]->rebuild();
			interTask.rebuild();
		}
		void print(const std::string& name)
		{
			interTask.Ass.print(name+"_ASS");
			for (int i = 0; i < domainCount; i++)
				subTasks[i]->App->print(name + "_APP_" + std::to_string(i));
		};
		int domainCount;///< Число подобластей
		ShurSubdomain **subTasks;///< Массив задач в подобластях
		ShurInterfaceTask interTask;///< Интерфейсная задача
	};
};

#endif //_LS_STAG_DOMAIN_DECOMPOSITOR_H_
