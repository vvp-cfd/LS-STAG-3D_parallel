/*
* \copyright (C) 2015 Puzikova V.V.
* \file      parallel_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      30.07.2015 14:46:21
* \brief     parallel_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_PARALLEL_TASK_DESCRIPTION_H_
#define _LS_STAG_PARALLEL_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include <iostream>
#include <fstream>
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ParallelTaskDescription

Данный файл содержит в себе реализацию класса ParallelTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Структура с параметрами запуска расчетов

Данная структура содержит параметры запуска расчетов, в т.ч. для параллельных режимов
*/
	struct ParallelTaskDescription
	{
		/// Режимы запуска расчетов
		enum TaskModes
		{
			SERIAL_TASK_MODE = 0,///< 0 --- последовательный режим работы
			SCALABILITY_TEST,    ///< 1 --- тестирование масштабируемости: 1 задача, несколько запусков с разным числом ядер
			TASKS_SERIA,         ///< 2 --- серия задач
			DDM,                 ///< 3 --- декомпозиция области
			DDM_SCALABILITY_TEST ///< 4 --- тестирование масштабируемости с декомпозицией области: 1 задача, несколько запусков с разным числом ядер

		};
		ParallelTaskDescription():// по умолчанию запускается одна задача в последовательном режиме
		tasksNumber(1),
		nodesNumberPerTask(1),
		mode(SERIAL_TASK_MODE)
		{
			;
		}
		void load(const std::string &name)
		{
			std::ifstream f;
			f.open(name.c_str());
			if(f.is_open())
			{
				int modeI;
				f >> modeI;
				mode = (TaskModes)modeI;
				f >> tasksNumber;
				f >> nodesNumberPerTask;
				f.close();
			}// если файла с параметрами паралелльного запуска нет, работаем в с 1-й задачей в последовательном режиме
		};
		TaskModes mode;        ///< режим запуска
		int tasksNumber;       ///< количество запускаемых задач --- для DDM - число подобластей (для DDM_SCALABILITY_TEST --- максимальное)
		int nodesNumberPerTask;///< количество ядер, выделяемых на решение одной задачи (для SCALABILITY_TEST --- максимальное)
	};
}

#endif //_LS_STAG_PARALLEL_TASK_DESCRIPTION_H_