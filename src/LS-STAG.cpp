// LS-STAG.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "TimeSteppingMethod/time_stepping_method_interface.h"
#include "Base/TaskDescription/parallel_task_description.h"

using namespace LS_STAG;

int main(int argc, char* argv[]) 
{
#if LS_STAG_USE_CUDA
	if (findCudaDevice(argc, (const char**)argv) < 0) 
	{
		std::cout << "CUDA device not found!\n";
		system("PAUSE");
		return -1;
	}
#endif
	setlocale(LC_ALL,"Russian");
	ParallelTaskDescription parallelParameters;
	parallelParameters.load(std::string(CONFIG_PATH) + "/ParallelTaskDescription.dat");
	switch(parallelParameters.mode)
	{
	  case ParallelTaskDescription::SERIAL_TASK_MODE:
		{
		TimeSteppingMethodInterface LSSTAG;
		LSSTAG.computation(1);
		system ("PAUSE");
		}
		break;
	  case ParallelTaskDescription::SCALABILITY_TEST:
		{
		for(int i = 1, m_end = parallelParameters.nodesNumberPerTask + 1; i < m_end; i++)
		{
			TimeSteppingMethodInterface LSSTAG;
			LSSTAG.computation(i);
		}
		}
		break;
	  case ParallelTaskDescription::TASKS_SERIA:
		{
		TimeSteppingMethodInterface LSSTAG;// пока поддерживает только 1 задачу
		LSSTAG.computation(parallelParameters.nodesNumberPerTask);
		}
		break;
	  case ParallelTaskDescription::DDM:// потом отключить LS_STAG_USE_CILK
		{// 03.04.19 --- декомпозиция области
		TimeSteppingMethodInterface LSSTAG;// пока поддерживает только 1 задачу
		LSSTAG.computation(parallelParameters.nodesNumberPerTask, parallelParameters.tasksNumber);
		system ("PAUSE");
		}
		break;
	  case ParallelTaskDescription::DDM_SCALABILITY_TEST:
	  {
		  for (int i = parallelParameters.nodesNumberPerTask; i > 0; i--)
		  {
			  TimeSteppingMethodInterface LSSTAG;
			  LSSTAG.computation(i, i);
		  }
	  }
	  break;
	}
	return 0;
}

