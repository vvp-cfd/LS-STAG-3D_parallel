/*
* \copyright (C) 2013 Puzikova V.V.
* \file      set_solutions.h
* \author    Puzikova Valeria (vvp) 
* \date      31.12.2013 18:11:44
* \brief     set_solutions.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SET_SOLUTIONS_H_
#define _LS_STAG_SET_SOLUTIONS_H_

// ================================================================== INCLUDES
#include "solution.h"
#include "SparseMatrix/sparse_matrix.h"
#include "TaskDescription/Tasks/Flow/Circle/flow_circle_mov_task_description.h"
#include "TecplotPrinters/sol_tecplot_printer.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SetSolutions

Данный файл содержит в себе реализацию класса SetSolutions
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(corr, LS_STAG_R_PD, _p, LS_STAG_R_PD, _predP, LS_STAG_CR_PPI, _OnAirfoil,
		const double, _boundP, const double, _coefCorrectionPredP);
	LS_STAG_DECLARATION_R_KERNEL_2D(corr, LS_STAG_R_PD, _u, LS_STAG_R_PD, _prevU, LS_STAG_CR_PCT, _xTypes,
		const double, value, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(corr, LS_STAG_R_PD, _u, LS_STAG_R_PD, _prevU, LS_STAG_CR_PCT, _xTypes, LS_STAG_CR_PD, values, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_R_KERNEL_2D(corrOutFlow, LS_STAG_R_PD, _u, LS_STAG_R_PD, _v, LS_STAG_R_PD, _w, const double, vBase, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_R_KERNEL_2D(corrOutFlowTB, LS_STAG_R_PD, _w, const int, _N, const int, _NMK2);
	LS_STAG_DECLARATION_KERNEL_2D(copyOutflow, LS_STAG_R_PD, _u, LS_STAG_CR_PD, _prevU, const int, _N, const int, _NM);
/*!
\brief Класс, реализующий работу с наборами решений

Данный класс реализует работу с решениями со всех необходимых одномоментно временных слоев:
производит коррекцию прогнозов, сохранение состояния расчета, загрузку состояния расчета при его возобновлении,
вывод в Tecplot
*/
	struct SetSolutions
	{
		typedef LSSTAGUniformMesh::MeshKinds kind;
		typedef BaseTaskDescription::TaskKinds task;
		typedef FlowTaskDescription::AirfoilKind airfoil;
		SetSolutions():
		boundP(20.0)
		{
			;
		}
		LS_STAG_KERNEL(corr, LS_STAG_R_PD _p, LS_STAG_R_PD _predP, LS_STAG_CR_PPI _OnAirfoil,
			const double _boundP, const double _coefCorrectionPredP)
		{
			const int num = _OnAirfoil[i].z;
			if (fabs(_p[num]) > _boundP) _p[num] = 0.0;
			if ((_predP[num] > (_boundP - _p[num]) * _coefCorrectionPredP) ||
				(_predP[num] < (-_boundP - _p[num]) * _coefCorrectionPredP))// сглаживаем скользящим средним
				_predP[num] = 0.0;// 4.06.13
		}
		LS_STAG_KERNEL_2D(corr, LS_STAG_R_PD _u, LS_STAG_R_PD _prevU, LS_STAG_CR_PCT _xTypes, const double value, const int _N, const int _NM)
		{
			for (int num = i * _N + j * _NM, m_end = num + _N; num < m_end; num++)
			{
				if (_xTypes[num] != SolidCell)
				{
					_u[num] = value;
					_prevU[num] = value;
				}
			}
		}
		LS_STAG_KERNEL_2D(corr, LS_STAG_R_PD _u, LS_STAG_R_PD _prevU, LS_STAG_CR_PCT _xTypes, LS_STAG_CR_PD values, const int _N, const int _M, 
			const int _NM)
		{
			const double value = values[i + j * _M];
			for (int num = i * _N + j * _NM, m_end = num + _N; num < m_end; num++)
			{
				if (_xTypes[num] != SolidCell)
				{
					_u[num] = value;
					_prevU[num] = value;
				}
			}
		}
		LS_STAG_KERNEL_2D(corrOutFlow, LS_STAG_R_PD _u, LS_STAG_R_PD _v, LS_STAG_R_PD _w, const double vBase, const int _N, const int _M, const int _NM)
		{
			int num = (i + 1) * _N - 1 + j * _NM;
			int num2 = num - i - 1 - j * _M;
			int num3 = num - j * _N;
			if (_u[num2] < 0.0)//minUxOutFlow) 
			{
				_u[num2] = vBase;
				_v[num3] = 0.0;
				_w[num] = 0.0;
			}
			num2--;
			if (_u[num2] < 0.0)//minUxOutFlow) 
			{
				_u[num2] = vBase;
				_v[num3 - 1] = 0.0;
				_w[num - 1] = 0.0;
			}
		}
		LS_STAG_KERNEL_2D(copyOutflow, LS_STAG_R_PD _u, LS_STAG_CR_PD _prevU, const int _N, const int _NM)
		{
			const int num = (i + 1) * _N - 1 + j * _NM;
			_u[num] = _prevU[num - 1];
		}
		LS_STAG_KERNEL_2D(corrOutFlowTB, LS_STAG_R_PD _w, const int _N, const int _NMK2)
		{// анти-вток через B/T (аналог corrOutFlow для выходных плоскостей z): зажим w;
		// z-сетка: K-1 слоёв (0..K-2) --- B: слой (0), T: слой (K-2) --- симметрично
			int num = i + j * _N;
			if (_w[num] > 0.0) _w[num] = 0.0;
			num += _NMK2;
			if (_w[num] < 0.0) _w[num] = 0.0;
		}
		void set(SparseMatrix* gradX, SparseMatrix* gradY, SparseMatrix* gradZ, const double _coef) 
		{
        coef = _coef;
        current.u.setMultipliers(gradX, m_predP, m_mesh->xMesh.inv_VOF);
        current.v.setMultipliers(gradY, m_predP, m_mesh->yMesh.inv_VOF);
        current.w.setMultipliers(gradZ, m_predP, m_mesh->zMesh.inv_VOF);
    }

    void init(BaseTaskDescription* _task, int* _it) 
	{
        m_it = _it;
        m_task = _task;
        m_mesh = m_task->m_mesh;
        m_baseMesh = &(m_mesh->baseMesh);
        m_xMesh   = &(m_mesh->xMesh);
        m_yMesh   = &(m_mesh->yMesh);
        m_zMesh   = &(m_mesh->zMesh);
        m_curTime = m_task->m_curTime;

        m_N = m_baseMesh->N;   // число ячеек по x
        m_M = m_baseMesh->M;   // число ячеек по y
        m_K = m_baseMesh->K;   // число ячеек по z
        m_NMK = m_N * m_M * m_K;

        m_type  = m_baseMesh->typeOfCell;
        m_theta = m_mesh->theta;
        m_delta = m_mesh->delta;

        // Инициализация решений
        current.init(m_xMesh->NMK, m_yMesh->NMK, m_zMesh->NMK, m_NMK);
        prediction.init(m_xMesh->NMK, m_yMesh->NMK, m_zMesh->NMK, m_NMK);
        previous.init(m_xMesh->NMK, m_yMesh->NMK, m_zMesh->NMK, m_NMK);

        current.set(0.0);
        prediction.set(0.0);
        previous.set(0.0);

        current.getCopy(m_u, m_v, m_w, m_p);
        previous.getCopy(m_prevU, m_prevV, m_prevW, m_prevP);
        prediction.getCopy(m_predU, m_predV, m_predW, m_predP);

        needScalePress = false;
        startCorrectionPredP = 5.0 * m_task->T_end;
        startCorrectionPress = startCorrectionPredP;

        // Границы выходного потока: диапазоны по y и z
        vStartOutFlow = 1;    vFinishOutFlow = m_M - 1;   // M – число ячеек по y
        wStartOutFlow = 1;    wFinishOutFlow = m_K - 1;   // K – число ячеек по z
			switch(m_task->id)
			{
				case task::TASK_FLOW:
					{
						existOutFlow = true;
						existOutFlowTB = true;
						FlowTaskDescription* aux = (FlowTaskDescription*)m_task;
						if(aux->id_airfoil != airfoil::AIRFOIL_ROTOR) 
						{
							if(m_task->id_motion == BaseTaskDescription::MotionKinds::MOTION_MOV)
							{
								FlowCircleMovTaskDescription* auxMovTask = (FlowCircleMovTaskDescription*)m_task;
								startCorrectionPress = auxMovTask->t_0;
								coefCorrectionPredP = m_task->t;//07.12.14
								startCorrectionPredP = startCorrectionPress;//07.12.14
							}
							needScalePress = true;
						}
						else
						{
							FlowRotorTaskDescription* aux2 = (FlowRotorTaskDescription*)m_task;
							startCorrectionPredP = aux2->rotor_time;
							startCorrectionPress = startCorrectionPredP;
							coefCorrectionPredP = 2.0 * m_task->t / 3.0;
						}
						if(m_task->useTM) coefCorrectionPredP = 2.0 * m_task->t / 3.0;
						if(m_task->useNNM && m_task->needXYMesh) coefCorrectionPredP = m_task->t;// 24.10.16
					}
					break;
				default:
					existOutFlow = false;
			}
			if(existOutFlow) initialVelocity();
			printerSol.init(m_task, m_it, &(current.u), &(current.v), &(current.w), &(current.p), m_task->PeriodPrintSol, false);
			ItSaveState = (int)floor(m_task->PeriodSaveState / m_task->t);
			ItCounterSaver = 1;
		};
		void init()
		{
			LS_STAG_PROFILE("SolInit");
			prediction.u.copy(m_u);
			prediction.v.copy(m_v);
			prediction.w.copy(m_w);
			prediction.p.set(0.0);
		};
		/*!
		Корректирует прогнозы скоостей после расчета поправки давления
		*/
		void correction()
		{
			LS_STAG_PROFILE("SolCorrection");
			previous.copy(&current);//shift
			if((*m_curTime) > startCorrectionPredP) correctionPredP();//for rotors
			// velocities
			current.u.setProduct();
			current.v.setProduct();
			current.w.setProduct();
			current.u.componentwiseProduct();
			current.v.componentwiseProduct();
			current.w.componentwiseProduct();
			current.u += m_predU;
			current.v += m_predV;
			current.w += m_predW;
			if(existOutFlow) correctionOutFlow();// 9.02.13 // test
			if(existOutFlowTB) correctionOutFlowTB();// анти-вток через B/T
			// pressure
			if(needScalePress && (((*m_curTime) < m_task->tScaleP))) prediction.p *= m_task->flowScaleP;
			prediction.p *= coef;
			current.p += m_predP;
			if((*m_curTime) > startCorrectionPress) checkPressureOnMovMesh();
		};
		void print()
		{
			LS_STAG_PROFILE("printSol");
			printerSol.print();
			checkNeedSaveState();
		};
		/*!
		Загружает состояние счета из файла при возобновлении расчета после перерыва
		*/
		void load(const std::string &name)
		{
			std::ifstream f;
			f.open(name.c_str(), std::ios::binary);
			f.read((char*)(m_it), sizeof(int));	
			f.read((char*)(m_curTime), sizeof(double));
			for(int i = 0; i < m_task->m_mesh->ibCount; i++)
			{
				f.read((char*)(&(m_mesh->pCur_C[i].x)), sizeof(double));
				f.read((char*)(&(m_mesh->pCur_C[i].y)), sizeof(double));
            f.read((char*)(&(m_mesh->pCur_C[i].z)), sizeof(double));
            f.read((char*)(&(m_mesh->pCurV_C[i].x)), sizeof(double));
            f.read((char*)(&(m_mesh->pCurV_C[i].y)), sizeof(double));
            f.read((char*)(&(m_mesh->pCurV_C[i].z)), sizeof(double));
			}
			m_mesh->m_cur_C.sync(m_mesh->pCur_C);
			m_mesh->m_curV_C.sync(m_mesh->pCurV_C);
			m_task->loadState(f);
			current.loadState(f);
			previous.loadState(f);
			f.close();
		};
	  private:
		void checkPressureOnMovMesh()
		{
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, corr, m_NMK, m_p, m_type)
		};
		/*!
		Учитывает начальное условие 
		*/
		void initialVelocity()
		{
			int num = 0;
			bool setValue = true;
			switch(m_task->id)
			{
			case task::TASK_FLOW:
				if (m_task->discontin)
				{
					FlowTaskDescription* aux = (FlowTaskDescription*)m_task;
					if ((aux->id_airfoil != airfoil::AIRFOIL_CIRCLES) && (aux->id_airfoil != airfoil::AIRFOIL_ELLIPTIC))
					{
						if (m_task->Re > RE_CRIT) // > criticalRe => //разрывное начальное условие
							num = (int)floor(m_xMesh->M * 0.5);
					}
				}
				break;
			}
			if(setValue) LS_STAG_PAR_FOR_2D_RANGE(SetSolutions, corr, num, m_xMesh->M, 0, m_xMesh->K, m_u, m_prevU, m_xMesh->typeOfCell, m_task->V_BASE, m_xMesh->N, m_xMesh->N * m_xMesh->M)
		};
		void correctionOutFlow()
		{
			LS_STAG_PAR_FOR_2D_RANGE(SetSolutions, corrOutFlow, vStartOutFlow, vFinishOutFlow, wStartOutFlow, wFinishOutFlow, m_u, m_v, m_w, m_task->V_BASE, m_baseMesh->N, m_baseMesh->M, m_baseMesh->N * m_baseMesh->M)
		};
		void correctionOutFlowTB()
		{// анти-вток через B/T (flow): зажим w на приграничных плоскостях
			LS_STAG_PAR_FOR_2D_RANGE(SetSolutions, corrOutFlowTB, 1, m_baseMesh->N, 1, m_baseMesh->M, m_w, m_baseMesh->N, (m_baseMesh->K - 2) * m_baseMesh->N * m_baseMesh->M)
		};
		void correctionPredP()
		{
			for(int j = 0; j < m_task->ibCount; j++)
			{
				LS_STAG_PAR_FOR(SetSolutions, corr, m_task->m_loads[j].countCutCells, m_p, m_predP, m_task->m_loads[j].OnAirfoil,
					boundP, coefCorrectionPredP)
			}
		};
	  public:
		/*!
		Cохраняет состояние расчета в файл, в имени которого указан текущий момент времени
		*/
		void saveState()
		{
			filename = m_task->outputDirectoryFull + _T("continueSolveTask__it") + strConv.toTString(*m_it) + _T("time") + strConv.toTString(*m_curTime) + _T("_") + m_task->label + _T(".dat");
			std::ofstream f;
			f.open(filename.c_str(), std::ios::binary);
			f.write((char*)(m_it), sizeof(int));	
			f.write((char*)(m_curTime), sizeof(double));
			for(int i = 0; i < m_task->m_mesh->ibCount; i++)
			{
				f.write((char*)(&(m_mesh->pCur_C[i].x)), sizeof(double));
				f.write((char*)(&(m_mesh->pCur_C[i].y)), sizeof(double));
            f.write((char*)(&(m_mesh->pCur_C[i].z)), sizeof(double));
            f.write((char*)(&(m_mesh->pCurV_C[i].x)), sizeof(double));
            f.write((char*)(&(m_mesh->pCurV_C[i].y)), sizeof(double));
            f.write((char*)(&(m_mesh->pCurV_C[i].z)), sizeof(double));
			}
			m_task->saveState(f);
			current.saveState(f);
			previous.saveState(f);
			f.close();
		};
		/*!
		Проверяет, надо ли сохранять на текущем шаге состояние расчета в файл и, если надо, вызывает функцию дял сохранения
		*/
		void checkNeedSaveState()
		{
			if(ItSaveState)
			{
				if((ItCounterSaver == ItSaveState) || ((*m_it) == 1))
				{
					saveState();
					ItCounterSaver = 1;
				}
				else ItCounterSaver++;
			}
		};
		SolTecplotPrinter printerSol;   ///< Служит для печати решений с LS-STAG-сетки в Tecplot Data Format
		StringConvertor strConv;
		std::tstring filename;
	  public: 
	  Solution current, prediction, previous;

private:
    double minUxOutFlow;
    double startCorrectionPress;
    double coef, startCorrectionPredP, coefCorrectionPredP, boundP;
    BaseTaskDescription* m_task;
    LSSTAGUniformMesh* m_mesh;
    LSSTAGBaseMeshParams *m_baseMesh, *m_xMesh, *m_yMesh, *m_zMesh; // добавлен m_zMesh
    coordinate *m_theta, *m_delta;
    double *m_u, *m_v, *m_w;             // сырые указатели для текущего решения
    double *m_prevU, *m_prevV, *m_prevW; // предыдущее решение
    double *m_predU, *m_predV, *m_predW; // предсказанное решение
    double *m_predP, *m_p, *m_prevP;
    double *m_curTime;
    CellTypes* m_type;
    int *m_it;
    int num, j, num2, i;
	int m_N, m_M, m_K;            // число ячеек по осям
    int m_NMK;                    // общее число ячеек = N * M * K
    int vStartOutFlow, vFinishOutFlow;
    int wStartOutFlow, wFinishOutFlow;
    int ItSaveState, ItCounterSaver;
    bool existOutFlow, existOutFlowTB, needScalePress;
	};
}

#endif //_LS_STAG_SET_SOLUTIONS_H_
