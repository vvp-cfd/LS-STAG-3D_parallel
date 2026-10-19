/*
* \copyright (C) 2013 Puzikova V.V.
* \file      base_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      02.10.2013 0:32:50
* \brief     base_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BASE_TASK_DESCRIPTION_H_
#define _LS_STAG_BASE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif
#include "../LS_STAG_Mesh/ls_stag_base_mesh_params.h"
#include "../LS_STAG_Mesh/ls_stag_uniform_mesh.h"
#include "../loads_calculator.h"
#include "../conditions_wnstb.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BaseTaskDescription

Данный файл содержит в себе реализацию класса BaseTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий работу с описанием задачи

Данный базовый класс реализует работу с описанием решаемой задачи
*/
	struct BaseTaskDescription
	{
		static const int DEFAULT_N_PRE = 0;
		static const int DEFAULT_N_POST = 3;// 4;// 2;// 8.11.14 1;
		static const int DEFAULT_N_MIN = 10000;
		static const int DEFAULT_MAX_ITER = 200;
		static const int DEFAULT_FGMRES_RESTART_NUM = 13;//20;	
		static const int DEFAULT_FGMRES_RESTART_NUM_P_FINEST = 64;//20;	
		/// Типы решаемых задач (различаются геометрией расчетной области и группами граничных условий)
		enum TaskKinds
		{
			TASK_UNKNOWN = -1,///< Тип задачи неизвестен
			TASK_FLOW = 1     ///< Обтекание тел
		};
		/// Виды движений погруженных тел
		enum MotionKinds
		{
			MOTION_UNKNOWN = -1,///< Вид неизвестен
			MOTION_NO,        ///< Непожвижные погруженные границы
			MOTION_OSC,       ///< Вынужденные крутильные колебания
			MOTION_MOV        ///< Подвижные погруженные границы
		};
		/// Виды колебаний погруженных тел
		enum OscKinds
		{
		   OSC_UNKNOWN = -1,     ///< Вид неизвестен
		   OSC_NO,               ///< Колебаний нет
		   OSC_VERTICAL,         ///< Вертикальные колебания
		   OSC_HORIZONTAL,       ///< Горизонтальные колебания
		   OSC_VISCO_ELASTIC     ///< Колебания при закреплении вязко-упругой связью
		};		
		/// Группы моделей ненютоновских жидкостей
		enum  NonNewtonianModelGroups
		{
			NNM_UNKNOWN = 0,            ///< Группа неизвестна
			NNM_VISCO_ELATIC_FLUIDS,    ///< Вязкоупругие жидкости --- 1
			NNM_GNF,                    ///< Обобщенные ньютоновские жидкости --- 2
			NNM_TIME_DEPENDENT_FLUIDS   ///< Нереостабильные жидкости --- 3
		};
		explicit BaseTaskDescription(TaskKinds _id = TASK_UNKNOWN, MotionKinds _id_motion = MOTION_UNKNOWN, LSSTAGUniformMesh* _mesh = (LSSTAGUniformMesh*)0, int _size = 0):
		id(_id),
		id_motion(_id_motion),
		m_mesh(_mesh),
		size(_size),
		MAX_ITER(DEFAULT_MAX_ITER),
		FGMRES_RESTART_NUM(DEFAULT_FGMRES_RESTART_NUM),
		FGMRES_RESTART_NUM_P(DEFAULT_FGMRES_RESTART_NUM),
		eps(DEFAULT_EPS),
		theta_eps(DEFAULT_THETA_EPS),
		flowScaleP(DEFAULT_FLOW_SCALE_P),//default
		tScaleP(DEFAULT_T_SCALE_P),// default
		rho(DIMENSIONLESS_RHO),// уже обезразмерили
		V_BASE(DEFAULT_V_BASE),// default
		L_BASE(DEFAULT_L_BASE),// default
		H_BASE(DEFAULT_H_BASE),// default
		n_pre(DEFAULT_N_PRE),// default
		n_post(DEFAULT_N_POST),// default
		N_MIN(DEFAULT_N_MIN),//default
		NUM_SUBDOMAINS(4),
		NUM_OVERLAP(8),
		THETA_PS(DEFAULT_THETA_PS),
		t(DEFAULT_T),
		Re(DEFAULT_RE),
		T_end(DEFAULT_T_END),
		We(-1.0),
		beta(0.0),
		ibCount(0),
		continueSolveTask(0),// если continue, то обновим при загрузке параметров расчёта
		needXYMesh(0),
		useTM(0),
		useNNM(0),
		existMov(false),
		PeriodPrintSol(DEFAULT_PERIOD_PRINT_SOL), 
		PeriodPrintRotor(DEFAULT_PERIOD_PRINT_ROTOR), 
		PeriodSaveState(DEFAULT_PERIOD_SAVE_STATE)
		{
			nu = V_BASE * L_BASE / Re;
			fcoef = 2.0 / (rho * V_BASE * V_BASE * L_BASE * H_BASE);// обновим после initMeshParam (там м.б. переопределен H_BASE)
			outputDirectory = _T("");
			outputDirectoryFull = _T("");
			label = _T(""); 
			#define CLASS_NAME_AIRFOIL_NUM BaseTaskDescription
		}
		~BaseTaskDescription()
		{
			if (ibCount > 0) delete[] m_loads; 
		};
		LS_STAG_KERNEL_INLINE int getAirfoilNum(const int _num, LS_STAG_CR_PI _cellOnAirfoil) { return _cellOnAirfoil ? _cellOnAirfoil[_num] : 0; };
		LS_STAG_KERNEL_INLINE void addCutCell(const int Airfoil, LS_STAG_R_PI _cutcellCounts, LS_STAG_R_P(tripletInt*) _airfoilNums,
			const int curNum, const int jAdd, const int kAdd, const int N, const int M)
		{
			_airfoilNums[Airfoil][_cutcellCounts[Airfoil]].set(curNum - jAdd - kAdd * M, curNum - kAdd * N, curNum);
			_cutcellCounts[Airfoil]++;
		};
		virtual const int* getPCellOnAirfoill() { return NULL;  }
		/*!
		Перестраивает сетку
		*/
		virtual void rebuildMesh(){};
		virtual void recalcBC(){};// пересчет г.у. (на напряжения --- вызываем после расчета напряжений)
		void init(double* _curTime)
		{// все параметры задачи к этому моменту уже загружены
		    m_curTime = _curTime;// указатель на текущее время
			if (ibCount > 0)
			{
				m_loads = new LoadsCalculator[ibCount];
				airfoilNums.init(ibCount, m_airfoilNums);
				cutcellCounts.init(ibCount, m_cutcellCounts);
			}
			buildMesh();
			if(id_motion == MOTION_MOV) m_mesh->m_U_g.set(0.0);// 10.11.14
		};
		virtual void saveState(std::ofstream &out)
		{
			for(int i = 0; i < ibCount; i++)
			{
				out.write((char*)&(m_loads[i].F_aer_x), sizeof(double));
				out.write((char*)&(m_loads[i].F_aer_y), sizeof(double));
				out.write((char*)&(m_loads[i].F_aer_z), sizeof(double));
				out.write((char*)&(m_loads[i].moment_aer), sizeof(double));
			}
		};//8.01.14
		virtual void loadState(std::ifstream &in)
		{
			for(int i = 0; i < ibCount; i++)
			{
				in.read((char*)&(m_loads[i].F_aer_x), sizeof(double));
				in.read((char*)&(m_loads[i].F_aer_y), sizeof(double));
				in.read((char*)&(m_loads[i].F_aer_z), sizeof(double));
				in.read((char*)&(m_loads[i].moment_aer), sizeof(double));
			}
		};
		/*!
		Выводит в файлы значения аэродинамических нагрузок для каждой из погруженных границ
		*/
		virtual void printLoads()
		{
			LS_STAG_PROFILE("printLoads");
			for(int i = 0; i < ibCount; i++)
			{
				m_loads[i].calculate();
				m_loads[i].print();
			}
		};
		/*!
		Загружает из файла параметры задачи
		*/
		virtual void load(std::ifstream &in)
		{
			int output_size;// число символов в outputDirectory
			in >> output_size;
			std::string auxOutputDirectory;
			getline(in, auxOutputDirectory, '\n');//д.б. записано без пробела после output_size, после него --- '\n'
			outputDirectory.resize(output_size);
			#if defined UNICODE && _WIN32 
				MultiByteToWideChar(1251, 0, &auxOutputDirectory[0], -1,	 &outputDirectory[0], output_size);
			#else
				outputDirectory = auxOutputDirectory;
			#endif
			in >> useNNM;
			if(useNNM) needXYMesh = 1;// !! 05.11.16
			in >> Re;
			in >> t;
			in >> T_end;
			in >> continueSolveTask;
			in >> PeriodPrintSol;
			in >> PeriodPrintRotor;
			in >> PeriodSaveState;
		};
		/*!
		Выводит в файл параметры задачи
		*/
		void print(std::ofstream &out)
		{ 
			std::string auxOutputDirectory;
			std::tstring name = outputDirectoryFull + _T("TaskDescription.dat");
			out.open(name.c_str());
			out << (int)id << " " << (int)id_motion << " " << getOscKind() << " " << getAirfoilKind() << " ";
			out << getRotorKind() << "\n";
			int output_size = (int)outputDirectory.size();
			out << output_size;
			auxOutputDirectory.resize(output_size);
			#if defined UNICODE && _WIN32 
				WideCharToMultiByte(1251, 0, &outputDirectory[0], -1, &auxOutputDirectory[0], output_size, NULL,	 NULL);
			#else
				auxOutputDirectory = outputDirectory;
			#endif
			out << auxOutputDirectory << "\n";
			out << useNNM << "\n";
			out << Re << "\n";
			out << t << "\n";
			out << T_end << "\n";
			out << 1 << "\n";//out << continueSolveTask << "\n";// т.к. этот файл формируется на тот случай, если захотим продолжить этот расчёт
			out << PeriodPrintSol << "\n";
			out << PeriodPrintRotor << "\n";
			out << PeriodSaveState << "\n";
			printExtParams(out);
		};
		void initSol(Solution* cur_sol, Solution* pred_sol, Solution* prev_sol)
		{
			m_mesh->initSol(cur_sol, pred_sol, prev_sol);
			for(int i = 0; i < ibCount; i++) m_loads[i].initSol();
		};
		void initGNF(double *_Nu, double *_NuXY, double *_NuXYS, double *_NuXYN, double *_NuXYW, double *_NuXZ, 
			double *_NuXZB, double *_NuXZT, double *_NuXZW, double *_NuYZ, double *_NuYZS, double *_NuYZN, double *_NuYZB, 
			double *_NuYZT)
		{
			pNu = _Nu;
			pNuXY = _NuXY; pNuXYS = _NuXYS; pNuXYN = _NuXYN; pNuXYW = _NuXYW;
			pNuXZ = _NuXZ; pNuXZB = _NuXZB; pNuXZT = _NuXZT; pNuXZW = _NuXZW;
			pNuYZ = _NuYZ; pNuYZS = _NuYZS; pNuYZN = _NuYZN; pNuYZB = _NuYZB; pNuYZT = _NuYZT;
			for(int i = 0; i < ibCount; i++) m_loads[i].initGNF(pNu);
		};
		void initTDF()
		{
			for(int i = 0; i < ibCount; i++) m_loads[i].initTDF(&nu);
		};
		virtual void checkExistMov()
		{
			if(!existMov)
			{
				if(m_mesh->movBody()) existMov = true;
			}
		};
		/*!
		\return расстояние от точки с координатами (х,у) до стенки
		*/
		virtual double getDistance(const double x, const double y, const double z){ return 1.0; };// для моделей турублентности (пока делаем только для flow)
	  protected:
		virtual void printExtParams(std::ofstream &out){};
		virtual void calcNodesCount(){};
		/*!
		Вычисляет значения коэффициентов заполнения ребер и корректирует их при необходимости
		*/
		virtual void buildAlpha(){ m_mesh->buildAlphaBaseBox(); };
		/*!
		Вычисляет значения коэффициентов заполнения ячеек
		*/
		virtual void buildTheta(){ m_mesh->buildThetaBaseBox(); };
		/*!
		Кодирует типы ячеек и граней
		*/
		virtual void buildTypeOfCells(){ m_mesh->buildTypeOfCellBaseBox(); };
		/*!
		Вычисляет скорости движения погруженной границы
		*/
		virtual void buildVib(){};
		/*!
		Вычисляет значения функции уровня в углах ячеек основной сетки
		\param zs номер ряда по z, с которого рассматриваем ячейки
		\param zf номер ряда по z, до которого рассматриваем ячейки
		\param vs номер ряда по вертикали, с которого рассматриваем ячейки
		\param vf номер ряда по вертикали, до которого рассматриваем ячейки
		\param hs номер ряда по горизонтали, с которого рассматриваем ячейки
		\param hf номер ряда по горизонтали, до которого рассматриваем ячейки
		*/
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf){};
		virtual void initMeshParams(){};
		virtual void buildBaseMesh()
		{	
			buildAlpha();
			buildTypeOfCells();
			buildTheta();
		};
		virtual void initCenter(){};
		virtual int getOscKind(){ return -1; };
		virtual int getAirfoilKind(){ return -1; };
		virtual int getRotorKind(){ return -1; };
		virtual void recalcViscoElasticParameters() 
		{
			if(We < 0.0) { nu = V_BASE * L_BASE / Re; }// в Ta изменим в initMeshParams
		}
		virtual void initBC()
		{
			m_U_W.init(m_mesh->baseMesh.M * m_mesh->baseMesh.K, pU_W);
			m_U_W.set(U_W);
		}
		void loadsInit()
		{
			LS_STAG_PAR_FOR(MeshFunction<int>, set, ibCount, 0, m_cutcellCounts)
			for (int i = 0; i < ibCount; i++) m_loads[i].init();
		}
		void loadsResize()
		{
			for (int i = 0; i < ibCount; i++)
			{
				m_loads[i].countCutCells = cutcellCounts.getElem(i);
				m_loads[i].resize();
			}
		}
		/*!
		Строит LS-STAG-сетку
		*/
		void buildMesh()
		{
			recalcViscoElasticParameters();
			m_mesh->time_coef = 1.0 / t;//25.12.13
			if((id_motion == MOTION_MOV) || (useNNM && needXYMesh))
				COEF_INTEGR_TIME = m_mesh->time_coef;
			else
				COEF_INTEGR_TIME = 1.5 / t;
			calcNodesCount();//base.N, base.M (для некоторых можем определить сразу в load()) + P_C.y
			label = strConv.toTString(m_mesh->baseMesh.N) + _T("x") + strConv.toTString(m_mesh->baseMesh.M) + _T("x") + strConv.toTString(m_mesh->baseMesh.K) + _T("Re") + strConv.toTString(Re);
			m_mesh->existXYMesh = needXYMesh;//7.02.14
			m_mesh->init();
			initCenter();//1.06.14
			initMeshParams(); 
			fcoef = 2.0 / (rho * V_BASE * V_BASE * L_BASE * H_BASE);
			for (int i = 0; i < ibCount; i++)
			{
				m_loads[i].set(m_mesh, nu, t, fcoef, outputDirectoryFull + strConv.toTString(i), label);
				airfoilNums.setElem(i, m_loads[i].OnAirfoil);
			}
			loadsInit();
			buildBaseMesh();
			loadsResize();
			initBC();// 24/10/2016
				m_mesh->buildCV_U();
				m_mesh->buildCV_V();
				m_mesh->buildCV_W();
			if(needXYMesh)
			{
			 m_mesh->buildCV_Txy();
			 m_mesh->buildCV_Txz();
			 m_mesh->buildCV_Tyz();
			if(useTM) m_mesh->buildInvDelta_ib();
			}
			m_mesh->buildDelta_ib();
			m_mesh->buildInvBaseDelta_ib();
			m_mesh->buildMeshFunctions();
			m_mesh->rebuildPrepare();//1.06.14
		};
		
		StringConvertor strConv;
	  public:
		std::tstring outputDirectoryFull;// c новой папкой для текущего расчёта
		std::tstring outputDirectory;// куда сохраняем результаты  // load
		std::tstring label; // NxMRe
		ConditionWNSTB TAU_XX, TAU_YY, TAU_ZZ, TAU_XY, TAU_XZ, TAU_YZ;// г.у. в в/у жидкостях
	  protected:
		extDouble m_DUDY_W, m_DUDZ_W;
		Array<tripletInt*> airfoilNums;
		tripletInt** m_airfoilNums;
		Array<int> cutcellCounts;
		int* m_cutcellCounts;
	  public:
		double Re, nu, t, T_end, L_BASE, H_BASE, V_BASE, W_N, W_S, W_W, W_E, W_B, W_T, U_N, U_S, U_W, U_E, U_B, U_T, V_N, V_S, V_B, V_T, V_W, V_E, rho; // load : Re, t, T_end (пока)
		double We, nu_e, lambda, beta, lambda_r;// для вязкоупругих жидкостей
		double COEF_INTEGR_TIME;// 1.5 / t или 1/t для подвижных границ, турб-ти и реологии
		double flowScaleP;//для task = "flow" несколько первых шагов умножаем прогноз Р на этот коэффициент
		double tScaleP;// до какого времени умножаем на flowScaleP   // load
		double INIT_T, t_n;// для схем с адаптивным шагом по времени // пока вообще не используем
		double eps;    // ограничение на невязку
		double theta_eps; //ограничение на коэффициент заполнения
		double THETA_PS;// (0.55) для интерфейсных гран.условий // 19.09.13 //AS
		double PeriodPrintSol, PeriodPrintRotor, PeriodSaveState; // load
	  protected:
		double fcoef;// для расчёта нагрузок
	  public:
		LoadsCalculator* m_loads;///< массив калькуляторов расчета нагрузок на погруженных границах
		LSSTAGUniformMesh* m_mesh;///< указатель на LS-STAG-сетку
		double* m_curTime;///< указатель на переменную, к которой хранится значение текущего времени
		double *pNu, *pNuXY;// массивы вязкостей
		double *pNuXYS, *pNuXYN, *pNuXYW;// массивы вязкостей на N и S границах расчетной области
    double *pNuXZ, *pNuXZB, *pNuXZT, *pNuXZW;
    double *pNuYZ, *pNuYZS, *pNuYZN, *pNuYZB, *pNuYZT;// for GNF
		double *pDUDY_W, *pDUDZ_W, *pU_W;// на входе задан профиль скорости
		extDouble m_U_W;// на входе задан профиль скорости (owner: pU_W; public для принтера --- pr-слой принтер тянет через setCopy)
		TaskKinds id;  ///< идентификатор типа задачи
		MotionKinds id_motion;///< идентификатор типа движения
		int discontin;///< 1 --- разрывное н.у.  // load
		int adaptiveT;///< = 1 --- автовыбор шага по времени // пока 0
		int adSolver;// 1 -- Additive Shwartz     // пока не используем
		int n_pre;///< число предсглаживаний  // MG
		int n_post;///< число постсглаживаний (полагаем, что предсглаживаний нет)  //MG
		int MAX_ITER;
		int FGMRES_RESTART_NUM;
		int FGMRES_RESTART_NUM_P;///< окно рестарта FGMRES для уравнения Пуассона (finest-сетки требуют шире)
		int N_MIN;    // число элементов, с которого начинается распараллеливание   // MG
		int NUM_SUBDOMAINS; // число подобластей (для распараллеливания решателя)  // 30.08.13  //AS
		int NUM_OVERLAP;// число рядов в перекрытии подобластей (всегда + 1)// 11.09.13 //AS
		int continueSolveTask;  ///< 1 --- загружаем параметры задачи из файла "continueSolveTask" // load
		int needXYMesh;///< 1 --- требуется построение сетки для касательных напряжений
		int useTM;///< 1  --- используется модель турбулентности
		int useNNM;///< 1 --- используется модель неньютоновской жидкости
		int ibCount;///< число погруженных границ в задаче  // load
		int size;// размер структуры (для каждой задачи свой)
		bool existMov;
	};
}
#endif //_LS_STAG_BASE_TASK_DESCRIPTION_H_
