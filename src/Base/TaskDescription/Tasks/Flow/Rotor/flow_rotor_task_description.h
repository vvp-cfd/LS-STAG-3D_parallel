/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_rotor_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 18:52:54
* \brief     flow_rotor_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ROTOR_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_ROTOR_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "../flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowRotorTaskDescription

Данный файл содержит в себе реализацию класса FlowRotorTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий работу с задачей об обтекании ротора

Данный базовый класс реализует работу с описанием решаемой задачи об обтекании ротора. Ротор может авторотировать или 
совершать автоколебания. Функции уровня роторов строятся в производных классах
*/
	struct FlowRotorTaskDescription : FlowTaskDescription
	{
		/// Тип ротора
		enum RotorKinds
		{
		   ROTOR_UNKNOWN = -1,  ///< Тип неизвестен
		   ROTOR_PLATE,         ///< Пластина
		   ROTOR_PROPELLER,     ///< Пропеллер (оперённый цилиндр)
		   ROTOR_SAVONIUS,      ///< Ротор Савониуса
		   ROTOR_DARRIEUS       ///< Ротор Дарье
		};
	  public:
		explicit FlowRotorTaskDescription(RotorKinds _id_rotor = ROTOR_UNKNOWN, int _size = sizeof(FlowRotorTaskDescription)):
		  FlowTaskDescription(AIRFOIL_ROTOR, MOTION_MOV, (LSSTAGUniformMesh*)&meshMov, _size),
		  id_rotor(_id_rotor),
		  TipSpeedRatio(DEFAULT_ROTOR_TIP_SPEED_RATIO),//default
		  omega_start(DEFAULT_ROTOR_OMEGA_START),// default
		  INERT(DEFAULT_ROTOR_INERT),
		  COEF_FRICTION(DEFAULT_ROTOR_COEF_FRICTION),// default
		  rotor_time(DEFAULT_ROTOR_TIME),//default
		  rotor_time2(DEFAULT_ROTOR_TIME2),//default
		  start(false),// всегда
		  Omega(0.0),// для первой секунды, чтоб расчёт не развалился
		  Omega_0(0.0),
		  filesExist(false),
		  withoutFinestMesh(false)// для грубосеточной отладки заменить на true
		  {
			  zone_y_3 = 0.75;// !
			  discontin = 0;// 09.12.14 existMov = true;
		  }
		  ~FlowRotorTaskDescription()
		  {
			  if(filesExist)
			  {
				  angleFile.close();
				  omegaFile.close();
			  }
		  };
		  LS_STAG_KERNEL_INLINE void checkPhiBlade(double& phiBlade)
		  {
			  if (phiBlade > NEAR_TWO_PI) phiBlade -= TWO_PI;
			  if (phiBlade < MINUS_EQ_DOUBLE) phiBlade = TWO_PI - phiBlade;
		  };
		  /*!
		  Определяет номер лопасти, ближайшей к точке, для которой вычисляем значение функции уровня
		  */
		  LS_STAG_KERNEL_INLINE double CurBlade(const coordinate& pToC, double& phi, const int numBlade, const double step2,
			  const double step, const double ANGLE, const double level, const RotorKinds id_rotor)
		  {
			  if (fabs(pToC.y) < EQ_DOUBLE)
			  {// Определяет полярный угол, соответствующий точке, для которой вычисляем значение функции уровня
				  phi = (pToC.x > MINUS_EQ_DOUBLE) ? 0.0 : PI;
			  }
			  else
			  {
				  if (fabs(pToC.x) < EQ_DOUBLE)
				  {
					  phi = (pToC.y > MINUS_EQ_DOUBLE) ? PI_2 : PI_32;
				  }
				  else
				  {
					  double phi_0 = atan(fabs(pToC.y / pToC.x));
					  phi = (pToC.x > MINUS_EQ_DOUBLE) ? ((pToC.y > MINUS_EQ_DOUBLE) ? phi_0 : (TWO_PI - phi_0))
						  : ((pToC.y > MINUS_EQ_DOUBLE) ? (PI - phi_0) : (PI + phi_0));
				  }
			  }
			  int numCurBlade = 0;
			  for (int num_aux = 0; num_aux < numBlade; num_aux++)
			  {
				  double phiBlade = ANGLE + step2 * num_aux;
				  FlowRotorTaskDescription::checkPhiBlade(phiBlade);
				  const double edge = phiBlade - level * step;
				  const double edge2 = phiBlade + (2.0 - level) * step;
				  if ((phi > edge) && (phi < (edge2 + EQ_DOUBLE)) && (edge > MINUS_EQ_DOUBLE) && (edge2 < NEAR_TWO_PI)) numCurBlade = num_aux;
				  if (((phi > (TWO_PI - edge)) || (phi < (edge2 + EQ_DOUBLE))) && (edge < MINUS_EQ_DOUBLE) && (edge2 < NEAR_TWO_PI)) numCurBlade = num_aux;
				  if (((phi > edge) || (phi < (edge2 - NEAR_TWO_PI))) && (edge > MINUS_EQ_DOUBLE) && (edge2 > TWO_PI)) numCurBlade = num_aux;
				  if ((id_rotor == ROTOR_DARRIEUS) && ((phi > (TWO_PI + edge)) || (phi < (edge2 + EQ_DOUBLE))) && (edge < MINUS_EQ_DOUBLE) && (edge2 < NEAR_TWO_PI))
					  numCurBlade = num_aux;// 12.06.13 --- проверить для остальных
			  }
			  return ANGLE + step2 * numCurBlade;
		  };
		void printLoads()
		{
			BaseTaskDescription::printLoads();
			angleFile << ANGLE << "\n";
			omegaFile << Omega << "\n";
		};
		void saveState(std::ofstream &out)
		{
			BaseTaskDescription::saveState(out);
			out.write((char*)(&ANGLE), sizeof(double));
			out.write((char*)(&Omega), sizeof(double));
		};
		void loadState(std::ifstream &in)
		{
			BaseTaskDescription::loadState(in);
			in.read((char*)(&ANGLE), sizeof(double));
			in.read((char*)(&Omega), sizeof(double));
		};
		virtual void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> ANGLE;
			in >> Omega;
			in >> COEF_FRICTION;
			in >> INERT;
			in >> TipSpeedRatio;
			in >> omega_start;
			in >> rotor_time;
			in >> rotor_time2;
			loadNumBlade(in);
		};
		void checkExistMov()// 1.06.14
		{
			if(!existMov)
			{
				if((*m_curTime) > rotor_time) existMov = true;
			}
		};
		void buildVib()
		{
			meshMov.rebuildAngleVelocity(Omega, P_C);
		};
	  protected:
		int getRotorKind(){ return (int)id_rotor; };
		void loadDiscontin(std::ifstream &in){  };
		void printDiscontin(std::ofstream &out){  };
		virtual void loadNumBlade(std::ifstream &in){ in >> numBlade; };
		virtual void printNumBlade(std::ofstream &out){ out << numBlade << "\n"; };
		void calcNodesCount()
		{
			FlowTaskDescription::calcNodesCount();
			if(!withoutFinestMesh)
			{// сетка измельчается (initFinestMesh) --- компенсируем расширение зон
				meshMov.baseMesh.N += (int)floor(9.5 * nn);//номер последнего узла по горизонтали
				meshMov.baseMesh.M += (int)floor(9.0 * nn); //номер последнего узла по вертикали 
			}
		};
		void calcPosition()
		{
			if((*m_curTime) > rotor_time) // if --- 28.01.13
			{
				if(!start)// 2.06.13 --- добавим начальную скорость для всех роторов
				{// для пропеллеров начальный угол 0, скорость != 0
					start = true;
					if(Omega < omega_start)// if continueSolveTask
					{
						Omega_0 = omega_start;
						Omega = omega_start;  
					}
				}
				calcOmega();
			}
		};
		virtual void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << ANGLE << "\n";
			out << Omega << "\n";
			out << COEF_FRICTION << "\n";
			out << INERT << "\n";
			out << TipSpeedRatio << "\n";
			out << omega_start << "\n";
			out << rotor_time << "\n";
			out << rotor_time2 << "\n";
			printNumBlade(out);
		};
		virtual void initMeshParams()
		{
			initRotorParams();
			FlowTaskDescription::initMeshParams();
			initMovMesh();
			if(!withoutFinestMesh) initFinestMesh();
			angleFilename = outputDirectoryFull + _T("ANGLE") + label + _T(".dat");
			omegaFilename = outputDirectoryFull + _T("OMEGA") + label + _T(".dat");
			angleFile.open(angleFilename.c_str());
			omegaFile.open(omegaFilename.c_str());
			filesExist = true;
		};
		virtual void initRotorParams()
		{
			step = PI / numBlade;
			step2 = step + step;
			level = 1.0;
		};
		virtual void calcOmega()
		{
			Omega_0 = Omega;
			if((*m_curTime) > rotor_time2) Omega += (m_loads->moment_aer - COEF_FRICTION * Omega_0) * t / INERT; // rotor_time2 --- test 4.06.13
			if(Omega > TipSpeedRatio) Omega = TipSpeedRatio; // 2.06.13
			if(Omega < -TipSpeedRatio) Omega = -TipSpeedRatio; // 2.06.13
			ANGLE += Omega * t;
			if(ANGLE > NEAR_TWO_PI) ANGLE -= TWO_PI;
			if(ANGLE < EQ_DOUBLE) ANGLE = TWO_PI - ANGLE;
		};
		LSSTAGMeshMov meshMov;///< неравномерная LS-STAG-сетка c подвижными погруженными границами
		std::tstring angleFilename, omegaFilename;
		std::ofstream angleFile, omegaFile;
		double step, step2, level, res;
	  public:
		double COEF_FRICTION;///< коэффициент трения в оси
		double INERT;
		double TipSpeedRatio; ///< быстроходность (накладывает ограничение на max угловую скорость)
		double omega_start; ///< начальная угловая скорость 
		double Omega, Omega_0;
		double rotor_time; ///< момент времени, с которого ротор начинает вращаться
		double rotor_time2; ///< момент времени, с которого ротор начинает вращаться под действием аэродинамической нагрузки
		RotorKinds id_rotor; ///< идентификатор типа ротора
		int numBlade; ///< число лопастей
		bool withoutFinestMesh; ///< true --- обычная сетка (без измельчения initFinestMesh)
	  private:
		bool start, filesExist;
	};
}

#endif //_LS_STAG_FLOW_ROTOR_TASK_DESCRIPTION_H_
