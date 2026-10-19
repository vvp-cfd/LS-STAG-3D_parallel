/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ls_stag_real_constants.h
* \author    Puzikova Valeria (vvp) 
* \date      23.10.2013 15:36:20
* \brief     ls_stag_real_constants.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_REAL_CONSTANTS_H_
#define _LS_STAG_LS_STAG_REAL_CONSTANTS_H_

// ================================================================== INCLUDES
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с используемыми в расчетах константами

Данный файл содержит в себе определения используемых в расчетах констант, а также значений некоторых переменных,
используемых по умолчанию
*/
// ===========================================================================

namespace LS_STAG
{
	//-------------------------------------
	static constexpr double PLUS_TINY = 1e-16;
	static constexpr double EQ_DOUBLE = 1e-5;
	static constexpr double THETA_EQ_ZERO = 0.05;
	static constexpr double EQ_TINY2_DOUBLE = 1e-5;// 1e-5;//-10
	static constexpr double EQ_TINY_DOUBLE = 1e-12;
	static constexpr double EQ_TINY3_DOUBLE = 1e-13;
	static constexpr double MINUS_EQ_DOUBLE = -EQ_DOUBLE;
	static constexpr double ONE_THIRD = 1.0 / 3.0;
	static const double TWO_ROOT = 1.4142135623730951;
	static const double FILTER_THRESH = 10.0;
	//-------------------------------------
	static constexpr double PI =  3.14159265358979323846;
	static constexpr double PI_2 = PI * 0.5;
	static constexpr double PI_32 = PI * 1.5;
	static constexpr double TWO_PI = 2.0 * PI;
	static constexpr double NEAR_TWO_PI = TWO_PI - EQ_DOUBLE;
	//-------------------------------------
	// default values for task description
	//-------------------------------------
	// BaseTaskDescription
	static const double DEFAULT_EPS = 1e-5;
static const double DEFAULT_EPS_TA = 1e-8;// TA: медленный вязкий транзиент — жёстче допуск (абсолютный критерий FGMRES)
	static constexpr double DEFAULT_THETA_EPS = 1e-4;
	static constexpr double DEFAULT_THETA_ONE_EPS = 1.0 + DEFAULT_THETA_EPS;
	static const double DEFAULT_FLOW_SCALE_P = 1e-5;
	static const double DEFAULT_T_SCALE_P = 1.0;
	static const double DIMENSIONLESS_RHO = 1.0;
	static const double DEFAULT_V_BASE = 1.0;
	static const double DEFAULT_L_BASE = 1.0;
	static const double DEFAULT_H_BASE = 1.0;
	static const double DEFAULT_T = 0.01;
	static const double DEFAULT_RE = 100.0;
	static const double DEFAULT_T_END = 100.0;
	static const double DEFAULT_THETA_PS = 0.0;// 0.55
	static const double DEFAULT_OMEGA_2D = 0.8559;//OMEGA_2D = 0.7085;
	static const double DEFAULT_OMEGA_3D = 0.805; // 0.6582;// vanKan 3D: 2*w-1=(1-w/2)^3; // 0.6528;
	static const double DEFAULT_PERIOD_PRINT_SOL = 1.0;
	static const double DEFAULT_PERIOD_PRINT_ROTOR = 0.0;
	static const double DEFAULT_PERIOD_SAVE_STATE = 1.0;
	static const int FRAME_STEP_COUT_LOG = 100;
	static const int FRAME_STEP_PROFILE_REPORT = 1000;
	//-------------------------------------
	// FlowTaskDescription
	static const double RE_CRIT = 46.0;
	//-------------------------------------
	// FlowEllipticTaskDescription
	static const double DEFAULT_ELLIPTIC_TH_RATIO = 0.2;
	//-------------------------------------
	// FlowRotorTaskDescription
	static const double DEFAULT_ROTOR_OMEGA_START = 0.4;
	static const double DEFAULT_ROTOR_TIP_SPEED_RATIO = 2.0;
	static const double DEFAULT_ROTOR_INERT = 10.0;
	static const double DEFAULT_ROTOR_COEF_FRICTION = 0.0;
	static const double DEFAULT_ROTOR_TIME = 0.01;
	static const double DEFAULT_ROTOR_TIME2 = 15.0;
	static const double DEFAULT_CUT_OFF_FREQUENCY = 5.0;
	static const double DEFAULT_SUPPRESSION_LEVEL = 3.0;
	//-------------------------------------
	// FlowRotorPropellerTaskDescription
	static const double DEFAULT_BLADE_LENGTH = 0.583;// длина лопасти
	//-------------------------------------
	// FlowRotorSavoniusTaskDescription
	static const double DEFAULT_SAVONIUS_ROTOR_TIME_AUX = 15.0;
	//-------------------------------------
	// FlowCircleVEOscTaskDescription
	static const double DEFAULT_VE_OSC_AIRFOIL_MASS = 39.75;
	static const double DEFAULT_VE_OSC_COEF_DEMP = 0.731;
	//-------------------------------------
	// FlowCircleVerOscTaskDescription
	static const double DEFAULT_VER_OSC_A_MAX = 0.2;
	static const double DEFAULT_VER_OSC_N_F = 1.2;
	static const double DEFAULT_VER_OSC_S_0 = 0.201;
	static const double DEFAULT_VER_OSC_MOV_START = 10.0;
	//-------------------------------------
	// FlowCircleHorOscTaskDescription
	static const double DEFAULT_HOR_OSC_A_MAX = 0.8;
	static const double DEFAULT_HOR_OSC_SE = 0.2;
	static const double DEFAULT_HOR_OSC_MOV_START = 10.0;
	//-------------------------------------
	// FlowCircleRotaryOscTaskDescription
	static const double DEFAULT_TANEDA_N_F = 20.0;
	static const double DEFAULT_TANEDA_S_0 = 0.17;
	static const double DEFAULT_ROT_OSC_OMEGA_0 = 6.0;
	static const double DEFAULT_ROT_OSC_SE = 0.74;
	//-------------------------------------
	// FlowCircleSMovTaskDescription
	static const double DEFAULT_CIRCLES_MOV_MASS = 4.7273;
	static const double DEFAULT_CIRCLES_MOV_COEF_DEMP = 0.00033;
	//-------------------------------------
	// default values for turbulent approaches
	//-------------------------------------
	static const double DEFAULT_LES_SMAGORINSKY_CONST = 0.2;
	//-------------------------------------
	// default values for turbulent models
	//-------------------------------------
	// SA-FT
	static const double DEFAULT_CHI_BASE_FT = 3.0;
	//-------------------------------------
	// SA-TL
	static const double DEFAULT_CHI_BASE_TL = 0.001;
	static const double DEFAULT_TIME_TL = 5.0;
	static const double DEFAULT_SA_SMAGORINSKY_CONST = 0.65;
	//-------------------------------------
	// K-W
	static const double DEFAULT_KW_SMAGORINSKY_CONST = 0.78;
	//-------------------------------------
	// K-E
	static const double DEFAULT_KE_SMAGORINSKY_CONST = 0.61;
}

#endif //_LS_STAG_LS_STAG_REAL_CONSTANTS_H_
