/*
* \copyright (C) 2016 Puzikova V.V.
* \file      sol_tecplot_printer.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2016 15:50:05
* \brief     sol_tecplot_printer.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SOL_TECPLOT_PRINTER_H_
#define _LS_STAG_SOL_TECPLOT_PRINTER_H_

// ================================================================== INCLUDES
#include "base_tecplot_printer.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SolTecplotPrinter

Данный файл содержит в себе реализацию класса SolTecplotPrinter
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий печать решений с LS-STAG-сетки в Tecplot Data Format

Данный класс реализует работу по выводу рассчитанных на LS-STAG-сетке полей неизвестных в Tecplot Data Format.
*/
	struct SolTecplotPrinter : BaseTecplotPrinter
	{
	  private:
		void defineMeshStr() override
    {
        auxStr  = strConv.toString(m_N + 2);
        auxStr2 = strConv.toString(m_M + 2);
        auxStr3 = strConv.toString(m_K + 2);
        filename = m_task->outputDirectoryFull + m_task->label + _T("Sol");
        Node = &(m_task->m_mesh->m_nodes);
        m_pUW = &(m_task->m_U_W);
        existOutFlow = false;
        if (m_task->id == BaseTaskDescription::TaskKinds::TASK_FLOW)
            existOutFlow = true;
    }
	inline void setPrintBaseMeshVars()
	{
			t_P = m_p[num];
	}
		inline void PrintNode(std::ofstream &out)
		{
			t_U = sqrt(t_Ux*t_Ux + t_Uy*t_Uy + t_Uz*t_Uz);
        out << t_x << " " << t_y << " " << t_z << " "
            << t_P  << " " << t_Ux << " " << t_Uy << " " << t_Uz << " " << t_U;
			out << std::endl;
		};
		inline void PrintWall(std::ofstream &out, double V, const double U, const double W)
		{ 
			t_x = t_xStart;
			setPrintBaseMeshVars();
			t_Ux = U;
			t_Uy = V;
			t_Uz = W;
			PrintNode(out);
			for(int i = 0; i < m_N; i++, num++, curD++)
			{
				t_x = curD->x;
				setPrintBaseMeshVars();
				PrintNode(out);
			}
			t_x += 0.5 * m_delta[num-1].x;
			PrintNode(out);
		};
		inline void PrintWallPlane(std::ofstream &out, double V, const double U, const double W)
		{
			t_y = t_yStart;
			PrintWall(out, V, U, W); 
			num -= m_N; curD -= m_N;
			for(int ii = 0; ii < m_M; ii++)
			{
				t_x = t_xStart;
				t_y = curD->y;
				PrintWall(out, V, U, W);
		    } 
			num -= m_N; curD -= m_N;
			t_y += 0.5 * m_delta[num].y;
			PrintWall(out, V, U, W);
		};
		inline void PrintWest(std::ofstream &out, const double V, const double W, const double W2)
		{
			t_y = curD->y; t_x = t_xStart;
			setPrintBaseMeshVars();
			t_Ux = (m_pU_W == NULL) ? 0.0 : m_pU_W[numW];
			t_Uy = m_task->V_W;
			t_Uz = m_task->W_W;
			PrintNode(out);
			t_x = curD->x;
			t_Ux = (m_pU_W == NULL) ? 0.0 : 0.5 * (m_u[num2] + m_pU_W[numW]);
			t_Uy = 0.5 * (m_v[num3] + V);
			t_Uz = 0.5 * (W2 + W);
			PrintNode(out);
			num++; numW++; num3++; curD++; num2++;
		};
		inline void PrintEast(std::ofstream &out, const double V, const double W, const double W2)
		{
			t_x = curD->x;
			setPrintBaseMeshVars();
			t_Ux = (existOutFlow) ? m_u[num2-1] : (0.5 * (m_u[num2-1] + m_task->U_E));
			t_Uy = 0.5 * (m_v[num3] + V);
			t_Uz = 0.5 * (W2 + W);
			PrintNode(out);
			t_x += 0.5 * m_delta[num].x;
			if(!existOutFlow)
			{
				t_Ux = m_task->U_E;
				t_Uy = m_task->V_E;
				t_Uz = m_task->W_E;
			};
			PrintNode(out);
			num++; num3++;
			curD++;
		};
		inline void PrintWest(std::ofstream &out, const double V, const double W)
		{
			PrintWest(out, V, W, m_w[num]);
		};
		inline void PrintEast(std::ofstream &out, const double V, const double W)
		{
			PrintEast(out, V, W, m_w[num]);
		};
		inline void PrintEdge(std::ofstream &out, const double V, const double W, const double W2)
		{
			PrintWest(out, V, W, W2);
			for(int i = 1; i < m_xN; i++, num++, curD++, num2++, num3++)
			{
				t_x = curD->x;
			setPrintBaseMeshVars();
				t_Ux = 0.5 * (m_u[num2] + m_u[num2-1]);
				t_Uy = 0.5 * (m_v[num3] + V);
				t_Uz = 0.5 * (W2 + W);
				PrintNode(out);
			}
			PrintEast(out, V, W, W2);
		};
		inline void PrintEdge(std::ofstream &out, const double V, const double W)
		{
			PrintEdge(out, V, W, m_w[num]);
		};
		inline void printWallLayerB(std::ofstream &out, const bool asB)
		{// вывод расчётного слоя k=0; для flow B-плоскость = его точная копия (только z=0)
			k = 0; num = 0; num2 = 0, num3 = 0, numW = 0;
			curD = m_Node; t_y = t_yStart; t_z = asB ? t_zStart : curD->z;
			PrintWall(out, m_task->V_S, m_task->U_S, m_task->W_S); 
			num = 0; curD = m_Node; j = 0;
			PrintEdge(out, m_task->V_S, m_task->W_B);//1-й ряд
			for(j = 1; j < m_yM; j++)
			{
				PrintWest(out, m_v[num3-m_N], m_task->W_B);
				num2 = num - j;        // k=0: correct x‑mesh stride (N‑1)
				num3 = num;            // k=0: correct y‑mesh stride (= base for k=0)
				for(int i = 1; i < m_xN; i++, num++, num2++, num3++, curD++)
				{
					t_x = curD->x;
			        setPrintBaseMeshVars();
					switch(m_type[num])
					{
						case SolidCell:
							t_Ux = 0.0;
							t_Uy = 0.0;
							t_Uz = 0.0;
							break;
						case FluidCell:
							t_Ux = 0.5 * (m_u[num2] + m_u[num2-1]);
							t_Uy = 0.5 * (m_v[num3] + m_v[num3-m_N]);
							t_Uz = 0.5 * (m_w[num] + m_task->W_B);
							break;
						default:
						{// 29.09.2026, вариант (б): простое θ-взвешенное (θa·va+θb·vb)/(θa+θb) — см. cellValue
								curT = &m_theta[num];
								t_Ux = cellValue(curT->x, m_u[num2], (curT - 1)->x, m_u[num2-1]);
								t_Uy = cellValue(curT->y, m_v[num3], (curT - m_N)->y, m_v[num3-m_N]);
								// B-плоскость: нижняя физическая стенка — значение W_B (апертура 1)
								t_Uz = cellValue(curT->z, m_w[num], curT->z, m_task->W_B);
						}
					}
					PrintNode(out);
				}
				PrintEast(out, m_v[num3-m_N], m_task->W_B); 
			}
			num3 -= m_N;
			PrintEdge(out, m_task->V_N, m_task->W_B);//последний ряд
			//верхняя стенка 
			num -= m_N; curD -= m_N;//21.01.14
			t_x = t_xStart;
			t_y += 0.5 * m_delta[num].y;
			PrintWall(out, m_task->V_N, m_task->U_N, m_task->W_N);
		};
		inline void printWallLayerT(std::ofstream &out, const bool asT)
		{// вывод расчётного слоя k=K-1; для flow T-плоскость = его точная копия (только z на границе)
			if (asT)
			{
				const int nT = m_NMK - m_NM;
				t_y = t_yStart;  t_z = m_Node[nT].z + 0.5 * m_delta[nT].z;
			}
			else { t_y = t_yStart;  t_z = curD->z; }
			PrintWall(out, m_task->V_S, m_task->U_S, m_task->W_S); 
			num -= m_N; curD -= m_N; j = 0;
			PrintEdge(out, m_task->V_S, m_w[num - m_NM], m_task->W_T);//1-й ряд
			for(j = 1; j < m_yM; j++)
			{
				PrintWest(out, m_v[num3-m_N], m_w[num - m_NM], m_task->W_T);
				num2 = num - (j + m_zK * m_M);
				num3 = num - m_zK * m_N;
				for(int i = 1; i < m_xN; i++, num++, num2++, num3++, curD++)
				{
					t_x = curD->x;
			        setPrintBaseMeshVars();
					switch(m_type[num])
					{
						case SolidCell:
							t_Ux = 0.0;
							t_Uy = 0.0;
							t_Uz = 0.0;
							break;
						case FluidCell:
							t_Ux = 0.5 * (m_u[num2] + m_u[num2-1]);
							t_Uy = 0.5 * (m_v[num3] + m_v[num3-m_N]);
							t_Uz = 0.5 * (m_w[num-m_NM] + m_task->W_T);
							break;
						default:
						{// 29.09.2026, вариант (б): простое θ-взвешенное (θa·va+θb·vb)/(θa+θb) — см. cellValue
								curT = &m_theta[num];
								t_Ux = cellValue(curT->x, m_u[num2], (curT - 1)->x, m_u[num2-1]);
								t_Uy = cellValue(curT->y, m_v[num3], (curT - m_N)->y, m_v[num3-m_N]);
								// T-плоскость: верхняя физическая стенка — значение W_T (апертура 1)
								t_Uz = cellValue(curT->z, m_w[num-m_NM], curT->z, m_task->W_T);
						}
					}
					PrintNode(out);
				}
				PrintEast(out, m_v[num3-m_N], m_w[num - m_NM], m_task->W_T); 
			}
			num3 -= m_N;
			PrintEdge(out, m_task->V_N, m_w[num - m_NM], m_task->W_T);//последний ряд
			//верхняя стенка 
			num -= m_N; curD -= m_N;//21.01.14
			t_x = t_xStart;
			t_y += 0.5 * m_delta[num].y;
			PrintWall(out, m_task->V_N, m_task->U_N, m_task->W_N);
		};
		/*!
		Выводит в Tecplot Data Format найденные на LS-STAG-сетке поля скоростей, давления и напряжений
		*/
		void printField()
		{
			std::ofstream f;
			if(m_task->pU_W) m_pUW->setCopy(m_pU_W);
			u->setCopy(m_u);
			v->setCopy(m_v);
			w->setCopy(m_w);
			p->setCopy(m_p);
			Node->setCopy(m_Node);
			delta->setCopy(m_delta);
			theta->setCopy(m_theta);
			type->setCopy(m_type);
			f.open(fullFileName.c_str());
			f << "TITLE     = \"Velocity Vectors\" " << std::endl;
			f << "VARIABLES = \"x\" " << std::endl;
			f << "\"y\"" << std::endl;
			f << "\"z\"" << std::endl;
			f << "\"P\"" << std::endl;
			f << "\"SPD_x\"" << std::endl;
			f << "\"SPD_y\"" << std::endl;
			f << "\"SPD_z\"" << std::endl;
			f << "\"SPD\"" << std::endl;
			f << "ZONE T=\"Data\" " << std::endl;       
			f << title << std::endl;
			//нижняя стенка
			t_xStart = 0.0, t_yStart = 0.0, t_zStart = 0.0;
			t_x = 0.0;
			t_y = 0.0;
			t_z = 0.0;
			num = 0; curD = m_Node; t_z = t_zStart; 
			if (existOutFlow) printWallLayerB(f, true);
			else PrintWallPlane(f, m_task->V_B, m_task->U_B, m_task->W_B);
			printWallLayerB(f, false);
			for(k = 1; k < m_zK; k++)
			{
				t_y = t_yStart;
				t_z = curD->z;
			PrintWall(f, m_task->V_S, m_task->U_S, m_task->W_S);
			num -= m_N; curD -= m_N; j = 0;
			PrintEdge(f, m_task->V_S, m_w[num-m_NM]);//1-й ряд
			for(j = 1; j < m_yM; j++)
			{
				PrintWest(f, m_v[num3-m_N], m_w[num-m_NM]);
				num2 = num - (j + k * m_M);
				num3 = num - k * m_N;
				for(int i = 1; i < m_xN; i++, num++, num2++, num3++, curD++)
				{
					t_x = curD->x;
			        setPrintBaseMeshVars();
					switch(m_type[num])
					{
						case SolidCell:
							t_Ux = 0.0;
							t_Uy = 0.0;
							t_Uz = 0.0;
							break;
						case FluidCell:
							t_Ux = 0.5 * (m_u[num2] + m_u[num2-1]);
							t_Uy = 0.5 * (m_v[num3] + m_v[num3-m_N]);
							t_Uz = 0.5 * (m_w[num] + m_w[num-m_NM]);
							break;
						default:
							{// 29.09.2026, вариант (б): простое θ-взвешенное (θa·va+θb·vb)/(θa+θb) — см. cellValue
								curT = &m_theta[num];
								t_Ux = cellValue(curT->x, m_u[num2], (curT - 1)->x, m_u[num2-1]);
								t_Uy = cellValue(curT->y, m_v[num3], (curT - m_N)->y, m_v[num3-m_N]);
								t_Uz = cellValue(curT->z, m_w[num], (curT - m_NM)->z, m_w[num-m_NM]);
							}	   
					}//switch
					PrintNode(f);
				}
				PrintEast(f, m_v[num3-m_N], m_w[num-m_NM]); 
			}
			num3 -= m_N;
			PrintEdge(f, m_task->V_N, m_w[num-m_NM]);//последний ряд
			//верхняя стенка 
			num -= m_N; curD -= m_N;//21.01.14
			t_x = t_xStart;
			t_y += 0.5 * m_delta[num].y;
			PrintWall(f, m_task->V_N, m_task->U_N, m_task->W_N);
		}// K - 1
			const int s_num = num, s_num2 = num2, s_num3 = num3, s_numW = numW;
			coordinate* s_curD = curD;
			printWallLayerT(f, false);
			if (existOutFlow)
			{
				num = s_num; num2 = s_num2; num3 = s_num3; numW = s_numW; curD = s_curD;
				printWallLayerT(f, true);
			}
			else
			{// без existOutFlow: T-плоскость рисуется от того же состояния, что и расчётный слой (иначе нумерация уходит за массивы)
				num = s_num; num2 = s_num2; num3 = s_num3; numW = s_numW; curD = s_curD;
				const int nT = m_NMK - m_NM;
				t_z = m_Node[nT].z + 0.5 * m_delta[nT].z;
				PrintWallPlane(f, m_task->V_T, m_task->U_T, m_task->W_T);
			}
			f.close();
		};
		inline double cellValue(const double ta, const double va, const double tb, const double vb)
		{// 29.09.2026, вариант (б): простое θ-взвешенное восстановление в центре усечённой ячейки
			const double sum = ta + tb;
			return (sum > EQ_DOUBLE) ? ((ta * va + tb * vb) / sum) : 0.0;
		}
		double t_U, t_Ux, t_Uy, t_Uz, t_P, tx, txx, ty, tyy, tz, tzz; 
		int num, j, k, num2, num3, numW, t_xStart, t_yStart, t_zStart;
		bool existOutFlow;
		double* m_pU_W = NULL;
		extDouble* m_pUW;// профиль входа на западной границе (wrapper на m_task->m_U_W); host-копия --- в m_pU_W
		Array<coordinate>* Node;
		coordinate* m_Node;
	};
};

#endif //_LS_STAG_SOL_TECPLOT_PRINTER_H_
