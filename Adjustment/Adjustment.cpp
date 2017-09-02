// Adjustment.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <math.h>

extern "C" int _declspec(dllexport)ParAdjComplete(unsigned int AdjType, unsigned int CMeas, unsigned int CUnk, double **A, double **P, double *f, double pffn, double X[], double V[], double pvv, double Me, double MX[], double MM[], double MV[]);








__declspec(dllexport) int __cdecl ParAdjComplete(unsigned int AdjType, unsigned int CMeas, unsigned int CUnk, double **A, double **P, double *f, double pffn, double X[], double V[], double pvv, double Me, double MX[], double MM[], double MV[])
{
//unsigned int AdjType - тип на изравнението/тежестите (без корелация/с корелация + вид на корелацията)
//unsigned int CMeas - брой на измерванията
//unsigned int CUnk - брой на неизвестните
//double **A - матрица А[CMeas][CUnk] (коефициенти)
//double **P - матрица P[CMeas][CMeas] (тежести)
//double *f - вектор f[CMeas] (свободни членове)
//double pffn - сума pff, CMeas пъти редуцирана
//double X[] - вектор X[CUnk] (неизвестни)
//double V[] - вектор V[CMeas] (поправки)
//double pvv - сума pvv
//double Me - СКГ за тежест единица
//double MX[] - MX[CUnk] СКГ на неизвестните
//double MM[] - MM[CMeas] СКГ на изравнените измервания
//double MV[] - MV[CMeas] СКГ на поправките

	int result = 0; //резултат от изпълнението на функцията
	unsigned int i, j, k; //броячи
	double **TEMP; //TEMP[CMeas][CUnk] - временна матрица
	double **QMEAS; //Q[CMeas][CMeas] на измерванията
	double **N; // N[CUnk][CUnk+1] - нормална матрица
	double **Q; //Q[CUnk][CUnk]
	double **QV; //Q[CMeas][CMeas] на поправките
	TEMP = new double *[CMeas];
	QMEAS = new double *[CMeas];
	N = new double *[CUnk];
	Q = new double *[CUnk];
	QV = new double *[CMeas];

	for (i=1;i<=CUnk;i++)
	{
		N[i] = new double [CUnk+1];
		Q[i] = new double [CUnk];
	}

	for(i=1;i<=CMeas;i++)
	{
		TEMP[i] = new double[CUnk];
		QMEAS[i] = new double[CMeas];
		QV[i] = new double[CMeas];
	}
	//TEMP[][] = A[][] * P [][]
	for(i=1;i<=CUnk;i++)
	{
		for(j=1;j<=CMeas;j++)
		{
			TEMP[i][j]=0.0;
			for(k=1;k<=CMeas;k++)
			{
				TEMP[i][j]=TEMP[i][j]+A[k][i]*P[k][j];
			}
		}
	}
	//N[][] = TEMP[][] * At [][]; total: N[][] = A[][] * P[][] * At[][]
	for(i=1;i<=CUnk;i++)
	{
		for(j=1;j<=CUnk;j++)
		{
			N[i][j]=0.0;
			for(k=1;k<=CMeas;k++)
			{
				N[i][j]=N[i][j]+TEMP[i][k]*A[k][j];
			}
		}
	}
	//F[] = TEMP[][] * f[]; total: F[] = A[][] * P [][] * f[]
	for(i=1;i<=CUnk;i++)
	{
		N[i][CUnk+1]=0.0;
		for(k=1;k<=CMeas;k++)
		{
			N[i][CUnk+1]=N[i][CUnk+1]+TEMP[i][k]*f[k];
		}
	}

	pffn=0.0;
	//pffn = P[][] * f[] * f[]
	for(i=1;i<=CMeas;i++)
	{
		for(j=1;j<CMeas;j++)
		{
			pffn=pffn+P[i][j]*f[i]*f[j];
		}
	}
	//Редуциране на нормалната матрица
	for(k=2;k<=CUnk+1;k++)
	{
		for(i=k;i<=CUnk;i++)
		{
			for(j=i;j<=CUnk+1;j++)
			{
				N[i][j]=N[i][j]-N[k-1][i]*N[k-1][j]/N[k-1][k-1];
			}
		pffn=pffn-N[k-1][CUnk+1]*N[k-1][CUnk+1]/N[k-1][k-1];
		}
	}
	//след оптимизиране ще отпадне
	for(i=1;i<=CUnk;i++)
	{
		for(j=1;j<=CUnk;j++)
		{
			if(i>j)
				N[i][j]=N[j][i];
		}
	}
	//Решаване на СНУ и намиране на неизвестните
	X[CUnk+1]=1.0;
	for(i=CUnk;i>0;i--)
	{
		X[i]=0.0;
		for(j=CUnk+1;j>i;j--)
		{
			X[i]=X[i]-(N[i][j]/N[i][i])*X[j];
		}
	}
	//V[] = A[][] * X[] + f[]
	for(i=1;i<=CMeas;i++)
	{
		V[i]=0.0;
		for(j=1;j<=CUnk;j++)
		{
			V[i]=V[i]+A[i][j]*X[j];
		}
		V[i]=V[i]+f[i];
	}

/*	for(i=1;i<=CMeas;i++)
	{
		V[i]=V[i]+f[i];
	}*/
	//pvv = P[][] * V[] * V[]
	pvv=0.0;
	for(i=1;i<=CMeas;i++)
	{
		for(j=1;j<=CMeas;j++)
		{
			pvv=pvv+P[i][j]*V[i]*V[j];
		}
	}

	Me=sqrt(pvv/(CMeas-CUnk));
	//подготовка за определяне на Q[][]
	for(i=1;i<=CUnk;i++)
	{
		for(j=1;j<CUnk;j++)
		{
			if(i==j)
				Q[i][j]=1/N[i][j];
			else
				Q[i][j]=0.0;
		}
	}
	//Определяне на Q[][]
	for(i=CUnk;i>0;i--)
	{
		for(j=i;j>0;j--)
		{
			for(k=CUnk;k>j;k--)
			{
				Q[i][j]=Q[i][j]-(N[j][k]/N[j][j])*Q[k][i];
			}
			if(i!=j)
				Q[j][i]=Q[i][j];
		}
	}
	//СКГ на неизвестните
	for(i=1;i<=CUnk;i++)
	{
		MX[i]=Me*sqrt(Q[i][i]);
	}
	//TEMP[][] = A[][] * Q[][]
	for(i=1;i<=CMeas;i++)
	{
		for(j=1;j<=CUnk;j++)
		{
			TEMP[i][j]=0.0;
			for(k=1;k<=CUnk;k++)
			{
				TEMP[i][j]=TEMP[i][j]+A[i][k]*Q[k][j];
			}
		}
	}
	//QMEAS[][] = TEMP[][] * At[][]; total: QMEAS[][] = A[][] * Q[][] * At[][]
	for(i=1;i<=CMeas;i++)
	{
		for(j=1;j<=CMeas;j++)
		{
			QMEAS[i][j]=0.0;
			for(k=1;k<=CUnk;k++)
			{
				QMEAS[i][j]=QMEAS[i][j]+TEMP[i][k]*A[j][k];
			}
		}
	}
	//СКГ на изравнените стойности на измерванията
	for(i=1;i<=CMeas;i++)
	{
		MM[i]=Me*sqrt(QMEAS[i][i]);

		
/*		if(radioKorelAdj->Checked==True)
		{
			//BCB, sdimain.cpp, rows 1800 - 1807
		}
*/
	}
	//QV (поправки)
	for(i=1;i<=CMeas;i++)
	{
		for(j=1;j<=CMeas;j++)
		{
			QV[i][j] = QMEAS[i][j] - 1/P[i][j];
		}
	}
	//СКГ на поправките
	for(i=1;i<=CMeas;i++)
	{
		MV[i] = Me*sqrt(QV[i][i]);
	}

	result = 1;


	return result;
}

