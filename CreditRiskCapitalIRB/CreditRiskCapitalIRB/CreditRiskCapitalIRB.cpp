#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <crtdbg.h>
#include "CalcDate.h"
#include "Structure.h"
#include "Util.h"

#ifndef DLLEXPORT(A) 
#ifdef WIN32
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
#include <crtdbg.h>
#elif _WIN64
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
#include <crtdbg.h>
#elif __linux__
#define DLLEXPORT(A) extern "C" A
#elif __hpux
#define DLLEXPORT(A) extern "C" A
#elif __unix__
#define DLLEXPORT(A) extern "C" A
#else
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
#endif
#endif


double Calc_Z(double CopulaCorr, double M, double Epsilon)
{
	double Z = sqrt(CopulaCorr) * M + sqrt(1 - CopulaCorr) * Epsilon;
	return Z;
}

double Calc_Tau(double Z, double PD)
{
	double Tau = -log(1 - CDF_N(Z)) / PD;
	return Tau;
}

double Sum(double* Array, long LenArray)
{
	long i;
	double s = 0.0;
	for (i = 0; i < LenArray; i++)
	{
		s += Array[i];
	}
	return s;
}

double Mean(double* Array, long LenArray)
{
	long i;
	double s = 0.0;
	for (i = 0; i < LenArray; i++)
	{
		s += Array[i] / (double)LenArray;
	}
	return s;
}

double Std(double* Array, long LenArray)
{
	long i;
	double mu = Mean(Array, LenArray);
	double var = 0.0;
	for (i = 0; i < LenArray; i++)
	{
		var += (Array[i] - mu) * (Array[i] - mu) / LenArray;
	}
	double stdev = sqrt(var);
	return stdev;
}

double WCDR(double PD, double Corr, double alpha)
{
	double WorstCaseDefaultRate = CDF_N((INV_CDF_N(PD) + sqrt(Corr) * INV_CDF_N(1 - alpha)) / (sqrt(1.0 - Corr)));
	return WorstCaseDefaultRate;
}

double EDR(
	double PD,
	double CopulaCorr
)
{
	long i;
	long NodeNumber = 10000;
	double dP = 1.0/(double)NodeNumber;
	double P;
	double* DefaultRate= (double*)malloc(sizeof(double) * NodeNumber);
	double MeanValue;

	P = 0.0;
	DefaultRate[0] = 0.0;
	for (i = 1; i < NodeNumber; i++)
	{
		P += dP;
		DefaultRate[i] = CDF_N((INV_CDF_N(PD) + sqrt(CopulaCorr) * INV_CDF_N(P)) / sqrt(1.0 - CopulaCorr));
	}

	MeanValue = Mean(DefaultRate, NodeNumber);
	free(DefaultRate);
	return MeanValue;
}

double StdevDR(
	double PD,
	double CopulaCorr
)
{
	long i;
	long NodeNumber = 10000;
	double dP = 1.0 / (double)NodeNumber;
	double P;
	double* DefaultRate = (double*)malloc(sizeof(double) * NodeNumber);
	double Stdev;

	P = 0.0;
	DefaultRate[0] = 0.0;
	for (i = 1; i < NodeNumber; i++)
	{
		P += dP;
		DefaultRate[i] = CDF_N((INV_CDF_N(PD) + sqrt(CopulaCorr) * INV_CDF_N(P)) / sqrt(1.0 - CopulaCorr));
	}

	Stdev = Std(DefaultRate, NodeNumber);
	free(DefaultRate);
	return Stdev;
}

long CreitRiskCapital_Formula(
	long Number_Loan,
	double *Loan_Port_Principal,
	double Recovery,
	double PD,
	long CopulaCorrAutoFlag,
	double CopulaCorr,
	double Alpha,
	double M,

	long MaturityFlag,
	double &T,
	long IsRetail,

	// Result Output
	double* CopulaCorrel,
	double* Expected_Loss,
	double* Credit_VaR,
	double* RiskCapital,
	double* ExpectedDefaultRate,
	double* Default_Rate_Std,
	double* WorstCaseDefaultRate,
	double* Beta,
	double* MA
)
{
	long i, j;
	double copula_corr;
	double TotalPrincipal;
	double LGD;
	double UnExpectedLoss;
	double EL;
	double* LoanPortCopy = (double*)malloc(sizeof(double) * Number_Loan);
	for (i = 0; i < Number_Loan; i++){
		LoanPortCopy[i] = Loan_Port_Principal[i];
	}

	bubble_sort(LoanPortCopy, Number_Loan);

	double* Cumm_Loss_Case = (double*)malloc(sizeof(double) * Number_Loan);
	Cumm_Loss_Case[0] = LoanPortCopy[0];
	for (i = 1; i < Number_Loan; i++){
		Cumm_Loss_Case[i] = Cumm_Loss_Case[i - 1] + LoanPortCopy[i];
	}

	LGD = 1.0 - Recovery;
	TotalPrincipal = Sum(Loan_Port_Principal, Number_Loan);

	if (CopulaCorrAutoFlag == 1) 
	{
		if (IsRetail == 0)
			copula_corr = 0.12 * (1.0 - exp(-50.0 * PD)) / (1.0 - exp(-50.0)) + 0.24 * (1.0 - (1.0 - exp(-50 * PD)) / (1.0 - exp(-50.0)));
		else if (IsRetail == 1)
			copula_corr = 0.03 * (1.0 - exp(-35.0 * PD)) / (1.0 - exp(-35.0)) + 0.16 * (1.0 - (1.0 - exp(-35.0 * PD)) / (1.0 - exp(-35.0)));
		else 
			copula_corr = 0.12 * (1.0 - exp(-50.0 * PD)) / (1.0 - exp(-50.0)) + 0.3 * (1.0 - (1.0 - exp(-50.0 * PD)) / (1.0 - exp(-50.0)));
	}
	else{
		copula_corr = CopulaCorr;
	}

	*CopulaCorrel = copula_corr;	
	*ExpectedDefaultRate = PD;

	
	*Default_Rate_Std = StdevDR(PD, copula_corr);
	*Expected_Loss = PD* TotalPrincipal * LGD;
	*Beta = (0.11852 - 0.05478 * log(PD)) * (0.11852 - 0.05478 * log(PD));
	*MA = (1.0 + (M - 2.5) * Beta[0]) / (1.0 - 1.5 * Beta[0]);
	*WorstCaseDefaultRate = WCDR(PD, copula_corr, Alpha);
	*Credit_VaR = TotalPrincipal * WorstCaseDefaultRate[0] * LGD;

	for (i = 0; i < Number_Loan-1; i++)
	{
		if (*Credit_VaR >= Cumm_Loss_Case[i] && *Credit_VaR < Cumm_Loss_Case[i + 1])
		{
			*Credit_VaR = Cumm_Loss_Case[i + 1];
			break;
		}
	}
	if (i == Number_Loan - 1)
	{
		*Credit_VaR = Cumm_Loss_Case[Number_Loan - 1];
	}

	if (MaturityFlag == 0)
		*RiskCapital = (Credit_VaR[0] - Expected_Loss[0]) * MA[0];
	else
		*RiskCapital = (Credit_VaR[0] - Expected_Loss[0]) * T;

	free(LoanPortCopy);
	free(Cumm_Loss_Case);
	return 1;
}


DLLEXPORT(long) Calc_CreditRiskCapital(
	long Number_Of_Portfolios,
	double* Port_Correlation_Reshaped,

	long* Number_Loan_Array,
	long* MaturityFlag_Array,
	double* Maturity_Array,
	double* RecoveryRate_Array,
	double* PD_Array,

	long* CopulaCorrAutoFlag_Array,
	double* Copula_Corr_Array,
	double* SignificanceLevel_Array,
	long* IsRetail,

	double *ReshapedLoanAmount,
	// Result Output,
	double* Output_Copula_Corr,
	double* Expected_Loss,
	double* Credit_VaR,
	double* RiskCapital,
	double* ExpectedDefaultRate,
	double* Default_Rate_Std,
	double* WorstCaseDefaultRate,
	double* Beta,
	double* MA,

	long Method,
	double* ResultData
)
{
	long i = 0;
	long j = 0;
	long k = 0;

	double v = 0;
	long sum_nbond;

	long Number_Loan;
	double* Loan_Port_Principal; 
	long MaturityFlag;
	double *T;
	double RecoveryRate;
	double PD;

	long CopulaCorrAutoFlag;
	double Alpha;
	double M;
	long ResultCode = 0;

	sum_nbond = 0;
	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		Loan_Port_Principal = ReshapedLoanAmount + sum_nbond;
		Number_Loan = Number_Loan_Array[i];
		MaturityFlag = MaturityFlag_Array[i];
		if (MaturityFlag == 1)
			T = &Maturity_Array[i];
		else
			T = &v;

		RecoveryRate = RecoveryRate_Array[i];
		PD = PD_Array[i];
		if (CopulaCorrAutoFlag_Array[i] == 0)
			CopulaCorrAutoFlag = 1;
		else
			CopulaCorrAutoFlag = 0;
		Alpha = SignificanceLevel_Array[i];
		if (IsRetail[i] == 1)
			M = 1.0;
		else
			M = 2.5;

		ResultCode = CreitRiskCapital_Formula(
			Number_Loan,
			Loan_Port_Principal,
			RecoveryRate,
			PD,
			CopulaCorrAutoFlag,
			Copula_Corr_Array[i],
			Alpha,
			M,
			MaturityFlag,
			*T,
			IsRetail[i],

			// Result Output
			Copula_Corr_Array + i,
			Expected_Loss + i,
			Credit_VaR + i,
			RiskCapital + i,
			ExpectedDefaultRate + i,
			Default_Rate_Std + i,
			WorstCaseDefaultRate + i,
			Beta + i,
			MA + i
		);
		

		sum_nbond += Number_Loan_Array[i];
	}
	
	// Correlation Matrix
	double** Corr;
	long shape_Corr[2] = { Number_Of_Portfolios, Number_Of_Portfolios };
	k = 0;
	Corr = make_array(Number_Of_Portfolios, Number_Of_Portfolios);
	for (i = 0 ; i < Number_Of_Portfolios; i++)
		for (j = 0; j < Number_Of_Portfolios; j++)
		{
			Corr[i][j] = Port_Correlation_Reshaped[k];
			k = k + 1;
		}

	long shape_VaR[2] = { 1,Number_Of_Portfolios };
	long shape_VaRT[2] = { Number_Of_Portfolios, 1 };
	double** VaR_1_N = make_array(1, Number_Of_Portfolios);
	double** VaR_N_1 = make_array(Number_Of_Portfolios, 1);
	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		VaR_1_N[0][i] = Credit_VaR[i];
		VaR_N_1[i][0] = Credit_VaR[i];
	}

	double** wC = Dot2dArray(VaR_1_N, shape_VaR, Corr, shape_Corr);
	double** wCw = Dot2dArray(wC, shape_VaR, VaR_N_1, shape_VaRT);

	double VaR = sqrt(wCw[0][0]);


	double EL = 0.0;
	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		EL += Expected_Loss[i];
	}

	ResultData[0] = VaR;
	ResultData[1] = EL;
	ResultData[2] = VaR - EL;

	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		Output_Copula_Corr[i] = Copula_Corr_Array[i];
	}
	
	free(VaR_1_N[0]);
	free(VaR_1_N);
	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		free(VaR_N_1[i]);
	}
	free(VaR_N_1);

	free(wC[0]);
	free(wC);
	free(wCw[0]);
	free(wCw);
	
	for (i = 0; i < Number_Of_Portfolios; i++)
	{
		free(Corr[i]);
	}
	free(Corr);
	_CrtDumpMemoryLeaks();
	return ResultCode;
}