#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef UTILITY
#include "Util.h"
#endif

#include "Structure.h"

double SABRIV(double alpha, double beta, double v, double rho, double Fut, double K, double T)
{
	long i;
	double z, x;
	double F = Fut;
	double numerator, denominator;
	double FKTerm = pow(F * K, (1.0 - beta) * 0.5);
	double lnFoverK = log(F/K);
	double SABRIV;

	numerator = alpha * (1.0 + ( (1.0-beta) * (1.0-beta) * alpha * alpha / (24.0 * FKTerm) + 0.25 * (rho * beta * v * alpha)/FKTerm + (2.0 - 3.0 * rho * rho )/24.0 * v * v) *T);
	denominator = FKTerm * (1.0 + (1.0 - beta) * (1.0 - beta) / 24.0 * lnFoverK * lnFoverK + (1.0 - beta) * (1.0 - beta) * (1.0 - beta) * (1.0 - beta) / 1920.0 * lnFoverK * lnFoverK * lnFoverK * lnFoverK);

	if (F != K)
	{
		z = v / alpha * FKTerm * lnFoverK;
		x = log((sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho));
		SABRIV = z / x * numerator / denominator;
	}
	else
	{
		SABRIV = numerator / denominator;
	}
	
	return SABRIV;
}

double SABRIV(double beta, double* params, double Fut, double K, double T)
{
	double v, rho, alpha;
	alpha = params[0];
	v = params[1];
	rho = params[2];
	return SABRIV(alpha, beta, v, rho, Fut, K, T);
}

void make_Jacov_SABR(
	long NParams,
	double* Params,
	long NResidual,
	double** TempJacovMatrix,
	double* ParamsUp,
	double* ParamsDn,
	double* TermVolNew,
	double* ParityVolNew,
	double* VolNew,
	double Beta,
	double Futures
)
{
	long i;
	long j;
	long n;

	double dparam_up = 0.0001;
	double Pup = 0.0;
	double Pdn = 0.0;
	double ErrorUp = 0.0;
	double ErrorDn = 0.0;

	long temp;

	for (i = 0; i < NResidual; i++)
	{
		for (j = 0; j < NParams; j++)
		{
			if (i == 80)
			{
				temp = 0;
			}

			for (n = 0; n < NParams; n++)
			{
				ParamsUp[n] = Params[n];
				ParamsDn[n] = Params[n];
			}

			if (j == 0 && Params[j] <= 0.0)
			{
				// Alpha 조건
				ParamsUp[j] = 0.0001 + 2.0 * dparam_up;
				ParamsDn[j] = 0.0001;
			}
			else if (j == 1 && Params[j] <= 0.0)
			{
				// Nu 조건
				ParamsUp[j] = 0.0001 + 2.0 * dparam_up;
				ParamsDn[j] = 0.0001;
			}
			else if (j == 2 && Params[j] >= 1.0 )
			{
				// Rho 조건
				ParamsUp[j] = 1.0;
				ParamsDn[j] = 1.0 - 2.0 * dparam_up;
			}
			else if (j == 2 && Params[j] <= -1.0)
			{
				// Rho 조건
				ParamsUp[j] = -1.0 + 2.0 * dparam_up;
				ParamsDn[j] = -1.0;
			}
			else
			{
				ParamsUp[j] = Params[j] + dparam_up;
				ParamsDn[j] = Params[j] - dparam_up;
			}

			Pup = SABRIV(Beta, ParamsUp, Futures, ParityVolNew[i], TermVolNew[i]);
			ErrorUp = VolNew[i] - Pup;
			Pdn = SABRIV(Beta, ParamsDn, Futures, ParityVolNew[i], TermVolNew[i]);
			ErrorDn = VolNew[i] - Pdn;
			TempJacovMatrix[i][j] = (ErrorUp - ErrorDn) / (2.0 * dparam_up);
		}
	}

}

void make_Residual_SABR(
	long NParams,
	double* Params,
	long NResidual,
	double* ResidualArray,
	double* TermVolNew,
	double* ParityVolNew,
	double* VolNew,
	double Futures,
	double Beta,
	double& absErrorSum,
	double* SABRVolNew
)
{
	long i;
	double s = 0.0;
	for (i = 0; i < NResidual; i++)
	{
		SABRVolNew[i] = SABRIV( Beta, Params, Futures, ParityVolNew[i], TermVolNew[i]);
	}
	for (i = 0; i < NResidual; i++) ResidualArray[i] = VolNew[i] - SABRVolNew[i];
	for (i = 0; i < NResidual; i++) s += fabs(ResidualArray[i]);
	absErrorSum = s;
}

void NextLambda(double ErrorSquareSum, double PrevErrorSquareSum, double* lambda, long& BreakFlag)
{
	double LambdaMax = 1000000;
	double LambdaMin = 0.0000001;

	if (ErrorSquareSum < PrevErrorSquareSum) *lambda *= 0.1;
	else *lambda *= 10.0;

	if (*lambda > LambdaMax) *lambda = LambdaMax;
	if (*lambda < LambdaMin) BreakFlag = 1;
}

void Levenberg_Marquardt_SABR(long NParams, long NResidual, double* NextParams, double* CurrentParams, double* lambda, double** Jacov, double* Residual, double& ParamSum, double** JT_J, double** Inverse_JT_J, double** JT_Res, double** ResultMatrix)
{
	long i;
	double mu = *lambda;
	double alphamax = 3.0, alphamin = 0.00001;
	double vmax = 3.0, vmin = 0.00001;
	double rhomax = 0.999, rhomin = -0.999;

	long n = NResidual, m = NParams;
	long Shape_J[2] = { n,m };
	long Shape_JT[2] = { m,n };
	long Shape_Residual[2] = { n,1 };

	// J' dot J                 Shape = m * m
	XprimeDotX(Jacov, Shape_J, JT_J);

	// J'J + mu * diag(J'J)     Shape = m * m
	for (i = 0; i < m; i++) JT_J[i][i] = JT_J[i][i] + mu;// *JT_J[i][i];

	// inv(J'J + mu * diag(J'J))
	long Shape_Inv_JT_J[2] = { m,m };

	MatrixInversion(JT_J, m, Inverse_JT_J);

	// J' dot Res               Shape = m * n dot n * 1 = m * 1
	long Shape_JT_Res[2] = { m,1 };
	XprimeY(Jacov, Shape_J, Residual, n, JT_Res);

	Dot2dArray(Inverse_JT_J, Shape_Inv_JT_J, JT_Res, Shape_JT_Res, ResultMatrix);

	for (i = 0; i < NParams; i++)
	{
		NextParams[i] = CurrentParams[i] - ResultMatrix[i][0];
	}

	NextParams[0] = min(alphamax, max(alphamin, NextParams[0]));
	NextParams[1] = min(vmax, max(vmin, NextParams[1]));
	NextParams[2] = min(rhomax, max(rhomin, NextParams[2]));

	double s = 0.0;
	for (i = 0; i < NParams; i++) s += fabs(ResultMatrix[i][0]);
	ParamSum = s;
}

void Levenberg_Marquardt_SABR(
	long NParams,
	double* Params,
	long NResidual,
	double* ResidualArray,
	double** TempJacovMatrix,
	double* ParamsUp,
	double* ParamsDn,

	double* TermVolNew,
	double* ParityVolNew,
	double* VolNew,
	double Futures,
	double* SABRVolNew,
	double Beta

	)
{
	long i;
	long n;
	long BreakFlag = 0;
	long Levenberg = 1;
	double StopCondition = 0.0001;

	double absErrorSum = 100000.0;
	double PrevAbsErrorSum = 0.0;
	double ParamSum = 10000.0;
	double lambda[1] = { 1.00 };
	double* NextParams = make_array(NParams);
	double** JT_J = make_array(NParams, NParams);
	double** Inverse_JT_J = make_array(NParams, NParams);
	double** JT_Res = make_array(NParams, 1);
	double** ResultMatrix = make_array(NParams, 1);

	if (Levenberg == 0)
	{
		for (n = 0; n < 50; n++)
		{
			make_Jacov_SABR(NParams, Params, NResidual, TempJacovMatrix, ParamsUp, ParamsDn, TermVolNew, ParityVolNew, VolNew,   Beta, Futures);

			Print_Array(TempJacovMatrix, NResidual, NParams);
			printf("\n");

			make_Residual_SABR(NParams, Params, NResidual, ResidualArray, TermVolNew, ParityVolNew, VolNew, Futures,  Beta, absErrorSum, SABRVolNew);

			if (n >= 1) NextLambda(absErrorSum, PrevAbsErrorSum, lambda, BreakFlag);

			Levenberg_Marquardt_SABR(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_J, Inverse_JT_J, JT_Res, ResultMatrix);
			for (i = 0; i < NParams; i++) Params[i] = NextParams[i];

			if (ParamSum < StopCondition && n > 10) break;
			if (BreakFlag == 1) break;

			PrevAbsErrorSum = absErrorSum;
		}
	}
	else
	{
		lambda[0] = 0.0;
		for (n = 0; n < 50; n++)
		{
			make_Jacov_SABR(NParams, Params, NResidual, TempJacovMatrix, ParamsUp, ParamsDn, TermVolNew, ParityVolNew, VolNew, Beta, Futures);

			//Print_Array(TempJacovMatrix, NResidual, NParams);
			//printf("\n");

			make_Residual_SABR(NParams, Params, NResidual, ResidualArray, TermVolNew, ParityVolNew, VolNew, Futures,  Beta, absErrorSum, SABRVolNew);

			//Print_Array(ResidualArray, NResidual);
			//printf("\n");

			if (n >= 1) NextLambda(absErrorSum, PrevAbsErrorSum, lambda, BreakFlag);

			Levenberg_Marquardt_SABR(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_J, Inverse_JT_J, JT_Res, ResultMatrix);
			for (i = 0; i < NParams; i++) Params[i] = NextParams[i];

			if (ParamSum < StopCondition && n > 10) break;

			PrevAbsErrorSum = absErrorSum;
		}
	}

	free(NextParams);
	for (i = 0; i < NParams; i++) free(JT_J[i]);
	for (i = 0; i < NParams; i++) free(Inverse_JT_J[i]);
	free(JT_J);
	free(Inverse_JT_J);
	for (i = 0; i < NParams; i++) free(JT_Res[i]);
	free(JT_Res);
	for (i = 0; i < NParams; i++) free(ResultMatrix[i]);
	free(ResultMatrix);

}

void main2(long NTermVol, double* TermVol, long NParityVol, double* ParityVol, double* Vol, double Futures, double Beta, double* Params, double* ResultLocVol)
{
	long i;
	long j;

	const long nparams = 3;
	long NResidual = NTermVol * NParityVol;

	double* TermVolNew = (double*)malloc(sizeof(double) * NTermVol * NParityVol);
	for (i = 0; i < NTermVol * NParityVol; i++) TermVolNew[i] = TermVol[i % NTermVol];
	double* ParityVolNew = (double*)malloc(sizeof(double) * NTermVol * NParityVol);
	for (i = 0; i < NTermVol * NParityVol; i++) ParityVolNew[i] = ParityVol[i / NTermVol];
	double* VolNew = (double*)malloc(sizeof(double) * NTermVol * NParityVol);
	for (i = 0; i < NTermVol * NParityVol; i++) VolNew[i] = Vol[i];
	double* SABRVolNew = (double*)malloc(sizeof(double) * NTermVol * NParityVol);
	for (i = 0; i < NTermVol * NParityVol; i++) SABRVolNew[i] = Vol[i];

	double* ResidualArray = (double*)malloc(sizeof(double) * NResidual);
	double** TempJacov = (double**)malloc(sizeof(double*) * NResidual);
	for (i = 0; i < NResidual; i++) TempJacov[i] = (double*)malloc(sizeof(double) * nparams);

	double TargetParams[nparams] = { 0.0, };
	double ParamsUp[nparams] = { 0.0, };
	double ParamsDn[nparams] = { 0.0, };
	for (i = 0; i < nparams; i++)
	{
		TargetParams[i] = Params[i];
		ParamsUp[i] = Params[i];
		ParamsDn[i] = Params[i];
	}

	Levenberg_Marquardt_SABR(nparams, TargetParams, NResidual, ResidualArray, TempJacov, ParamsUp, ParamsDn, TermVolNew, ParityVolNew, VolNew, Futures, SABRVolNew, Beta);
	for (i = 0; i < NTermVol * NParityVol; i++) ResultLocVol[i] = SABRVolNew[i];
	for (i = 0; i < nparams; i++) Params[i] = TargetParams[i];
	free(TermVolNew);
	free(ParityVolNew);
	free(ResidualArray);
	for (i = 0; i < NResidual; i++) free(TempJacov[i]);
	free(TempJacov);
	free(SABRVolNew);
	free(VolNew);
}

int main()
{
	long i;
	long j;
	long k;

	const long NTermVol = 10;
	double TermVol[NTermVol] = { 0.33, 0.67, 1.0, 1.33, 1.67, 2.0, 2.33, 2.67, 3.0 , 3.33 };

	const long NParityVol = 9;
	double ParityVol[NParityVol] = { 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3 };

	const long NVol = NTermVol * NParityVol;
	double Vol[NVol] = {    0.3995300,	0.3543100, 	0.3410300, 	0.3120800,  0.2920000,  0.2783500 	,0.2685600 	,0.2584700 	,0.2477000 	,0.2433400,
							0.3587200, 	0.3049600, 	0.2950800, 	0.2692700, 	0.2513900, 	0.2385400, 	0.2305500, 	0.2237300, 	0.2153000, 	0.2125000,
							0.2974900, 	0.2594500, 	0.2474600, 	0.2300200, 	0.2175900, 	0.2072700, 	0.2020400, 	0.1983300, 	0.1936200, 	0.1922400,
							0.2413000, 	0.2179900, 	0.2082200, 	0.1979000, 	0.1906900, 	0.1839800, 	0.1819600, 	0.1807300, 	0.1782800, 	0.1783200,
							0.1858300, 	0.1775800, 	0.1708200, 	0.1696900, 	0.1689100, 	0.1661900, 	0.1671200, 	0.1680300, 	0.1674100, 	0.1688500,
							0.1340200, 	0.1410400, 	0.1426700, 	0.1474900, 	0.1510600, 	0.1518200, 	0.1553700, 	0.1582300, 	0.1591600, 	0.1616100,
							0.1032700, 	0.1141300, 	0.1281000, 	0.1334600, 	0.1378600, 	0.1417300, 	0.1468200, 	0.1510900, 	0.1538600, 	0.1567300,
							0.1098400, 	0.1083000, 	0.1213600, 	0.1260800, 	0.1312800, 	0.1370000, 	0.1428800, 	0.1480200, 	0.1520600, 	0.1551400,
							0.1148200, 	0.1139000, 	0.1179300, 	0.1247500, 	0.1310300, 	0.1354900, 	0.1418000, 	0.1477500, 	0.1526000, 	0.1559700};
	double ResultImpliedVolReshaped[NVol] = { 0.0 };

	double alpha, v, rho, Fut, Beta;
	double ResultParams[3] = { 0., };
	
	double** TermParams = (double**)malloc(sizeof(double*) * NTermVol);
	for (i = 0; i < NTermVol; i++) TermParams[i] = (double*)malloc(sizeof(double) * 4);

	double* VolArray = (double*)malloc(sizeof(double) *NParityVol);
	double** ResultImVol = (double**)malloc(sizeof(double*) * NParityVol);
	for (i = 0; i < NParityVol; i++) ResultImVol[i] = (double*)malloc(sizeof(double) * NTermVol);
	double* TempVolResult = (double*)malloc(sizeof(double) * NParityVol);

	for (i = 0; i < NTermVol; i++) TermParams[i][3] = TermVol[i];

	for (i = 0; i < NTermVol; i++)
	{
		alpha = 0.15;
		v = 0.25;
		rho = -0.05;
		Fut = 1.0;
		TermParams[i][0] = alpha;
		TermParams[i][1] = v;
		TermParams[i][2] = rho;
		Beta = 1.0;

		k = 0;
		for (j = 0; j < NParityVol; j++)
		{
			VolArray[k] = Vol[i + j * NTermVol];
			k += 1;
		}
		main2(1, TermVol + i, NParityVol, ParityVol, VolArray, Fut, Beta, TermParams[i], TempVolResult);

		for (j = 0; j < NParityVol; j++) ResultImVol[j][i] = TempVolResult[j];
		
	}

	printf("\nSABR Implied Volatility\n");
	Print_Array(ResultImVol, NParityVol, NTermVol);
	printf("\n\n Alpha,      Vega,      Rho  ,  Term \n");
	Print_Array(TermParams, NTermVol, 4);

	k = 0;
	for (i = 0; i < NParityVol; i++)
	{
		for (j = 0; j < NTermVol; j++)
		{
			ResultImpliedVolReshaped[k] = ResultImVol[i][j];
			k += 1;
		}
	}

	long N_Rf = 2;
	double RfTerm[2] = { 1.0, 2.0 };
	double RfRate[2] = { 0.02, 0.02 };

	long N_Div = 2;
	double DivTerm[2] = { 1.0, 2.0 };
	double DivRate[2] = { 0.02, 0.02 };

	curveinfo RfCurve(N_Rf, RfTerm, RfRate);
	curveinfo DivCurve(N_Div, DivTerm, DivRate);

	volinfo VolMatrix(NParityVol, ParityVol, NTermVol, TermVol, ResultImpliedVolReshaped);
	
	VolMatrix.set_localvol(&RfCurve, &DivCurve);
	printf("\n\nSABR Local Volatility\n\n");
	Print_Array(VolMatrix.LocalVolMat, NParityVol, NTermVol);

	for (i = 0; i < NTermVol; i++) free(TermParams[i]);
	free(TermParams);
	free(VolArray);
	for (i = 0; i < NParityVol; i++) free(ResultImVol[i]);
	free(ResultImVol);
	free(TempVolResult);
}