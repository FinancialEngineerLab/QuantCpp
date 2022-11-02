#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef UTILITY
#include "Util.h"
#endif

#ifndef STRUCTURE
#include "Structure.h"
#endif

double SABRIV(double alpha, double beta, double v, double rho, double Fut, double K, double T)
{
	long i;
	double z, x;
	double F = Fut;
	double numerator, denominator;
	double FKTerm = pow(F * K, (1.0 - beta) * 0.5);
	double lnFoverK = log(F / K);
	double SABRIV;

	numerator = alpha * (1.0 + ((1.0 - beta) * (1.0 - beta) * alpha * alpha / (24.0 * FKTerm) + 0.25 * (rho * beta * v * alpha) / FKTerm + (2.0 - 3.0 * rho * rho) / 24.0 * v * v) * T);
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
	long NTermFutures,
	double* TermFuturesArray,
	double* FuturesArray
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
	double Futures = 1.0;

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
			else if (j == 2 && Params[j] >= 1.0)
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
			Futures = Interpolate_Linear(TermFuturesArray, FuturesArray, NTermFutures, TermVolNew[i]);
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
	long NTermFutures,
	double* TermFuturesArray,
	double* FuturesArray,
	double Beta,
	double& ErrorSquareSum,
	double* SABRVolNew
)
{
	long i;
	double s = 0.0;
	double Futures = 1.0;
	for (i = 0; i < NResidual; i++)
	{
		Futures = Interpolate_Linear(TermFuturesArray, FuturesArray, NTermFutures, TermVolNew[i]);
		SABRVolNew[i] = SABRIV(Beta, Params, Futures, ParityVolNew[i], TermVolNew[i]);
	}
	for (i = 0; i < NResidual; i++) ResidualArray[i] = VolNew[i] - SABRVolNew[i];
	for (i = 0; i < NResidual; i++) s += (ResidualArray[i] * ResidualArray[i]);
	ErrorSquareSum = s;
}

void NextLambda(double ErrorSquareSum, double PrevErrorSquareSum, double* lambda, long& BreakFlag)
{
	double LambdaMax = 1000000;
	double LambdaMin = 0.00001;

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
	rounding(JT_J, m, m, 5);

	// J'J + mu * diag(J'J)     Shape = m * m
	for (i = 0; i < m; i++) JT_J[i][i] = JT_J[i][i] + mu;//*JT_J[i][i];

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
	long NTermFutures,
	double* TermFuturesArray,
	double* FuturesArray,
	double* SABRVolNew,
	double Beta

)
{
	long i;
	long n;
	long BreakFlag = 0;
	long Levenberg = 1;
	double StopCondition = 0.00001;

	double ErrorSquareSum = 100000.0;
	double PrevErrorSquareSum = 0.0;
	double ParamSum = 0.0;
	double minerror = 1000000.0;
	double lambda[1] = { 1.00 };
	double* NextParams = make_array(NParams);
	double** JT_J = make_array(NParams, NParams);
	double** Inverse_JT_J = make_array(NParams, NParams);
	double** JT_Res = make_array(NParams, 1);
	double** ResultMatrix = make_array(NParams, 1);
	double FirstErrorSquare = 1.0;
	double* argminparam = make_array(NParams);
	double* argminvol = make_array(NResidual);
	for (n = 0; n < 50; n++)
	{
		make_Jacov_SABR(NParams, Params, NResidual, TempJacovMatrix, ParamsUp, ParamsDn, TermVolNew, ParityVolNew, VolNew, Beta, NTermFutures, TermFuturesArray, FuturesArray);

		make_Residual_SABR(NParams, Params, NResidual, ResidualArray, TermVolNew, ParityVolNew, VolNew, NTermFutures, TermFuturesArray, FuturesArray, Beta, ErrorSquareSum, SABRVolNew);

		if (n == 0) FirstErrorSquare = ErrorSquareSum + 0.0;

		if (n >= 1)
		{
			NextLambda(ErrorSquareSum, PrevErrorSquareSum, lambda, BreakFlag);
			if (ErrorSquareSum < minerror || n == 1)
			{
				minerror = ErrorSquareSum;
				for (i = 0; i < NParams; i++) argminparam[i] = Params[i];
				for (i = 0; i < NResidual; i++) argminvol[i] = SABRVolNew[i];
			}
		}

		Levenberg_Marquardt_SABR(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_J, Inverse_JT_J, JT_Res, ResultMatrix);
		for (i = 0; i < NParams; i++) Params[i] = NextParams[i];

		if ((ParamSum < StopCondition || lambda[0] < 1.0e-5) && n > 10) break;
		if (ErrorSquareSum / FirstErrorSquare < 0.001) break;

		PrevErrorSquareSum = ErrorSquareSum;
	}


	for (i = 0; i < NParams; i++) Params[i] = argminparam[i];
	for (i = 0; i < NResidual; i++) SABRVolNew[i] = argminvol[i];

	free(NextParams);
	for (i = 0; i < NParams; i++) free(JT_J[i]);
	for (i = 0; i < NParams; i++) free(Inverse_JT_J[i]);
	free(JT_J);
	free(Inverse_JT_J);
	for (i = 0; i < NParams; i++) free(JT_Res[i]);
	free(JT_Res);
	for (i = 0; i < NParams; i++) free(ResultMatrix[i]);
	free(ResultMatrix);
	free(argminparam);
	free(argminvol);
}

void SABRCalibration(long NTermVol, double* TermVol, long NParityVol, double* ParityVol, double* Vol, long NTermFutures, double* TermFuturesArray, double* FuturesArray, double Beta, double* Params, double* ResultLocVol)
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

	Levenberg_Marquardt_SABR(nparams, TargetParams, NResidual, ResidualArray, TempJacov, ParamsUp, ParamsDn, TermVolNew, ParityVolNew, VolNew, NTermFutures, TermFuturesArray, FuturesArray, SABRVolNew, Beta);
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

long SABR_Vol(
	long N_Rf,
	double* RfTerm,
	double* RfRate,
	long N_Div,
	double* DivTerm,

	double* DivRate,
	long NTermVol,
	double* TermVol,
	long NParityVol,
	double* ParityVol,

	double* Vol,
	long CalcLocalVolFlag,
	double Beta,
	double* ResultImpliedVolReshaped,
	double* ResultLocalVolReshaped,

	double* ResultParams,
	double* Futures
)
{
	long i;
	long j;
	long k;

	double r, d;

	curveinfo RfCurve(N_Rf, RfTerm, RfRate);
	curveinfo DivCurve(N_Div, DivTerm, DivRate);

	long NVol = NTermVol * NParityVol;

	long NTermFutures = NTermVol;
	double* TermFuturesArray = TermVol;
	double* FuturesArray = (double*)malloc(sizeof(double) * NTermFutures);
	for (i = 0; i < NTermFutures; i++)
	{
		r = Interpolate_Linear(RfTerm, RfRate, N_Rf, TermFuturesArray[i]);
		d = Interpolate_Linear(DivTerm, DivRate, N_Div, TermFuturesArray[i]);
		FuturesArray[i] = exp(r - d);
		Futures[i] = FuturesArray[i];
	}

	double alpha = 0.15;
	double v = 0.25;
	double rho = -0.05;

	double** TermParams = make_array(NTermVol, 4);
	double** ResultImVol = make_array(NParityVol, NTermVol);

	double* VolArray = (double*)malloc(sizeof(double) * NParityVol);
	double* TempVolResult = (double*)malloc(sizeof(double) * NParityVol);

	for (i = 0; i < NTermVol; i++)
	{

		TermParams[i][0] = alpha;
		TermParams[i][1] = v;
		TermParams[i][2] = rho;
		TermParams[i][3] = TermVol[i];
		k = 0;
		for (j = 0; j < NParityVol; j++)
		{
			VolArray[k] = Vol[i + j * NTermVol];
			k += 1;
		}
		SABRCalibration(1, TermVol + i, NParityVol, ParityVol, VolArray, NTermFutures, TermFuturesArray, FuturesArray, Beta, TermParams[i], TempVolResult);

		for (j = 0; j < NParityVol; j++) ResultImVol[j][i] = TempVolResult[j];
	}

	k = 0;
	for (i = 0; i < NParityVol; i++)
	{
		for (j = 0; j < NTermVol; j++)
		{
			ResultImpliedVolReshaped[k] = ResultImVol[i][j];
			k += 1;
		}
	}

	volinfo VolMatrix(NParityVol, ParityVol, NTermVol, TermVol, ResultImpliedVolReshaped);
	if (CalcLocalVolFlag != 0)
	{
		VolMatrix.set_localvol(&RfCurve, &DivCurve);
		k = 0;
		for (i = 0; i < NParityVol; i++)
		{
			for (j = 0; j < NTermVol; j++)
			{
				ResultLocalVolReshaped[k] = VolMatrix.LocalVolMat[i][j];
				k += 1;
			}
		}
	}

	k = 0;
	for (i = 0; i < NTermVol; i++)
	{
		for (j = 0; j < 3; j++)
		{
			ResultParams[k] = TermParams[i][j];
			k += 1;
		}
	}

	for (i = 0; i < NTermVol; i++) free(TermParams[i]);
	free(TermParams);
	free(VolArray);
	for (i = 0; i < NParityVol; i++) free(ResultImVol[i]);
	free(ResultImVol);
	free(TempVolResult);
	free(FuturesArray);

	return 1;
}

DLLEXPORT(long) Excel_SABR_Vol(
	long N_Rf,							// Risk Free Term 개수
	double* RfTerm,						// Risk Free Term
	double* RfRate,						// Risk Free Rate
	long N_Div,							// Dividend Term 개수
	double* DivTerm,					// Dividend Term 

	double* DivRate,					// Dividend Rate
	long NTermVol,						// Implied Vol Surface의 Tenor 개수
	double* TermVol,					// Implied Vol Surface의 Tenor
	long NParityVol,					// Implied Vol Surface의 Moneyness 개수
	double* ParityVol,					// Implied Vol Surface의 Moneyness

	double* Vol,						// Implied Vol Surface.Reshape(0)
	long CalcLocalVolFlag,				// Local Volatility 산출여부
	double Beta,						// SABR Beta Parameter
	double* ResultImpliedVolReshaped,	// Output : SABR Implied Vol
	double* ResultLocalVolReshaped,		// Output : SABR LocalVol

	double* ResultParams,				// Term별 SABR 파라미터
	double* Futures						// Futures 가격
)
{
	long i, j;
	long ResultCode = 0;

	ResultCode = SABR_Vol(
		N_Rf, RfTerm, RfRate, N_Div, DivTerm,
		DivRate, NTermVol, TermVol, NParityVol, ParityVol,
		Vol, CalcLocalVolFlag, Beta, ResultImpliedVolReshaped, ResultLocalVolReshaped,
		ResultParams, Futures
	);

	//_CrtDumpMemoryLeaks();
	return ResultCode;
}