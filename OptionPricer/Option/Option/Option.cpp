#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "CalcDate.h"
#include "Structure.h"
#ifndef UTILITY
#include "Util.h"
#endif
#include <crtdbg.h>

double Sumation(long NArray, double* Array)
{
	long i;
	double s = 0.0;
	for (i = 0; i < NArray; i++)
	{
		s += Array[i];
	}
	return s;
}

long SaveErrorName(char* Error, char ErrorName[100])
{
	long k;
	long i;
	long Ascii;
	k = 0;
	for (i = 0; i < 100; i++)
	{
		Ascii = (long)ErrorName[i];
		if ((Ascii >= 48 && Ascii <= 57) || (Ascii >= 65 && Ascii <= 90) || (Ascii >= 97 && Ascii <= 122) || (Ascii == 32))
		{
			Error[k] = ErrorName[i];
			k++;
		}
		else
		{
			break;
		}
	}
	return -1;
}
/*
double Call_Down_In(
	double S,
	double X,
	double H,
	double T,
	double r_disc,
	double r_ref,
	double q,
	double sig
)
{
	long i, j;
	double lambda = (r_ref - q + sig * sig * 0.5) / (sig * sig);
	double y = log( H*H/ )
}*/

long BSOption(
	long TypeFlag,			// Option Type: 1이면 콜, 2이면 풋
	double S0,				// 기초자산 현재가
	double X,				// 행사가격
	curveinfo* DiscCurve,	// 할인커브
	curveinfo* RefCurve,	// 레퍼런스커브
	curveinfo* DivCurve,	// 배당커브
	double QuantoCorr,		// Quanto COrrelation
	curveinfo* FXVolCurve,	// FXVol 커브
	volinfo* VolMat,		// Vol 정보
	long DivTypeFlag,		// Div 배당타입
	double T_Ref,
	double T_Pay,
	double* ResultPrice
)
{
	long i;
	long N;
	double dvd;
	double vol;
	double fxvol;
	double r, r_disc, Rf_Quanto, frate;
	double Price, Delta, Gamma, Vega, Theta, Rho;

	if (DivTypeFlag == 0)
	{
		dvd = DivCurve->Rate[0];
	}
	else if (DivTypeFlag == 1)
	{
		dvd = DivCurve->Interpolated_Rate(T_Ref);
	}
	else if (DivTypeFlag == 2)
	{
		N = 0;
		for (i = 0; i < DivCurve->nterm(); i++)
		{
			if (DivCurve->Term[i] > T_Ref) break;
			N += 1;
		}
		dvd = 1.0/T_Ref * Sumation(N, DivCurve->Rate) / S0;
	}
	vol = VolMat->Calc_Implied_Volatility(T_Ref, S0 / X);
	fxvol = FXVolCurve->Interpolated_Rate(T_Ref);
	r = RefCurve->Interpolated_Rate(T_Ref);
	r_disc = DiscCurve->Interpolated_Rate(T_Pay);
	Rf_Quanto = r - fxvol * vol * QuantoCorr;

	double d1 = (log(S0 / X) + (r - dvd + 0.5 * vol * vol) * T_Ref) / (vol * sqrt(T_Ref));
	double d2 = d1 - vol * sqrt(T_Ref);

	if (T_Ref <= 0.0)
	{
		if (TypeFlag == 1)
		{
			Price = max(0, S0 - X) * exp(-r_disc * T_Pay);
			Delta = 1.0;
			Gamma = 0.0;
			Vega = 0.0;
			Theta = 0.0;
			Rho = 0.0;
		}
		else
		{
			Price = max(0, X - S0) * exp(-r_disc * T_Pay);
			Delta = -1.0;
			Gamma = 0.0;
			Vega = 0.0;
			Theta = 0.0;
			Rho = 0.0;
		}
	}
	else
	{
		if (TypeFlag == 1) 
		{ //call
			if (X <= 0.0)
			{
				Price = S0 * exp(-dvd * T_Ref) - X * exp(-Rf_Quanto * T_Ref);
				Delta = exp(-dvd * T_Ref);
				Gamma = 0.0;
				Vega = 0.0;
				Theta = dvd * S0 * exp(-dvd * T_Ref) - r * X * exp(-r * T_Ref);
				Rho = X * T_Ref * exp(-r * T_Ref);
			}
			else
			{
				Price = S0 * exp(-dvd * T_Ref) * CDF_N(d1) - X * exp(-Rf_Quanto * T_Ref) * CDF_N(d2);
				Delta = exp(-dvd * T_Ref) * CDF_N(d1);
				Gamma = exp(-dvd * T_Ref) * PDF_N(d1) / (S0 * vol * sqrt(T_Ref));
				Vega = exp(-dvd * T_Ref) * S0 * sqrt(T_Ref) * PDF_N(d1);
				Theta = -0.5 * S0 * PDF_N(d1) * vol * exp(-dvd * T_Ref) / sqrt(T_Ref) + dvd * S0 * CDF_N(d1) * exp(-dvd * T_Ref) - r * X * exp(-r * T_Ref) * CDF_N(d2);
				Rho = X * T_Ref * exp(-r * T_Ref) * CDF_N(d2);
			}
		}
		else
		{	// put
			if (X <= 0.0)
			{
				Price = 0.0;
				Delta = 0.0;
				Gamma = 0.0;
				Vega = 0.0;
				Theta = 0.0;
				Rho = 0.0;
			}
			else
			{
				Price = -S0 * exp(-dvd * T_Ref) * CDF_N(-d1) + X * exp(-Rf_Quanto * T_Ref) * CDF_N(-d2);
				Delta = exp(-dvd * T_Ref) * (CDF_N(d1) - 1.0);
				Gamma = exp(-dvd * T_Ref) * PDF_N(d1) / (S0 * vol * sqrt(T_Ref));
				Vega = exp(-dvd * T_Ref) * S0 * sqrt(T_Ref) * PDF_N(d1);
				Theta = 0.5 * S0 * PDF_N(d1) * vol * exp(-dvd * T_Ref) / sqrt(T_Ref) - dvd * S0 * CDF_N(-d1) * exp(-dvd * T_Ref) + r * X * exp(-r * T_Ref) * CDF_N(-d2);
				Rho = -X * T_Ref * exp(-r * T_Ref) * CDF_N(-d2);
			}
		}
		if (T_Pay - T_Ref <= 0.0) frate = 0.0;
		else frate = Calc_Forward_Rate(DiscCurve->Term, DiscCurve->Rate, DiscCurve->nterm(), T_Ref, T_Pay);
		Price = Price / (1.0 + frate * (T_Pay - T_Ref));
		Delta = Delta / (1.0 + frate * (T_Pay - T_Ref));
		Gamma = Gamma / (1.0 + frate * (T_Pay - T_Ref));
		Vega = Vega / (1.0 + frate * (T_Pay - T_Ref))/100.0;
		Theta = Theta / (1.0 + frate * (T_Pay - T_Ref)) * 1.0/ 365.0;
		Rho = Rho / (1.0 + frate * (T_Pay - T_Ref))/100.0;
	}
	ResultPrice[0] = Price;
	ResultPrice[1] = Delta;
	ResultPrice[2] = Gamma;
	ResultPrice[3] = Vega;
	ResultPrice[4] = Theta;
	ResultPrice[5] = Rho;
	return 1;
}

long BSDigitalOption(
	long TypeFlag,			// Option Type: 1이면 콜, 2이면 풋
	double S0,				// 기초자산 현재가
	double X,				// 행사가격
	curveinfo* DiscCurve,	// 할인커브
	curveinfo* RefCurve,	// 레퍼런스커브
	curveinfo* DivCurve,	// 배당커브
	double QuantoCorr,		// Quanto COrrelation
	curveinfo* FXVolCurve,	// FXVol 커브
	volinfo* VolMat,		// Vol 정보
	long DivTypeFlag,		// Div 배당타입
	double T_Ref,
	double T_Pay,
	double* ResultPrice
)
{
	long i;
	long N;
	double dvd;
	double vol;
	double fxvol;
	double r, r_disc, Rf_Quanto, frate;
	double Price, Delta, Gamma, Vega, Theta, Rho;
	double fix_amt = 1.0;

	if (DivTypeFlag == 0)
	{
		dvd = DivCurve->Rate[0];
	}
	else if (DivTypeFlag == 1)
	{
		dvd = DivCurve->Interpolated_Rate(T_Ref);
	}
	else if (DivTypeFlag == 2)
	{
		N = 0;
		for (i = 0; i < DivCurve->nterm(); i++)
		{
			if (DivCurve->Term[i] > T_Ref) break;
			N += 1;
		}
		dvd = 1.0 / T_Ref * Sumation(N, DivCurve->Rate) / S0;
	}
	vol = VolMat->Calc_Implied_Volatility(T_Ref, S0 / X);
	fxvol = FXVolCurve->Interpolated_Rate(T_Ref);
	r = RefCurve->Interpolated_Rate(T_Ref);
	r_disc = DiscCurve->Interpolated_Rate(T_Pay);
	Rf_Quanto = r - fxvol * vol * QuantoCorr;

	double d1 = (log(S0 / X) + (r - dvd + 0.5 * vol * vol) * T_Ref) / (vol * sqrt(T_Ref));
	double d2 = d1 - vol * sqrt(T_Ref);
	double theta_d2 = -(r - dvd - 0.5 * vol * vol) / (vol * sqrt(T_Ref)) + 0.5 * d2 / T_Ref;
	if (T_Ref <= 0.0)
	{
		if (TypeFlag == 1)
		{
			if (S0 >= X) Price = fix_amt * exp(-r_disc * T_Pay);
			else Price = 0.0;

			Delta = 0.0;
			Gamma = 0.0;
			Vega = 0.0;
			Theta = 0.0;
			Rho = 0.0;
		}
		else
		{
			if (S0 <= X) Price = fix_amt * exp(-r_disc * T_Pay);
			else Price = 0.0;

			Delta = 0.0;
			Gamma = 0.0;
			Vega = 0.0;
			Theta = 0.0;
			Rho = 0.0;
		}
	}
	else
	{
		if (TypeFlag == 1)
		{ //call
			if (X <= 0.0)
			{
				Price = fix_amt * exp(-Rf_Quanto * T_Ref);
				Delta = 0.0;
				Gamma = 0.0;
				Vega = 0.0;
				Theta = dvd * S0 * exp(-dvd * T_Ref) - r * X * exp(-r * T_Ref);
				Rho = X * T_Ref * exp(-r * T_Ref);
			}
			else
			{
				Price = fix_amt * exp(-Rf_Quanto * T_Ref) * CDF_N(d2) ;
				Delta = fix_amt * exp(-Rf_Quanto * T_Ref) * PDF_N(d2) / (vol * sqrt(T_Ref) * S0);
				Gamma = -fix_amt * exp(-Rf_Quanto * T_Ref) * PDF_N(d2) / (vol * vol * (T_Ref)*S0 * S0) * (d2 + (vol * sqrt(T_Ref)));
				Vega = -fix_amt * exp(-Rf_Quanto * T_Ref) * PDF_N(d2) * (sqrt(T_Ref) + d2 / vol);
				Theta = fix_amt * exp(-Rf_Quanto * T_Ref) * (Rf_Quanto * CDF_N(d2) + PDF_N(d2) * theta_d2);
				Rho = -fix_amt * exp(-Rf_Quanto * T_Ref) * (T_Ref * CDF_N(d2) - PDF_N(d2) * sqrt(T_Ref) / vol);
			}
		}
		else
		{	// put
			if (X <= 0.0)
			{
				Price = 0.0;
				Delta = 0.0;
				Gamma = 0.0;
				Vega = 0.0;
				Theta = 0.0;
				Rho = 0.0;
			}
			else
			{
				Price = fix_amt * exp(-Rf_Quanto * T_Ref) * CDF_N(-d2);
				Delta = -fix_amt * exp(-Rf_Quanto * T_Ref) * PDF_N(d2) / (vol * sqrt(T_Ref) * S0);
				Gamma = fix_amt*exp(-Rf_Quanto * T_Ref) * PDF_N(d2) / (vol * vol * (T_Ref)*S0 * S0) * (d2 + (vol * sqrt(T_Ref)));
				Vega = fix_amt* exp(-Rf_Quanto * T_Ref)* PDF_N(d2)* (sqrt(T_Ref) + d2 / vol);
				Theta = fix_amt * exp(-Rf_Quanto * T_Ref) * (Rf_Quanto * CDF_N(-d2) - PDF_N(d2) * theta_d2);
				Rho = -(fix_amt * (exp(-Rf_Quanto * T_Ref) * (PDF_N(d2) * sqrt(T_Ref) / vol - T_Ref * CDF_N(-1.0 * d2))));
			}
		}
		if (T_Pay - T_Ref <= 0.0) frate = 0.0;
		else frate = Calc_Forward_Rate(DiscCurve->Term, DiscCurve->Rate, DiscCurve->nterm(), T_Ref, T_Pay);
		Price = Price / (1.0 + frate * (T_Pay - T_Ref));
		Delta = Delta / (1.0 + frate * (T_Pay - T_Ref));
		Gamma = Gamma / (1.0 + frate * (T_Pay - T_Ref));
		Vega = Vega / (1.0 + frate * (T_Pay - T_Ref)) / 100.0;
		Theta = Theta / (1.0 + frate * (T_Pay - T_Ref)) * 1.0 / 365.0;
		Rho = Rho / (1.0 + frate * (T_Pay - T_Ref)) / 100.0;
	}
	ResultPrice[0] = Price;
	ResultPrice[1] = Delta;
	ResultPrice[2] = Gamma;
	ResultPrice[3] = Vega;
	ResultPrice[4] = Theta;
	ResultPrice[5] = Rho;
	return 1;
}

double AsianOptionSimul(
	long TypeFlag,			// Option Type: 1이면 콜, 2이면 풋
	double S0,				// 기초자산 현재가
	double X,				// 행사가격
	double MeanPrice,		// 평균가격
	curveinfo* DiscCurve,	// 할인커브
	curveinfo* RefCurve,	// 레퍼런스커브
	curveinfo* DivCurve,	// 배당커브
	double QuantoCorr,		// Quanto COrrelation
	curveinfo* FXVolCurve,	// FXVol 커브
	volinfo* VolMat,		// Vol 정보
	long DivTypeFlag,		// Div 배당타입
	long VolTypeFlag,
	long Days_CalcToMeanStart,
	long Days_CalcToMaturity,
	double T_CalcToPay,
	double** FixedRandn,
	long S_UpDnFlag
)
{
	long i, j;
	long N;
	long idx;
	long NSimul = 10000;
	long idx_T, idx_S;
	long Days_MeanStartToMaturity = Days_CalcToMaturity - Days_CalcToMeanStart;

	double OptPrice = 0.0;
	double eps;
	double vol, simvol = 0.0;
	double S_Start, StockPrice;
	double dt = 1.0 / 365.0;
	double T_MeanStart_To_Maturity;
	double w_MeanStart_To_CalcDate, T_MeanStart_To_CalcDate;
	double w_Simul_To_Maturity, T_Simul_To_Maturity;
	double T, T2;
	double Drift, Diffusion;
	if (Days_CalcToMeanStart < 0)
	{
		T_MeanStart_To_CalcDate = ((double)(-Days_CalcToMeanStart)) / 365.0;
		T_MeanStart_To_Maturity = ((double)(Days_CalcToMaturity - Days_CalcToMeanStart)) / 365.0;
		w_MeanStart_To_CalcDate = T_MeanStart_To_CalcDate / T_MeanStart_To_Maturity;
		w_Simul_To_Maturity = 1.0 - w_MeanStart_To_CalcDate;
		N = Days_CalcToMaturity;
	}
	else if (Days_CalcToMeanStart == 0)
	{
		T_MeanStart_To_CalcDate = 1.0 / 365.0;
		T_MeanStart_To_Maturity = ((double)(Days_CalcToMaturity - Days_CalcToMeanStart)) / 365.0;
		w_MeanStart_To_CalcDate = T_MeanStart_To_CalcDate / T_MeanStart_To_Maturity;
		w_Simul_To_Maturity = 1.0 - w_MeanStart_To_CalcDate;
		N = Days_CalcToMaturity;
		MeanPrice = S0;
	}
	else
	{
		w_MeanStart_To_CalcDate = 0.0;
		w_Simul_To_Maturity = 1.0;
		N = Days_CalcToMaturity - Days_CalcToMeanStart;
	}

	if (S_UpDnFlag == 0) S_Start = 1.0;
	else if (S_UpDnFlag == 1) S_Start = 1.01;
	else S_Start = 0.99;

	double* DailyForwardRate = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		DailyForwardRate[i] = Calc_Forward_Rate(RefCurve->Term, RefCurve->Rate, RefCurve->nterm(), ((double)(i)) * dt, ((double)(i + 1)) * dt);
	}
	double* DailyDiv = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	if (DivTypeFlag == 0)
	{
		for (i = 0; i < Days_CalcToMaturity; i++)
		{
			DailyDiv[i] = DivCurve->Rate[0];
		}
	}
	else if (DivTypeFlag == 1)
	{
		for (i = 0; i < Days_CalcToMaturity; i++)
		{
			DailyDiv[i] = Calc_Forward_Rate(DivCurve->Term, DivCurve->Rate, DivCurve->nterm(), ((double)(i)) * dt, ((double)(i + 1)) * dt);
		}
	}
	else if (DivTypeFlag == 2)
	{
		for (i = 0; i < Days_CalcToMaturity; i++) DailyDiv[i] = 0.0;

		idx = 0;
		for (j = 0; j < DivCurve->nterm(); j++)
		{
			T2 = DivCurve->Term[j];
			for (i = idx; i < Days_CalcToMaturity; i++)
			{
				T = ((double)(i + 1)) * dt;
				if (T >= T2)
				{
					DailyDiv[i] = DivCurve->Rate[j];
					idx = i;
					break;
				}
			}

		}
	}
	double* DailyFXVol = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		DailyFXVol[i] = Interpolate_Linear(FXVolCurve->Term, FXVolCurve->Rate, FXVolCurve->nterm(), ((double)(i + 1)) * dt);
	}

	if (VolTypeFlag == 2) simvol = VolMat->Calc_Implied_Volatility(((double)Days_CalcToMaturity)/365.0, S0/X);
	double DF = Calc_Discount_Factor(DiscCurve->Term, DiscCurve->Rate, DiscCurve->nterm(), T_CalcToPay);
	double AvgStock;
	if (Days_CalcToMaturity <= 0)
	{
		if (TypeFlag == 1)
		{
			return DF * max(MeanPrice - X, 0.0);
		}
		else
		{
			return DF* max(X - MeanPrice, 0.0) ;
		}
	}

	for (i = 0; i < NSimul; i++)
	{
		StockPrice = S_Start;
		idx_T = 0;
		idx_S = 0;
		T = 0.0;
		AvgStock = 0.0;
		for (j = 0; j < Days_CalcToMaturity; j++)
		{
			T = T + dt;
			eps = FixedRandn[i][j];
			if (VolTypeFlag == 2) vol = simvol;
			else vol = VolMat->Calc_Volatility_for_Simulation(T, StockPrice, &idx_T, &idx_S);
			
			if (DivTypeFlag == 2)
			{
				Drift = (DailyForwardRate[j] - 0.5 * vol * vol- QuantoCorr * vol * DailyFXVol[j]) * dt;
				Diffusion = vol * sqrt(dt) * eps;
				StockPrice = (StockPrice - DailyDiv[j]) * exp(Drift + Diffusion);
			}
			else
			{
				Drift = (DailyForwardRate[j] - DailyDiv[j] - 0.5 * vol * vol - QuantoCorr * vol * DailyFXVol[j]) * dt;
				Diffusion = vol * sqrt(dt) * eps;
				StockPrice = (StockPrice ) * exp(Drift + Diffusion);
			}
			
			if (j >= Days_CalcToMeanStart)
			{
				AvgStock += StockPrice/(double)N;
			}
		}
		StockPrice = S0 * StockPrice;
		if (Days_MeanStartToMaturity == 0) AvgStock = StockPrice;
		else
		{
			AvgStock = S0 * AvgStock;
			AvgStock = AvgStock * w_Simul_To_Maturity + w_MeanStart_To_CalcDate * MeanPrice;
		}

		if (TypeFlag == 1) OptPrice += DF * max(AvgStock - X, 0.0) / (double)NSimul;
		else OptPrice += DF * max(X - AvgStock, 0.0) / (double)NSimul;
	}

	free(DailyForwardRate);
	free(DailyDiv);
	free(DailyFXVol);
	return OptPrice;
}

long AsianOptionSimul(
	long TypeFlag,			// Option Type: 1이면 콜, 2이면 풋
	double S0,				// 기초자산 현재가
	double X,				// 행사가격
	double MeanPrice,		// 평균가격
	curveinfo* DiscCurve,	// 할인커브
	curveinfo* RefCurve,	// 레퍼런스커브
	curveinfo* DivCurve,	// 배당커브
	double QuantoCorr,		// Quanto COrrelation
	curveinfo* FXVolCurve,	// FXVol 커브
	volinfo* VolMat,		// Vol 정보
	long DivTypeFlag,		// Div 배당타입
	long VolTypeFlag,
	long Days_CalcToMeanStart,
	long Days_CalcToMaturity,
	double T_CalcToPay,
	long GreekFlag,
	double* ResultPrice
)
{
	long i, j, k;
	long NSimul = 10000;
	double dt = 1.0 / 365.0;
	double** FixedRandn = make_array_randn(NSimul, Days_CalcToMaturity + 10);
	double OptPrice, OptPriceU, OptPriceD, Delta, Gamma, Vega, Theta, Rho;
	randnorm(0);

	if (GreekFlag == 0)
	{
		OptPrice = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
									RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
									DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
									0);
		ResultPrice[0] = OptPrice;
	}
	else
	{
		OptPrice = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			0);
		ResultPrice[0] = OptPrice;
		OptPriceU = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			1);
		OptPriceD = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			-1);
		Delta = (OptPriceU - OptPriceD) / (2.0 * 0.01 * S0);
		Gamma = (OptPriceU + OptPriceD - 2.0 * OptPrice) / (0.01 * 0.01 * S0 * S0);
		ResultPrice[1] = Delta;
		ResultPrice[2] = Gamma;
		double* VolReshaped_Up = (double*)malloc(sizeof(double) * (VolMat->size(0) * VolMat->size(1)));
		double* VolReshaped_Dn = (double*)malloc(sizeof(double) * (VolMat->size(0) * VolMat->size(1)));
		k = 0;
		for (i = 0 ; i < VolMat->size(0); i++)
			for (j = 0; j < VolMat->size(1); j++)
			{
				VolReshaped_Up[k] = VolMat->Vol_Matrix[i][j] + 0.01;
				VolReshaped_Dn[k] = max(0.0, VolMat->Vol_Matrix[i][j] - 0.01);
				k += 1;
			}
		volinfo VolUp(VolMat->size(0), VolMat->Parity, VolMat->size(1), VolMat->Term, VolReshaped_Up);
		volinfo VolDn(VolMat->size(0), VolMat->Parity, VolMat->size(1), VolMat->Term, VolReshaped_Dn);
		if (VolTypeFlag == 0)
		{
			long NTempRate = 1;
			double TempRefRateTerm[1] = { 1.0 };
			double TempRefRate[1] = { 0.0 };
			if (DivTypeFlag == 2)
			{
				curveinfo TempCurve(NTempRate, TempRefRateTerm, TempRefRate);
				VolUp.set_localvol(&TempCurve, &TempCurve, 2.0, 0.01);
				fillna_Interpolate(VolUp.LocalVolMat, VolUp.size(0), VolUp.size(1));
				VolDn.set_localvol(&TempCurve, &TempCurve, 2.0, 0.01);
				fillna_Interpolate(VolDn.LocalVolMat, VolDn.size(0), VolDn.size(1));
			}
			else
			{
				VolUp.set_localvol(RefCurve, DivCurve, 2.0, 0.01);
				fillna_Interpolate(VolUp.LocalVolMat, VolUp.size(0), VolUp.size(1));
				VolDn.set_localvol(RefCurve, DivCurve, 2.0, 0.01);
				fillna_Interpolate(VolDn.LocalVolMat, VolDn.size(0), VolDn.size(1));
			}
		}
		else
		{
			VolUp.LocalVolMat = VolUp.Vol_Matrix;
			VolUp.Parity_Locvol = VolUp.Parity;
			VolUp.Term_Locvol = VolUp.Term;
			VolDn.LocalVolMat = VolDn.Vol_Matrix;
			VolDn.Parity_Locvol = VolDn.Parity;
			VolDn.Term_Locvol = VolDn.Term;
		}
		OptPriceU = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, &VolUp,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			0);
		OptPriceD = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, &VolDn,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			0);
		Vega = (OptPriceU - OptPriceD) / 2.0;
		ResultPrice[3] = Vega;

		if (Days_CalcToMaturity > 0)
		{
			long Days_CalcToMeanStart_Up = Days_CalcToMeanStart + 1;
			long Days_CalcToMaturity_Up = Days_CalcToMaturity + 1;
			long Days_CalcToMeanStart_Dn = Days_CalcToMeanStart - 1;
			long Days_CalcToMaturity_Dn = Days_CalcToMaturity - 1;
			double T_CalcToPay_Up = T_CalcToPay + dt;
			double T_CalcToPay_Dn = T_CalcToPay - dt;
			OptPriceU = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
				RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
				DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart_Up, Days_CalcToMaturity_Up, T_CalcToPay_Up, FixedRandn,
				0);
			OptPriceD = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
				RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
				DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart_Dn, Days_CalcToMaturity_Dn, T_CalcToPay_Dn, FixedRandn,
				0);
			Theta = -(OptPriceU - OptPriceD) / 2.0;
		}
		else
		{
			Theta = 0.0;
		}
		ResultPrice[4] = Theta;

		for (i = 0; i < DiscCurve->nterm(); i++)
		{
			DiscCurve->Rate[i] = DiscCurve->Rate[i] + 0.0001;
		}
		for (i = 0; i < RefCurve->nterm(); i++)
		{
			RefCurve->Rate[i] = RefCurve->Rate[i] + 0.0001;
		}
		OptPriceU = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			0);

		for (i = 0; i < DiscCurve->nterm(); i++)
		{
			DiscCurve->Rate[i] = DiscCurve->Rate[i] - 0.0002;
		}
		for (i = 0; i < RefCurve->nterm(); i++)
		{
			RefCurve->Rate[i] = RefCurve->Rate[i] - 0.0002;
		}
		OptPriceD = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, DiscCurve,
			RefCurve, DivCurve, QuantoCorr, FXVolCurve, VolMat,
			DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, FixedRandn,
			0);

		Rho = (OptPriceU - OptPriceD) / 2.0;
		for (i = 0; i < DiscCurve->nterm(); i++)
		{
			DiscCurve->Rate[i] = DiscCurve->Rate[i] + 0.0001;
		}
		for (i = 0; i < RefCurve->nterm(); i++)
		{
			RefCurve->Rate[i] = RefCurve->Rate[i] + 0.0001;
		}
		ResultPrice[5] = Rho *100.0;

		free(VolReshaped_Up);
		free(VolReshaped_Dn);
	}
	for (i = 0; i < NSimul; i++) free(FixedRandn[i]);
	free(FixedRandn);
	return 1;
}
double an(
	double dt,
	double Rf,
	double Div,
	double QuantoCorr,
	double Vol,
	double FXVol,
	long GreedIdx
)
{
	double value;
	double n = (double)GreedIdx;
	value = 0.5 * dt * ((Rf - Div - QuantoCorr * FXVol * Vol) * n - Vol * Vol * n * n);
	return value;
}

double bn(
	double dt,
	double Rf_Disc,
	double Vol,
	long GreedIdx
)
{
	double value;
	double n = (double)GreedIdx;
	value = 1.0 + dt * (Rf_Disc + Vol * Vol * n * n);
	return value;
}

double cn(
	double dt,
	double Rf,
	double Div,
	double QuantoCorr,
	double Vol,
	double FXVol,
	long GreedIdx
)
{
	double value;
	double n = (double)GreedIdx;
	value = 0.5 * dt * (-(Rf - Div - QuantoCorr * Vol * FXVol) * n - Vol * Vol * n * n);
	return value;
}

long AmericanOptionFDM(
	long TypeFlag,				// Option Type: 1이면 콜, 2이면 풋
	long AmericanOptionType,	// 0:아메리칸옵션 1:아메리칸디지털옵션
	double S0,					// 기초자산 현재가
	double X,					// 행사가격
	double MeanPrice,			// 평균가격
	curveinfo* DiscCurve,		// 할인커브
	curveinfo* RefCurve,		// 레퍼런스커브
	curveinfo* DivCurve,		// 배당커브
	double QuantoCorr,			// Quanto COrrelation
	curveinfo* FXVolCurve,		// FXVol 커브
	volinfo* VolMat,			// Vol 정보
	long DivTypeFlag,			// Div 배당타입
	long VolTypeFlag,
	long Days_CalcToMeanStart,
	long Days_CalcToMaturity,
	double T_CalcToPay,
	double* ResultPrice
)
{
	long i, j, k,N;
	long idx;
	double d;
	double dt = 1.0 / 365.0, T,T2;
	double dS = 0.01;
	double T_Ref = ((double)Days_CalcToMaturity) / 365.0;
	long idxT, idxS;
	double S, fxvol, r_ref, r_disc;
	long NGreed = 200;
	double vol, Rf_Quanto, r;
	double value, valueU, valueD;
	double* DailyForwardRate = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	double* DailyDiscRate = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		DailyForwardRate[i] = Calc_Forward_Rate(RefCurve->Term, RefCurve->Rate, RefCurve->nterm(), ((double)(i)) * dt, ((double)(i + 1)) * dt);
		DailyDiscRate[i] = Calc_Forward_Rate(DiscCurve->Term, DiscCurve->Rate, DiscCurve->nterm(), ((double)(i)) * dt, ((double)(i + 1)) * dt);
	}
	double* DailyDiv = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	if (DivTypeFlag == 0)
	{
		for (i = 0; i < Days_CalcToMaturity; i++)
		{
			DailyDiv[i] = DivCurve->Rate[0];
		}
	}
	else if (DivTypeFlag == 1)
	{
		for (i = 0; i < Days_CalcToMaturity; i++)
		{
			DailyDiv[i] = Calc_Forward_Rate(DivCurve->Term, DivCurve->Rate, DivCurve->nterm(), ((double)(i)) * dt, ((double)(i + 1)) * dt);
		}
	}
	else if (DivTypeFlag == 2)
	{
		N = 0;
		for (i = 0; i < DivCurve->nterm(); i++)
		{
			if (DivCurve->Term[i] > T_Ref) break;
			N += 1;
		}
		d = 1.0 / T_Ref * Sumation(N, DivCurve->Rate) / S0;
		for (i = 0 ; i < Days_CalcToMaturity; i++) DailyDiv[i] = d;
	}

	double* DailyFXVol = (double*)malloc(sizeof(double) * Days_CalcToMaturity);
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		DailyFXVol[i] = Interpolate_Linear(FXVolCurve->Term, FXVolCurve->Rate, FXVolCurve->nterm(), ((double)(i + 1)) * dt);
	}

	double a0, b0, c0, aN, bN, cN;
	double** An = make_array(Days_CalcToMaturity, NGreed);
	double** Bn = make_array(Days_CalcToMaturity, NGreed);
	double** Cn = make_array(Days_CalcToMaturity, NGreed);
	double* S_Array = make_array(NGreed);
	double* Temp_An = make_array(NGreed);
	double* Temp_Bn = make_array(NGreed);
	double* Temp_Cn = make_array(NGreed);
	for (i = 0; i < NGreed; i++) S_Array[i] = (double)i * dS;
	T = 0.0;
	idxT = 0;
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		T = (double)(i+1) * dt;
		idxS = 0;
		r_ref = DailyForwardRate[i];
		r_disc = DailyDiscRate[i];
		fxvol = DailyDiscRate[i];
		d = DailyDiv[i];

		for (j = 0; j < NGreed; j++)
		{
			S = (double)j * dS;
			if (VolTypeFlag == 2) vol = VolMat->Calc_Implied_Volatility(((double)Days_CalcToMaturity) / 365.0, S0 / X);
			else vol = VolMat->Calc_Volatility_for_Simulation(T, S, &idxT, &idxS);
			An[i][j] = an(dt, r_ref, d, QuantoCorr, vol, fxvol, j);
			Bn[i][j] = bn(dt, r_disc, vol, j);
			Cn[i][j] = cn(dt, r_ref, d, QuantoCorr, vol, fxvol, j);
		}
	}
	
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		a0 = An[i][0];
		b0 = Bn[i][0];
		c0 = Cn[i][0];
		aN = An[i][NGreed - 1];
		bN = Bn[i][NGreed - 1];
		cN = Cn[i][NGreed - 1];
		Bn[i][0] = 2.0 * a0 + b0;
		Cn[i][0] = c0 - a0;
		Bn[i][NGreed - 1] = bN + 2.0 * cN;
		An[i][NGreed - 1] = aN - cN;
	}
	
	double* Payoff = make_array(NGreed);
	if (AmericanOptionType == 0)
	{
		for (i = 0; i < NGreed; i++)
		{
			if (TypeFlag == 1) Payoff[i] = max(S_Array[i] * S0 - X, 0.0);
			else Payoff[i] = max(X - S_Array[i] * S0, 0.0);
		}
	}
	else
	{
		for (i = 0; i < NGreed; i++)
		{
			if (TypeFlag == 1)
			{
				Payoff[i] = 1.0 * (double)(S_Array[i] * S0 - X >= 0.0);
			}
			else
			{
				Payoff[i] = 1.0 * (double)(X - S_Array[i] * S0 >= 0.0);
			}
		}
	}

	for (i = Days_CalcToMaturity - 1; i >= 0; i--)
	{
		Tri_diagonal_Fast(An[i],    // Tridiagonal Matrix에서 왼쪽 대각선 행렬
			Bn[i],     // Tridiagonal Matrix에서 중앙 대각선 행렬
			Cn[i],    // Tridiagonal Matrix에서 오른쪽 대각선 행렬
			Payoff,     // V(T1)
			NGreed,           // 대각선 행렬들의 길이
			Temp_An,// 임시 행렬1 
			Temp_Bn, // 임시 행렬2
			Temp_Cn // 임시 행렬3
		);
		if (i >= Days_CalcToMeanStart)
		{
			if (AmericanOptionType == 0)
			{
				for (j = 0; j < NGreed; j++)
				{
					if (TypeFlag == 1) Payoff[j] = max(Payoff[j], max(S_Array[j] * S0 - X, 0.0));
					else Payoff[j] = max(Payoff[j], max(X - S_Array[j] * S0, 0.0));
				}
			}
			else
			{
				for (j = 0; j < NGreed; j++)
				{
					if (TypeFlag == 1) Payoff[j] = max(Payoff[j], 1.0 * (double)(S_Array[i] * S0 - X >= 0.0));
					else Payoff[j] = max(Payoff[j],1.0 * (double)(X - S_Array[i] * S0 >= 0.0));
				}
			}
		}
	}
	double frate;
	if (T_CalcToPay - T_Ref > 0.0) frate = Calc_Forward_Rate(DiscCurve->Term, DiscCurve->Rate, DiscCurve->nterm(), T_Ref, T_CalcToPay);
	else frate = 0.0;
	double DF_t_T = 1.0 / (1.0 + frate * (T_CalcToPay - T_Ref));
	value = Interpolate_Linear(S_Array, Payoff, NGreed, S0 / X) * DF_t_T;
	valueU = Interpolate_Linear(S_Array, Payoff, NGreed, S0 * 1.01 / X) * DF_t_T;
	valueD = Interpolate_Linear(S_Array, Payoff, NGreed, S0 * 0.99 / X) * DF_t_T;

	ResultPrice[0] = value;
	ResultPrice[1] = (valueU - valueD) / (2.0 * 0.01 * S0);
	ResultPrice[2] = (valueU + valueD - 2.0 * value) / (S0 * S0 * 0.01 * 0.01);
	free(DailyForwardRate);
	free(DailyDiscRate);
	free(DailyDiv);
	free(DailyFXVol);
	for (i = 0; i < Days_CalcToMaturity; i++)
	{
		free(An[i]);
		free(Bn[i]);
		free(Cn[i]);
	}
	free(An);
	free(Bn);
	free(Cn);
	free(Temp_An);
	free(Temp_Bn);
	free(Temp_Cn);
	free(S_Array);
	free(Payoff);
	return 1;
}

long AmericanOptionFDM(
	long TypeFlag,				// Option Type: 1이면 콜, 2이면 풋
	long AmericanOptionType,	// 0:아메리칸옵션 1:아메리칸디지털옵션
	double S0,					// 기초자산 현재가
	double X,					// 행사가격
	double MeanPrice,			// 평균가격
	curveinfo* DiscCurve,		// 할인커브
	curveinfo* RefCurve,		// 레퍼런스커브
	curveinfo* DivCurve,		// 배당커브
	double QuantoCorr,			// Quanto COrrelation
	curveinfo* FXVolCurve,		// FXVol 커브
	volinfo* VolMat,			// Vol 정보
	long DivTypeFlag,			// Div 배당타입
	long VolTypeFlag,
	long Days_CalcToMeanStart,
	long Days_CalcToMaturity,
	double T_CalcToPay,
	long GreekFlag,
	double* ResultPrice
)
{
	long i, j,k;
	long ResultCode;
	double valueU, valueD;
	double Vega, Theta, Rho;
	double TempResultPrice[10] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
	ResultCode = AmericanOptionFDM(
		TypeFlag, AmericanOptionType, S0, X, MeanPrice,
		DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
		VolMat, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
		T_CalcToPay, ResultPrice);

	if (GreekFlag > 0)
	{
		double* VolReshaped_Up = (double*)malloc(sizeof(double) * (VolMat->size(0) * VolMat->size(1)));
		double* VolReshaped_Dn = (double*)malloc(sizeof(double) * (VolMat->size(0) * VolMat->size(1)));
		k = 0;
		for (i = 0; i < VolMat->size(0); i++)
			for (j = 0; j < VolMat->size(1); j++)
			{
				VolReshaped_Up[k] = VolMat->Vol_Matrix[i][j] + 0.01;
				VolReshaped_Dn[k] = max(0.001, VolMat->Vol_Matrix[i][j] - 0.01);
				k = k + 1;
			}

		volinfo VolUp(VolMat->size(0), VolMat->Parity, VolMat->size(1), VolMat->Term, VolReshaped_Up);
		volinfo VolDn(VolMat->size(0), VolMat->Parity, VolMat->size(1), VolMat->Term, VolReshaped_Dn);
		if (VolTypeFlag == 0)
		{
			long NTempRate = 1;
			double TempRefRateTerm[1] = { 1.0 };
			double TempRefRate[1] = { 0.0 };
			if (DivTypeFlag == 2)
			{
				curveinfo TempCurve(NTempRate, TempRefRateTerm, TempRefRate);
				VolUp.set_localvol(&TempCurve, &TempCurve, 2.0, 0.01);
				fillna_Interpolate(VolUp.LocalVolMat, VolUp.size(0), VolUp.size(1));
				VolDn.set_localvol(&TempCurve, &TempCurve, 2.0, 0.01);
				fillna_Interpolate(VolDn.LocalVolMat, VolDn.size(0), VolDn.size(1));
			}
			else
			{
				VolUp.set_localvol(RefCurve, DivCurve, 2.0, 0.01);
				fillna_Interpolate(VolUp.LocalVolMat, VolUp.size(0), VolUp.size(1));
				VolDn.set_localvol(RefCurve, DivCurve, 2.0, 0.01);
				fillna_Interpolate(VolDn.LocalVolMat, VolDn.size(0), VolDn.size(1));
			}
		}
		else
		{
			VolUp.LocalVolMat = VolUp.Vol_Matrix;
			VolUp.Parity_Locvol = VolUp.Parity;
			VolUp.Term_Locvol = VolUp.Term;
			VolDn.LocalVolMat = VolDn.Vol_Matrix;
			VolDn.Parity_Locvol = VolDn.Parity;
			VolDn.Term_Locvol = VolDn.Term;
		}

		ResultCode = AmericanOptionFDM(
			TypeFlag, AmericanOptionType, S0, X, MeanPrice,
			DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
			&VolUp, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
			T_CalcToPay, TempResultPrice);
		valueU = TempResultPrice[0];
		ResultCode = AmericanOptionFDM(
			TypeFlag, AmericanOptionType, S0, X, MeanPrice,
			DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
			&VolDn, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
			T_CalcToPay, TempResultPrice);
		valueD = TempResultPrice[0];
		Vega = (valueU - valueD) / 2.0;
		ResultPrice[3] = Vega;
		if (Days_CalcToMaturity > 0)
		{
			long Days_CalcToMeanStart_Up = Days_CalcToMeanStart + 1;
			long Days_CalcToMaturity_Up = Days_CalcToMaturity + 1;
			long Days_CalcToMeanStart_Dn = Days_CalcToMeanStart - 1;
			long Days_CalcToMaturity_Dn = Days_CalcToMaturity - 1;
			double T_CalcToPay_Up = T_CalcToPay + 1.0 / 365.0;
			double T_CalcToPay_Dn = T_CalcToPay - 1.0 / 365.0;
			ResultCode = AmericanOptionFDM(
				TypeFlag, AmericanOptionType, S0, X, MeanPrice,
				DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
				VolMat, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart_Up, Days_CalcToMaturity_Up,
				T_CalcToPay_Up, TempResultPrice);
			valueU = TempResultPrice[0];
			ResultCode = AmericanOptionFDM(
				TypeFlag, AmericanOptionType, S0, X, MeanPrice,
				DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
				VolMat, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart_Dn, Days_CalcToMaturity_Dn,
				T_CalcToPay_Dn, TempResultPrice);
			valueD = TempResultPrice[0];
			Theta = (-(valueU - valueD)) / 2.0;
			for (i = 0; i < DiscCurve->nterm(); i++)
			{
				DiscCurve->Rate[i] = DiscCurve->Rate[i] + 0.0001;
			}
			for (i = 0; i < RefCurve->nterm(); i++)
			{
				RefCurve->Rate[i] = RefCurve->Rate[i] + 0.0001;
			}

			ResultCode = AmericanOptionFDM(
				TypeFlag, AmericanOptionType, S0, X, MeanPrice,
				DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
				VolMat, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
				T_CalcToPay, TempResultPrice);
			valueU = TempResultPrice[0];

			for (i = 0; i < DiscCurve->nterm(); i++)
			{
				DiscCurve->Rate[i] = DiscCurve->Rate[i] - 0.0002;
			}
			for (i = 0; i < RefCurve->nterm(); i++)
			{
				RefCurve->Rate[i] = RefCurve->Rate[i] - 0.0002;
			}
			ResultCode = AmericanOptionFDM(
				TypeFlag, AmericanOptionType, S0, X, MeanPrice,
				DiscCurve, RefCurve, DivCurve, QuantoCorr, FXVolCurve,
				VolMat, DivTypeFlag, VolTypeFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
				T_CalcToPay, TempResultPrice);
			valueD = TempResultPrice[0];
			Rho = (valueU - valueD) / 2.0 *100.0;
			for (i = 0; i < DiscCurve->nterm(); i++)
			{
				DiscCurve->Rate[i] = DiscCurve->Rate[i] + 0.0001;
			}
			for (i = 0; i < RefCurve->nterm(); i++)
			{
				RefCurve->Rate[i] = RefCurve->Rate[i] + 0.0001;
			}
		}
		else
		{
			Theta = 0.0;
		}
		ResultPrice[4] = Theta;
		ResultPrice[5] = Rho;
		free(VolReshaped_Up);
		free(VolReshaped_Dn);
	}
	return ResultCode;
}

long PricingOption(
	long LongShort,
	long CallPut,
	long OptionType,
	long CalcDate,
	long MeanStartDate,

	long MaturityDate,
	long PayDate,
	long NRateDisc,
	double* RateTermDisc,
	double* RateDisc,

	long NRateRef,
	double* RateTermRef,
	double* RateRef,
	double StockPrice,
	double StrikePrice,

	double MeanPrice,
	long DivType,
	long NDiv,
	double* DivTerm,
	double* DivRate,

	double QuantoCorr,
	long NFXVol,
	double* FXVolTerm,
	double* FXVol,
	long ImvolLocalVolFlag,

	long NParity,
	double* Parity,
	long NTermVol,
	double* TermVol,
	double* VolReshaped,

	long GreekFlag,
	double* ResultPrice
)
{
	long ResultCode = 0;
	long i, j,k;
	long TypeFlag;
	double S0 = StockPrice;
	double X = StrikePrice;

	curveinfo DiscCurve(NRateDisc, RateTermDisc, RateDisc);
	curveinfo RefCurve(NRateRef, RateTermRef, RateRef);
	curveinfo DivCurve(NDiv, DivTerm, DivRate);
	curveinfo FXVolCurve(NFXVol, FXVolTerm, FXVol);
	volinfo VolMat(NParity, Parity, NTermVol, TermVol, VolReshaped);
	if (ImvolLocalVolFlag == 0)
	{
		long NTempRate = 1;
		double TempRefRateTerm[1] = { 1.0 };
		double TempRefRate[1] = { 0.0 };
		if (DivType == 2)
		{
			curveinfo TempCurve(NTempRate, TempRefRateTerm, TempRefRate);
			VolMat.set_localvol(&TempCurve, &TempCurve, 2.0, 0.01);
		}
		else
		{
			VolMat.set_localvol(&RefCurve, &DivCurve, 2.0, 0.01);
		}
		fillna_Interpolate(VolMat.LocalVolMat, VolMat.size(0), VolMat.size(1));
	}
	else
	{
		VolMat.LocalVolMat = VolMat.Vol_Matrix;
		VolMat.Term_Locvol = VolMat.Term;
		VolMat.Parity_Locvol = VolMat.Parity;
	}

	long Days_CalcToMeanStart = DayCountAtoB(CalcDate, MeanStartDate);
	long Days_CalcToMaturity = DayCountAtoB(CalcDate, MaturityDate);
	long Days_CalcToPay = DayCountAtoB(CalcDate, PayDate);

	double T_CalcToMeanStart = ((double)Days_CalcToMeanStart) / 365.0;
	double T_CalcToMaturity = ((double)Days_CalcToMaturity) / 365.0;
	double T_CalcToPay = ((double)Days_CalcToPay) / 365.0;

	if (CallPut == 0) TypeFlag = 1;
	else TypeFlag = 2;

	if (OptionType == 0)
	{
		ResultCode = BSOption(TypeFlag,S0,X,&DiscCurve,	&RefCurve,
			&DivCurve,	QuantoCorr,	&FXVolCurve,&VolMat,	DivType,	
			T_CalcToMaturity,T_CalcToPay,ResultPrice);
	}
	else if (OptionType == 1)
	{
		ResultCode = BSDigitalOption(TypeFlag, S0, X, &DiscCurve, &RefCurve,
			&DivCurve, QuantoCorr, &FXVolCurve, &VolMat, DivType,
			T_CalcToMaturity, T_CalcToPay, ResultPrice);
	}
	else if (OptionType == 2)
	{
		ResultCode = AsianOptionSimul(TypeFlag, S0, X, MeanPrice, &DiscCurve,
			&RefCurve, &DivCurve, QuantoCorr, &FXVolCurve, &VolMat,
			DivType, ImvolLocalVolFlag, Days_CalcToMeanStart, Days_CalcToMaturity, T_CalcToPay, GreekFlag, ResultPrice);
	}
	else if (OptionType == 3 || OptionType == 4)
	{
		long AmericanOptionType = 0;
		if (OptionType == 3) AmericanOptionType = 0;
		else AmericanOptionType = 1;

		ResultCode = AmericanOptionFDM(TypeFlag, AmericanOptionType, S0, X, MeanPrice,
			&DiscCurve, &RefCurve, &DivCurve, QuantoCorr, &FXVolCurve,
			&VolMat, DivType, ImvolLocalVolFlag, Days_CalcToMeanStart, Days_CalcToMaturity,
			T_CalcToPay, GreekFlag, ResultPrice);
	}
	k = 0;
	for (i = 0 ; i < NParity; i++)
		for (j = 0; j < NTermVol; j++)
		{
			VolReshaped[k] = VolMat.LocalVolMat[i][j];
			k += 1;
		}
	if (LongShort != 0)
	{
		for (i = 0; i < 6; i++) ResultPrice[i] *= -1.0;
	}
	return 1;
}

long ErrorCheck_OPTIONPRICING_Excel(
	long LongShort,
	long CallPut,
	long OptionType,
	long CalcDateExl,
	long MeanStartDateExl,

	long MaturityDateExl,
	long PayDateExl,
	long NRateDisc,
	double* RateTermDisc,
	double* RateDisc,

	long NRateRef,
	double* RateTermRef,
	double* RateRef,
	double StockPrice,
	double StrikePrice,

	double MeanPrice,
	long DivType,
	long NDiv,
	double* DivTerm,
	double* DivRate,

	double QuantoCorr,
	long NFXVol,
	double* FXVolTerm,
	double* FXVol,
	long ImvolLocalVolFlag,

	long NParity,
	double* Parity,
	long NTermVol,
	double* TermVol,
	double* VolReshaped,

	long GreekFlag,
	double* ResultPrice,
	char* Error
)
{
	long i, j;
	char ErrorName[100];

	if (LongShort != 0 && LongShort != 1)
	{
		strcpy_s(ErrorName, "Check LongShortFlag");
		return SaveErrorName(Error, ErrorName);
	}

	if (CallPut != 0 && CallPut != 1)
	{
		strcpy_s(ErrorName, "Check CallPutFlag");
		return SaveErrorName(Error, ErrorName);
	}

	if (OptionType < 0 || OptionType > 4)
	{
		strcpy_s(ErrorName, "Check OptionType");
		return SaveErrorName(Error, ErrorName);
	}

	if (CalcDateExl < 0)
	{
		strcpy_s(ErrorName, "PricingDate is too early");
		return SaveErrorName(Error, ErrorName);
	}


	if (CalcDateExl > PayDateExl)
	{
		strcpy_s(ErrorName, "PricingDate is too late");
		return SaveErrorName(Error, ErrorName);
	}

	if (MaturityDateExl > PayDateExl)
	{
		strcpy_s(ErrorName, "MaturityDate is too late");
		return SaveErrorName(Error, ErrorName);
	}

	if (MeanStartDateExl > MaturityDateExl)
	{
		strcpy_s(ErrorName, "MeanStartDate is too late");
		return SaveErrorName(Error, ErrorName);
	}

	if (NRateDisc < 1)
	{
		strcpy_s(ErrorName, "Check Number Disc Rate");
		return SaveErrorName(Error, ErrorName);
	}
	for (i = 0; i < NRateDisc; i++)
	{
		if (i > 0)
		{
			if (RateTermDisc[i] < RateTermDisc[i - 1])
			{
				strcpy_s(ErrorName, "Sort Disc Rate Term");
				return SaveErrorName(Error, ErrorName);
			}
		}
	}

	if (NRateRef < 1)
	{
		strcpy_s(ErrorName, "Check Number Ref Rate");
		return SaveErrorName(Error, ErrorName);
	}

	for (i = 0; i < NRateRef; i++)
	{
		if (i > 0)
		{
			if (RateTermRef[i] < RateTermRef[i - 1])
			{
				strcpy_s(ErrorName, "Sort Ref Rate Term");
				return SaveErrorName(Error, ErrorName);
			}
		}
	}

	if (NDiv < 1)
	{
		strcpy_s(ErrorName, "Check Number Div Rate");
		return SaveErrorName(Error, ErrorName);
	}
	for (i = 0; i < NDiv; i++)
	{
		if (i > 0)
		{
			if (DivTerm[i] < DivTerm[i - 1])
			{
				strcpy_s(ErrorName, "Sort Div Term");
				return SaveErrorName(Error, ErrorName);
			}
		}	
	}

	if (DivType < 0 || DivType > 2)
	{
		strcpy_s(ErrorName, "Check Dividend Type");
		return SaveErrorName(Error, ErrorName);
	}

	if (NFXVol < 1)
	{
		strcpy_s(ErrorName, "Check Number FXVol Rate");
		return SaveErrorName(Error, ErrorName);
	}
	for (i = 0; i < NFXVol; i++)
	{
		if (i > 0)
		{
			if (FXVolTerm[i] < FXVolTerm[i - 1])
			{
				strcpy_s(ErrorName, "Sort FXVol Term");
				return SaveErrorName(Error, ErrorName);
			}
		}
	}

	if (ImvolLocalVolFlag < 0 || ImvolLocalVolFlag > 2)
	{
		strcpy_s(ErrorName, "Check Vol Type");
		return SaveErrorName(Error, ErrorName);
	}

	if (NParity < 1)
	{
		strcpy_s(ErrorName, "Check Vol Parity Number");
		return SaveErrorName(Error, ErrorName);
	}

	if (NTermVol < 1)
	{
		strcpy_s(ErrorName, "Check Vol Term Number");
		return SaveErrorName(Error, ErrorName);
	}

	for (i = 0; i < NParity; i++)
	{
		if (i > 0)
		{
			if (Parity[i] < Parity[i - 1])
			{
				strcpy_s(ErrorName, "Sort Parity");
				return SaveErrorName(Error, ErrorName);
			}
		}
	}
	for (i = 0; i < NTermVol; i++)
	{
		if (i > 0)
		{
			if (TermVol[i] < TermVol[i - 1])
			{
				strcpy_s(ErrorName, "Sort Vol Term");
				return SaveErrorName(Error, ErrorName);
			}
		}
	}

	long k = 0;
	for (i = 0 ; i < NParity; i++)
		for (j = 0; j < NTermVol; j++)
		{
			if (VolReshaped[k] < 0.0)
			{
				strcpy_s(ErrorName, "Vol Must be Positive");
				return SaveErrorName(Error, ErrorName);
			}
			k += 1;
		}
	return 1;
}

DLLEXPORT(long) OPTIONPRICING_Excel(
	long LongShort,
	long CallPut,
	long OptionType,
	long CalcDateExl,
	long MeanStartDateExl,

	long MaturityDateExl,
	long PayDateExl,
	long NRateDisc,
	double *RateTermDisc,
	double *RateDisc,

	long NRateRef,
	double* RateTermRef,
	double* RateRef,
	double StockPrice,
	double StrikePrice,

	double MeanPrice,
	long DivType,
	long NDiv,
	double* DivTerm,
	double* DivRate,
	
	double QuantoCorr,
	long NFXVol,
	double *FXVolTerm,
	double *FXVol,
	long ImvolLocalVolFlag,
	
	long NParity,
	double* Parity,
	long NTermVol,
	double* TermVol,
	double* VolReshaped,
	
	long GreekFlag,
	double* ResultPrice,
	char* Error
)
{
	long i, j;
	long ResultCode;
	long ErrorCode;

	ErrorCode = ErrorCheck_OPTIONPRICING_Excel(
		LongShort,CallPut,OptionType,CalcDateExl,MeanStartDateExl,
		MaturityDateExl,PayDateExl,NRateDisc,RateTermDisc,RateDisc,
		NRateRef,RateTermRef,RateRef,StockPrice,StrikePrice,
		MeanPrice,DivType,NDiv,DivTerm,DivRate,
		QuantoCorr,NFXVol,FXVolTerm,FXVol,ImvolLocalVolFlag,
		NParity,Parity,NTermVol,TermVol,VolReshaped,
		GreekFlag,ResultPrice,Error);
	if (ErrorCode < 0) return ErrorCode;

	long CalcDate = ExcelDateToCDate(CalcDateExl);
	long MeanStartDate = ExcelDateToCDate(MeanStartDateExl);
	long MaturityDate = ExcelDateToCDate(MaturityDateExl);
	long PayDate = ExcelDateToCDate(PayDateExl);

	ResultCode = PricingOption(LongShort, CallPut, OptionType, CalcDate, MeanStartDate,
		MaturityDate, PayDate, NRateDisc, RateTermDisc, RateDisc,
		NRateRef, RateTermRef, RateRef, StockPrice, StrikePrice,
		MeanPrice, DivType, NDiv, DivTerm, DivRate,
		QuantoCorr, NFXVol, FXVolTerm, FXVol, ImvolLocalVolFlag,
		NParity, Parity, NTermVol, TermVol, VolReshaped,
		GreekFlag, ResultPrice);

	return 1;
}