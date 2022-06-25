#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <crtdbg.h>

#include "IRStructure.h"

#ifndef DLLEXPORT(A)
#ifdef WIN32
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
#elif _WIN64
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
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

#ifndef strcpy_s
#include <string.h>
#endif


long SaveErrorName2(char* Error, char ErrorName[100])
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

double Day_to_T(long Day, long DaysOneYear)
{
	return ((double)Day) / ((double)DaysOneYear);
}

double CalcZeroRateFromSwapRate(double Spread, long NResetDate, long* ResetDays, double* DeltaT, long Convention1Y, long NRate, double* RateTerm, double* Rate, double *Time)
{
	long i;
	double a, b;
	double T;
	double ZeroRate;
	double PV;
	double r;

	a = 1.0;
	for (i = 0; i < NResetDate - 1; i++)
	{
		T = Day_to_T(ResetDays[i], Convention1Y);
		ZeroRate = CubicSpline(NRate, RateTerm, Rate, T);
		a -= Spread * DeltaT[i] * exp(-ZeroRate * T);
	}

	b = 1.0 + Spread * DeltaT[NResetDate - 1];
	PV = a / b;
	T = Day_to_T(ResetDays[NResetDate - 1], Convention1Y);

	r = -1.0 / T * log(PV);
	*Time = T;
	return r;
}

double CalcZeroRateFromCouponBond(double NA, double BondPrice, double CPNRate, double Maturity, long NCPNDate, long N_Annual_Cpn,long NRate, double* RateTerm, double* Rate, double* Time)
{
	long i;
	double r;
	double a = BondPrice ;
	double ZeroRate;
	double T;
	
	for (i = 0; i < NCPNDate - 1; i++)
	{
		T = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
		ZeroRate = CubicSpline(NRate, RateTerm, Rate, T);
		a -= NA * CPNRate * 1.0/(double)N_Annual_Cpn * exp(-ZeroRate * T);
	}
	double b = NA + NA * CPNRate * 1.0 / (double)N_Annual_Cpn;
	double DF = a / b;

	r = -1.0 / Maturity * log(DF);
	*Time = Maturity;
	return r;
}

long ZeroCurveGenerate_Bond(
	long PriceDate_C,

	long NShortTerm,
	long* ShortTermType,
	long* ShortTerm_Maturity_C,
	double* ShortTerm_Rate,
	long* ShortTerm_OneYConvention,

	long NSwap,
	long* SwapMaturity_C,
	long* NCPN_Ann,
	double* SwapRate,
	long* Swap_OneYConvention,

	long NResultCurve,
	double* ResultCurveTerm,
	double* ResultCurve
)
{
	long i;
	long N;

	double r; 
	
	long NCPNBond =NSwap;
	long* BondMaturityDate = SwapMaturity_C;
	double* YTM = SwapRate;
	long* Bond_OneYConvention = Swap_OneYConvention;

	long NZero = (NShortTerm + NSwap);
	double* ZeroTerm = (double*)malloc(sizeof(double) * NZero);
	double* ZeroRate = (double*)malloc(sizeof(double) * NZero);
	N = NShortTerm;
	for (i = 0; i < N; i++)
	{
		ZeroTerm[i] = ((double)DayCountAtoB(PriceDate_C,ShortTerm_Maturity_C[i]))/(double)ShortTerm_OneYConvention[i];
		ZeroRate[i] = ShortTerm_Rate[i];
	}

	CPN_Bond* Bond_Infos = new CPN_Bond[NCPNBond];
	for (i = 0; i < NCPNBond; i++)
	{
		(Bond_Infos + i)->initialize(PriceDate_C, BondMaturityDate[i], NCPN_Ann[i], YTM[i], YTM[i], Bond_OneYConvention[i]);
	}

	for (i = 0; i < NCPNBond; i++)
	{
		r = CalcZeroRateFromCouponBond((Bond_Infos + i)->NA, (Bond_Infos + i)->PV, (Bond_Infos + i)->YTM, (Bond_Infos + i)->Maturity, (Bond_Infos + i)->NCPNTime,
										(Bond_Infos + i)->N_Annual_Cpn, N, ZeroTerm, ZeroRate, ZeroTerm + i + NShortTerm);
		ZeroRate[i + NShortTerm] = r;
		N += 1;
	}

	for (i = 0; i < N; i++)
	{
		ResultCurveTerm[i] = ZeroTerm[i];
		ResultCurve[i] = ZeroRate[i];
	}
	
	delete[] Bond_Infos;
	free(ZeroTerm);
	free(ZeroRate);
	return 1;
}

long ZeroCurveGenerate_IRS(
	long PriceDate_C,

	long NShortTerm,
	long* ShortTermType,
	long* ShortTerm_Maturity_C,
	double* ShortTerm_Rate,
	long* ShortTerm_OneYConvention,

	long NForward,
	long* ForwardStart_C,
	long* ForwardEnd_C,
	double* ForwardRate,
	long* Forward_OneYConvention,

	long NSwap,
	long* SwapMaturity_C,
	long* NCPN_Ann,
	double* SwapRate,
	long* Swap_OneYConvention,

	long NResultCurve,
	double* ResultCurveTerm,
	double* ResultCurve
)
{
	long i;
	long j;
	long k;
	long N;
	long EndFlag = 0;

	long PricingDate = 20211228;
	long D1;
	long D2;
	double r;
	double T;
	double T1;
	double T2;
	double DF_0t, DF_tT, DF_0T;
	long ForwardUsingFlag = 0;

	MoneyMarket* MMF_Infos = new MoneyMarket[NShortTerm];
	for (i = 0; i < NShortTerm; i++)
	{
		(MMF_Infos + i)->initialize(ShortTermType[i], PriceDate_C, ShortTerm_Maturity_C[i], ShortTerm_Rate[i], ShortTerm_OneYConvention[i]);
	}

	IR_Forward* IRForward_Infos;
	if (NForward > 1) IRForward_Infos = new IR_Forward[NForward + 1];
	else IRForward_Infos = new IR_Forward[NForward];

	for (i = 0; i < NForward; i++) (IRForward_Infos + i)->initialize(PricingDate, ForwardStart_C[i], ForwardEnd_C[i], ForwardRate[i], Forward_OneYConvention[i]);

	if (NForward > 1)
	{
		DF_tT = 1.0;
		for (i = 0; i < NForward; i++)
		{
			if (i < NForward - 1)
			{
				T1 = ((double)DayCountAtoB(PricingDate, ForwardStart_C[i]))/365.0;
				T2 = ((double)DayCountAtoB(PricingDate, ForwardStart_C[i+1]))/365.0;
				DF_tT *= (1.0 / (1.0 + ForwardRate[i] * (T2 - T1)));
			}
			else
			{
				T1 = ((double)DayCountAtoB(PricingDate, ForwardStart_C[i])) / 365.0;
				T2 = ((double)DayCountAtoB(PricingDate, ForwardEnd_C[i])) / 365.0;
				DF_tT *= (1.0 / (1.0 + ForwardRate[i] * (T2 - T1)));
			}
		}
		T =( (double)DayCountAtoB(ForwardStart_C[0], ForwardEnd_C[NForward - 1]))/365.0;
		r = 1.0 / T * (1.0 / DF_tT - 1.0);

		(IRForward_Infos + NForward)->initialize(PricingDate, ForwardStart_C[0], ForwardEnd_C[NForward -1], r, 365);
		NForward += 1;
	}

	for (i = 0 ; i < NShortTerm; i++)
		for (j = 0; j < NForward; j++)
		{
			D1 = (MMF_Infos + i)->Maturity_Days;
			D2 = (IRForward_Infos + j)->forward_start_day;
			if (D1 >= D2)
			{
				ForwardUsingFlag = 1;
				break;
			}
		}
	long nCurve1Under = NShortTerm;
	double* CurveTerm1Under;
	double* Curve1Under;

	if (ForwardUsingFlag == 1)
	{
		k = 0;
		for (i = 0; i < NShortTerm; i++)
		{
			for (j = 0; j < NForward; j++)
			{
				D1 = (MMF_Infos + i)->Maturity_Days;
				D2 = (IRForward_Infos + j)->forward_start_day;
				if (D1 >= D2)
				{
					k += 1;
					if (j == NForward - 1)
					{
						EndFlag = 1;
						break;
					}
				}
			}

			if (EndFlag == 1) break;
		}

		nCurve1Under += k;
		CurveTerm1Under = (double*)malloc(sizeof(double) * nCurve1Under);
		Curve1Under = (double*)malloc(sizeof(double) * nCurve1Under);

		EndFlag = 0;
		k = 0;
		for (i = 0; i < NShortTerm; i++)
		{
			CurveTerm1Under[i] = Day_to_T((MMF_Infos + i)->Maturity_Days, (MMF_Infos + i)->OneYearConvention);
			Curve1Under[i] = (MMF_Infos + i)->ZeroRate;
			for (j = 0; j < NForward; j++)
			{
				D1 = (MMF_Infos + i)->Maturity_Days;
				D2 = (IRForward_Infos + j)->forward_start_day;
				if (D1 >= D2)
				{
					ForwardUsingFlag = max(1, ForwardUsingFlag);
					T = Day_to_T((IRForward_Infos + j)->forward_end_day, (IRForward_Infos + j)->OneYearConvention);
					DF_0t = (MMF_Infos + i)->Discount_Factor;
					DF_tT = (1.0 / (1.0 + (double)( (IRForward_Infos + j)->forward_end_day - D1) / 365.0 *(IRForward_Infos + j)->forward_rate));
					DF_0T = DF_0t * DF_tT;
					r = -1.0 / T * log(DF_0T);
					CurveTerm1Under[k + NShortTerm] = T;
					Curve1Under[k + NShortTerm] = r;
					k += 1;
					if (j == NForward - 1)
					{
						EndFlag = 1;
						break;
					}
				}
			}

			if (EndFlag == 1) break;
		}
	}
	else
	{
		nCurve1Under = NShortTerm + NForward;
		CurveTerm1Under = (double*)malloc(sizeof(double) * nCurve1Under);
		Curve1Under = (double*)malloc(sizeof(double) * nCurve1Under);
		for (i = 0; i < NShortTerm; i++)
		{
			CurveTerm1Under[i] = Day_to_T((MMF_Infos + i)->Maturity_Days, (MMF_Infos + i)->OneYearConvention);
			Curve1Under[i] = (MMF_Infos + i)->ZeroRate;
		}
		r = (MMF_Infos + NShortTerm - 1)->ZeroRate;
		for (i = 0; i < NForward; i++)
		{
			T = Day_to_T((IRForward_Infos + i)->forward_end_day, (IRForward_Infos + i)->OneYearConvention);
			T1 = Day_to_T((IRForward_Infos + i)->forward_start_day, (IRForward_Infos + i)->OneYearConvention);
			DF_0t = exp(-r * T1);
			DF_tT = (1.0 / (1.0 + (double)((IRForward_Infos + i)->forward_end_day - (IRForward_Infos + i)->forward_start_day) / (double)((IRForward_Infos + i)->OneYearConvention) * (IRForward_Infos + i)->forward_rate));
			DF_0T = DF_0t * DF_tT;
			r = -1.0 / T * log(DF_0T);
			CurveTerm1Under[NShortTerm + i] = T;
			Curve1Under[NShortTerm + i] = r;
		}
	}
	
	IR_Swap* Swap_Infos = new IR_Swap[NSwap];
	for (i = 0; i < NSwap; i++)
	{
		(Swap_Infos + i)->initialize(PriceDate_C, SwapMaturity_C[i], NCPN_Ann[i], SwapRate[i], Swap_OneYConvention[i]);
	}

	long NZero = (nCurve1Under + NSwap);
	double* ZeroTerm = (double*)malloc(sizeof(double) * NZero);
	double* ZeroRate = (double*)malloc(sizeof(double) * NZero);
	N = nCurve1Under;
	for (i = 0; i < N; i++)
	{
		ZeroTerm[i] = CurveTerm1Under[i];
		ZeroRate[i] = Curve1Under[i];
	}

	for (i = 0; i < NSwap; i++)
	{
		r = CalcZeroRateFromSwapRate((Swap_Infos + i)->SwapRate, (Swap_Infos + i)->NResetTime, (Swap_Infos + i)->ResetDays,  (Swap_Infos + i)->DeltaT, (Swap_Infos + i)->OneYearConvention, N, ZeroTerm, ZeroRate, ZeroTerm + i + nCurve1Under);
		ZeroRate[i + nCurve1Under] = r;
		N += 1;
	}

	for (i = 0; i < N; i++)
	{
		ResultCurveTerm[i] = ZeroTerm[i];
		ResultCurve[i] = ZeroRate[i];
	}

	delete[] MMF_Infos;
	delete[] IRForward_Infos;
	delete[] Swap_Infos;
	free(CurveTerm1Under);
	free(Curve1Under);
	return 1;
}

long ErrorCheckZeroMakerExcel(
	long PriceDate_Excel,
	long MarketDataFlag,

	long NShortTerm,
	long* ShortTermType,
	long* ShortTerm_Maturity_Excel,
	double* ShortTerm_Rate,
	long* ShortTerm_OneYConvention,

	long NForward,
	long* ForwardStart_Excel,
	long* ForwardEnd_Excel,
	double* ForwardRate,
	long* Forward_OneYConvention,

	long NSwap,
	long* SwapMaturity_Excel,
	long* NCPN_Ann,
	double* SwapRate,
	long* Swap_OneYConvention,

	char* Error
)
{
	long i;
	long j;
	char ErrorName[100];

	if (PriceDate_Excel < 0 || PriceDate_Excel >= ShortTerm_Maturity_Excel[0] || PriceDate_Excel >= ForwardStart_Excel[0] || PriceDate_Excel >= SwapMaturity_Excel[0])
	{
		strcpy_s(ErrorName, "Check PriceDate");
		return SaveErrorName2(Error, ErrorName);
	}

	if (NShortTerm < 1)
	{
		strcpy_s(ErrorName, "Check Short Term Rate Number");
		return SaveErrorName2(Error, ErrorName);
	}

	for (i = 0; i < NShortTerm; i++)
	{
		if (i > 0)
		{
			if (ShortTerm_Maturity_Excel[i] < ShortTerm_Maturity_Excel[i - 1])
			{
				strcpy_s(ErrorName, "Sort Short Term Maturity");
				return SaveErrorName2(Error, ErrorName);
			}
		}

		if (ShortTermType[i] < 0 || ShortTermType[i] > 1)
		{
			strcpy_s(ErrorName, "Check Short Term Product Type");
			return SaveErrorName2(Error, ErrorName);
		}

		if (ShortTerm_Maturity_Excel[i] - PriceDate_Excel >= 365)
		{
			strcpy_s(ErrorName, "Check Short Term Rate Maturity");
			return SaveErrorName2(Error, ErrorName);
		}

		if (ShortTerm_OneYConvention[i] != 360 && ShortTerm_OneYConvention[i] != 365)
		{
			strcpy_s(ErrorName, "Check Short Term Rate Y Convention");
			return SaveErrorName2(Error, ErrorName);
		}
	}

	if (NForward < 1)
	{
		strcpy_s(ErrorName, "Check Forward Number");
		return SaveErrorName2(Error, ErrorName);
	}

	for (i = 0; i < NForward; i++)
	{
		if (i > 0)
		{
			if (ForwardStart_Excel[i] < ForwardStart_Excel[i - 1])
			{
				strcpy_s(ErrorName, "Sort Forward Start Date");
				return SaveErrorName2(Error, ErrorName);
			}
		}
		if (ForwardStart_Excel[i] - PriceDate_Excel > 366 || ForwardEnd_Excel[i] - PriceDate_Excel > 366)
		{
			strcpy_s(ErrorName, "Check IR Forward Maturity");
			return SaveErrorName2(Error, ErrorName);
		}

		if (Forward_OneYConvention[i] != 360 && Forward_OneYConvention[i] != 365)
		{
			strcpy_s(ErrorName, "Check Forward Y Convention");
			return SaveErrorName2(Error, ErrorName);
		}
	}

	if (NSwap < 1)
	{
		strcpy_s(ErrorName, "Check Swap Number");
		return SaveErrorName2(Error, ErrorName);
	}

	for (i = 0; i < NSwap; i++)
	{
		if (i > 0)
		{
			if (SwapMaturity_Excel[i] < SwapMaturity_Excel[i - 1])
			{
				strcpy_s(ErrorName, "Sort Swap Maturity Date");
				return SaveErrorName2(Error, ErrorName);
			}
		}

		if (NCPN_Ann[i] != 1 && NCPN_Ann[i] != 2 && NCPN_Ann[i] != 3 && NCPN_Ann[i] != 4 && NCPN_Ann[i] != 6 && NCPN_Ann[i] != 12)
		{
			strcpy_s(ErrorName, "Check Swap NCPN Ann");
			return SaveErrorName2(Error, ErrorName);
		}

		if (Swap_OneYConvention[i] != 360 && Swap_OneYConvention[i] != 365)
		{
			strcpy_s(ErrorName, "Check Forward Y Convention");
			return SaveErrorName2(Error, ErrorName);
		}

	}

	return 1;
}

DLLEXPORT(long) ZeroMakerExcel(
	long PriceDate_Excel,			// Pricing Date ExcelType
	long MarketDataFlag,			// 0: MMF, IRForward, Swap 사용 1: 단기채YTM, 쿠폰채YTM사용

	long NShortTerm,				// Short Term 이자율 개수
	long* ShortTermType,			// 상품타입 0:OverNight, 1: CD, Libor 등
	long* ShortTerm_Maturity_Excel,	// Short Term 이자율 만기일 ExcelType Array
	double* ShortTerm_Rate,			// Short Term 이자율 Rate Array
	long* ShortTerm_OneYConvention,	// DayCount 1Y (360 or 365)

	long NForward,					// Forward 이자율 개수
	long* ForwardStart_Excel,		// Forward금리 추정시작일 ExcelType Array
	long* ForwardEnd_Excel,			// Forward금리 추정종료일 ExcelType Array
	double* ForwardRate,			// Forward Rate Array
	long* Forward_OneYConvention,	// Forward Rate DayCount 1Y (360 or 365) Array

	long NSwap,						// 입력된 Swap 개수
	long* SwapMaturity_Excel,		// Swap 만기 ExcelType Array
	long* NCPN_Ann,					// Swap 연 이자지급 회수 Array
	double* SwapRate,				// Swap Rate Array
	long* Swap_OneYConvention,		// Swap DayCount 1Y (360 or 365)

	long NResultCurve,				// Out: 결과 커브 Term 개수
	double* ResultCurveTerm,		// Out: 결과 커브 Term Array
	double* ResultCurve,			// Out: 결과 커브 Rate Array
	char* Error						// Out: 에러메시지
)
{
	long i;
	long ResultCode;

	ResultCode = ErrorCheckZeroMakerExcel(PriceDate_Excel, MarketDataFlag, NShortTerm, ShortTermType, ShortTerm_Maturity_Excel,
		ShortTerm_Rate, ShortTerm_OneYConvention, NForward, ForwardStart_Excel, ForwardEnd_Excel, ForwardRate, Forward_OneYConvention,
		NSwap, SwapMaturity_Excel, NCPN_Ann, SwapRate, Swap_OneYConvention, Error);

	if (ResultCode < 0) return ResultCode;

	long PriceDate_C = ExcelDateToCDate(PriceDate_Excel);
	long* ShortTerm_Maturity_C = (long*)malloc(sizeof(long) * NShortTerm);
	long* ForwardStart_C = (long*)malloc(sizeof(long) * NForward);
	long* ForwardEnd_C = (long*)malloc(sizeof(long) * NForward);
	long* SwapMaturity_C = (long*)malloc(sizeof(long) * NSwap);
	for (i = 0; i < NShortTerm; i++)
	{
		ShortTerm_Maturity_C[i] = ExcelDateToCDate(ShortTerm_Maturity_Excel[i]);
	}

	for (i = 0; i < NForward; i++)
	{
		ForwardStart_C[i] = ExcelDateToCDate(ForwardStart_Excel[i]);
		ForwardEnd_C[i] = ExcelDateToCDate(ForwardEnd_Excel[i]);
	}

	for (i = 0; i < NSwap; i++)
	{
		SwapMaturity_C[i] = ExcelDateToCDate(SwapMaturity_Excel[i]);
	}

	for (i = 0; i < NForward; i++)
	{
		if (ForwardEnd_C[i] >= SwapMaturity_C[0]) NForward -= 1;
	}
	if (MarketDataFlag == 0)
	{
		ResultCode = ZeroCurveGenerate_IRS(
			PriceDate_C, NShortTerm, ShortTermType, ShortTerm_Maturity_C, ShortTerm_Rate, ShortTerm_OneYConvention,
			NForward, ForwardStart_C, ForwardEnd_C, ForwardRate, Forward_OneYConvention,
			NSwap, SwapMaturity_C, NCPN_Ann, SwapRate, Swap_OneYConvention,
			NResultCurve, ResultCurveTerm, ResultCurve);
	}
	else
	{
		ResultCode = ZeroCurveGenerate_Bond(
			PriceDate_C, NShortTerm, ShortTermType, ShortTerm_Maturity_C, ShortTerm_Rate, ShortTerm_OneYConvention,
			NSwap, SwapMaturity_C, NCPN_Ann, SwapRate, Swap_OneYConvention,
			NResultCurve, ResultCurveTerm, ResultCurve);
	}

	free(ShortTerm_Maturity_C);
	free(ForwardStart_C);
	free(ForwardEnd_C);
	free(SwapMaturity_C);
	return ResultCode;
}

DLLEXPORT(double) Calc_IRS(long PricingDate_Excel, long SwapMaturityDate_Excel, long Number_Ann_CPN, long Convention1Y, long NRate, double* RateTerm, double* Rate)
{
	long i;
	long NResetTime;

	double a, b;
	double t;
	double dt = 1.0 / (double)Number_Ann_CPN;
	double LeftTime;
	double Maturity;
	double Accrued_dt;

	Maturity = (SwapMaturityDate_Excel - PricingDate_Excel) / (double)Convention1Y;

	for (i = 0; i < (long)(Maturity + 10.0) * Number_Ann_CPN; i++)
	{
		LeftTime = Maturity - ((double)i) * 1.0 / ((double)Number_Ann_CPN);
		if (LeftTime < 0.00000000001)
		{
			NResetTime = i;
			Accrued_dt = -LeftTime;
			break;
		}
	}

	long* ResetDays = (long*)malloc(sizeof(long) * NResetTime);
	double* DeltaT = (double*)malloc(sizeof(double) * NResetTime);

	for (i = 0; i < NResetTime; i++)
	{
		ResetDays[NResetTime - 1 - i] = (long)((Maturity - ((double)i) * 1.0 / ((double)Number_Ann_CPN)) * 365.0);
	}

	DeltaT[0] = ((double)ResetDays[0]) / (double)Convention1Y;
	for (i = 1; i < NResetTime; i++)
	{
		DeltaT[i] = ((double)ResetDays[i] - (double)ResetDays[i - 1]) / (double)Convention1Y;
	}

	double T;
	double ZeroRate;
	b = 0.0;
	for (i = 0; i < NResetTime; i++)
	{
		T = Day_to_T(ResetDays[i], Convention1Y);
		ZeroRate = CubicSpline(NRate, RateTerm, Rate, T);
		b += DeltaT[i] * exp(-ZeroRate * T);
	}
	ZeroRate = CubicSpline(NRate, RateTerm, Rate, T);
	a = 1.0 - exp(-ZeroRate * T);

	free(ResetDays);
	free(DeltaT);

	return a / b;
}
