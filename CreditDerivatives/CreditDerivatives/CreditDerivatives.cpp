#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalcDate.h"
#include "Util.h"


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

DLLEXPORT(long) Calc_Hazard_Rate_From_CDSCurve(
	long NPremiumCurve,				// Premium Leg 커브 Term 개수
	double* PremiumCurveTerm,		// Premium Leg 커브 Term Array
	double* PremiumCurve,			// Premium Leg 커브 Rate Array

	long NProtectionCurve,			// Protection Leg 커브 Term 개수
	double* ProtectionCurveTerm,	// Protection Leg 커브 Term
	double* ProtectionCurve,		// Protection Leg 커브 Rate Array

	long NCDSCurve,					// CDS Term 개수
	double* CDSCurveTerm,			// CDS Term Array
	double* CDSCurve,				// CDS Spread Array

	double Recovery,				// Recovery Rate
	long HazardRateCalcFlag,		// Hazard Rate Calc 방법 0: Continuous Annual Hazard Rate 1: 단리 파산확률
	long NCPN_Ann,					// 연 이자지급 수			

	double* ResultHazardTerm,		// OutPut Hazard Term
	double* ResultHazard			// OutPut Hazard Rate
);

double Calc_RiskyZeroBond(
	double Lambda,
	double rf_rate,
	double T,
	double RR,
	long CalcQMethod
)
{
	long i;
	// 10 + 1개의 Node로 나눠서 적분식 계산
	double* P0t = (double*)malloc(sizeof(double) * 11);
	double* Q0t = (double*)malloc(sizeof(double) * 11);
	P0t[0] = 1.0;
	Q0t[0] = 1.0;

	for (i = 0; i < 10; i++)
	{
		P0t[i + 1] = exp(-rf_rate * (double)(i + 1) / 10.0 * T);
		if (CalcQMethod == 0)
		{
			Q0t[i + 1] = exp(-Lambda * (double)(i + 1) / 10.0 * T);
		}
		else
		{
			Q0t[i + 1] = 1.0 - Lambda * (double)(i + 1) * T;
		}
	}

	double zerobondprice = exp(-(Lambda + rf_rate) * T);
	for (i = 1; i < 11; i++)
		zerobondprice += RR * P0t[i] * (Q0t[i - 1] - Q0t[i]); // integral_0_to_T Z(0,t) d(-Q(0,t))

	if (P0t) free(P0t);
	if (Q0t) free(Q0t);

	return zerobondprice;
}

double Calc_RiskyCouponBondCleanPrice(
	long NT_Coupon,
	double* T_Coupon,
	double* T_Coupon_Pay,
	double* CouponRate,
	double* dt,
	long NHazard,
	double* HazardTerm,
	double* HazardRate,
	long NRiskFree,
	double* RiskFreeTerm,
	double* RiskFreeRate,
	double T_Mat,
	double NotionalAmount,
	double Recovery,
	long CalcQMethod
)
{
	long i;
	double Discounted_Coupon = 0.0;
	double Discounted_Notional = 0.0;
	double Clean_Price;
	double lambda;
	double Rf;
	double Rf_Pay;
	// 쿠폰 Date가 쿠폰 PayDate와 차이나는경우
	double Fwd_Rate;

	if (NT_Coupon > 0)
	{
		for (i = 0; i < NT_Coupon; i++)
		{
			lambda = Calc_Zero_Rate(HazardTerm, HazardRate, NHazard, T_Coupon[i]);
			Rf = Calc_Zero_Rate(RiskFreeTerm, RiskFreeRate, NRiskFree, T_Coupon[i]);
			Discounted_Coupon += NotionalAmount * CouponRate[i] * dt[i] * Calc_RiskyZeroBond(lambda, Rf, T_Coupon[i], Recovery, CalcQMethod);

			if (T_Coupon[i] < T_Coupon_Pay[i])
			{
				Rf_Pay = Calc_Zero_Rate(RiskFreeTerm, RiskFreeRate, NRiskFree, T_Coupon_Pay[i]);
				Fwd_Rate = (Rf_Pay * T_Coupon_Pay[i] - Rf * T_Coupon[i]) / (T_Coupon_Pay[i] - T_Coupon[i]);
				Discounted_Coupon = Discounted_Coupon * 1.0 / (1.0 + Fwd_Rate * (T_Coupon_Pay[i] - T_Coupon[i]));
			}

		}
	}
	lambda = Calc_Zero_Rate(HazardTerm, HazardRate, NHazard, T_Mat);
	Rf = Calc_Zero_Rate(RiskFreeTerm, RiskFreeRate, NRiskFree, T_Mat);
	Discounted_Notional = NotionalAmount * Calc_RiskyZeroBond(lambda, Rf, T_Mat, Recovery, CalcQMethod);

	Clean_Price = Discounted_Coupon + Discounted_Notional;
	return Clean_Price;
}

long RiskyCouponBond_With_Schedule(
	long EffectiveDate,
	long CPNInputFlag,
	double Recovery,
	long CalcQMethod,
	// CPNInputFlag가 1일 때의 변수들
	long NCoupon,
	long* CouponDate,
	long* CouponPayDate,
	double* CouponRate,
	long PrincipalDate,
	double Principal,

	long NHazard,
	double* HazardTerm,
	double* HazardRate,

	long NRiskFree,
	double* RiskFreeTerm,
	double* RiskFreeRate,

	long PricingDate,
	long TextFlag,
	long GreekFlag,
	double* ResultPrice,
	double* KeyRate_PV01_Convex,
	double* KeyRate_HazardDelta_Gamma
)
{
	long i, j, k;
	long NT;
	long DateDiff;
	long DateDiffPay;
	double* T;
	double* T_Pay;
	double* C;
	double DeltaT;
	double* dt;
	long NumberOfEffectiveCoupon;
	Recovery = min(1.0, Recovery);
	double T_Mat = ((double)DayCountAtoB(PricingDate, PrincipalDate)) / 365.0;
	if (T_Mat < 0)
	{
		T_Mat = 0.0;
		Principal = 0.0;
	}
	double CleanPrice;

	NT = 0;
	for (i = 0; i < NCoupon; i++)
	{
		DateDiff = DayCountAtoB(EffectiveDate, CouponDate[i]);
		if (DateDiff > 0)
		{
			NT += 1;
		}
	}

	T = (double*)malloc(sizeof(double) * NT);
	T_Pay = (double*)malloc(sizeof(double) * NT);
	C = (double*)malloc(sizeof(double) * NT);
	dt = (double*)malloc(sizeof(double) * NT);
	k = 0;
	for (i = 0; i < NCoupon; i++)
	{
		DateDiff = DayCountAtoB(PricingDate, CouponDate[i]);
		DateDiffPay = DayCountAtoB(PricingDate, CouponPayDate[i]);
		if (DateDiff > 0)
		{
			T[k] = ((double)DateDiff) / 365.0;
			T_Pay[k] = ((double)DateDiffPay) / 365.0;
			C[k] = CouponRate[i];
			k += 1;
		}

	}

	NumberOfEffectiveCoupon = k;

	k = 0;
	for (i = 0; i < NCoupon; i++)
	{
		DateDiff = DayCountAtoB(PricingDate, CouponDate[i]);
		if (i == 0){
			DeltaT = ((double)DayCountAtoB(EffectiveDate, CouponDate[i]))/365.0 ;
		}
		else{
			DeltaT = ((double)DayCountAtoB(CouponDate[i - 1], CouponDate[i]))/365.0;
		}

		if (DateDiff > 0)
		{
			dt[k] = DeltaT;
			k += 1;
		}

	}
	CleanPrice = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
		NHazard, HazardTerm, HazardRate,
		NRiskFree, RiskFreeTerm, RiskFreeRate,
		T_Mat, Principal, Recovery, CalcQMethod);

	ResultPrice[0] = CleanPrice;

	long Accrued_Flag = 0;
	double Accrued_dt;
	double Accrued_CPN_Rate;
	for (i = 0; i < NCoupon - 1; i++)
		if (PricingDate > CouponDate[i - 1] && PricingDate < CouponDate[i])
		{
			Accrued_Flag = 1;
			Accrued_dt = ((double)DayCountAtoB(CouponDate[i - 1], PricingDate))/365.0;
			Accrued_CPN_Rate = CouponRate[i - 1];
			break;
		}
	if (Accrued_Flag == 1)
	{
		ResultPrice[1] = Principal * Accrued_dt * Accrued_CPN_Rate;
		ResultPrice[2] = ResultPrice[0] + ResultPrice[1];
	}
	else
	{
		ResultPrice[1] = 0.0;
		ResultPrice[2] = CleanPrice;
	}


	double P_Ru;
	double P_Rd;
	double* CurveRUp;
	double* CurveRDn;
	double* CreditCurveLambdaUp;
	double* CreditCurveLambdaDn;
	double PV01;
	double Convexity;
	double HazardDelta;
	double HazardGamma;
	double RecoveryDelta;
	if (GreekFlag >= 1)
	{
		CurveRUp = (double*)malloc(sizeof(double) * NRiskFree);
		CurveRDn = (double*)malloc(sizeof(double) * NRiskFree);
		for (i = 0; i < NRiskFree; i++)
		{
			CurveRUp[i] = RiskFreeRate[i] + 0.0001;
			CurveRDn[i] = RiskFreeRate[i] - 0.0001;
		}

		P_Ru = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, CurveRUp,
			T_Mat, Principal, Recovery, CalcQMethod);

		P_Rd = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, CurveRDn,
			T_Mat, Principal, Recovery, CalcQMethod);

		PV01 = (P_Rd - P_Ru) / (2.0);
		Convexity = (P_Rd + P_Ru - 2.0 * CleanPrice);
		ResultPrice[3] = PV01;
		ResultPrice[4] = Convexity;

		CreditCurveLambdaUp = (double*)malloc(sizeof(double) * NHazard);
		CreditCurveLambdaDn = (double*)malloc(sizeof(double) * NHazard);

		for (i = 0; i < NHazard; i++)
		{
			CreditCurveLambdaUp[i] = HazardRate[i] + 0.0001;
			CreditCurveLambdaDn[i] = HazardRate[i] - 0.0001;
		}

		P_Ru = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, CreditCurveLambdaUp,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, Principal, Recovery, CalcQMethod);

		P_Rd = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, CreditCurveLambdaDn,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, Principal, Recovery, CalcQMethod);

		HazardDelta = (P_Rd - P_Ru) / (2.0);
		HazardGamma = (P_Rd + P_Ru - 2.0 * CleanPrice);

		ResultPrice[5] = HazardDelta;
		ResultPrice[6] = HazardGamma;

		P_Ru = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, Principal, min(Recovery + 0.01, 1.0), CalcQMethod);
		P_Rd = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, Principal, Recovery - 0.01, CalcQMethod);


		RecoveryDelta = (P_Ru - P_Rd) / (2.0);
		ResultPrice[7] = RecoveryDelta;

		if (GreekFlag >= 2)
		{
			for (i = 0; i < NRiskFree; i++)
			{

				for (j = 0; j < NRiskFree; j++)
				{
					if (i == j)
					{
						CurveRUp[j] = RiskFreeRate[j] + 0.0001;
						CurveRDn[j] = RiskFreeRate[j] - 0.0001;
					}
					else
					{
						CurveRUp[j] = RiskFreeRate[j];
						CurveRDn[j] = RiskFreeRate[j];
					}
				}

				P_Ru = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
					NHazard, HazardTerm, HazardRate,
					NRiskFree, RiskFreeTerm, CurveRUp,
					T_Mat, Principal, Recovery, CalcQMethod);

				P_Rd = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
					NHazard, HazardTerm, HazardRate,
					NRiskFree, RiskFreeTerm, CurveRDn,
					T_Mat, Principal, Recovery, CalcQMethod);

				PV01 = (P_Rd - P_Ru) / (2.0);
				Convexity = (P_Rd + P_Ru - 2.0 * CleanPrice);
				KeyRate_PV01_Convex[i] = PV01;
				KeyRate_PV01_Convex[i + NRiskFree] = Convexity;
			}

			for (i = 0; i < NHazard; i++)
			{

				for (j = 0; j < NHazard; j++)
				{
					if (i == j)
					{
						CreditCurveLambdaUp[j] = HazardRate[j] + 0.0001;
						CreditCurveLambdaDn[j] = HazardRate[j] - 0.0001;
					}
					else
					{
						CreditCurveLambdaUp[j] = HazardRate[j];
						CreditCurveLambdaDn[j] = HazardRate[j];
					}
				}

				P_Ru = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
					NHazard, HazardTerm, CreditCurveLambdaUp,
					NRiskFree, RiskFreeTerm, RiskFreeRate,
					T_Mat, Principal, Recovery, CalcQMethod);

				P_Rd = Calc_RiskyCouponBondCleanPrice(NumberOfEffectiveCoupon, T, T_Pay, C, dt,
					NHazard, HazardTerm, CreditCurveLambdaDn,
					NRiskFree, RiskFreeTerm, RiskFreeRate,
					T_Mat, Principal, Recovery, CalcQMethod);

				HazardDelta = (P_Rd - P_Ru) / (2.0);
				HazardGamma = (P_Rd + P_Ru - 2.0 * CleanPrice);
				KeyRate_HazardDelta_Gamma[i] = HazardDelta;
				KeyRate_HazardDelta_Gamma[i + NHazard] = HazardGamma;
			}
		}

		if (CurveRUp) free(CurveRUp);
		if (CurveRDn) free(CurveRDn);
		if (CreditCurveLambdaUp) free(CreditCurveLambdaUp);
		if (CreditCurveLambdaDn) free(CreditCurveLambdaDn);
	}

	if (T) free(T);
	if (T_Pay) free(T_Pay);
	if (C) free(C);
	if (dt) free(dt);

	return 1;
}

long RiskyCouponBond(
	long EffectiveDate,
	long CPNInputFlag,
	double Recovery,
	long CalcQMethod,
	// CPNInputFlag가 0일 때의 추가변수들
	double CPNRate_0,
	long NCPN_Ann_0,
	double NotionalAmount_0,
	double T_Mat_0,
	long NHazard,
	double* HazardTerm,
	double* HazardRate,

	long NRiskFree,
	double* RiskFreeTerm,
	double* RiskFreeRate,

	long PricingDate,
	long TextFlag,
	long GreekFlag,
	double* ResultPrice,
	double* KeyRate_PV01_Convex,
	double* KeyRate_HazardDelta_Gamma
)
{
	long i, j;
	long NT;

	double* T;
	double* C;

	double* dt;

	double T_Mat = T_Mat_0;

	double Accrued_dt;
	double LeftTime;
	for (i = 0; i < (long)(T_Mat_0 + 10.0) * NCPN_Ann_0; i++)
	{
		LeftTime = T_Mat - ((double)i) * 1.0 / ((double)NCPN_Ann_0);
		if (LeftTime < 0.00000001) // 부동소수점 문제로 인해
		{
			NT = i;
			Accrued_dt = -LeftTime;
			break;
		}
	}

	double CleanPrice;

	//NT = (long)(T_Mat_0 * NCPN_Ann_0);
	T = (double*)malloc(sizeof(double) * NT);
	C = (double*)malloc(sizeof(double) * NT);
	dt = (double*)malloc(sizeof(double) * NT);

	for (i = 0; i < NT; i++)
	{
		dt[i] = 1.0 / ((double)NCPN_Ann_0);
		T[NT - 1 - i] = T_Mat - ((double)i) * 1.0 / ((double)NCPN_Ann_0);;
		C[i] = CPNRate_0;
	}

	double* T_Pay = T;

	CleanPrice = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
		NHazard, HazardTerm, HazardRate,
		NRiskFree, RiskFreeTerm, RiskFreeRate,
		T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

	ResultPrice[0] = CleanPrice;
	ResultPrice[1] = Accrued_dt * CPNRate_0 * NotionalAmount_0;
	ResultPrice[2] = CleanPrice + ResultPrice[1];

	double P_Ru;
	double P_Rd;
	double* CurveRUp;
	double* CurveRDn;
	double* CreditCurveLambdaUp;
	double* CreditCurveLambdaDn;
	double PV01;
	double Convexity;
	double HazardDelta;
	double HazardGamma;
	double RecoveryDelta;
	if (GreekFlag >= 1)
	{
		CurveRUp = (double*)malloc(sizeof(double) * NRiskFree);
		CurveRDn = (double*)malloc(sizeof(double) * NRiskFree);
		for (i = 0; i < NRiskFree; i++)
		{
			CurveRUp[i] = RiskFreeRate[i] + 0.0001;
			CurveRDn[i] = RiskFreeRate[i] - 0.0001;
		}

		P_Ru = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, CurveRUp,
			T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

		P_Rd = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, CurveRDn,
			T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

		PV01 = (P_Rd - P_Ru) / (2.0);
		Convexity = (P_Rd + P_Ru - 2.0 * CleanPrice);
		ResultPrice[3] = PV01;
		ResultPrice[4] = Convexity;

		CreditCurveLambdaUp = (double*)malloc(sizeof(double) * NHazard);
		CreditCurveLambdaDn = (double*)malloc(sizeof(double) * NHazard);

		for (i = 0; i < NHazard; i++)
		{
			CreditCurveLambdaUp[i] = HazardRate[i] + 0.0001;
			CreditCurveLambdaDn[i] = HazardRate[i] - 0.0001;
		}

		P_Ru = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, CreditCurveLambdaUp,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

		P_Rd = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, CreditCurveLambdaDn,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

		HazardDelta = (P_Rd - P_Ru) / (2.0);
		HazardGamma = (P_Rd + P_Ru - 2.0 * CleanPrice);

		ResultPrice[5] = HazardDelta;
		ResultPrice[6] = HazardGamma;

		P_Ru = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, NotionalAmount_0, min(Recovery + 0.01, 1.0), CalcQMethod);
		P_Rd = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
			NHazard, HazardTerm, HazardRate,
			NRiskFree, RiskFreeTerm, RiskFreeRate,
			T_Mat, NotionalAmount_0, Recovery - 0.01, CalcQMethod);


		RecoveryDelta = (P_Ru - P_Rd) / (2.0);
		ResultPrice[7] = RecoveryDelta;

		if (GreekFlag >= 2)
		{
			for (i = 0; i < NRiskFree; i++)
			{

				for (j = 0; j < NRiskFree; j++)
				{
					if (i == j)
					{
						CurveRUp[j] = RiskFreeRate[j] + 0.0001;
						CurveRDn[j] = RiskFreeRate[j] - 0.0001;
					}
					else
					{
						CurveRUp[j] = RiskFreeRate[j];
						CurveRDn[j] = RiskFreeRate[j];
					}
				}

				P_Ru = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
					NHazard, HazardTerm, HazardRate,
					NRiskFree, RiskFreeTerm, CurveRUp,
					T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

				P_Rd = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
					NHazard, HazardTerm, HazardRate,
					NRiskFree, RiskFreeTerm, CurveRDn,
					T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

				PV01 = (P_Rd - P_Ru) / (2.0);
				Convexity = (P_Rd + P_Ru - 2.0 * CleanPrice);
				KeyRate_PV01_Convex[i] = PV01;
				KeyRate_PV01_Convex[i + NRiskFree] = Convexity;
			}

			for (i = 0; i < NHazard; i++)
			{

				for (j = 0; j < NHazard; j++)
				{
					if (i == j)
					{
						CreditCurveLambdaUp[j] = HazardRate[j] + 0.0001;
						CreditCurveLambdaDn[j] = HazardRate[j] - 0.0001;
					}
					else
					{
						CreditCurveLambdaUp[j] = HazardRate[j];
						CreditCurveLambdaDn[j] = HazardRate[j];
					}
				}

				P_Ru = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
					NHazard, HazardTerm, CreditCurveLambdaUp,
					NRiskFree, RiskFreeTerm, RiskFreeRate,
					T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

				P_Rd = Calc_RiskyCouponBondCleanPrice(NT, T, T_Pay, C, dt,
					NHazard, HazardTerm, CreditCurveLambdaDn,
					NRiskFree, RiskFreeTerm, RiskFreeRate,
					T_Mat, NotionalAmount_0, Recovery, CalcQMethod);

				HazardDelta = (P_Rd - P_Ru) / (2.0);
				HazardGamma = (P_Rd + P_Ru - 2.0 * CleanPrice);
				KeyRate_HazardDelta_Gamma[i] = HazardDelta;
				KeyRate_HazardDelta_Gamma[i + NHazard] = HazardGamma;
			}
		}

		if (CurveRUp) free(CurveRUp);
		if (CurveRDn) free(CurveRDn);
		if (CreditCurveLambdaUp) free(CreditCurveLambdaUp);
		if (CreditCurveLambdaDn) free(CreditCurveLambdaDn);
	}

	if (T) free(T);
	if (C) free(C);
	if (dt) free(dt);

	return 1;
}

long Inputcheck_RiskyCouponBond(
	long EffectiveDate,
	long CPNInputFlag,
	double Recovery,
	long CalcQMethod,
	// CPNInputFlag가 0일 때의 추가변수들
	double CPNRate_0,
	long NCPN_Ann_0,
	double NotionalAmount_0,
	double T_Mat_0,
	// CPNInputFlag가 1일 때의 변수들
	long NCoupon,
	long* CouponDate,
	long* CouponPayDate,
	double* CouponRate,
	long PrincipalDate,
	double Principal,

	long NHazard,
	double* HazardTerm,
	double* HazardRate,

	long NRiskFree,
	double* RiskFreeTerm,
	double* RiskFreeRate,

	long PricingDate,
	long TextFlag,
	long GreekFlag)
{
	long i;

	long EFFECTIVE_DATE_ERROR = -101;
	long CPN_INPUT_FLAG_ERROR = -103;
	long RECOVERY_ERROR = -105;
	long CALC_Q_METHOD_ERROR = -107;
	long CPN_RATE_ERROR = -109;
	long NCPN_ANN_ERROR = -111;
	long NOTIONAL_ERROR = -113;
	long MATURITY_ERROR = -115;
	long NCOUPON_ERROR = -117;
	long COUPONDATE_ERROR = -119;
	long COUPONRATE_ERROR = -121;
	long NOTIONAL_DATE_ERROR = -123;
	long NHAZARD_ERROR = -125;
	long HAZARD_TERM_ERROR = -127;
	long HAZARD_RATE_ERROR = -129;
	long NRISKFREE_ERROR = -131;
	long RISKFREETERM_ERROR = -133;
	long RISKFREERATE_ERROR = -135;
	long PRICINGDATE_ERROR = -137;
	long TEXTFLAGERROR = -139;
	long GREEKFLAGERROR = -141;

	if (EffectiveDate < 19000101 || EffectiveDate > 99991231) return EFFECTIVE_DATE_ERROR;

	if (CPNInputFlag < 0 || CPNInputFlag >= 2) return CPN_INPUT_FLAG_ERROR;

	if (Recovery < 0.0 || Recovery >= 1.0) return RECOVERY_ERROR;

	if (CalcQMethod < 0 || CalcQMethod >= 3) return CALC_Q_METHOD_ERROR;

	if (CPNInputFlag == 0)
	{
		if (NCPN_Ann_0 < 1) return NCPN_ANN_ERROR;

		if (NotionalAmount_0 < 0.0) return NOTIONAL_ERROR;

		if (T_Mat_0 < 0.0) return MATURITY_ERROR;

	}
	else
	{
		if (NCoupon < 0) return NCOUPON_ERROR;

		for (i = 1; i < NCoupon; i++)
		{
			if (CouponDate[i] < CouponDate[i - 1]) return COUPONDATE_ERROR;
			if (CouponPayDate[i] < CouponPayDate[i - 1]) return COUPONDATE_ERROR;
		}

		if (PrincipalDate < 19000101 || PrincipalDate > 99991231) return NOTIONAL_DATE_ERROR;

		if (Principal < 0) return NOTIONAL_ERROR;

		//Paydate기준으로 수정 예정
		if (PricingDate > CouponDate[NCoupon - 1]) return PRICINGDATE_ERROR;
	}

	// Hazard 관련 에러
	if (NHazard < 1) return NHAZARD_ERROR;

	if (HazardTerm[0] < 0.0) return HAZARD_TERM_ERROR;

	for (i = 1; i < NHazard; i++)
	{
		if (HazardTerm[i] < HazardTerm[i - 1] || HazardTerm[i] < 0.0) return HAZARD_TERM_ERROR;

		//이 조건을 만족해야함 e^(-lambda1 x T1) < e^(-lambda2 x T2)
		if (HazardRate[i] < HazardRate[i - 1] * HazardTerm[i - 1] / HazardTerm[i]) return HAZARD_RATE_ERROR;
	}

	if (NRiskFree < 1) return NRISKFREE_ERROR;

	return 0;
}

DLLEXPORT(long) Calc_RiskyCouponBond(
	long EffectiveDate,
	long CPNInputFlag,
	double Recovery,
	long CalcQMethod,
	// CPNInputFlag가 0일 때의 추가변수들
	double CPNRate_0,
	long NCPN_Ann_0,
	double NotionalAmount_0,
	double T_Mat_0,
	// CPNInputFlag가 1일 때의 변수들
	long NCoupon,
	long* CouponDate,
	long* CouponPayDate,
	double* CouponRate,
	long PrincipalDate,
	double Principal,

	long InputType,
	long NHazardCDS,
	double* HazardCDSTerm,
	double* HazardCDSRate,

	long NRiskFree,
	double* RiskFreeTerm,
	double* RiskFreeRate,

	long PricingDate,
	long TextFlag,
	long GreekFlag,
	double* ResultPrice,
	double* KeyRate_PV01_Convex,
	double* KeyRate_HazardDelta_Gamma
)
{


	// CDS Curve 기반으로 자체 Calibration할 경우의 변수
	double* TempHazardRate = (double*)malloc(sizeof(double) * NHazardCDS);
	double* TempHazardTerm = (double*)malloc(sizeof(double) * NHazardCDS);
	long ResultCode = 0;

	ResultCode = Inputcheck_RiskyCouponBond(EffectiveDate, CPNInputFlag, Recovery, CalcQMethod,
		CPNRate_0, NCPN_Ann_0, NotionalAmount_0, T_Mat_0,
		NCoupon, CouponDate, CouponPayDate, CouponRate, PrincipalDate,
		Principal, NHazardCDS, HazardCDSTerm, HazardCDSRate,
		NRiskFree, RiskFreeTerm, RiskFreeRate,
		PricingDate, TextFlag, GreekFlag);

	if (ResultCode < 0)
		return ResultCode;

	if (CPNInputFlag == 1)
	{
		if (InputType == 0)
		{
			ResultCode = RiskyCouponBond_With_Schedule(EffectiveDate, CPNInputFlag, Recovery, CalcQMethod, NCoupon,
				CouponDate, CouponPayDate, CouponRate, PrincipalDate, Principal, NHazardCDS,
				HazardCDSTerm, HazardCDSRate, NRiskFree, RiskFreeTerm, RiskFreeRate,
				PricingDate, TextFlag, GreekFlag, ResultPrice, KeyRate_PV01_Convex, KeyRate_HazardDelta_Gamma);
		}
		else
		{
			ResultCode = Calc_Hazard_Rate_From_CDSCurve(NRiskFree, RiskFreeTerm, RiskFreeRate,
				NRiskFree, RiskFreeTerm, RiskFreeRate,
				NHazardCDS, HazardCDSTerm, HazardCDSRate,
				Recovery, 0, 2, TempHazardTerm, TempHazardRate);

			ResultCode = RiskyCouponBond_With_Schedule(EffectiveDate, CPNInputFlag, Recovery, CalcQMethod, NCoupon,
				CouponDate, CouponPayDate, CouponRate, PrincipalDate, Principal, NHazardCDS,
				TempHazardTerm, TempHazardRate, NRiskFree, RiskFreeTerm, RiskFreeRate,
				PricingDate, TextFlag, GreekFlag, ResultPrice, KeyRate_PV01_Convex, KeyRate_HazardDelta_Gamma);
		}
	}
	else
	{
		if (InputType == 0)
		{
			ResultCode = RiskyCouponBond(EffectiveDate, CPNInputFlag, Recovery, CalcQMethod, CPNRate_0, NCPN_Ann_0, NotionalAmount_0, T_Mat_0, NHazardCDS, HazardCDSTerm, HazardCDSRate, NRiskFree, RiskFreeTerm, RiskFreeRate, PricingDate, TextFlag, GreekFlag, ResultPrice, KeyRate_PV01_Convex, KeyRate_HazardDelta_Gamma);
		}
		else
		{
			ResultCode = Calc_Hazard_Rate_From_CDSCurve(NRiskFree, RiskFreeTerm, RiskFreeRate,
				NRiskFree, RiskFreeTerm, RiskFreeRate,
				NHazardCDS, HazardCDSTerm, HazardCDSRate,
				Recovery, 0, 2, TempHazardTerm, TempHazardRate);

			ResultCode = RiskyCouponBond(EffectiveDate, CPNInputFlag, Recovery, CalcQMethod, CPNRate_0, NCPN_Ann_0, NotionalAmount_0, T_Mat_0, NHazardCDS, HazardCDSTerm, HazardCDSRate, NRiskFree, RiskFreeTerm, RiskFreeRate, PricingDate, TextFlag, GreekFlag, ResultPrice, KeyRate_PV01_Convex, KeyRate_HazardDelta_Gamma);
		}
	}


	if (TempHazardRate) free(TempHazardRate);
	if (TempHazardTerm) free(TempHazardTerm);
	return 1;
}




double Calibrate_HazardRate_CreditSpread(
	double riskless_rate,
	double risky_rate,
	double T,
	double RR
)
{

	double ObjValue = exp(-risky_rate * T);

	double value_a = 1.0;
	double value_b = 0.0;
	double value;
	double dblErrorRange = 0.0001;
	double dblCalcPrice = Calc_RiskyZeroBond(value_a, riskless_rate, T, RR, 0);
	if (dblCalcPrice <= 0)  return 0.0;

	double dblPriceGab = ObjValue - dblCalcPrice;
	if (fabs(dblPriceGab) < dblErrorRange)
		return value_a;

	while (dblPriceGab < 0)
	{
		value_b = value_a;
		value_a *= 2.0;

		dblCalcPrice = Calc_RiskyZeroBond(value_a, riskless_rate, T, RR, 0);
		dblPriceGab = ObjValue - dblCalcPrice;

		if (fabs(dblPriceGab) < dblErrorRange)
			return value_a;
	}

	value = (value_a + value_b) / 2.0;
	for (int i = 0; i < 1000; i++)
	{
		dblCalcPrice = Calc_RiskyZeroBond(value, riskless_rate, T, RR, 0);
		dblPriceGab = ObjValue - dblCalcPrice;

		// 절대값 < 설정 Ytm 오차
		if (fabs(dblPriceGab) < dblErrorRange)
			break;

		// Ytm < 설정 Ytm 오차 || Ytm == 0
		if (value < dblErrorRange || value == 0.0)
			return 0.0;

		if (dblPriceGab > 0.0)
			value_a = value;
		else
			value_b = value;

		value = (value_a + value_b) / 2.0;
	}
	return value;
}

double Calibrate_HazardRate_CreditSpread2(
	double riskless_rate,
	double risky_rate,
	double T,
	double RR
)
{
	double Lambda = (risky_rate - riskless_rate) / (1 - RR);
	return Lambda;
}

long Calibrate_HazardRateStructure_CreditSpread(
	long NRefZero,
	double* RefZeroTerm,
	double* ZeroRiskFree,
	double* ZeroRiskyBond,
	double Recovery,
	double* ResultHazardStructure,
	long DefaultRateCalcMethod
)
{
	long i;
	for (i = 0; i < NRefZero; i++)
	{
		if (DefaultRateCalcMethod == 0)
		{
			ResultHazardStructure[i] = Calibrate_HazardRate_CreditSpread(ZeroRiskFree[i], ZeroRiskyBond[i], RefZeroTerm[i], Recovery);
		}
		else
		{
			ResultHazardStructure[i] = Calibrate_HazardRate_CreditSpread2(ZeroRiskFree[i], ZeroRiskyBond[i], RefZeroTerm[i], Recovery);
		}
	}
	return 1;
}

double Calc_CDS_From_Hazard(
	long NPremiumCurve,
	double* PremiumCurveTerm,
	double* PremiumCurve,

	long NProtectionCurve,
	double* ProtectionCurveTerm,
	double* ProtectionCurve,

	long DataFlag,
	long NHazardRate,
	double* HazardRateTerm,
	double* HazardRate,
	long CalcQMethod,

	double Recovery,
	long ScheduleInputFlag,
	long NCPN_Ann,
	double Maturity,
	long NSchedule,
	long* ResetDateExcelDate,
	long* PayDateExcelDate,
	long PricingDateExcelDate,

	double* Premium_Leg,
	double* Protection_Leg,
	double* Premium_Schedule,
	double NotionalAmount
)
{
	long i;
	long k;
	double *dt_Premium; //= 1.0 / ((double)NCPN_Ann);
	double dt_Protection = 1.0 / 12.0;
	double t;

	double Result_Spread;
	double Premium_Value_At_T;
	// 이자지급횟수  + 1 만큼 할당
	long N_PremiumNode;
	double* Z_Premium;
	double* Q_Premium;
	if (ScheduleInputFlag == 0)
	{
		N_PremiumNode = (long)((double)NCPN_Ann * Maturity + 0.5) + 1;
		Z_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);
		dt_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);
		for (i = 0; i < N_PremiumNode; i++)
		{
			dt_Premium[i] = 1.0 / ((double)NCPN_Ann);
		}
		Z_Premium[0] = 1.0;
		t = dt_Premium[0];
		for (i = 1; i < N_PremiumNode; i++)
		{
			Z_Premium[i] = Calc_Discount_Factor(PremiumCurveTerm, PremiumCurve, NPremiumCurve, t);
			t += dt_Premium[i];
		}

		Q_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);
		Q_Premium[0] = 1.0;
		t = dt_Premium[0];
		for (i = 1; i < N_PremiumNode; i++)
		{
			if (CalcQMethod == 0)
			{
				Q_Premium[i] = exp(-Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t);
			}
			else
			{
				Q_Premium[i] = 1.0 - Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t;
			}
			t += dt_Premium[i];
		}
	}
	else
	{
		Maturity = (PayDateExcelDate[NSchedule - 1] - PricingDateExcelDate) / 365.0;
		N_PremiumNode = 1;
		for (i = 0; i < NSchedule; i++)
		{
			if (PricingDateExcelDate < PayDateExcelDate[i])
			{
				N_PremiumNode += 1;
			}
		}
		Z_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);
		dt_Premium = (double*)malloc(sizeof(double) * N_PremiumNode -1);
		Q_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);

		Z_Premium[0] = 1.0;
		Q_Premium[0] = 1.0;
		k = 1;
		for (i = 0; i < NSchedule; i++)
		{
			if (PricingDateExcelDate < PayDateExcelDate[i])
			{
				dt_Premium[k-1] = (PayDateExcelDate[i] - ResetDateExcelDate[i])/365.0;
				t = (PayDateExcelDate[i] - PricingDateExcelDate) / 365.0;
				Z_Premium[k] = Calc_Discount_Factor(PremiumCurveTerm, PremiumCurve, NPremiumCurve, t);

				if (CalcQMethod == 0) {
					Q_Premium[k] = exp(-Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t);
				}
				else{
					Q_Premium[k] = 1.0 - Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t;
				}
				k+=1;

			}
		}

	}
	// 12 * 만기(년) 수 + 1만큼 할당 
	long N_ProtectionNode = (long)(12.0 * Maturity + 0.5) + 1;
	double* Z_Protection = (double*)malloc(sizeof(double) * N_ProtectionNode);
	Z_Protection[0] = 1.0;
	t = dt_Protection;
	for (i = 1; i < N_ProtectionNode; i++)
	{
		Z_Protection[i] = Calc_Discount_Factor(ProtectionCurveTerm, ProtectionCurve, NProtectionCurve, t);
		t += dt_Protection;
	}

	double* Q_Protection = (double*)malloc(sizeof(double) * N_ProtectionNode);
	Q_Protection[0] = 1.0;
	t = dt_Protection;
	for (i = 1; i < N_ProtectionNode; i++)
	{
		if (CalcQMethod == 0)
		{
			Q_Protection[i] = exp(-Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t);
		}
		else
		{
			Q_Protection[i] = 1.0 - Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t;
		}
		t += dt_Protection;
	}

	double Premium_Value, Protection_Value;

	Protection_Value = 0.0;
	for (i = 1; i < N_ProtectionNode; i++)
		Protection_Value += 0.5 * (1.0 - Recovery) * (Z_Protection[i - 1] + Z_Protection[i]) * max(0.0, (Q_Protection[i - 1] - Q_Protection[i]));

	Premium_Value = 0.0;
	for (i = 1; i < N_PremiumNode; i++)
	{
		Premium_Value_At_T = 0.5 * dt_Premium[i-1] * Z_Premium[i] * (Q_Premium[i - 1] + Q_Premium[i]);
		Premium_Value += Premium_Value_At_T;
		Premium_Schedule[i - 1] = Premium_Value_At_T;
	}

	Result_Spread = Protection_Value / Premium_Value;
	*Protection_Leg = Protection_Value * NotionalAmount;
	*Premium_Leg = Premium_Value * Result_Spread * NotionalAmount;

	for (i = 0; i < N_PremiumNode - 1; i++)
	{
		Premium_Schedule[i] *= Result_Spread * NotionalAmount;
	}

	if (Z_Premium) free(Z_Premium);
	if (Z_Protection) free(Z_Protection);
	if (Q_Premium) free(Q_Premium);
	if (Q_Protection) free(Q_Protection);
	if (dt_Premium) free(dt_Premium);
	return Result_Spread;
}

double Calc_CDS_Greeks_From_Hazard(
	long NPremiumCurve,
	double* PremiumCurveTerm,
	double* PremiumCurve,

	long NProtectionCurve,
	double* ProtectionCurveTerm,
	double* ProtectionCurve,

	long DataFlag,
	long NHazardRate,
	double* HazardRateTerm,
	double* HazardRate,
	long CalcQMethod,

	double Recovery,
	long NCPN_Ann,
	double Maturity,

	double* Premium_Leg,
	double* Protection_Leg,
	double*
)
{
	long i;
	double dt_Premium = 1.0 / ((double)NCPN_Ann);
	double dt_Protection = 1.0 / 12.0;
	double t;

	double Result_Spread;

	// 이자지급횟수  + 1 만큼 할당
	long N_PremiumNode = (long)((double)NCPN_Ann * Maturity + 0.5) + 1;
	double* Z_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);

	Z_Premium[0] = 1.0;
	t = dt_Premium;
	for (i = 1; i < N_PremiumNode; i++)
	{
		Z_Premium[i] = Calc_Discount_Factor(PremiumCurveTerm, PremiumCurve, NPremiumCurve, t);
		t += dt_Premium;
	}

	double* Q_Premium = (double*)malloc(sizeof(double) * N_PremiumNode);
	Q_Premium[0] = 1.0;
	t = dt_Premium;
	for (i = 1; i < N_PremiumNode; i++)
	{
		if (CalcQMethod == 0)
		{
			Q_Premium[i] = exp(-Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t);
		}
		else
		{
			Q_Premium[i] = 1.0 - Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t;
		}
		t += dt_Premium;
	}

	// 12 * 만기(년) 수 + 1만큼 할당 
	long N_ProtectionNode = (long)(12.0 * Maturity) + 1;
	double* Z_Protection = (double*)malloc(sizeof(double) * N_ProtectionNode);
	Z_Protection[0] = 1.0;
	t = dt_Protection;
	for (i = 1; i < N_ProtectionNode; i++)
	{
		Z_Protection[i] = Calc_Discount_Factor(ProtectionCurveTerm, ProtectionCurve, NProtectionCurve, t);
		t += dt_Protection;
	}

	double* Q_Protection = (double*)malloc(sizeof(double) * N_ProtectionNode);
	Q_Protection[0] = 1.0;
	t = dt_Protection;
	for (i = 1; i < N_ProtectionNode; i++)
	{
		if (CalcQMethod == 0)
		{
			Q_Protection[i] = exp(-Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t);
		}
		else
		{
			Q_Protection[i] = 1.0 - Calc_Zero_Rate(HazardRateTerm, HazardRate, NHazardRate, t) * t;
		}
		t += dt_Protection;
	}

	double Premium_Value, Protection_Value;

	Protection_Value = 0.0;
	for (i = 1; i < N_ProtectionNode; i++)
		Protection_Value += 0.5 * (1.0 - Recovery) * (Z_Protection[i - 1] + Z_Protection[i]) * max(0.0, (Q_Protection[i - 1] - Q_Protection[i]));

	Premium_Value = 0.0;
	for (i = 1; i < N_PremiumNode; i++)
	{
		Premium_Value += 0.5 * dt_Premium * Z_Premium[i] * (Q_Premium[i - 1] + Q_Premium[i]);
	}

	Result_Spread = Protection_Value / Premium_Value;
	*Protection_Leg = Protection_Value;
	*Premium_Leg = Premium_Value * Result_Spread;

	if (Z_Premium) free(Z_Premium);
	if (Z_Protection) free(Z_Protection);
	if (Q_Premium) free(Q_Premium);
	if (Q_Protection) free(Q_Protection);
	return Result_Spread;
}

double Calc_CDS_From_CreditSpread(
	long NPremiumCurve,
	double* PremiumCurveTerm,
	double* PremiumCurve,

	long NProtectionCurve,
	double* ProtectionCurveTerm,
	double* ProtectionCurve,

	long DataFlag,
	long NRefZero,
	double* RefZeroTerm,
	double* ZeroRiskFree,
	double* ZeroRiskyBond,
	long DefaultRateCalcMethod,

	double Recovery,
	long ScheduleInputFlag,
	long NCPN_Ann,
	double Maturity,
	long NSchedule,
	long* ResetDateExcelDate,
	long* PayDateExcelDate,
	long PricingDateExcelDate,

	double* Premium_Leg,
	double* Protection_Leg,
	double* ResultHazardRate,
	double* Premium_Schedule,
	long GreekFlag,
	double* ResultGreeks,
	double OLD_CDS_Spread,
	double NotionalAmount
)
{
	long i;
	long temp;
	double ResultSpread;
	double Temp_Protection_Leg = 0.0;
	double Temp_Premium_Leg = 0.0;
	double P_Ru, P_Rd;

	double* HazardRate = (double*)malloc(sizeof(double) * NRefZero);
	double* PremiumCurveRup = (double*)malloc(sizeof(double) * NPremiumCurve);
	double* PremiumCurveRdn = (double*)malloc(sizeof(double) * NPremiumCurve);
	double* ProtectionCurveRup = (double*)malloc(sizeof(double) * NProtectionCurve);
	double* ProtectionCurveRdn = (double*)malloc(sizeof(double) * NProtectionCurve);
	double* HazardUp = (double*)malloc(sizeof(double) * NRefZero);
	double* HazardDn = (double*)malloc(sizeof(double) * NRefZero);
	double* Temp_Premium_Schedule = (double*)malloc(sizeof(double) * (long)(NCPN_Ann * Maturity + 2.0));

	temp = Calibrate_HazardRateStructure_CreditSpread(NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, Recovery, HazardRate, DefaultRateCalcMethod);

	long NHazardRate = NRefZero;
	double* HazardRateTerm = RefZeroTerm;

	ResultSpread = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
		NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
		DataFlag, NHazardRate, HazardRateTerm, HazardRate, 0,
		Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, Premium_Leg, Protection_Leg, Premium_Schedule, NotionalAmount);


	for (i = 0; i < NRefZero; i++)
	{
		ResultHazardRate[i] = 1.0 - exp(-HazardRate[i] * HazardRateTerm[i]);
	}

	double Protection_Value_U;
	double RPV01_U;
	double Premium_Value_U;
	double Value_ProtectionBuyer_U;
	double Protection_Value_D;
	double RPV01_D;
	double Premium_Value_D;
	double Value_ProtectionBuyer_D;
	if (GreekFlag >= 1)
	{
		for (i = 0; i < NPremiumCurve; i++)
		{
			PremiumCurveRup[i] = PremiumCurve[i] + 0.0001;
			PremiumCurveRdn[i] = PremiumCurve[i] - 0.0001;
		}

		for (i = 0; i < NProtectionCurve; i++)
		{
			ProtectionCurveRup[i] = ProtectionCurve[i] + 0.0001;
			ProtectionCurveRdn[i] = ProtectionCurve[i] - 0.0001;
		}

		P_Ru = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurveRup,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurveRup,
			DataFlag, NHazardRate, HazardRateTerm, HazardRate, 0,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

		Protection_Value_U = Temp_Protection_Leg;
		RPV01_U = Temp_Premium_Leg / P_Ru;
		Premium_Value_U = OLD_CDS_Spread * RPV01_U;
		Value_ProtectionBuyer_U = Protection_Value_U - Premium_Value_U;

		P_Rd = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurveRdn,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurveRdn,
			DataFlag, NHazardRate, HazardRateTerm, HazardRate, 0,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

		Protection_Value_D = Temp_Protection_Leg;
		RPV01_D = Temp_Premium_Leg / P_Rd;
		Premium_Value_D = OLD_CDS_Spread * RPV01_D;
		Value_ProtectionBuyer_D = Protection_Value_D - Premium_Value_D;

		ResultGreeks[0] = (Value_ProtectionBuyer_U - Value_ProtectionBuyer_D) / 2.0;
		ResultGreeks[1] = (Premium_Value_U - Premium_Value_D) / 2.0;
		ResultGreeks[2] = (Protection_Value_U - Protection_Value_D) / 2.0;

		for (i = 0; i < NRefZero; i++)
		{
			HazardUp[i] = HazardRate[i] + 0.0001;
			HazardDn[i] = HazardRate[i] - 0.0001;
		}

		P_Ru = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NHazardRate, HazardRateTerm, HazardUp, 0,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

		Protection_Value_U = Temp_Protection_Leg;
		RPV01_U = Temp_Premium_Leg / P_Ru;
		Premium_Value_U = OLD_CDS_Spread * RPV01_U;
		Value_ProtectionBuyer_U = Protection_Value_U - Premium_Value_U;

		P_Rd = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NHazardRate, HazardRateTerm, HazardDn, 0,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

		Protection_Value_D = Temp_Protection_Leg;
		RPV01_D = Temp_Premium_Leg / P_Rd;
		Premium_Value_D = OLD_CDS_Spread * RPV01_D;
		Value_ProtectionBuyer_D = Protection_Value_D - Premium_Value_D;

		ResultGreeks[3] = (Value_ProtectionBuyer_U - Value_ProtectionBuyer_D) / 2.0;
		ResultGreeks[4] = (Premium_Value_U - Premium_Value_D) / 2.0;
		ResultGreeks[5] = (Protection_Value_U - Protection_Value_D) / 2.0;

	}
	if (HazardRate) free(HazardRate);
	if (PremiumCurveRup) free(PremiumCurveRup);
	if (PremiumCurveRdn) free(PremiumCurveRdn);
	if (ProtectionCurveRup) free(ProtectionCurveRup);
	if (ProtectionCurveRdn) free(ProtectionCurveRdn);
	if (HazardUp) free(HazardUp);
	if (HazardDn) free(HazardDn);
	if (Temp_Premium_Schedule) free(Temp_Premium_Schedule);

	return ResultSpread;
}

long Inputcheck_Calc_CDS_Spread(
	long NPremiumCurve,
	double* PremiumCurveTerm,
	double* PremiumCurve,

	long NProtectionCurve,
	double* ProtectionCurveTerm,
	double* ProtectionCurve,

	long DataFlag,
	long NHazardRate,
	double* HazardRateTerm,
	double* HazardRate,
	long CalcQMethod,

	long NRefZero,
	double* RefZeroTerm,
	double* ZeroRiskFree,
	double* ZeroRiskyBond,
	long DefaultRateCalcMethod,

	double Recovery,
	long NCPN_Ann,
	double Maturity,
	double PrevCDSSpread,
	long CDSPosition,
	long TextFlag,
	long GreekFlag)
{
	long i;
	long NPremiumCurve_Error = -1, PremiumCurveTerm_Error = -2, PremiumCurve_Error = -3;
	long NProtectionCurve_Error = -4, ProtectionCurveTerm_Error = -5, ProtectionCurve_Error = -6;
	long DataFlag_Error = -7, NHazardRate_Error = -8, HazardRateTerm_Error = -9, HazardRate_Error = -10, CalcQMethod_Error = -11;
	long NRefZero_Error = -12, RefZeroTerm_Error = -13, ZeroRiskFree_Error = -14, ZeroRiskyBond_Error = -15, DefaultRateCalcMethod_Error = -16;
	long Recovery_Error = -17, NCPN_Ann_Error = -18, Maturity_Error = -19;

	if (NPremiumCurve < 1) return NPremiumCurve_Error;
	if (PremiumCurveTerm[0] < 0.0) return PremiumCurveTerm_Error;
	for (i = 1; i < NPremiumCurve; i++)
	{
		if (PremiumCurveTerm[i] < PremiumCurveTerm[i - 1]) return PremiumCurveTerm_Error;
	}

	if (NProtectionCurve < 1) return NProtectionCurve_Error;
	if (ProtectionCurveTerm[0] < 0.0) return ProtectionCurveTerm_Error;
	for (i = 1; i < NProtectionCurve; i++)
	{
		if (ProtectionCurveTerm[i] < ProtectionCurveTerm[i - 1]) return ProtectionCurveTerm_Error;
	}

	if (DataFlag == 0)
	{
		if (NHazardRate < 1) return NHazardRate_Error;

		if (HazardRateTerm[0] < 0.0) return HazardRateTerm_Error;

		for (i = 0; i < NHazardRate; i++)
		{
			if (HazardRateTerm[i] < HazardRateTerm[i - 1]) return HazardRateTerm_Error;
		}

	}
	else
	{
		if (NRefZero < 1) return (NRefZero_Error);

		if (RefZeroTerm[0] < 0.0) return RefZeroTerm_Error;

		for (i = 0; i < NRefZero; i++)
		{
			if (RefZeroTerm[i] < RefZeroTerm[i - 1]) return RefZeroTerm_Error;
			if (ZeroRiskFree[i] > ZeroRiskyBond[i]) return RefZeroTerm_Error;
		}
	}
	if (Recovery < 0.0 || Recovery > 1.0) return Recovery_Error;
	if (PrevCDSSpread < 0.0) return -20;

	return 0;
}

DLLEXPORT(long) Calc_CDS_Spread(
	long NPremiumCurve,				// Premium Leg 커브 Term 개수
	double* PremiumCurveTerm,		// Premium Leg 커브 Term Array
	double* PremiumCurve,			// Premium Leg 커브 Rate Array

	long NProtectionCurve,			// Protection Leg 커브 Term 개수
	double* ProtectionCurveTerm,	// Protection Leg 커브 Term
	double* ProtectionCurve,		// Protection Leg 커브 Rate Array

	long DataFlag,					// Credit 인풋 유형 0: Hazard Rate, 1: 위험,무위험금리
	long NHazardRate,				// Hazard Rate Term 개수
	double* HazardRateTerm,			// Hazard Rate Term Array
	double* HazardRate,				// Hazard Rate Array
	long CalcQMethod,				// 0:생존율 = e^(-lambda * t) 계산, 1:생존률 = 1 - lambda * t 계산

	long NRefZero,					// CalcQMethod == 1일 때 제로금리개수
	double* RefZeroTerm,			// ZeroRate Term Array
	double* ZeroRiskFree,			// RiskFree ZeroRate Array
	double* ZeroRiskyBond,			// Risky ZeroRate Array
	long DefaultRateCalcMethod,		// HazardRate 산출방식 0:RiskyBond Pricer 1: 간단산출

	double Recovery,				// 회수율
	long ScheduleInputFlag,			// 스케줄입력방식 0: 간단입력, 1: 직접입력
	long NCPN_Ann,					// ScheduleInputFlag0: 연 이자지급 수
	double Maturity,				// ScheduleInputFlag0: 만기
	long NSchedule,					// ScheduleInputFlag1: 스케줄개수
	long* ResetDateExcelDate,		// ScheduleInputFlag1: 리셋일 엑셀타입 Array
	long* PayDateExcelDate,			// ScheduleInputFlag1: 지급일 엑셀타입 Array
	long PricingDateExcelDate,		// ScheduleInputFlag1: 가격계산일 엑셀타입

	long TextFlag,					// 텍스트DumpFlag (미완성)
	long GreekFlag,					// Greek산출Flag
	double* Result_Value,			// 결과 0: Result_Spread 1:Premium_Leg 2:Protection Leg
	double* ResultHazardRate,		// 산출된 Hazard Rate Array
	double* Premium_Schedule,		// 프리미엄스케줄 Array
	double NotionalAmount			// 액면가액
)
{
	double Premium_Leg = 0.0;
	double Protection_Leg = 0.0;
	// CDS Spread 계산하는 함수는 GreekFlag 0으로 고정
	GreekFlag = 0;
	double ResultGreeks[10] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };

	long i;

	long ResultCode = Inputcheck_Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve, NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
		DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
		NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, DefaultRateCalcMethod,
		Recovery, NCPN_Ann, Maturity, 0.01, 1,
		TextFlag, GreekFlag);

	if (ResultCode < 0) return ResultCode;

	double Result_Spread;

	if (DataFlag == 0)
	{
		Result_Spread = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Premium_Leg, &Protection_Leg, Premium_Schedule, NotionalAmount);
		if (CalcQMethod == 0)
		{
			for (i = 0; i < NHazardRate; i++)
			{
				ResultHazardRate[i] = 1.0 - exp(-HazardRate[i] * HazardRateTerm[i]);
			}
		}
		else
		{
			for (i = 0; i < NHazardRate; i++)
			{
				ResultHazardRate[i] = HazardRate[i] * HazardRateTerm[i];
			}
		}
	}
	else
	{
		Result_Spread = Calc_CDS_From_CreditSpread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, DefaultRateCalcMethod,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Premium_Leg, &Protection_Leg, ResultHazardRate, Premium_Schedule, GreekFlag, ResultGreeks, 0.0, NotionalAmount);
	}

	Result_Value[0] = Result_Spread;
	Result_Value[1] = Premium_Leg;
	Result_Value[2] = Protection_Leg;
	Result_Value[3] = 0.0;
	Result_Value[4] = 0.0;
	Result_Value[5] = 0.0;

	return 1;
}



DLLEXPORT(long) Calc_OLD_CDS_Value(
	long NPremiumCurve,				// Premium Leg 커브 Term 개수
	double* PremiumCurveTerm,		// Premium Leg 커브 Term Array
	double* PremiumCurve,			// Premium Leg 커브 Rate Array

	long NProtectionCurve,			// Protection Leg 커브 Term 개수
	double* ProtectionCurveTerm,	// Protection Leg 커브 Term
	double* ProtectionCurve,		// Protection Leg 커브 Rate Array

	long DataFlag,					// Credit 인풋 유형 0: Hazard Rate,1: 위험,무위험금리
	long NHazardRate,				// Hazard Rate Term 개수
	double* HazardRateTerm,			// Hazard Rate Term Array
	double* HazardRate,				// Hazard Rate Array
	long CalcQMethod,				// 0:생존율 e^(-H * t) 계산, 1:생존률 1 - H * t 계산

	long NRefZero,					// CalcQMethod == 1일 때 제로금리개수
	double* RefZeroTerm,			// ZeroRate Term Array
	double* ZeroRiskFree,			// RiskFree ZeroRate Array
	double* ZeroRiskyBond,			// Risky ZeroRate Array
	long DefaultRateCalcMethod,		// HazardRate 산출방식 0:RiskyBond Pricer 1: 간단산출

	double Recovery,				// 회수율
	long ScheduleInputFlag,			// 스케줄입력방식 0: 간단입력, 1: 직접입력
	long NCPN_Ann,					// ScheduleInputFlag0: 연 이자지급 수
	double Maturity,				// ScheduleInputFlag0: 만기
	long NSchedule,					// ScheduleInputFlag1: 스케줄개수
	long* ResetDateExcelDate,		// ScheduleInputFlag1: 리셋일 엑셀타입 Array
	long* PayDateExcelDate,			// ScheduleInputFlag1: 지급일 엑셀타입 Array
	long PricingDateExcelDate,		// ScheduleInputFlag1: 가격계산일 엑셀타입


	double OLD_CDS_Spread,			// 예전 계약 CDS Spread
	long Protection_Position,		// 포지션1: 롱 -1: 숏
	long TextFlag,					// TextFlag (미완성)
	long GreekFlag,					// Greek산출Flag
	double* Result_Value,			// 결과 [0]: Result_Spread [1]:Premium_Leg [2]:Protection Leg
	double* ResultHazardRate,		// 산출된 Hazard Rate Array
	double* ResultGreeks,			// [0~2]IR PV01 Net, Pre, Pro [3~5]부도율 PV01 Net, Pre, Pro 
	double* Premium_Schedule,		// 프리미엄스케줄 Array
	double NotionalAmount			// 액면가액
)
{

	double Premium_Leg = 0.0;
	double Protection_Leg = 0.0;
	// Greeks용 임시변수
	double Temp_Premium_Leg;
	double Temp_Protection_Leg;

	double* Temp_Premium_Schedule;
	double* PremiumCurveRup;
	double* PremiumCurveRdn;
	double* ProtectionCurveRup;
	double* ProtectionCurveRdn;
	double* HazardUp;
	double* HazardDn;
	//
	double RPV01;
	double Protection_Value;
	double Premium_Value;
	double Value_ProtectionBuyer;
	long N_PremiumNode = (long)((double)NCPN_Ann * Maturity + 0.5) + 1;

	double P_Rup;
	double P_Rdn;

	long i;


	long ResultCode = Inputcheck_Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve, NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
		DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
		NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, DefaultRateCalcMethod,
		Recovery, NCPN_Ann, Maturity, OLD_CDS_Spread, Protection_Position,
		TextFlag, GreekFlag);

	if (ResultCode < 0) return ResultCode;

	double Result_Spread;
	double Protection_Value_U;
	double Protection_Value_D;
	// Greeks용 임시변수
	double Premium_Value_U;
	double Premium_Value_D;
	double Value_ProtectionBuyer_U;
	double Value_ProtectionBuyer_D;
	double RPV01_U;
	double RPV01_D;
	
	//////////////////////////////////////////////////////////////////////
	// DataFlag0 : Hazard Rate 입력되어있음, DataFlag1: ZeroRate Risky - Rf 입력되어있음
	//////////////////////////////////////////////////////////////////////

	if (DataFlag == 0)
	{
		///////////////////////////////////
		// Value 계산을 위해 현재기준 Premium_Leg, Protection_Leg, CDS Spread 계산
		///////////////////////////////////

		Result_Spread = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Premium_Leg, &Protection_Leg, Premium_Schedule, NotionalAmount);
		if (CalcQMethod == 0)
		{
			for (i = 0; i < NHazardRate; i++)
			{
				ResultHazardRate[i] = 1.0 - exp(-HazardRate[i] * HazardRateTerm[i]);
			}
		}
		else
		{
			for (i = 0; i < NHazardRate; i++)
			{
				ResultHazardRate[i] = HazardRate[i] * HazardRateTerm[i];
			}
		}


		if (GreekFlag >= 1)
		{
			//////////////////////////////////////////////////////////////////////
			// OLD CDS 의 그릭스 계산 -> Value = ProtectionValue - 구스프레드/현스프레드 * 현PremiumValue
			//////////////////////////////////////////////////////////////////////

			Temp_Premium_Schedule = (double*)malloc(sizeof(double) * (long)(NCPN_Ann * Maturity + 2.0));
			PremiumCurveRup = (double*)malloc(sizeof(double) * NPremiumCurve);
			PremiumCurveRdn = (double*)malloc(sizeof(double) * NPremiumCurve);
			ProtectionCurveRup = (double*)malloc(sizeof(double) * NProtectionCurve);
			ProtectionCurveRdn = (double*)malloc(sizeof(double) * NProtectionCurve);
			HazardUp = (double*)malloc(sizeof(double) * NHazardRate);
			HazardDn = (double*)malloc(sizeof(double) * NHazardRate);

			for (i = 0; i < NPremiumCurve; i++)
			{
				PremiumCurveRup[i] = PremiumCurve[i] + 0.0001;
				PremiumCurveRdn[i] = PremiumCurve[i] - 0.0001;
			}

			for (i = 0; i < NProtectionCurve; i++)
			{
				ProtectionCurveRup[i] = ProtectionCurve[i] + 0.0001;
				ProtectionCurveRdn[i] = ProtectionCurve[i] - 0.0001;
			}

			for (i = 0; i < NHazardRate; i++)
			{
				HazardUp[i] = HazardRate[i] + 0.0001;
				HazardDn[i] = HazardRate[i] - 0.0001;
			}

			P_Rup = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurveRup,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurveRup,
				DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

			Protection_Value_U = Temp_Protection_Leg;
			RPV01_U = Temp_Premium_Leg / P_Rup;
			Premium_Value_U = OLD_CDS_Spread * RPV01_U;
			Value_ProtectionBuyer_U = Protection_Value_U - Premium_Value_U;

			P_Rdn = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurveRdn,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurveRdn,
				DataFlag, NHazardRate, HazardRateTerm, HazardRate, CalcQMethod,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);


			Protection_Value_D = Temp_Protection_Leg;
			RPV01_D = Temp_Premium_Leg / P_Rdn;
			Premium_Value_D = OLD_CDS_Spread * RPV01_D;
			Value_ProtectionBuyer_D = Protection_Value_D - Premium_Value_D;

			ResultGreeks[0] = (Value_ProtectionBuyer_U - Value_ProtectionBuyer_D) / 2.0;
			ResultGreeks[1] = (Premium_Value_U - Premium_Value_D) / 2.0;
			ResultGreeks[2] = (Protection_Value_U - Protection_Value_D) / 2.0;

			P_Rup = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				DataFlag, NHazardRate, HazardRateTerm, HazardUp, CalcQMethod,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

			Protection_Value_U = Temp_Protection_Leg;
			RPV01_U = Temp_Premium_Leg / P_Rup;
			Premium_Value_U = OLD_CDS_Spread * RPV01_U;
			Value_ProtectionBuyer_U = Protection_Value_U - Premium_Value_U;

			P_Rdn = Calc_CDS_From_Hazard(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				DataFlag, NHazardRate, HazardRateTerm, HazardDn, CalcQMethod,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Temp_Premium_Leg, &Temp_Protection_Leg, Temp_Premium_Schedule, NotionalAmount);

			Protection_Value_D = Temp_Protection_Leg;
			RPV01_D = Temp_Premium_Leg / P_Rdn;
			Premium_Value_D = OLD_CDS_Spread * RPV01_D;
			Value_ProtectionBuyer_D = Protection_Value_D - Premium_Value_D;

			ResultGreeks[3] = (Value_ProtectionBuyer_U - Value_ProtectionBuyer_D) / 2.0;
			ResultGreeks[4] = (Premium_Value_U - Premium_Value_D) / 2.0;
			ResultGreeks[5] = (Protection_Value_U - Protection_Value_D) / 2.0;

			if (Temp_Premium_Schedule) free(Temp_Premium_Schedule);
			if (PremiumCurveRup) free(PremiumCurveRup);
			if (PremiumCurveRdn) free(PremiumCurveRdn);
			if (ProtectionCurveRup) free(ProtectionCurveRup);
			if (ProtectionCurveRdn) free(ProtectionCurveRdn);

		}
	}
	else
	{
		Result_Spread = Calc_CDS_From_CreditSpread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
			NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
			DataFlag, NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, DefaultRateCalcMethod,
			Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, &Premium_Leg, &Protection_Leg, ResultHazardRate, Premium_Schedule, GreekFlag, ResultGreeks, OLD_CDS_Spread, NotionalAmount);
	}

	Result_Value[0] = Result_Spread;
	Result_Value[1] = Premium_Leg;
	Result_Value[2] = Protection_Leg;

	Protection_Value = Protection_Leg;
	RPV01 = Premium_Leg / Result_Spread;
	Premium_Value = OLD_CDS_Spread * RPV01;
	for (i = 0; i < N_PremiumNode - 1; i++)
	{
		Premium_Schedule[i] = OLD_CDS_Spread / Result_Spread * Premium_Schedule[i];
	}
	Result_Value[3] = Premium_Value;
	Result_Value[4] = Protection_Value;
	Value_ProtectionBuyer = Protection_Value - Premium_Value;
	if (Protection_Position == 1)
		Result_Value[5] = Value_ProtectionBuyer;
	else
		Result_Value[5] = -Value_ProtectionBuyer;


	return 1;
}

long Inputcheck_Hazard_Rate_From_CDSCurve(
	long NPremiumCurve,
	double* PremiumCurveTerm,
	double* PremiumCurve,

	long NProtectionCurve,
	double* ProtectionCurveTerm,
	double* ProtectionCurve,

	long NCDSCurve,
	double* CDSCurveTerm,
	double* CDSCurve,

	double Recovery,
	long HazardRateCalcFlag,
	long NCPN_Ann)
{
	long i;
	long NPremiumCurve_Error = -1, PremiumCurveTerm_Error = -2, PremiumCurve_Error = -3;
	long NProtectionCurve_Error = -4, ProtectionCurveTerm_Error = -5, ProtectionCurve_Error = -6;
	long NCDSCurve_Error = -7, CDSCurveTerm_Error = -8, CDSCurve_Error = -9;
	long Recovery_Error = -10, HazardRateCalcFlag_Error = -11, NCPN_Ann_Error = -12, TextFlag_Error = -13, GreekFlag_Error = -14;

	if (NPremiumCurve < 1) return NPremiumCurve_Error;
	if (PremiumCurveTerm[0] < 0.0) return PremiumCurveTerm_Error;
	for (i = 1; i < NPremiumCurve; i++)
	{
		if (PremiumCurveTerm[i] < PremiumCurveTerm[i - 1]) return PremiumCurveTerm_Error;
	}

	if (NProtectionCurve < 1) return NProtectionCurve_Error;
	if (ProtectionCurveTerm[0] < 0.0) return ProtectionCurveTerm_Error;
	for (i = 1; i < NProtectionCurve; i++)
	{
		if (ProtectionCurveTerm[i] < ProtectionCurveTerm[i - 1]) return ProtectionCurveTerm_Error;
	}

	if (NCDSCurve < 1) return NCDSCurve_Error;
	if (CDSCurveTerm[0] < 0.0) return ProtectionCurveTerm_Error;
	for (i = 1; i < NCDSCurve; i++)
	{
		if (CDSCurveTerm[i] < CDSCurveTerm[i - 1]) return CDSCurveTerm_Error;
	}

	if (Recovery < 0.0 || Recovery > 1.0) return Recovery_Error;

	return 1;
}

DLLEXPORT(long) Calc_Hazard_Rate_From_CDSCurve(
	long NPremiumCurve,				// Premium Leg 커브 Term 개수
	double* PremiumCurveTerm,		// Premium Leg 커브 Term Array
	double* PremiumCurve,			// Premium Leg 커브 Rate Array

	long NProtectionCurve,			// Protection Leg 커브 Term 개수
	double* ProtectionCurveTerm,	// Protection Leg 커브 Term
	double* ProtectionCurve,		// Protection Leg 커브 Rate Array

	long NCDSCurve,					// CDS Term 개수
	double* CDSCurveTerm,			// CDS Term Array
	double* CDSCurve,				// CDS Spread Array

	double Recovery,				// Recovery Rate
	long HazardRateCalcFlag,		// Hazard Rate Calc 방법 0: Continuous Annual Hazard Rate 1: 단리 파산확률
	long NCPN_Ann,					// 연 이자지급 수			

	double* ResultHazardTerm,		// OutPut Hazard Term
	double* ResultHazard			// OutPut Hazard Rate
)
{
	long i, j;

	////////////////
	// 에러체크
	////////////////

	long ResultCode = Inputcheck_Hazard_Rate_From_CDSCurve(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
		NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
		NCDSCurve, CDSCurveTerm, CDSCurve,
		Recovery, HazardRateCalcFlag, NCPN_Ann);

	if (ResultCode < 0) return ResultCode;

	double ObjValue;
	double value_a;
	double value_b;
	double value;
	double dblErrorRange = 0.00001;

	long NHazardRate;
	double* HazardRate = (double*)malloc(sizeof(double) * NCDSCurve);
	double* HazardRateTerm = CDSCurveTerm;
	double* ResultHazardRate = (double*)malloc(sizeof(double) * NCDSCurve);

	double Price1, Price2, Error1, Error2, Price, Error;
	double Maturity;
	double* Premium_Schedule = (double*)malloc(sizeof(double) * ((long)(NCPN_Ann * CDSCurveTerm[NCDSCurve - 1] + 0.5) + 2));
	long TempNumber = 0;
	long EndFlag = 0;

	////////////////
	// 더미 변수들
	////////////////

	long NRefZero = 1;
	double RefZeroTerm[1] = { 0.0 };
	double ZeroRiskFree[1] = { 0.0 };
	double ZeroRiskyBond[1] = { 0.0 };

	double TempResultPrice[10] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
	double MinimumValue = 0.0;
	long ScheduleInputFlag = 0;
	long NSchedule = 0;
	long* ResetDateExcelDate = &NSchedule;
	long* PayDateExcelDate = &NSchedule;
	long PricingDateExcelDate = 20000101;

	////////////////

	for (i = 0; i < NCDSCurve; i++)
		ResultHazardTerm[i] = CDSCurveTerm[i];

	for (i = 0; i < NCDSCurve; i++)
	{
		value_a = 0.00001;
		value_b = 0.99;

		if (i > 0) MinimumValue = value * Maturity / CDSCurveTerm[i]; // e^(-lambda2 * T2) < e^(-lambda1 * T1) 생존확률은 우하향해야함
		
		EndFlag = 0;
		NHazardRate = i + 1;
		Maturity = CDSCurveTerm[i];
		ObjValue = CDSCurve[i];
		HazardRate[i] = value_a;
		if (HazardRateCalcFlag == 0)
		{
			TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				0, NHazardRate, HazardRateTerm, HazardRate, 0,
				NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
		}
		else
		{
			TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				0, NHazardRate, HazardRateTerm, HazardRate, 1,
				NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
		}

		Price1 = TempResultPrice[0];
		Error1 = ObjValue - Price1;
		if (fabs(Error1) < dblErrorRange)
		{
			ResultHazard[i] = value_a;
			EndFlag = 1;
		}

		HazardRate[i] = value_b;
		if (HazardRateCalcFlag == 0)
		{
			TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				0, NHazardRate, HazardRateTerm, HazardRate, 0,
				NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
		}
		else
		{
			TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
				NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
				0, NHazardRate, HazardRateTerm, HazardRate, 1,
				NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
				Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
		}
		Price2 = TempResultPrice[0];
		Error2 = ObjValue - Price2;
		if (fabs(Error2) < dblErrorRange)
		{
			ResultHazard[i] = value_b;
			EndFlag = 1;
		}

		value = (value_a + value_b) / 2.0;
		for (j = 0; j < 2000; j++)
		{
			HazardRate[i] = value;
			if (HazardRateCalcFlag == 0)
			{
				TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
					NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
					0, NHazardRate, HazardRateTerm, HazardRate, 0,
					NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
					Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
			}
			else
			{
				TempNumber = Calc_CDS_Spread(NPremiumCurve, PremiumCurveTerm, PremiumCurve,
					NProtectionCurve, ProtectionCurveTerm, ProtectionCurve,
					0, NHazardRate, HazardRateTerm, HazardRate, 1,
					NRefZero, RefZeroTerm, ZeroRiskFree, ZeroRiskyBond, 0,
					Recovery, ScheduleInputFlag, NCPN_Ann, Maturity, NSchedule, ResetDateExcelDate, PayDateExcelDate, PricingDateExcelDate, 0, 0, TempResultPrice, ResultHazardRate, Premium_Schedule, 1.0);
			}

			Price = TempResultPrice[0];
			Error = ObjValue - Price;
			if (fabs(Error) < dblErrorRange)
			{
				ResultHazard[i] = value;
				break;
			}
			else
			{
				if (Error > 0)
				{
					value_a = max(value, MinimumValue);
					value = (value_a + value_b) / 2.0;
				}
				else
				{
					value_b = value;
					value = (value_a + value_b) / 2.0;
				}
			}
			if (TempNumber < 0 && i > 0)
			{
				ResultHazard[i] = ResultHazard[i - 1];
				break;
			}
		}
		if (j == 2000 && i > 0) ResultHazard[i] = ResultHazard[i - 1];
	}

	if (HazardRate) free(HazardRate);
	if (ResultHazardRate) free(ResultHazardRate);
	if (Premium_Schedule) free(Premium_Schedule);

	return 1;
}