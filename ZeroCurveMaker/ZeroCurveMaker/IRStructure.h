#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CalcDate.h"
#include "Structure.h"
#include "Util.h"

class MoneyMarket {
private: 
	
public: 
	long Product_Type;      // 상품타입 0: Overnight 이자율, 1: 만기 30일 미만 이자율 2: 만기 30일 이상 이자율
	long Maturity_Days;     // 만기까지 일수
	double Interest_Rate;   // 이자율
	long OneYearConvention; // 360일이면 360입력 365면 365입력

	double Discount_Factor;
	double ZeroRate;

	MoneyMarket(){}

	MoneyMarket(long type, long PricingDate, long MaturityDate, double rate)
	{
		Product_Type = type;
		Maturity_Days = DayCountAtoB(PricingDate, MaturityDate);
		Interest_Rate = rate;
		OneYearConvention = 365;
		double T = ((double)(Maturity_Days)) / ((double)(OneYearConvention));
		Discount_Factor = 1.0 / (1.0 + Interest_Rate * T);
		ZeroRate = -1.0 / T * log(Discount_Factor);
	}


	MoneyMarket(long type, long PricingDate, long MaturityDate, double rate, long YearConvention)
	{
		Product_Type = type;
		Maturity_Days = DayCountAtoB(PricingDate, MaturityDate);
		Interest_Rate = rate;
		OneYearConvention = YearConvention;
		double T = (double)(Maturity_Days) / (double)(OneYearConvention);
		Discount_Factor = 1.0 / (1.0 + Interest_Rate * T);
		ZeroRate = -1.0 / T * log(Discount_Factor);
	}

	void initialize(long type, long PricingDate, long MaturityDate, double rate, long YearConvention)
	{
		Product_Type = type;
		Maturity_Days = DayCountAtoB(PricingDate, MaturityDate);
		Interest_Rate = rate;
		OneYearConvention = YearConvention;
		double T = (double)(Maturity_Days) / (double)(OneYearConvention);
		Discount_Factor = 1.0 / (1.0 + Interest_Rate * T);
		ZeroRate = -1.0 / T * log(Discount_Factor);
	}
};

class IR_Forward {
private:

public:
	double forward_rate;
	long forward_start_day;
	long forward_end_day;
	long OneYearConvention; // 360일이면 360입력 365면 365입력
	double forward_discount_factor;

	IR_Forward(){}

	IR_Forward(long PricingDate, long ForwardStartDate, long ForwardEndDate, double F_Rate, long YearConvention)
	{
		double T;

		forward_rate = F_Rate;
		forward_start_day = DayCountAtoB(PricingDate, ForwardStartDate);
		forward_end_day = DayCountAtoB(PricingDate, ForwardEndDate); 
		OneYearConvention = YearConvention;
		
		T = ((double)(forward_end_day - forward_start_day)) / (double)(OneYearConvention);
		
		if (T <= 1.0) forward_discount_factor = 1.0 / (1.0 + forward_rate * T);
		else forward_discount_factor = exp(-forward_rate * T);
	}

	IR_Forward(long PricingDate, long ForwardStartDate, long ForwardEndDate, double F_Rate)
	{
		double T;

		forward_rate = F_Rate;
		forward_start_day = DayCountAtoB(PricingDate, ForwardStartDate);
		forward_end_day = DayCountAtoB(PricingDate, ForwardEndDate);
		OneYearConvention = 365;
		
		T = ((double)(forward_end_day - forward_start_day)) / (double)(OneYearConvention);
		
		if (T <= 1.0) forward_discount_factor = 1.0 / (1.0 + forward_rate * T);
		else forward_discount_factor = exp(-forward_rate * T);
	}

	void initialize(long PricingDate, long ForwardStartDate, long ForwardEndDate, double F_Rate, long YearConvention)
	{
		double T;

		forward_rate = F_Rate;
		forward_start_day = DayCountAtoB(PricingDate, ForwardStartDate);
		forward_end_day = DayCountAtoB(PricingDate, ForwardEndDate);
		OneYearConvention = YearConvention;

		T = ((double)(forward_end_day - forward_start_day)) / (double)(OneYearConvention);

		if (T <= 1.0) forward_discount_factor = 1.0 / (1.0 + forward_rate * T);
		else forward_discount_factor = exp(-forward_rate * T);
	}

	IR_Forward(curveinfo& Curve, long PricingDate, long ForwardStartDate, long ForwardEndDate, long YearConvention)
	{
		double T;
		double T1;
		double T2;
		long LengthArray;

		LengthArray = Curve.nterm();

		forward_start_day = DayCountAtoB(PricingDate, ForwardStartDate);
		forward_end_day = DayCountAtoB(PricingDate, ForwardEndDate);
		OneYearConvention = YearConvention;
		
		T1 = (double)forward_start_day / (double)(OneYearConvention);
		T2 = (double)forward_end_day / (double)(OneYearConvention);
		T = T2 - T1;
		
		forward_rate = Calc_Forward_Rate(Curve.Term, Curve.Rate, LengthArray, T1, T2);

		if (T <= 1.0) forward_discount_factor = 1.0 / (1.0 + forward_rate * T);
		else forward_discount_factor = exp(-forward_rate * T);
	}

	void initialize(curveinfo& Curve, long PricingDate, long ForwardStartDate, long ForwardEndDate, long YearConvention)
	{
		double T;
		double T1;
		double T2;
		long LengthArray;

		LengthArray = Curve.nterm();

		forward_start_day = DayCountAtoB(PricingDate, ForwardStartDate);
		forward_end_day = DayCountAtoB(PricingDate, ForwardEndDate);
		OneYearConvention = YearConvention;

		T1 = (double)forward_start_day / (double)(OneYearConvention);
		T2 = (double)forward_end_day / (double)(OneYearConvention);
		T = T2 - T1;

		forward_rate = Calc_Forward_Rate(Curve.Term, Curve.Rate, LengthArray, T1, T2);

		if (T <= 1.0) forward_discount_factor = 1.0 / (1.0 + forward_rate * T);
		else forward_discount_factor = exp(-forward_rate * T);
	}
};

class IR_Swap {
private:
	long DynamicFlag = 0;
public:
	double SwapRate;
	long Maturity_Day;
	double Maturity;
	long NResetTime;
	long* ResetDays;
	double* DeltaT;
	long N_Annual_Cpn;
	long OneYearConvention;

	IR_Swap() {}

	IR_Swap(long PricingDate, long SwapMaturityDate, long Number_Ann_CPN, double Rate, long Convention1Y)
	{
		long i;
		long n;

		double Accrued_dt;
		double LeftTime;

		DynamicFlag = 1;

		SwapRate = Rate;
		N_Annual_Cpn = Number_Ann_CPN;
		OneYearConvention = Convention1Y;
		Maturity_Day = DayCountAtoB(PricingDate,SwapMaturityDate);

		Maturity = ((double)DayCountAtoB(PricingDate, SwapMaturityDate))/(double)Convention1Y;
		
		for (i = 0; i < (long)(Maturity + 10.0) * N_Annual_Cpn; i++)
		{
			LeftTime = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			if (LeftTime < 0.00000000001)
			{
				NResetTime = i;
				Accrued_dt = -LeftTime;
				break;
			}
		}
		
		ResetDays = (long*)malloc(sizeof(long) * NResetTime);
		DeltaT= (double*)malloc(sizeof(double) * NResetTime);

		for (i = 0; i < NResetTime; i++)
		{
			ResetDays[NResetTime - 1 - i] = (long)((Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn)) * 365.0);
		}

		DeltaT[0] = ((double)ResetDays[0]) / (double)OneYearConvention;
		for (i = 1; i < NResetTime; i++)
		{
			DeltaT[i] = ((double)ResetDays[i] - (double)ResetDays[i - 1]) / (double)OneYearConvention;
		}
	}

	void initialize(long PricingDate, long SwapMaturityDate, long Number_Ann_CPN, double Rate, long Convention1Y)
	{
		long i;
		long n;

		double Accrued_dt;
		double LeftTime;

		DynamicFlag = 1;

		SwapRate = Rate;
		N_Annual_Cpn = Number_Ann_CPN;
		OneYearConvention = Convention1Y;
		Maturity_Day = DayCountAtoB(PricingDate, SwapMaturityDate);

		Maturity = ((double)DayCountAtoB(PricingDate, SwapMaturityDate)) / (double)Convention1Y;

		for (i = 0; i < (long)(Maturity + 10.0) * N_Annual_Cpn; i++)
		{
			LeftTime = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			if (LeftTime < 0.00000000001)
			{
				NResetTime = i;
				Accrued_dt = -LeftTime;
				break;
			}
		}

		ResetDays = (long*)malloc(sizeof(long) * NResetTime);
		DeltaT = (double*)malloc(sizeof(double) * NResetTime);

		for (i = 0; i < NResetTime; i++)
		{
			ResetDays[NResetTime - 1 - i] = (long)((Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn)) * 365.0 );
		}

		DeltaT[0] = ((double)ResetDays[0]) / (double)OneYearConvention;
		for (i = 1; i < NResetTime; i++)
		{
			DeltaT[i] = ((double)ResetDays[i] - (double)ResetDays[i - 1]) / (double)OneYearConvention;
		}
	}

	~IR_Swap()
	{
		if (DynamicFlag != 0)
		{
			free(ResetDays);
			free(DeltaT);
		}
	}
};

double Accrued_Interest(double Notional, double CpnRate, double Accrued_dt)
{
	return Notional * CpnRate * Accrued_dt;
}

class CPN_Bond {
private:
	long DynamicFlag = 0;
public:
	double PV;
	double YTM;
	double NA;
	long Maturity_Day;
	double Maturity;
	double CPNRate;
	long NCPNTime;
	long N_Annual_Cpn;
	long OneYearConvention;
	double AccruedInterest;

	CPN_Bond() {}

	CPN_Bond(long PricingDate, long BondMaturityDate, long Number_Ann_CPN, double Yield_To_Maturity, double CouponRate, long Convention1Y)
	{
		long i;
		long n;

		double Accrued_dt;
		double LeftTime;
		NA = 100.0;
		DynamicFlag = 1;

		YTM = Yield_To_Maturity;
		N_Annual_Cpn = Number_Ann_CPN;
		OneYearConvention = Convention1Y;
		Maturity_Day = DayCountAtoB(PricingDate, BondMaturityDate);
		CPNRate = CouponRate;
		Maturity = ((double)DayCountAtoB(PricingDate, BondMaturityDate)) / (double)Convention1Y;

		for (i = 0; i < (long)(Maturity + 10.0) * N_Annual_Cpn; i++)
		{
			LeftTime = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			if (LeftTime < 0.0001)
			{
				NCPNTime = i;
				Accrued_dt = -LeftTime;
				break;
			}
		}

		double T = 0.0;
		double DF = 1.0;
		PV = 0.0;
		for (i = 0; i < NCPNTime; i++)
		{
			T = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			DF = 1.0 / pow((1.0 + YTM), T);
			PV += 1.0 / (double)N_Annual_Cpn * NA * CPNRate * DF;
		}
		PV += NA * DF;
		AccruedInterest = Accrued_Interest(NA, CPNRate, Accrued_dt);
	}

	void initialize(long PricingDate, long BondMaturityDate, long Number_Ann_CPN, double Yield_To_Maturity, double CouponRate, long Convention1Y)
	{
		long i;
		long n;

		double Accrued_dt;
		double LeftTime;
		NA = 100.0;
		DynamicFlag = 1;

		YTM = Yield_To_Maturity;
		N_Annual_Cpn = Number_Ann_CPN;
		OneYearConvention = Convention1Y;
		Maturity_Day = DayCountAtoB(PricingDate, BondMaturityDate);
		CPNRate = CouponRate;
		Maturity = ((double)DayCountAtoB(PricingDate, BondMaturityDate)) / (double)Convention1Y;

		for (i = 0; i < (long)(Maturity + 10.0) * N_Annual_Cpn; i++)
		{
			LeftTime = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			if (LeftTime < 0.0001)
			{
				NCPNTime = i;
				Accrued_dt = -LeftTime;
				break;
			}
		}

		double T = 0.0;
		double DF = 1.0;
		double dt = 1.0 / (double)N_Annual_Cpn;
		PV = 0.0;
		for (i = 0; i < NCPNTime; i++)
		{
			T = Maturity - ((double)i) * 1.0 / ((double)N_Annual_Cpn);
			DF = 1.0/pow(1.0 + CPNRate * dt, (double)N_Annual_Cpn * T);  //exp(-YTM * T);
			PV += dt * NA * CPNRate * DF;
		}
		PV += NA * 1.0 / pow(1.0 + CPNRate * dt, (double)N_Annual_Cpn * Maturity);
		AccruedInterest = Accrued_Interest(NA, CPNRate, Accrued_dt);
	}

	~CPN_Bond()
	{
		if (DynamicFlag != 0)
		{
		}
	}
};