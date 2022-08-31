#include <stdio.h>
#include <math.h>

#include "CalcDate.h"
#include "Structure.h"
#include "Util.h"

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

long sumation(long* array, long narray)
{
	long i;
	long s = 0;
	for (i = 0; i < narray; i++)
	{
		s += array[i];
	}
	return s;
}

double sumation(double* array, long narray)
{
	long i;
	double s = 0;
	for (i = 0; i < narray; i++)
	{
		s += array[i];
	}
	return s;
}

long isin(long x, long* array, long narray)
{
	long i;
	long s = 0;
	for (i = 0; i < narray; i++)
	{
		if (x == array[i])
		{
			s = 1;
			break;
		}
	}
	return s;
}

long isweekend(long ExlDate)
{
	// ������ 1�̸� �Ͽ���, 2�̸� ������, 3�̸� ȭ����, 4�̸� ������, 5�̸� �����, 6�̸� �ݿ���, 0�̸� �����
	long MOD7;
	if (ExlDate > 0)
	{
		MOD7 = ExlDate % 7;
		if (MOD7 == 1 || MOD7 == 0) return 1;
	}
	return 0;
}

double Continuous_ForwardRate(curveinfo& Curve, double T0, double T1)
{
	double r0, r1;
	double Result;
	if (T0 != T1)
	{
		r0 = Curve.Interpolated_Rate(T0);
		r1 = Curve.Interpolated_Rate(T1);
		Result = (r1 * T1 - r0 * T0) / (T1 - T0);
	}
	else
	{
		T1 = T0 + 1.0 / 365;
		r0 = Curve.Interpolated_Rate(T0);
		r1 = Curve.Interpolated_Rate(T1);
		Result = (r1 * T1 - r0 * T0) / (T1 - T0);
	}
	return Result;
}

long NPaymentSwap(double T_OptionMaturity, double T_SwapMaturity, double PayFreqOfMonth)
{
	double dT = (T_SwapMaturity - T_OptionMaturity + 0.00001);
	double Num_AnnPayment = 12.0 / PayFreqOfMonth;
	long N;
	N = (long)(dT * Num_AnnPayment + 0.5);
	return N;
}

long MappingPaymentDates(double T_SwapMaturity, double FreqMonth, long* TempDatesArray, long NDates)
{
	long i;
	long DayDiff = (long)(FreqMonth / 12.0 * 365.0);

	TempDatesArray[NDates - 1] = (long)(365.0 * T_SwapMaturity + 0.5);
	for (i = NDates - 2; i >= 0; i--)
	{
		TempDatesArray[i] = TempDatesArray[i + 1] - DayDiff;
	}
	return NDates;
}

/////////////////////////////
// Forward Swap Rate ��� ///
/////////////////////////////
double FSR(
	double* Term,			// �Ⱓ���� Term Array
	double* Rate,			// �Ⱓ���� Rate Array
	long NTerm,				// �Ⱓ���� Term ����
	double T_Option,		// Forward Swap ������
	double Tenor,			// ���� ����
	double FreqMonth		// Monthȯ�� ���� Frequency
)
{
	long i;
	double Swap_Rate;

	long ndates = NPaymentSwap(T_Option, T_Option + Tenor, FreqMonth);
	long* dates = (long*)malloc(sizeof(long) * (ndates));
	ndates = MappingPaymentDates(T_Option + Tenor, FreqMonth, dates, ndates);

	double* P0T = (double*)malloc(sizeof(double) * (ndates + 1));
	P0T[0] = Calc_Discount_Factor(Term, Rate, NTerm, T_Option); // �ɼǸ������
	for (i = 1; i < ndates + 1; i++)
	{
		P0T[i] = Calc_Discount_Factor(Term, Rate, NTerm, (double)dates[i - 1] / 365.0);
	}

	double a, b, dt;
	a = P0T[0] - P0T[ndates];
	b = 0.0;
	for (i = 0; i < ndates; i++)
	{
		if (i == 0) dt = (double)dates[0] / 365.0 - T_Option;
		else dt = ((double)(dates[i] - dates[i - 1])) / 365.0;
		b += dt * P0T[i + 1];
	}
	Swap_Rate = a / b;


	if (dates) free(dates);
	if (P0T) free(P0T);
	return Swap_Rate;
}

typedef struct schd_info {
	long PriceDate_C;			// PricingDate CType

	long ReferenceType;			// Reference Rate type
	long FixedFlotype;			// Fix or Flo Flag 0:Fix 1: Flo
	long DayCount;				// DayCountConvention 0:365 1:365
	double NotionalAmount;		// Notional Amount
	long NAFlag;				// Notional ���޿���

	long RefSwapFlag;			// ���۷����ݸ��� ���ұݸ���������
	long NSwapPayAnnual;		// ���۷����ݸ��� ���ұݸ���� �� �������� ��
	double RefSwapMaturity;		// ���۷��� �ݸ��� ���ұݸ���� ����

	long NCF;					// �����帧����
	long* ForwardStart_C;		// �ݸ�����������
	long* ForwardEnd_C;			// �ݸ�����������
	long* StartDate_C;			// Fraction ������(�����)
	long* EndDate_C;			// Fraction ������(�⸻��)
	long* PayDate_C;			// ������
	long NotionalPayDate_C;		// �׸�ݾ� ������
	long* Days_ForwardStart;	// ���� to ����������
	long* Days_ForwardEnd;		// ���� to ����������
	long* Days_StartDate;		// ���� to �����
	long* Days_EndDate;			// ���� to �⸻��
	long* Days_PayDate;			// ���� to ������

	long HolidayFlag_Ref;		// ���ʱݸ� Holiday Calc Flag
	long NHolidays_Ref;			// ���ʱݸ� Holiday ����
	long* Days_Holidays_Ref;	// ���ʱݸ� ���� to Holiday

	double* FixedRefRate;		// ���� Ȯ���ݸ� ������
	double* Slope;				// ���ʱݸ��� ���� ���̿��� ����
	double* CPN;				// ����������

	long Days_Notional;			// ���� to Notional ������

	long LockOutRef;			//  LockOut ��¥ N������
	long LookBackRef;			//  LookBack ��¥ 
	long ObservationShift;		//  Observation Shift �� ������

	long *NWeekendDate;			//  �ָ�����(���� = NCF)
	long** WeekendList;			//  �ָ� Array ����Ʈ

	long NRefHistory;
	long* RefHistoryDate;
	double* RefHistory;

} SCHD;

double Calc_Current_IRS(
	long CalcCRSFlag,
	curveinfo& Float_DiscCurve,
	curveinfo& Float_RefCurve,
	curveinfo& Fix_DiscCurve,
	curveinfo& FX_FloCurve,
	curveinfo& FX_FixCurve,
	SCHD* Float_Schedule,
	SCHD* Fix_Schedule
)
{
	long i;
	double IRSpread;
	double denominator;

	if (Float_Schedule->DayCount == 0) denominator = 365.0;
	else denominator = 360.0;

	double P0;
	double P1;
	double P_Pay;
	double T0;
	double T1;
	double TPay;
	double FloatValue, FixValue;
	double dt_Fix;
	double FX = 1.0;
	double NA = 1.0;

	P0 = 1.0;
	T0 = 0.0;
	FloatValue = 0.0;
	for (i = 0; i < Float_Schedule->NCF; i++)
	{
		T1 = ((double)(Float_Schedule->Days_ForwardEnd[i] - Float_Schedule->Days_ForwardStart[0]))/denominator;
		TPay = ((double)(Float_Schedule->Days_PayDate[i] - Float_Schedule->Days_ForwardStart[0])) / denominator;
		P1 = exp(-Float_RefCurve.Interpolated_Rate(T1) * T1);
		P_Pay = exp(-Float_DiscCurve.Interpolated_Rate(TPay) * TPay);
		
		if (CalcCRSFlag == 0)
		{
			NA = 1.0;
			FX = 1.0;
		}
		else
		{
			NA = Float_Schedule->NotionalAmount;
			FX = FX_FloCurve.Interpolated_Rate(T1);
		}		

		FloatValue += NA * FX * (P0 / P1 - 1.0) * P_Pay;
		P0 = P1;
	}
	
	if (Fix_Schedule->DayCount == 0) denominator = 365.0;
	else denominator = 360.0;

	FixValue = 0.0;
	for (i = 0; i < Fix_Schedule->NCF; i++)
	{
		dt_Fix = ((double)(Fix_Schedule->Days_EndDate[i] - Fix_Schedule->Days_StartDate[i])) / denominator;
		TPay = ((double)(Fix_Schedule->Days_PayDate[i] - Fix_Schedule->Days_StartDate[0])) / denominator;
		P_Pay = exp(-Fix_DiscCurve.Interpolated_Rate(TPay) * TPay);

		T1 = ((double)(Fix_Schedule->Days_EndDate[i] - Fix_Schedule->Days_StartDate[0])) / denominator;

		if (CalcCRSFlag == 0)
		{
			NA = 1.0;
			FX = 1.0;
		}
		else
		{
			NA = Fix_Schedule->NotionalAmount;
			FX = FX_FixCurve.Interpolated_Rate(T1);
		}

		FixValue += NA* FX * P_Pay * dt_Fix;
	}

	IRSpread = FloatValue / FixValue;
	return IRSpread;
}

double SOFR_ForwardRate_Compound(
	curveinfo& RefCurve,
	long Ref_Days_0,
	long Ref_Days_1,
	long FixedRateFlag,
	double FixedRate,
	long HolidayCalcType,
	long NHolidays,
	long* Days_Holidays,
	long LockOut,
	long LookBack,
	long LookBackObservationShift,
	double& AnnualizedOISRate
)
{
	long i;
	long j;
	long k;
	long n;



	long HolidayFlag = 0;
	long HolidayFlag2 = 0;
	long nHoliday = 0;
	long TimePos = 0;
	double dt = 1.0 / 365.0;
	double Prev_PI;
	double TodayRate;
	double ForwardRate;
	double PrevRate, PrevRate2;
	double t = 0.0;
	double PI_0;
	double T;

	long ObservShiftFlag = 0;
	if (LookBack > 0 && LookBackObservationShift > 0) 	ObservShiftFlag = 1;

	if (LockOut < 0) LockOut = 0;
	
	///////////////////////////
	// N������ �� LockOutDay ���
	///////////////////////////

	long LockOutDay = Ref_Days_1;
	if (LockOut > 0)
	{
		k = 0;
		for (i = 1; i < 30; i++)
		{
			LockOutDay -= 1;
			// N�����ϱ��� ��¥ �ڷΰ���
			if (isin(LockOutDay, Days_Holidays, NHolidays) == 0) k += 1;
			if (k == LockOut) break;
		}
	}

	double LockOutDayRate = 0.0;

	///////////////////////////
	// Holiday Rate�� Interpolate�� �� ����� ����
	///////////////////////////
	double TargetRate[2] = { 0.0,0.0 };
	double TargetT[2] = { 0.0,0.0 };

	long lookbackday;
	if (Ref_Days_0 < 0 && ObservShiftFlag == 0)
	{
		///////////////////////////
		// Forward Start ��¥�� ���� ����Ϻ��� �տ� ���� ���
		///////////////////////////

		if (FixedRateFlag == 1)
		{
			if (FixedRate > 0.0 || FixedRate < -0.00002) 
			{
				// �Էµ� FixedRate�� �����ϴ°��
				Prev_PI = (1.0 + (-(double)Ref_Days_0)/365.0 * FixedRate);
			}
			else
			{
				// �Էµ� FixedRate�� �������� ������ ���� OverNight�ݸ� ���
				Prev_PI = (1.0 + (-(double)Ref_Days_0) / 365.0 * Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(),0.0, &TimePos));
			}
		}
		else 
		{
			// Flag�� ��� ���Ÿ�ĥ�� �ݸ� ��� ���� ������ �ݸ� ���
			Prev_PI = (1.0 + (-(double)Ref_Days_0) / 365.0 * Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), 0.0, &TimePos));
		}
		// ������� ���űݸ� Compound
		TodayRate = Calc_Forward_Rate(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), 0.0, 1.0 / 365.0);
		PrevRate = TodayRate;

		// ���ñ��� �ݸ� Compound
		PI_0 = Prev_PI * (1.0 + TodayRate * dt);

		// ���Ϻ��� �ݸ� Compound
		if (HolidayCalcType == 0) 
		{
			//////////////////
			// HolidayRate Forward Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0; // ù ��¥�� ������ �����Ϸ� ó��
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);

					if (LockOutDay == i) LockOutDayRate = ForwardRate;
				}
			}
		}
		else if (HolidayCalcType == 1)
		{
			//////////////////
			// HolidayRate Backward Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}


					if (nHoliday == 0) ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					else
					{
						// ������ ���� ù �������� ��¥
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 120; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}
						ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(),t, &TimePos);
					}

					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);

					if (LockOutDay == i) LockOutDayRate = ForwardRate;
				}
			}
		}
		else if (HolidayCalcType == 2)
		{
			//////////////////
			// HolidayRate Interpolated Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == i) LockOutDayRate = ForwardRate;

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}

					if (nHoliday == 0)
					{
						PI_0 *= (1.0 + ForwardRate * dt);
						if (LockOutDay == i) LockOutDayRate = ForwardRate;
					}
					else
					{
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 120; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}

						PI_0 *= (1.0 + ForwardRate * dt);
						TargetT[0] = 0.0;
						TargetT[1] = (double)nHoliday;
						TargetRate[0] = ForwardRate;
						PrevRate2 = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
						TargetRate[1] = PrevRate2;

						for (j = 0; j < nHoliday; j++)
						{
							if (i + j + 1 > LockOutDay) PI_0 *= (1.0 + ForwardRate * dt);
							else PI_0 *= (1.0 + Interpolate_Linear(TargetT, TargetRate, 2, ((double)(j + 1))) * dt);
						}

					}
				}
			}
		}
		else
		{
			for (i = 1; i < Ref_Days_1; i++)
			{
				ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), ((double)i) / 365.0, &TimePos);
				PI_0 *= (1.0 + ForwardRate * dt);
			}
		}
	}
	else if (Ref_Days_0 >= 0 && ObservShiftFlag == 0)
	{
		PI_0 = 1.0;
		if (HolidayCalcType == 0)
		{
			//////////////////
			// HolidayRate Forward Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				if (i == Ref_Days_0) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						k = 0;
						for (n = 1; n < 120; n++)
						{
							///////////////////////////
							// N������ �� LookBackDay ���
							///////////////////////////

							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == i) { LockOutDayRate = ForwardRate;}

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
				}
			}
		}
		else if (HolidayCalcType == 1)
		{
			//////////////////
			// HolidayRate Backward Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				if (i == Ref_Days_0) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					if (nHoliday == 0) ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					else
					{
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 120; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}
						ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					}

					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == i) { LockOutDayRate = ForwardRate;}

					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
				}
			}
		}
		else if (HolidayCalcType == 2)
		{
			//////////////////
			// HolidayRate Interpolated Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				if (i == Ref_Days_0) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == i)
					{
						LockOutDayRate = ForwardRate;
					}

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}

					if (nHoliday == 0)
					{
						PI_0 *= (1.0 + ForwardRate *  dt);
					}
					else
					{
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 30; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}

						PI_0 *= (1.0 + ForwardRate * dt);
						TargetT[0] = 0.0;
						TargetT[1] = ((double)(nHoliday)) ;
						TargetRate[0] = ForwardRate;
						PrevRate2 = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
						TargetRate[1] = PrevRate2;
						for (j = 0; j < nHoliday; j++)
						{
							if (i + j + 1 > LockOutDay) PI_0 *= (1.0 + ForwardRate * dt);
							else PI_0 *= (1.0 + Interpolate_Linear(TargetT, TargetRate, 2, ((double)(j + 1))) * dt);
						}
					}
				}
			}
		}
		else
		{
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				HolidayFlag = isin(i, Days_Holidays, NHolidays);

				ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), ((double)i) / 365.0, &TimePos);
				PI_0 *= (1.0 + ForwardRate * dt);
			}
		}
	}
	else if (Ref_Days_0 < 0 && ObservShiftFlag != 0)
	{
		if (FixedRateFlag == 1)
		{
			if (FixedRate > 0.0 || FixedRate < -0.00002) // �Էµ� FixedRate�� �����ϴ°��
			{
				Prev_PI = (1.0 + (-(double)Ref_Days_0) / 365.0 * FixedRate);
			}
			else
			{
				Prev_PI = (1.0 + (-(double)Ref_Days_0) / 365.0 * Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), 0.0, &TimePos));
			}
		}
		else // Flag�� ��� ���Ÿ�ĥ�� �ݸ� ��� ���� ������ �ݸ� ���
		{
			Prev_PI = (1.0 + (-(double)Ref_Days_0) / 365.0 * Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), 0.0, &TimePos));
		}
		// ������� ���űݸ� Compound
		TodayRate = Calc_Forward_Rate(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), 0.0, 1.0 / 365.0);
		PrevRate = TodayRate;

		// ���ñݸ� Compound
		PI_0 = Prev_PI * (1.0 + TodayRate * dt);
		// ���Ϻ��� �ݸ� Compound
		if (HolidayCalcType == 0) 
		{
			//////////////////
			// HolidayRate Forward Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);

					if (LockOutDay == i) LockOutDayRate = ForwardRate;
				}
			}
		}
		else if (HolidayCalcType == 1)
		{
			//////////////////
			// HolidayRate Backward Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}


					if (nHoliday == 0) ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					else
					{
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 120; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}
						ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					}

					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);

					if (LockOutDay == i) LockOutDayRate = ForwardRate;
				}
			}
		}
		else if (HolidayCalcType == 2)
		{
			//////////////////
			// HolidayRate Interpolated Fill
			//////////////////
			for (i = 1; i < Ref_Days_1; i++)
			{
				if (i == 1) HolidayFlag = 0;
				else HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (i > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == i) LockOutDayRate = ForwardRate;

					nHoliday = 0;
					for (j = i + 1; j < Ref_Days_1; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}

					if (nHoliday == 0)
					{
						PI_0 *= (1.0 + ForwardRate * dt);
						if (LockOutDay == i) LockOutDayRate = ForwardRate;
					}
					else
					{
						lookbackday = i + nHoliday + 1;
						if (LookBack < 1) t = ((double)lookbackday) / 365.0;
						else
						{
							k = 0;
							for (n = 1; n < 120; n++)
							{
								lookbackday -= 1;
								// N�����ϱ��� ����
								if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
								if (k == LookBack) break;
							}
							t = ((double)lookbackday) / 365.0;
						}

						PI_0 *= (1.0 + ForwardRate * dt);
						TargetT[0] = 0.0;
						TargetT[1] = (double)nHoliday;
						TargetRate[0] = ForwardRate;
						PrevRate2 = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
						TargetRate[1] = PrevRate2;

						for (j = 0; j < nHoliday; j++)
						{
							if (i + j + 1 > LockOutDay) PI_0 *= (1.0 + ForwardRate * dt);
							else PI_0 *= (1.0 + Interpolate_Linear(TargetT, TargetRate, 2, ((double)(j + 1))) * dt);
						}

					}
				}
			}
		}
		else
		{
			for (i = 1; i < Ref_Days_1; i++)
			{
				ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), ((double)i) / 365.0, &TimePos);
				PI_0 *= (1.0 + ForwardRate * dt);
			}
		}
	}
	else
	{
		long LastResetDate = Ref_Days_1+0;
		k = 0;
		for (n = 1; n < 120; n++)
		{
			LastResetDate -= 1;
			if (isin(LastResetDate, Days_Holidays, NHolidays) == 0) k += 1;
			if (k == LookBack) break;
		}

		PI_0 = 1.0;
		if (HolidayCalcType == 0)
		{
			//////////////////
			// HolidayRate Forward Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (lookbackday > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == lookbackday)
					{
						LockOutDayRate = ForwardRate;
					}

					nHoliday = 0;
					for (j = lookbackday + 1; j < LastResetDate; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
				}
			}
		}
		else if (HolidayCalcType == 1)
		{
			//////////////////
			// HolidayRate Backward Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					nHoliday = 0;
					for (j = lookbackday + 1; j < LastResetDate; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}
					if (nHoliday == 0) ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					else
					{
						t = ((double)(lookbackday + nHoliday + 1))/365.0;
						ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					}

					if (lookbackday > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == lookbackday) LockOutDayRate = ForwardRate;

					PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
				}
			}
		}
		else if (HolidayCalcType == 2)
		{
			//////////////////
			// HolidayRate Interpolated Fill
			//////////////////
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				HolidayFlag = isin(i, Days_Holidays, NHolidays);

				if (HolidayFlag == 0)
				{
					lookbackday = i + 0;

					if (LookBack < 1) t = ((double)i) / 365.0;
					else
					{
						///////////////////////////
						// N������ �� LookBackDay ���
						///////////////////////////

						k = 0;
						for (n = 1; n < 120; n++)
						{
							lookbackday -= 1;
							// N�����ϱ��� ����
							if (isin(lookbackday, Days_Holidays, NHolidays) == 0) k += 1;
							if (k == LookBack) break;
						}
						t = ((double)lookbackday) / 365.0;
					}

					ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);
					if (lookbackday > LockOutDay && LockOutDayRate != 0.0) ForwardRate = LockOutDayRate;
					if (LockOutDay == lookbackday)
					{
						LockOutDayRate = ForwardRate;
					}

					nHoliday = 0;
					for (j = lookbackday + 1; j < LastResetDate; j++)
					{
						HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
						if (HolidayFlag2 == 0) break;
						else nHoliday += 1;
					}

					if (nHoliday == 0)
					{
						PI_0 *= (1.0 + ForwardRate * dt);
					}
					else
					{
						t = ((double)(lookbackday + nHoliday + 1)) / 365.0;
						PrevRate2 = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), t, &TimePos);

						PI_0 *= (1.0 + ForwardRate * dt);
						TargetT[0] = 0.0;
						TargetT[1] = ((double)(nHoliday));
						TargetRate[0] = ForwardRate;
						TargetRate[1] = PrevRate2;
						for (j = 0; j < nHoliday; j++)
						{
							if (i + j + 1 > LockOutDay) PI_0 *= (1.0 + ForwardRate * dt);
							else PI_0 *= (1.0 + Interpolate_Linear(TargetT, TargetRate, 2, ((double)(j + 1))) * dt);
						}
					}
				}
			}
		}
		else
		{
			for (i = Ref_Days_0; i < Ref_Days_1; i++)
			{
				HolidayFlag = isin(i, Days_Holidays, NHolidays);

				ForwardRate = Calc_Forward_Rate_Daily(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), ((double)i) / 365.0, &TimePos);
				PI_0 *= (1.0 + ForwardRate * dt);
			}
		}
	}


	T = (double)(Ref_Days_1 - Ref_Days_0) * dt;
	AnnualizedOISRate = (PI_0 -1.0) * 1.0/T;
	return PI_0 - 1.0;
}

double Calc_Current_SOFR_Swap(
	long CalcCRSFlag,
	curveinfo& Float_DiscCurve,
	curveinfo& Float_RefCurve,
	curveinfo& Fix_DiscCurve,
	curveinfo& Float_FXCurve,
	curveinfo& Fix_FXCurve,
	SCHD* Float_Schedule,
	SCHD* Fix_Schedule
)
{
	long i;
	long Day0, Day1;

	double SOFR_Swap_Spread;
	double denominator;

	if (Float_Schedule->DayCount == 0)denominator = 365.0;
	else denominator = 360.0;

	double P_Pay;
	double TPay;
	double FloatValue, FixValue;
	double dt_Fix;
	double SOFR_Compound;
	double Annual_OIS = 0.0;
	double FX = 1.0;
	double NA = 1.0;

	Day0 = 0;
	FloatValue = 0.0;
	for (i = 0; i < Float_Schedule->NCF; i++)
	{
		Day1 = Float_Schedule->Days_ForwardEnd[i] - Float_Schedule->Days_ForwardStart[0];
		TPay = ((double)(Float_Schedule->Days_PayDate[i] - Float_Schedule->Days_StartDate[0])) / denominator;
		SOFR_Compound = SOFR_ForwardRate_Compound(Float_RefCurve, Day0, Day1, 0, 0.0, 
			Float_Schedule->HolidayFlag_Ref, Float_Schedule->NHolidays_Ref, Float_Schedule->Days_Holidays_Ref, Float_Schedule->LockOutRef,Float_Schedule->LookBackRef, Float_Schedule->ObservationShift, Annual_OIS);
		P_Pay = exp(-Float_DiscCurve.Interpolated_Rate(TPay) * TPay) ;
		
		if (CalcCRSFlag == 0)
		{
			NA = 1.0;
			FX = 1.0;
		}
		else
		{
			NA = Float_Schedule->NotionalAmount;
			FX = Float_FXCurve.Interpolated_Rate(((double)(Float_Schedule->Days_ForwardEnd[i] - Float_Schedule->Days_StartDate[0])) / denominator);
		}
		
		FloatValue += NA * FX * SOFR_Compound * P_Pay;
		Day0 = Day1;
	}

	if (Fix_Schedule->DayCount == 0) denominator = 365.0;
	else denominator = 360.0;

	FixValue = 0.0;
	for (i = 0; i < Fix_Schedule->NCF; i++)
	{
		dt_Fix = ((double)(Fix_Schedule->Days_EndDate[i] - Fix_Schedule->Days_StartDate[i])) / denominator;
		TPay = ((double)(Fix_Schedule->Days_PayDate[i] - Fix_Schedule->Days_StartDate[0])) / denominator;
		P_Pay = exp(-Fix_DiscCurve.Interpolated_Rate(TPay) * TPay);

		if (CalcCRSFlag == 0)
		{
			NA = 1.0;
			FX = 1.0;
		}
		else
		{
			NA = Fix_Schedule->NotionalAmount;
			FX = Fix_FXCurve.Interpolated_Rate(((double)(Fix_Schedule->Days_EndDate[i] - Fix_Schedule->Days_StartDate[0])) / denominator);
		}
		FixValue += NA*FX * P_Pay * dt_Fix;
	}

	SOFR_Swap_Spread = FloatValue / FixValue;
	return SOFR_Swap_Spread;
}

double Calc_Forward_SOFR_Swap(
	curveinfo& Float_DiscCurve,
	curveinfo& Float_RefCurve,
	curveinfo& Fix_DiscCurve,
	double T_Option,
	double Tenor,
	double FreqMonth,
	long HolidayCalcFlag,
	long NHolidays,
	long* Days_Holidays,
	long Float_RefLockOut,
	long Float_RefLookBack,
	long ObservationShift
)
{
	long i;
	long k;
	
	long Day0, Day1;
	double dt;

	long ndates = NPaymentSwap(T_Option, T_Option + Tenor, FreqMonth);
	long* dates = (long*)malloc(sizeof(long) * (ndates));
	ndates = MappingPaymentDates(T_Option + Tenor, FreqMonth, dates, ndates);
	double* P0T = (double*)malloc(sizeof(double) * (ndates + 1));
	P0T[0] = Calc_Discount_Factor(Fix_DiscCurve.Term, Fix_DiscCurve.Rate, Fix_DiscCurve.nterm(), T_Option); // �ɼǸ������
	for (i = 1; i < ndates + 1; i++) P0T[i] = Calc_Discount_Factor(Fix_DiscCurve.Term, Fix_DiscCurve.Rate, Fix_DiscCurve.nterm(), ((double)dates[i - 1]) / 365.0);

	double* P0T_Float = (double*)malloc(sizeof(double) * (ndates + 1));
	P0T_Float[0] = Calc_Discount_Factor(Float_DiscCurve.Term, Float_DiscCurve.Rate, Float_DiscCurve.nterm(), T_Option);
	for (i = 1; i < ndates + 1; i++) P0T_Float[i] = Calc_Discount_Factor(Float_DiscCurve.Term, Float_DiscCurve.Rate, Float_DiscCurve.nterm(), ((double)dates[i - 1]) / 365.0);

	double SOFR_Swap_Spread;

	double P_Pay;

	double FloatValue, FixValue;

	double SOFR_Compound;
	double Annual_OIS = 0.0;	

	Day0 = (long)(T_Option * 365.0 + 0.001);

	long NHoliday_Adjust;
	NHoliday_Adjust = 0;
	for (i = 0; i < NHolidays; i++)
	{
		if (Days_Holidays[i] >= Day0){
			NHoliday_Adjust += 1;
		}
	}
	long* Days_Holidays_Adjust = (long*)malloc(sizeof(long) * max(NHoliday_Adjust, 1));
	k = 0;
	for (i = 0; i < NHolidays; i++)
	{
		if (Days_Holidays[i] >= Day0) {
			Days_Holidays_Adjust[k] = Days_Holidays[i] - Day0;
			k+= 1;
		}
	}

	FloatValue = 0.0;
	for (i = 0; i < ndates; i++)
	{
		Day1 = dates[i];
		SOFR_Compound = SOFR_ForwardRate_Compound(Float_RefCurve, Day0, Day1, 0, 0.0, HolidayCalcFlag, NHoliday_Adjust, Days_Holidays_Adjust, Float_RefLockOut, Float_RefLookBack, ObservationShift, Annual_OIS);
		P_Pay = P0T_Float[i + 1];
		FloatValue += SOFR_Compound * P_Pay;
		Day0 = Day1;
	}

	FixValue = 0.0;
	for (i = 0; i < ndates; i++)
	{
		if (i == 0) dt = (double)dates[0] / 365.0 - T_Option;
		else dt = ((double)(dates[i] - dates[i - 1])) / 365.0;
		FixValue += dt * P0T[i + 1];
	}


	SOFR_Swap_Spread = FloatValue / FixValue;

	free(dates);
	free(P0T_Float);
	free(Days_Holidays_Adjust);
	free(P0T);
	return SOFR_Swap_Spread;
}

void Floating_PartialValue(
	long CalcCRSFlag,
	curveinfo& DiscCurve,
	curveinfo& RefCurve,
	curveinfo& FXCurve,
	double Ref_T0,
	double Ref_T1,
	double Frac_T0,
	double Frac_T1,
	double Pay_T,
	long FixedRateFlag,
	double FixedRate,
	double Slope,
	double FixedAmount,
	double Notional,
	// ���
	double* ResultRefRate,
	double* ResultCF,
	double* ResultDisc,
	double* ResultDiscCF
)
{
	double dt_Forward;
	double dt_CPN;

	double Disc0, Disc1, DiscPay;
	double ForwardRate;
	double CF;
	double FX = 1.0;

	if (CalcCRSFlag == 0) FX = 1.0;
	else FX = FXCurve.Interpolated_Rate(Frac_T1);

	dt_Forward = (Ref_T1 - Ref_T0);
	dt_CPN = (Frac_T1 - Frac_T0);

	if (FixedRateFlag == 0)
	{
		Disc0 = exp(-RefCurve.Interpolated_Rate(Ref_T0) * Ref_T0);
		Disc1 = exp(-RefCurve.Interpolated_Rate(Ref_T1) * Ref_T1);
		if (dt_Forward < 1.0) ForwardRate = 1.0 / dt_Forward * (Disc0 / Disc1 - 1.0);
		else ForwardRate = Continuous_ForwardRate(RefCurve, Ref_T0, Ref_T1);
	}
	else{
		if (FixedRate > 0.00001 || FixedRate < -0.00001)  ForwardRate = FixedRate;
		else
		{
			Disc0 = exp(-RefCurve.Interpolated_Rate(Ref_T0) * Ref_T0);
			Disc1 = exp(-RefCurve.Interpolated_Rate(Ref_T1) * Ref_T1);
			if (dt_Forward < 1.0) ForwardRate = 1.0 / dt_Forward * (Disc0 / Disc1 - 1.0);
			else ForwardRate = Continuous_ForwardRate(RefCurve, Ref_T0, Ref_T1);
		}
	}
	CF = (Slope*ForwardRate + FixedAmount )* dt_CPN  ;
	DiscPay = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);

	*ResultRefRate = ForwardRate;
	*ResultCF = Notional * CF * FX;
	*ResultDisc = DiscPay;
	*ResultDiscCF = Notional * DiscPay * CF;
}

void FixedLeg_PartialValue(
	long CalcCRSFlag,
	curveinfo& DiscCurve,
	curveinfo& FXCurve,
	double Frac_T0,
	double Frac_T1,
	double Pay_T,
	double FixedAmount,
	double Notional,
	// ���
	double* ResultRefRate,
	double* ResultCF,
	double* ResultDisc,
	double* ResultDiscCF
)
{
	double dt_CPN;

	double DiscPay;
	double CF;
	double FX = 1.0;
	
	dt_CPN = (Frac_T1 - Frac_T0);
	CF = (FixedAmount) * dt_CPN;
	DiscPay = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);
	if (CalcCRSFlag == 0) FX = 1.0;
	else FX = FXCurve.Interpolated_Rate(Frac_T1);

	*ResultRefRate = 0.0;
	*ResultCF = FX*Notional * CF;
	*ResultDisc = DiscPay;
	*ResultDiscCF = Notional * DiscPay * CF;
}

double LegValue(
	long CRSFlag,
	SCHD* Schedule,
	curveinfo& DiscCurve,
	curveinfo& RefCurve,
	curveinfo& FXCurve,
	double* ResultRefRate,
	double* ResultCPN,
	double* ResultDF,
	double* DiscCFArray
)
{
	long i;
	long FixedRateFlag ;
	long PrevFlag;
	double denominator;

	if (Schedule->DayCount == 0) denominator = 365.0;
	else denominator = 360.0;

	///////////////
	// SOFR ���� //
	///////////////

	long UseHistorySOFR;
	long* NSaturSunDay_List = Schedule->NWeekendDate;
	long** SaturSunDay_List = Schedule->WeekendList;
	long NRefHistory = Schedule->NRefHistory;
	long* RefHistoryDate = Schedule->RefHistoryDate;
	double* RefHistory = Schedule->RefHistory;

	double Ref_T0, Ref_T1, SwapStartT, SwapEndT;
	double Frac_T0, Frac_T1;
	double Pay_T;
	double value =0.0;
	double FreqMonth;
	double FXRate = 1.0;
	double OIS_Annual = 0.0, OIS_Compound = 0.0;

	if (Schedule->FixedFlotype == 1)
	{
		//////////////
		// �����ݸ� //
		//////////////

		if (Schedule->ReferenceType == 0) 
		{
			///////////////////////////////////////////
			// �Ϲ����� Libor, CD, Zero Rate ����   ///
			///////////////////////////////////////////

			for (i = 0; i < Schedule->NCF; i++)
			{
				// �̹� ���� ������ �������� Ȯ��
				PrevFlag = 1;
				if (Schedule->Days_PayDate[i] > 0)
					PrevFlag = 0;

				if (PrevFlag == 0)
				{
					Ref_T0 = ((double)Schedule->Days_ForwardStart[i]) / denominator;
					Ref_T1 = ((double)Schedule->Days_ForwardEnd[i]) / denominator;
					Frac_T0 = ((double)Schedule->Days_StartDate[i]) / denominator;
					Frac_T1 = ((double)Schedule->Days_EndDate[i]) / denominator;
					Pay_T = ((double)Schedule->Days_PayDate[i]) / denominator;

					if (Schedule->Days_ForwardStart[i] < 0) FixedRateFlag = 1;
					else FixedRateFlag = 0;

					Floating_PartialValue(
						CRSFlag, DiscCurve, RefCurve, FXCurve, Ref_T0, Ref_T1, Frac_T0,
						Frac_T1, Pay_T, FixedRateFlag, Schedule->FixedRefRate[i], Schedule->Slope[i], Schedule->CPN[i], Schedule->NotionalAmount,
						ResultRefRate + i, ResultCPN + i, ResultDF + i, DiscCFArray + i);
				}
				else
				{
					ResultRefRate[i] = Schedule->FixedRefRate[i];
					ResultCPN[i] = 0.0;
					ResultDF[i] = 0.0;
					DiscCFArray[i] = 0.0;
				}
			}
		}
		else if(Schedule->ReferenceType ==1) 
		{
			//////////////////////////////////
			// �����ڻ��� Swap Rate ����   ///
			//////////////////////////////////

			for (i = 0; i < Schedule->NCF; i++)
			{
				// �̹� ���� ������ �������� Ȯ��
				PrevFlag = 1;
				if (Schedule->Days_PayDate[i] > 0) PrevFlag = 0;

				if (PrevFlag == 0)
				{
					SwapStartT = ((double)Schedule->Days_ForwardStart[i]) / denominator;
					SwapEndT = SwapStartT + Schedule->RefSwapMaturity;
					Pay_T = ((double)Schedule->Days_PayDate[i]) / denominator;

					if (Schedule->Days_ForwardStart[i] <= 0) FixedRateFlag = 1;
					else FixedRateFlag = 0;

					if (CRSFlag == 0) FXRate = 1.0;
					else FXRate = FXCurve.Interpolated_Rate(SwapEndT);

					FreqMonth = 12.0 / (double)Schedule->NSwapPayAnnual;
					if (FixedRateFlag == 0)
					{
						ResultRefRate[i] = FSR(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), SwapStartT, Schedule->RefSwapMaturity, FreqMonth);
						ResultCPN[i] = FXRate * Schedule->NotionalAmount * (ResultRefRate[i] * Schedule->Slope[i] + Schedule->CPN[i]) * ((double)(Schedule->Days_EndDate[i] - Schedule->Days_StartDate[i])) / denominator;
						ResultDF[i] = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);
						DiscCFArray[i] = ResultCPN[i] * ResultDF[i];
					}
					else
					{
						//////////////////////////////////////////////////////////////////////////
						// �Է��ص� ���� �ݸ��� �����ϸ� �װ� ����ϰ� �ƴϸ� ���� ���ұݸ���� //
						//////////////////////////////////////////////////////////////////////////

						if (Schedule->FixedRefRate[i] > 0.00001 || Schedule->FixedRefRate[i] < -0.00001) ResultRefRate[i] = Schedule->FixedRefRate[i];
						else ResultRefRate[i] = FSR(RefCurve.Term, RefCurve.Rate, RefCurve.nterm(), SwapStartT, Schedule->RefSwapMaturity, FreqMonth);
						
						ResultCPN[i] = FXRate* Schedule->NotionalAmount * (ResultRefRate[i] * Schedule->Slope[i] + Schedule->CPN[i]) * ((double)(Schedule->Days_EndDate[i] - Schedule->Days_StartDate[i])) / denominator;
						ResultDF[i] = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);
						DiscCFArray[i] = ResultCPN[i] * ResultDF[i];
					}
				}
				else
				{
					ResultRefRate[i] = Schedule->FixedRefRate[i];
					ResultCPN[i] = 0.0;
					ResultDF[i] = 0.0;
					DiscCFArray[i] = 0.0;
				}
			}
		}
		else if (Schedule->ReferenceType == 2) 
		{
			//////////////////////////////////
			// �����ڻ��� SOFR Rate ����   ///
			//////////////////////////////////

			for (i = 0; i < Schedule->NCF; i++)
			{
				// �̹� ���� ������ �������� Ȯ��
				PrevFlag = 1;
				if (Schedule->Days_PayDate[i] > 0) PrevFlag = 0;

				if (PrevFlag == 0)
				{
					SwapStartT = ((double)Schedule->Days_ForwardStart[i]) / denominator;
					SwapEndT = SwapStartT + Schedule->RefSwapMaturity;
					Pay_T = ((double)Schedule->Days_PayDate[i]) / denominator;

					if (Schedule->Days_ForwardStart[i] < 0) FixedRateFlag = 1;
					else FixedRateFlag = 0;

					if (CRSFlag == 0) FXRate = 1.0;
					else FXRate = Interpolate_Linear(FXCurve.Term, FXCurve.Rate, FXCurve.nterm(), SwapEndT);

					FreqMonth = 12.0 / (double)Schedule->NSwapPayAnnual;
					OIS_Compound = SOFR_ForwardRate_Compound(RefCurve, Schedule->Days_ForwardStart[i], Schedule->Days_ForwardEnd[i], FixedRateFlag, Schedule->FixedRefRate[i], Schedule->HolidayFlag_Ref, Schedule->NHolidays_Ref, Schedule->Days_Holidays_Ref,Schedule->LockOutRef, Schedule->LookBackRef, Schedule->ObservationShift,OIS_Annual);
					ResultRefRate[i] = OIS_Annual;
					ResultCPN[i] = FXRate* Schedule->NotionalAmount * (OIS_Compound * Schedule->Slope[i] + Schedule->CPN[i]);
					ResultDF[i] = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);
					DiscCFArray[i] = ResultCPN[i] * ResultDF[i];
				}
				else
				{
					ResultRefRate[i] = Schedule->FixedRefRate[i];
					ResultCPN[i] = 0.0;
					ResultDF[i] = 0.0;
					DiscCFArray[i] = 0.0;
				}
			}
		}	
		else if (Schedule->ReferenceType == 3) 
		{
			///////////////////////////////////////
			// �����ڻ��� SOFR Swap Rate ����   ///
			///////////////////////////////////////

			for (i = 0; i < Schedule->NCF; i++)
			{
				// �̹� ���� ������ �������� Ȯ��
				PrevFlag = 1;
				if (Schedule->Days_PayDate[i] > 0) PrevFlag = 0;

				if (PrevFlag == 0)
				{
					SwapStartT = ((double)Schedule->Days_ForwardStart[i]) / denominator;
					SwapEndT = SwapStartT + Schedule->RefSwapMaturity;
					Pay_T = ((double)Schedule->Days_PayDate[i]) / denominator;

					if (Schedule->Days_ForwardStart[i] <= 0) FixedRateFlag = 1;
					else FixedRateFlag = 0;

					if (CRSFlag == 0) FXRate = 1.0;
					else FXRate = FXCurve.Interpolated_Rate(SwapEndT);

					FreqMonth = 12.0 / (double)Schedule->NSwapPayAnnual;
					ResultRefRate[i] = Calc_Forward_SOFR_Swap(RefCurve, RefCurve, DiscCurve, SwapStartT, Schedule->RefSwapMaturity, FreqMonth, Schedule->HolidayFlag_Ref, Schedule->NHolidays_Ref, Schedule->Days_Holidays_Ref, Schedule->LockOutRef, Schedule->LookBackRef, Schedule->ObservationShift);
					ResultCPN[i] = FXRate* Schedule->NotionalAmount * (ResultRefRate[i] * Schedule->Slope[i] + Schedule->CPN[i]) * ((double)(Schedule->Days_EndDate[i] - Schedule->Days_StartDate[i])) / denominator;
					ResultDF[i] = exp(-DiscCurve.Interpolated_Rate(Pay_T) * Pay_T);
					DiscCFArray[i] = ResultCPN[i] * ResultDF[i];
				}
				else
				{
					ResultRefRate[i] = Schedule->FixedRefRate[i];
					ResultCPN[i] = 0.0;
					ResultDF[i] = 0.0;
					DiscCFArray[i] = 0.0;
				}
			}
		}
	}
	else
	{
		//////////////
		// �����ݸ� //
		//////////////

		for (i = 0; i < Schedule->NCF; i++)
		{
			// �̹� ���� ������ �������� Ȯ��
			PrevFlag = 1;
			if (Schedule->Days_PayDate[i] > 0) PrevFlag = 0;

			if (PrevFlag == 0)
			{
				Frac_T0 = ((double)Schedule->Days_StartDate[i]) / denominator;
				Frac_T1 = ((double)Schedule->Days_EndDate[i]) / denominator;
				Pay_T = ((double)Schedule->Days_PayDate[i]) / denominator;

				FixedLeg_PartialValue(
					CRSFlag, DiscCurve, FXCurve, Frac_T0,Frac_T1,Pay_T,Schedule->CPN[i], Schedule->NotionalAmount,
					ResultRefRate + i, ResultCPN + i, ResultDF + i, DiscCFArray + i);

				ResultRefRate[i] = Schedule->CPN[i];
			}
			else
			{
				ResultRefRate[i] = 0.0;
				ResultCPN[i] = 0.0;
				ResultDF[i] = 0.0;
				DiscCFArray[i] = 0.0;
			}
		}
	}
	value = sumation(DiscCFArray, Schedule->NCF);
	return value;
}
long SwapPricer(
	long PriceDate_C,
	long CalcCRSFlag,

	SCHD* RcvSchedule,
	long RcvDisc_NTerm,
	double* RcvDisc_Term,
	double* RcvDisc_Rate,
	long RcvRef_NTerm,
	double* RcvRef_Term,
	double* RcvRef_Rate,

	SCHD* PaySchedule,
	long PayDisc_NTerm,
	double* PayDisc_Term,
	double* PayDisc_Rate,
	long PayRef_NTerm,
	double* PayRef_Term,
	double* PayRef_Rate,

	long Rcv_NFX,
	double* Rcv_TermFX,
	double* Rcv_FX,

	long Pay_NFX,
	double* Pay_TermFX,
	double* Pay_FX,

	long GreekFlag,
	double* ResultPrice,
	double* ResultRefRate,
	double* ResultCPN,
	double* ResultDF,

	double* PV01,
	double* KeyRateRcvPV01,
	double* KeyRatePayPV01,
	long PricingOnly
)
{
	long i, j;
	long N, M;
	long s;
	long CalcCurrentSwapFlag;

	double RcvValueUp, RcvValueDn, PayValueUp, PayValueDn;

	double Rcv_T_Mat = ((double)DayCountAtoB(PriceDate_C, RcvSchedule->NotionalPayDate_C)) / 365.0;
	double Pay_T_Mat = ((double)DayCountAtoB(PriceDate_C, PaySchedule->NotionalPayDate_C)) / 365.0;

	double Rcv_FXMat;
	if (CalcCRSFlag == 0) Rcv_FXMat = 1.0;
	else Rcv_FXMat = Interpolate_Linear(Rcv_TermFX, Rcv_FX, Rcv_NFX, Rcv_T_Mat);
	
	double Pay_FXMat;
	if (CalcCRSFlag == 0) Pay_FXMat = 1.0;
	else Pay_FXMat = Interpolate_Linear(Pay_TermFX, Pay_FX, Pay_NFX, Pay_T_Mat);
	
	if (RcvSchedule->NAFlag != 0) RcvSchedule->NAFlag = 1;
	if (PaySchedule->NAFlag != 0) PaySchedule->NAFlag = 1;

	curveinfo Rcv_DiscCurve(RcvDisc_NTerm, RcvDisc_Term, RcvDisc_Rate);
	curveinfo Rcv_RefCurve(RcvRef_NTerm, RcvRef_Term, RcvRef_Rate);
	curveinfo Rcv_FXCurve(Rcv_NFX, Rcv_TermFX, Rcv_FX);

	curveinfo Pay_DiscCurve(PayDisc_NTerm, PayDisc_Term, PayDisc_Rate);
	curveinfo Pay_RefCurve(PayRef_NTerm, PayRef_Term, PayRef_Rate);
	curveinfo Pay_FXCurve(Pay_NFX, Pay_TermFX, Pay_FX);

	double RcvValue, PayValue, CurrentIRS;
	double *Rcv_DiscCFArray = (double*)malloc(sizeof(double) * RcvSchedule->NCF);
	double *Pay_DiscCFArray = (double*)malloc(sizeof(double) * PaySchedule->NCF);

	////////////////////
	// �� Leg Pricing //
	////////////////////

	RcvValue = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve, Rcv_RefCurve, Rcv_FXCurve, ResultRefRate, ResultCPN, ResultDF, Rcv_DiscCFArray);
	PayValue = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve, Pay_RefCurve, Pay_FXCurve, ResultRefRate + RcvSchedule->NCF, ResultCPN + RcvSchedule->NCF, ResultDF + RcvSchedule->NCF, Rcv_DiscCFArray);

	//////////////////
	// Notional �ݿ�
	//////////////////

	if (RcvSchedule->NAFlag == 1)
	{
		ResultDF[RcvSchedule->NCF + PaySchedule->NCF] = Calc_Discount_Factor(Rcv_DiscCurve.Term, Rcv_DiscCurve.Rate, Rcv_DiscCurve.nterm(), ((double)(RcvSchedule->Days_Notional)) / 365.0);// ResultDF[RcvSchedule->NCF - 1];
		ResultRefRate[RcvSchedule->NCF + PaySchedule->NCF] = 0.0;
		ResultCPN[RcvSchedule->NCF + PaySchedule->NCF] = RcvSchedule->NotionalAmount* Rcv_FXMat;
		if (CalcCRSFlag == 0) RcvValue += RcvSchedule->NotionalAmount * ResultDF[RcvSchedule->NCF + PaySchedule->NCF];
		else RcvValue += RcvSchedule->NotionalAmount * ResultDF[RcvSchedule->NCF + PaySchedule->NCF] * Rcv_FXMat;
	}

	if (PaySchedule->NAFlag == 1)
	{
		ResultDF[RcvSchedule->NCF + PaySchedule->NCF + 1] = Calc_Discount_Factor(Pay_DiscCurve.Term, Pay_DiscCurve.Rate, Pay_DiscCurve.nterm(), ((double)(PaySchedule->Days_Notional)) / 365.0);//ResultDF[RcvSchedule->NCF + PaySchedule->NCF - 1];
		ResultRefRate[RcvSchedule->NCF + PaySchedule->NCF + 1] = 0.0;
		ResultCPN[RcvSchedule->NCF + PaySchedule->NCF + 1] = PaySchedule->NotionalAmount * Pay_FXMat;
		if (CalcCRSFlag == 0) PayValue += PaySchedule->NotionalAmount * ResultDF[RcvSchedule->NCF + PaySchedule->NCF + 1];
		else PayValue += PaySchedule->NotionalAmount * ResultDF[RcvSchedule->NCF + PaySchedule->NCF + 1] * Pay_FXMat;
	}

	CalcCurrentSwapFlag = RcvSchedule->FixedFlotype + PaySchedule->FixedFlotype;
	if (CalcCurrentSwapFlag == 1 && PricingOnly == 1)
	{
		if (RcvSchedule->FixedFlotype == 0)
		{
			if (RcvSchedule->ReferenceType < 2 && PaySchedule->ReferenceType < 2) CurrentIRS = Calc_Current_IRS(CalcCRSFlag, Pay_DiscCurve, Pay_RefCurve, Rcv_DiscCurve, Pay_FXCurve, Rcv_FXCurve, PaySchedule, RcvSchedule);
			else CurrentIRS = Calc_Current_SOFR_Swap(CalcCRSFlag, Pay_DiscCurve, Pay_RefCurve, Rcv_DiscCurve, Pay_FXCurve, Rcv_FXCurve, PaySchedule, RcvSchedule);
		}
		else
		{
			if (RcvSchedule->ReferenceType < 2 && PaySchedule->ReferenceType < 2) CurrentIRS = Calc_Current_IRS(CalcCRSFlag, Rcv_DiscCurve, Rcv_RefCurve, Pay_DiscCurve, Rcv_FXCurve, Pay_FXCurve, RcvSchedule, PaySchedule);
			else CurrentIRS = Calc_Current_SOFR_Swap(CalcCRSFlag, Rcv_DiscCurve, Rcv_RefCurve, Pay_DiscCurve, Rcv_FXCurve, Pay_FXCurve, RcvSchedule, PaySchedule);
		}
	}
	else
	{
		CurrentIRS = 0.0;
	}

	
	if (GreekFlag > 0)
	{
		double T;

		double* Rcv_Temp_RefRate = (double*)malloc(sizeof(double) * RcvSchedule->NCF);
		double* Rcv_Temp_ResultCPN = (double*)malloc(sizeof(double) * RcvSchedule->NCF);
		double* Rcv_Temp_ResultDF = (double*)malloc(sizeof(double) * RcvSchedule->NCF);
		double* Rcv_Temp_DiscCFArray = (double*)malloc(sizeof(double) * RcvSchedule->NCF);

		double* Pay_Temp_RefRate = (double*)malloc(sizeof(double) * PaySchedule->NCF);
		double* Pay_Temp_ResultCPN = (double*)malloc(sizeof(double) * PaySchedule->NCF);
		double* Pay_Temp_ResultDF = (double*)malloc(sizeof(double) * PaySchedule->NCF);
		double* Pay_Temp_DiscCFArray = (double*)malloc(sizeof(double) * PaySchedule->NCF);

		// Receive Leg
		double* RcvDisc_Rate_Up = (double*)malloc(sizeof(double) * RcvDisc_NTerm);
		double* RcvDisc_Rate_Dn = (double*)malloc(sizeof(double) * RcvDisc_NTerm);
		for (i = 0; i < RcvDisc_NTerm; i++)
		{
			RcvDisc_Rate_Up[i] = RcvDisc_Rate[i] + 0.0001;
			RcvDisc_Rate_Dn[i] = RcvDisc_Rate[i] - 0.0001;
		}

		double* RcvRef_Rate_Up = (double*)malloc(sizeof(double) * RcvRef_NTerm);
		double* RcvRef_Rate_Dn = (double*)malloc(sizeof(double) * RcvRef_NTerm);
		for (i = 0; i < RcvRef_NTerm; i++)
		{
			RcvRef_Rate_Up[i] = RcvRef_Rate[i] + 0.0001;
			RcvRef_Rate_Dn[i] = RcvRef_Rate[i] - 0.0001;
		}


		// Pay Leg
		double* PayDisc_Rate_Up = (double*)malloc(sizeof(double) * PayDisc_NTerm);
		double* PayDisc_Rate_Dn = (double*)malloc(sizeof(double) * PayDisc_NTerm);
		for (i = 0; i < PayDisc_NTerm; i++)
		{
			PayDisc_Rate_Up[i] = PayDisc_Rate[i] + 0.0001;
			PayDisc_Rate_Dn[i] = PayDisc_Rate[i] - 0.0001;
		}

		double* PayRef_Rate_Up = (double*)malloc(sizeof(double) * PayRef_NTerm);
		double* PayRef_Rate_Dn = (double*)malloc(sizeof(double) * PayRef_NTerm);
		for (i = 0; i < PayRef_NTerm; i++)
		{
			PayRef_Rate_Up[i] = PayRef_Rate[i] + 0.0001;
			PayRef_Rate_Dn[i] = PayRef_Rate[i] - 0.0001;
		}

		curveinfo Rcv_DiscCurve_Up(RcvDisc_NTerm, RcvDisc_Term, RcvDisc_Rate_Up);
		curveinfo Rcv_DiscCurve_Dn(RcvDisc_NTerm, RcvDisc_Term, RcvDisc_Rate_Up);

		curveinfo Rcv_RefCurve_Up(RcvRef_NTerm, RcvRef_Term, RcvRef_Rate_Up);
		curveinfo Rcv_RefCurve_Dn(RcvRef_NTerm, RcvRef_Term, RcvRef_Rate_Dn);

		curveinfo Pay_DiscCurve_Up(PayDisc_NTerm, PayDisc_Term, PayDisc_Rate_Up);
		curveinfo Pay_DiscCurve_Dn(PayDisc_NTerm, PayDisc_Term, PayDisc_Rate_Dn);

		curveinfo Pay_RefCurve_Up(PayRef_NTerm, PayRef_Term, PayRef_Rate_Up);
		curveinfo Pay_RefCurve_Dn(PayRef_NTerm, PayRef_Term, PayRef_Rate_Dn);

		// RCV Leg

		RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Up, Rcv_RefCurve, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
		RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
		RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Dn, Rcv_RefCurve, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
		RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
		PV01[0] = (RcvValueUp - RcvValueDn) * 0.5;

		RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve, Rcv_RefCurve_Up, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
		RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
		RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve, Rcv_RefCurve_Dn, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
		RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
		PV01[1] = (RcvValueUp - RcvValueDn) * 0.5;

		N = Rcv_DiscCurve.nterm();
		M = Rcv_RefCurve.nterm();

		if (N == M)
		{
			s = 0;
			for (i = 0; i < N; i++)
			{
				if (Rcv_DiscCurve.Term[i] == Rcv_RefCurve.Term[i])
				{
					s += 1;
				}
			}

			if (s == N)
			{
				RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Up, Rcv_RefCurve_Up, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Dn, Rcv_RefCurve_Dn, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				PV01[2] = (RcvValueUp - RcvValueDn) * 0.5;
			}
		}

		for (i = 0; i < N; i++)
		{
			T = ((double)(RcvSchedule->Days_PayDate[RcvSchedule->NCF - 1])) / 365.0;
			if (RcvDisc_Term[i] <= T)
			{
				for (j = 0; j < N; j++)
				{
					if (i == j)
					{
						Rcv_DiscCurve_Up.Rate[j] = RcvDisc_Rate[j] + 0.0001;
						Rcv_DiscCurve_Dn.Rate[j] = RcvDisc_Rate[j] - 0.0001;
					}
					else
					{
						Rcv_DiscCurve_Up.Rate[j] = RcvDisc_Rate[j];
						Rcv_DiscCurve_Dn.Rate[j] = RcvDisc_Rate[j];
					}
				}
				RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Up, Rcv_RefCurve, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Dn, Rcv_RefCurve, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				KeyRateRcvPV01[i] = (RcvValueUp - RcvValueDn) * 0.5;
			}
			else
			{
				KeyRateRcvPV01[i] = 0.0;
			}
		}

		for (i = 0; i < M; i++)
		{
			T = ((double)(RcvSchedule->Days_PayDate[RcvSchedule->NCF - 1])) / 365.0;
			if (RcvRef_Term[i] <= T)
			{
				for (j = 0; j < M; j++)
				{
					if (i == j)
					{
						Rcv_RefCurve_Up.Rate[j] = RcvRef_Rate[j] + 0.0001;
						Rcv_RefCurve_Dn.Rate[j] = RcvRef_Rate[j] - 0.0001;
					}
					else
					{
						Rcv_RefCurve_Up.Rate[j] = RcvRef_Rate[j];
						Rcv_RefCurve_Dn.Rate[j] = RcvRef_Rate[j];
					}
				}
				RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve, Rcv_RefCurve_Up, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve, Rcv_RefCurve_Dn, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
				RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
				KeyRateRcvPV01[i + N] = (RcvValueUp - RcvValueDn) * 0.5;
			}
			else
			{
				KeyRateRcvPV01[i + N] = 0.0;
			}
		}

		if (N == M)
		{
			s = 0;
			for (i = 0; i < N; i++)
			{
				if (Rcv_DiscCurve.Term[i] == Rcv_RefCurve.Term[i])
				{
					s += 1;
				}
			}

			T = ((double)(RcvSchedule->Days_PayDate[RcvSchedule->NCF - 1])) / 365.0;

			if (s == N)
			{
				for (i = 0; i < N; i++)
				{
					if (RcvRef_Term[i] <= T)
					{
						for (j = 0; j < N; j++)
						{
							if (i == j)
							{
								Rcv_DiscCurve_Up.Rate[j] = RcvDisc_Rate[j] + 0.0001;
								Rcv_DiscCurve_Dn.Rate[j] = RcvDisc_Rate[j] - 0.0001;
								Rcv_RefCurve_Up.Rate[j] = RcvRef_Rate[j] + 0.0001;
								Rcv_RefCurve_Dn.Rate[j] = RcvRef_Rate[j] - 0.0001;
							}
							else
							{
								Rcv_DiscCurve_Up.Rate[j] = RcvDisc_Rate[j];
								Rcv_DiscCurve_Dn.Rate[j] = RcvDisc_Rate[j];
								Rcv_RefCurve_Up.Rate[j] = RcvRef_Rate[j];
								Rcv_RefCurve_Dn.Rate[j] = RcvRef_Rate[j];
							}
						}
						RcvValueUp = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Up, Rcv_RefCurve_Up, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
						RcvValueUp += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
						RcvValueDn = LegValue(CalcCRSFlag, RcvSchedule, Rcv_DiscCurve_Dn, Rcv_RefCurve_Dn, Rcv_FXCurve, Rcv_Temp_RefRate, Rcv_Temp_ResultCPN, Rcv_Temp_ResultDF, Rcv_Temp_DiscCFArray);
						RcvValueDn += (RcvSchedule->NotionalAmount * (double)RcvSchedule->NAFlag) * Rcv_Temp_ResultDF[RcvSchedule->NCF - 1] * Rcv_FXMat;
						KeyRateRcvPV01[i + N + M] = (RcvValueUp - RcvValueDn) * 0.5;
					}
					else
					{
						KeyRateRcvPV01[i + N + M] = 0.0;
					}
				}
			}
		}

		// PAY Leg
		PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Up, Pay_RefCurve, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
		PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
		PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Dn, Pay_RefCurve, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
		PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
		PV01[3] = (PayValueUp - PayValueDn) * 0.5;

		PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve, Pay_RefCurve_Up, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
		PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
		PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve, Pay_RefCurve_Dn, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
		PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
		PV01[4] = (PayValueUp - PayValueDn) * 0.5;

		N = Pay_DiscCurve.nterm();
		M = Pay_RefCurve.nterm();

		if (N == M)
		{
			s = 0;
			for (i = 0; i < N; i++)
			{
				if (Pay_DiscCurve.Term[i] == Pay_RefCurve.Term[i])
				{
					s += 1;
				}
			}

			if (s == N)
			{
				PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Up, Pay_RefCurve_Up, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Dn, Pay_RefCurve_Dn, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				PV01[5] =  (PayValueUp - PayValueDn) * 0.5;
			}
		}

		for (i = 0; i < N; i++)
		{
			T = ((double)(PaySchedule->Days_PayDate[PaySchedule->NCF - 1])) / 365.0;
			if (PayDisc_Term[i] <= T)
			{
				for (j = 0; j < N; j++)
				{
					if (i == j)
					{
						Pay_DiscCurve_Up.Rate[j] = PayDisc_Rate[j] + 0.0001;
						Pay_DiscCurve_Dn.Rate[j] = PayDisc_Rate[j] - 0.0001;
					}
					else
					{
						Pay_DiscCurve_Up.Rate[j] = PayDisc_Rate[j];
						Pay_DiscCurve_Dn.Rate[j] = PayDisc_Rate[j];
					}
				}
				PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Up, Pay_RefCurve, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Dn, Pay_RefCurve, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				KeyRatePayPV01[i] = (PayValueUp - PayValueDn) * 0.5;
			}
			else
			{
				KeyRatePayPV01[i] = 0.0;
			}

		}

		for (i = 0; i < M; i++)
		{
			T = ((double)(PaySchedule->Days_PayDate[PaySchedule->NCF - 1])) / 365.0;
			if (PayRef_Term[i] <= T)
			{
				for (j = 0; j < M; j++)
				{
					if (i == j)
					{
						Pay_RefCurve_Up.Rate[j] = PayRef_Rate[j] + 0.0001;
						Pay_RefCurve_Dn.Rate[j] = PayRef_Rate[j] - 0.0001;
					}
					else
					{
						Pay_RefCurve_Up.Rate[j] = PayRef_Rate[j];
						Pay_RefCurve_Dn.Rate[j] = PayRef_Rate[j];
					}
				}
				PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve, Pay_RefCurve_Up, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve, Pay_RefCurve_Dn, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
				PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
				KeyRatePayPV01[i + N] = (PayValueUp - PayValueDn) * 0.5;
			}
			else
			{
				KeyRatePayPV01[i + N] = 0.0;
			}
		}

		if (N == M)
		{
			s = 0;
			for (i = 0; i < N; i++)
			{
				if (Pay_DiscCurve.Term[i] == Pay_RefCurve.Term[i])
				{
					s += 1;
				}
			}


			T = ((double)(PaySchedule->Days_PayDate[PaySchedule->NCF - 1])) / 365.0;
			if (s == N)
			{
				for (i = 0; i < N; i++)
				{
					if (PayRef_Term[i] <= T)
					{
						for (j = 0; j < N; j++)
						{
							if (i == j)
							{
								Pay_DiscCurve_Up.Rate[j] = PayDisc_Rate[j] + 0.0001;
								Pay_DiscCurve_Dn.Rate[j] = PayDisc_Rate[j] - 0.0001;
								Pay_RefCurve_Up.Rate[j] = PayRef_Rate[j] + 0.0001;
								Pay_RefCurve_Dn.Rate[j] = PayRef_Rate[j] - 0.0001;
							}
							else
							{
								Pay_DiscCurve_Up.Rate[j] = PayDisc_Rate[j];
								Pay_DiscCurve_Dn.Rate[j] = PayDisc_Rate[j];
								Pay_RefCurve_Up.Rate[j] = PayRef_Rate[j];
								Pay_RefCurve_Dn.Rate[j] = PayRef_Rate[j];
							}
						}
						PayValueUp = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Up, Pay_RefCurve_Up, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
						PayValueUp += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
						PayValueDn = LegValue(CalcCRSFlag, PaySchedule, Pay_DiscCurve_Dn, Pay_RefCurve_Dn, Pay_FXCurve, Pay_Temp_RefRate, Pay_Temp_ResultCPN, Pay_Temp_ResultDF, Pay_Temp_DiscCFArray);
						PayValueDn += (PaySchedule->NotionalAmount * (double)PaySchedule->NAFlag) * Pay_Temp_ResultDF[PaySchedule->NCF - 1] * Pay_FXMat;
						KeyRatePayPV01[i + N + M] =  (PayValueUp - PayValueDn) * 0.5;
					}
					else
					{
						KeyRatePayPV01[i + N + M] = 0.0;
					}
				}
			}
		}

		free(Rcv_Temp_RefRate);
		free(Rcv_Temp_ResultCPN);
		free(Rcv_Temp_ResultDF);
		free(Rcv_Temp_DiscCFArray);

		free(Pay_Temp_RefRate);
		free(Pay_Temp_ResultCPN);
		free(Pay_Temp_ResultDF);
		free(Pay_Temp_DiscCFArray);


		free(RcvDisc_Rate_Up);
		free(RcvDisc_Rate_Dn);

		free(RcvRef_Rate_Up);
		free(RcvRef_Rate_Dn);

		free(PayDisc_Rate_Up);
		free(PayDisc_Rate_Dn);

		free(PayRef_Rate_Up);
		free(PayRef_Rate_Dn);
	}

	ResultPrice[0] = CurrentIRS;
	ResultPrice[1] = RcvValue;
	ResultPrice[2] = PayValue;

	free(Rcv_DiscCFArray);
	free(Pay_DiscCFArray);

	return 1;
}

long ErrorCheckIRSwap_Excel(
	long PriceDate_Exl,				// PricingDate ExcelType
	long GreekFlag,					// Greek���� Flag
	long NAFlag,					// Notional Amount ��� Flag
	long* CRS_Flag,					// [0]CRS Pricing Flag [1]FX���� Term ����
	double* CRS_Info,				// [0~FX����-1] FX Forward Term, [FX����-1~2*FX����-1] FX Forward

	long Rcv_RefRateType,			// Rcv ���ʱݸ� 0: Libor/CD 1: Swap 2: SOFR 3:SOFR Swap
	long Rcv_SwapYearlyNPayment,	// Rcv_RefRateType�� 1, 3�� �� ���� �� ����ȸ��
	double Rcv_SwapMaturity,		// Rcv_RefRateType�� 1, 3�� �� ���Ҹ���
	long Rcv_FixFloFlag,			// Rcv Fix/Flo Flag
	long Rcv_DayCount,				// Rcv DayCountConvention 0:Act365  1: Act360

	double Rcv_NotionalAMT,			// Rcv Leg Notional Amount
	long Rcv_NotionalPayDate,		// Rcv Leg Notional Payment Date
	long RcvDisc_NTerm,				// Rcv Leg ���� �ݸ� Term ����
	double* RcvDisc_Term,			// Rcv Leg ���� �ݸ� Term Array
	double* RcvDisc_Rate,			// Rcv Leg ���� �ݸ� Rate Array 

	long RcvRef_NTerm,				// Rcv Leg ���۷��� �ݸ� Term ����
	double* RcvRef_Term,			// Rcv Leg ���۷��� �ݸ� Term Array
	double* RcvRef_Rate,			// Rcv Leg ���۷��� �ݸ� Rate Array
	long NRcvCF,					// Rcv Leg CashFlow ����
	long* RcvCashFlowSchedule,		// Rcv Forward Start, End, ���, �⸻, ������
	double* Rcv_Slope,				// Rcv Leg �����ݸ� ���� Array

	double* Rcv_CPN,				// Rcv Leg �������� Array
	double* Rcv_FixedRefRate,		// Rcv Leg ���� Ȯ���ݸ� Array
	long Pay_RefRateType,			// Pay ���ʱݸ� 0: Libor/CD 1: Swap 2: SOFR 3:SOFR Swap
	long Pay_SwapYearlyNPayment,	// Pay_RefRateType�� 1, 3�� �� ���� �� ����ȸ��
	double Pay_SwapMaturity,		// Pay_RefRateType�� 1, 3�� �� ���Ҹ���

	long Pay_FixFloFlag,			// Pay Fix/Flo Flag
	long Pay_DayCount,				// Pay DayCountConvention 0:Act365  1: Act360
	double Pay_NotionalAMT,			// Pay Leg Notional Amount
	long Pay_NotionalPayDate,		// Pay Leg Notional Payment Date
	long PayDisc_NTerm,				// Pay Leg ���� �ݸ� Term ����

	double* PayDisc_Term,			// Pay Leg ���� �ݸ� Term Array
	double* PayDisc_Rate,			// Pay Leg ���� �ݸ� Rate Array 
	long PayRef_NTerm,				// Pay Leg ���۷��� �ݸ� Term ����
	double* PayRef_Term,			// Pay Leg ���۷��� �ݸ� Term Array
	double* PayRef_Rate,			// Pay Leg ���� �ݸ� Rate Array 

	long NPayCF,					// Pay Leg CashFlow ����
	long* PayCashFlowSchedule,		// Pay Forward Start, End, ���, �⸻, ������
	double* Pay_Slope,				// Pay Leg �����ݸ� ���� Array
	double* Pay_CPN,				// Pay Leg �������� Array
	double* Pay_FixedRefRate,		// Pay Leg ���� Ȯ���ݸ� Array

	long* HolidayCalcFlag,			// Holiday���� ��ǲ Flag [0]: Rcv [1]: Pay
	long* NHolidays,				// Holiday ���� [0] NRcvRef [1] NPayRef
	long* Holidays					// Holiday Exceltype
)
{
	long i;


	long* Rcv_ForwardStartExl = RcvCashFlowSchedule;
	long* Rcv_ForwardEndExl = RcvCashFlowSchedule + NRcvCF;
	long* Rcv_StartDateExl = RcvCashFlowSchedule + 2 * NRcvCF;
	long* Rcv_EndDateExl = RcvCashFlowSchedule + 3 * NRcvCF;
	long* Rcv_PayDateExl = RcvCashFlowSchedule + 4 * NRcvCF;

	long* Pay_ForwardStartExl = PayCashFlowSchedule;
	long* Pay_ForwardEndExl = PayCashFlowSchedule + NPayCF;
	long* Pay_StartDateExl = PayCashFlowSchedule + 2 * NPayCF;
	long* Pay_EndDateExl = PayCashFlowSchedule + 3 * NPayCF;
	long* Pay_PayDateExl = PayCashFlowSchedule + 4 * NPayCF;

	if (PriceDate_Exl <= 0 || PriceDate_Exl > Rcv_PayDateExl[NRcvCF -1] || PriceDate_Exl > Pay_PayDateExl[NPayCF - 1] || PriceDate_Exl < Rcv_ForwardStartExl[0] - 1000 || PriceDate_Exl < Pay_ForwardStartExl[0] - 1000) return -1;

	if (GreekFlag != 0 && GreekFlag != 1) return -2;

	if (NAFlag != 0 && NAFlag != 1) return -3;

	long CRSPricingFlag = CRS_Flag[0];
	long CRS_NFX_Term = CRS_Flag[1];
	double* CRS_FX_Term = CRS_Info;
	double* CRS_FX_Rate = CRS_Info + CRS_NFX_Term;

	if (CRSPricingFlag != 0 && CRSPricingFlag != 1) return -4;

	if (CRS_FX_Rate[0] < 0.0) return -5;

	for (i = 0; i < CRS_NFX_Term; i++)
	{
		if (i > 0)
		{
			if (CRS_FX_Term[i] < CRS_FX_Term[i - 1]) return -5;
			if (CRS_FX_Rate[i] < 0.0) return -5;
		}
	}

	if (Rcv_RefRateType < 0 || Rcv_RefRateType > 5) return -6;
	if (Rcv_SwapYearlyNPayment <= 0 || Rcv_SwapYearlyNPayment > 12) return -7;
	if (Rcv_SwapMaturity < 0.0) return -8;
	if (Rcv_FixFloFlag != 0 && Rcv_FixFloFlag != 1) return -9;
	if (Rcv_DayCount != 0 && Rcv_DayCount != 1) return -10;
	if (Rcv_NotionalAMT < 0.0) return -11;
	if (Rcv_NotionalPayDate < PriceDate_Exl) return -12;
	if (RcvDisc_NTerm < 1) return -13;
	
	for (i = 0; i < RcvDisc_NTerm; i++)
	{
		if (i > 0)
		{
			if (RcvDisc_Term[i] < RcvDisc_Term[i - 1]) return -14;
		}
	}

	if (RcvRef_NTerm < 1) return -16;
	for (i = 0; i < RcvRef_NTerm; i++)
	{
		if (i > 0)
		{
			if (RcvRef_Term[i] < RcvRef_Term[i - 1]) return -17;
		}
	}

	for (i = 0; i < NRcvCF; i++)
	{
		if (Rcv_ForwardStartExl[i] >= Rcv_ForwardEndExl[i] || Rcv_PayDateExl[i] < Rcv_ForwardEndExl[i]) return -20;
	}

	if (Rcv_RefRateType < 0 || Rcv_RefRateType > 5) return -6;
	if (Rcv_RefRateType == 1 || Rcv_RefRateType == 3)
	{
		if (Rcv_SwapYearlyNPayment <= 0 || Rcv_SwapYearlyNPayment > 12) return -7;
	}
	if (Rcv_SwapMaturity < 0.0) return -8;
	if (Rcv_FixFloFlag != 0 && Rcv_FixFloFlag != 1) return -9;
	if (Rcv_DayCount != 0 && Rcv_DayCount != 1) return -10;
	if (Rcv_NotionalAMT < 0.0) return -11;
	if (Rcv_NotionalPayDate < PriceDate_Exl) return -12;
	if (RcvDisc_NTerm < 1) return -13;
	
	for (i = 0; i < RcvDisc_NTerm; i++)
	{
		if (i > 0)
		{
			if (RcvDisc_Term[i] < RcvDisc_Term[i - 1]) return -14;
		}
	}

	if (RcvRef_NTerm < 1) return -16;
	for (i = 0; i < RcvRef_NTerm; i++)
	{
		if (i > 0)
		{
			if (RcvRef_Term[i] < RcvRef_Term[i - 1]) return -17;
		}
	}

	for (i = 0; i < NRcvCF; i++)
	{
		if (Rcv_ForwardStartExl[i] >= Rcv_ForwardEndExl[i] || Rcv_PayDateExl[i] < Rcv_ForwardEndExl[i]) return -20;
	}

	//
	if (Pay_RefRateType < 0 || Pay_RefRateType > 5) return -6;
	if (Pay_RefRateType == 1 || Pay_RefRateType == 3)
	{
		if (Pay_SwapYearlyNPayment <= 0 || Pay_SwapYearlyNPayment > 12) return -7;
	}
	if (Pay_SwapMaturity < 0.0) return -8;
	if (Pay_FixFloFlag != 0 && Pay_FixFloFlag != 1) return -9;
	if (Pay_DayCount != 0 && Pay_DayCount != 1) return -10;
	if (Pay_NotionalAMT < 0.0) return -11;
	if (Pay_NotionalPayDate < PriceDate_Exl) return -12;
	if (PayDisc_NTerm < 1) return -13;

	for (i = 0; i < PayDisc_NTerm; i++)
	{
		if (i > 0)
		{
			if (PayDisc_Term[i] < PayDisc_Term[i - 1]) return -14;
		}
	}

	if (PayRef_NTerm < 1) return -16;
	for (i = 0; i < PayRef_NTerm; i++)
	{
		if (i > 0)
		{
			if (PayRef_Term[i] < PayRef_Term[i - 1]) return -17;
		}
	}

	for (i = 0; i < NPayCF; i++)
	{
		if (Pay_ForwardStartExl[i] >= Pay_ForwardEndExl[i] || Pay_PayDateExl[i] < Pay_ForwardEndExl[i]) return -20;
	}

	return 0;
}

DLLEXPORT(long) IRSwap_Excel(
	long PriceDate_Exl,				// PricingDate ExcelType
	long GreekFlag,					// Greek���� Flag
	long NAFlag,					// Notional Amount ��� Flag
	long* CRS_Flag,					// [0]CRS Pricing Flag [1]FX���� Term ����
	double* CRS_Info,				// [0~FX����-1] FX Forward Term, [FX����-1~2*FX����-1] FX Forward

	long Rcv_RefRateType,			// Rcv ���ʱݸ� 0: Libor/CD 1: Swap 2: SOFR 3:SOFR Swap
	long Rcv_SwapYearlyNPayment,	// Rcv_RefRateType�� 1, 3�� �� ���� �� ����ȸ��
	double Rcv_SwapMaturity,		// Rcv_RefRateType�� 1, 3�� �� ���Ҹ���
	long Rcv_FixFloFlag,			// Rcv Fix/Flo Flag
	long Rcv_DayCount,				// Rcv DayCountConvention 0:Act365  1: Act360

	double Rcv_NotionalAMT,			// Rcv Leg Notional Amount
	long Rcv_NotionalPayDate,		// Rcv Leg Notional Payment Date
	long RcvDisc_NTerm,				// Rcv Leg ���� �ݸ� Term ����
	double* RcvDisc_Term,			// Rcv Leg ���� �ݸ� Term Array
	double* RcvDisc_Rate,			// Rcv Leg ���� �ݸ� Rate Array 

	long RcvRef_NTerm,				// Rcv Leg ���۷��� �ݸ� Term ����
	double* RcvRef_Term,			// Rcv Leg ���۷��� �ݸ� Term Array
	double* RcvRef_Rate,			// Rcv Leg ���۷��� �ݸ� Rate Array
	long NRcvCF,					// Rcv Leg CashFlow ����
	long* RcvCashFlowSchedule,		// Rcv Forward Start, End, ���, �⸻, ������

	double* Rcv_Slope,				// Rcv Leg �����ݸ� ���� Array
	double* Rcv_CPN,				// Rcv Leg �������� Array
	double* Rcv_FixedRefRate,		// Rcv Leg ���� Ȯ���ݸ� Array
	long Pay_RefRateType,			// Pay ���ʱݸ� 0: Libor/CD 1: Swap 2: SOFR 3:SOFR Swap
	long Pay_SwapYearlyNPayment,	// Pay_RefRateType�� 1, 3�� �� ���� �� ����ȸ��

	double Pay_SwapMaturity,		// Pay_RefRateType�� 1, 3�� �� ���Ҹ���
	long Pay_FixFloFlag,			// Pay Fix/Flo Flag
	long Pay_DayCount,				// Pay DayCountConvention 0:Act365  1: Act360
	double Pay_NotionalAMT,			// Pay Leg Notional Amount
	long Pay_NotionalPayDate,		// Pay Leg Notional Payment Date

	long PayDisc_NTerm,				// Pay Leg ���� �ݸ� Term ����
	double* PayDisc_Term,			// Pay Leg ���� �ݸ� Term Array
	double* PayDisc_Rate,			// Pay Leg ���� �ݸ� Rate Array 
	long PayRef_NTerm,				// Pay Leg ���۷��� �ݸ� Term ����
	double* PayRef_Term,			// Pay Leg ���۷��� �ݸ� Term Array

	double* PayRef_Rate,			// Pay Leg ���� �ݸ� Rate Array 
	long NPayCF,					// Pay Leg CashFlow ����
	long* PayCashFlowSchedule,		// Pay Forward Start, End, ���, �⸻, ������
	double* Pay_Slope,				// Pay Leg �����ݸ� ���� Array
	double* Pay_CPN,				// Pay Leg �������� Array

	double* Pay_FixedRefRate,		// Pay Leg ���� Ȯ���ݸ� Array
	double* ResultPrice,			// Output ����� [0] Current Swap Rate [1] Rcv Value [2] Payment Value
	double* ResultRefRate,			// Output ���ʱݸ� Array
	double* ResultCPN,				// Output ���� ���� Array
	double* ResultDF,				// Output Discount Factor Array

	double* PV01,					// Output PV01[0]RcvDisc [1]RcvRef [2]both [3]PayDisc [4]PayRef [5]both
	double* KeyRateRcvPV01,			// Output Rcv Key Rate PV01 .rehaped(-1)
	double* KeyRatePayPV01,			// Output Pay KeyRate PV01 .reshaped(-1)
	long* SOFRConv,					// [0~2] Rcv LockOut LookBackFlag [3~5] Pay LockOut LookBackFlag
	long* HolidayCalcFlag,			// Holiday���� ��ǲ Flag [0]: Rcv [1]: Pay

	long* NHolidays,				// Holiday ���� [0] NRcvRef [1] NPayRef
	long* Holidays,					// Holiday Exceltype
	long* NHistory,
	long* HistoryDateExl,
	double* HistoryRate
)
{
	long i;
	long j;
	long k;
	long ResultCode = 0;

	ResultCode = ErrorCheckIRSwap_Excel(
		PriceDate_Exl,				GreekFlag,					 NAFlag,					CRS_Flag,					CRS_Info,				
		Rcv_RefRateType,			Rcv_SwapYearlyNPayment,	Rcv_SwapMaturity,		Rcv_FixFloFlag,			Rcv_DayCount,				
		Rcv_NotionalAMT,			Rcv_NotionalPayDate,		RcvDisc_NTerm,				RcvDisc_Term,			 RcvDisc_Rate,			
		RcvRef_NTerm,				RcvRef_Term,			RcvRef_Rate,			 NRcvCF, RcvCashFlowSchedule,			 Rcv_Slope,
		Rcv_CPN,				Rcv_FixedRefRate,		Pay_RefRateType,			Pay_SwapYearlyNPayment,	Pay_SwapMaturity,		
		 Pay_FixFloFlag,			Pay_DayCount,			Pay_NotionalAMT,			 Pay_NotionalPayDate,	PayDisc_NTerm,				
		PayDisc_Term,			PayDisc_Rate,			PayRef_NTerm,				PayRef_Term,			 PayRef_Rate,			
		NPayCF, PayCashFlowSchedule,			Pay_Slope,				Pay_CPN,				Pay_FixedRefRate,		HolidayCalcFlag,
		NHolidays,				Holidays					);

	if (ResultCode < 0) return ResultCode;

	long* Rcv_ForwardStartExl = RcvCashFlowSchedule;
	long* Rcv_ForwardEndExl = RcvCashFlowSchedule + NRcvCF;
	long* Rcv_StartDateExl = RcvCashFlowSchedule + 2 * NRcvCF;
	long* Rcv_EndDateExl = RcvCashFlowSchedule + 3 * NRcvCF;
	long* Rcv_PayDateExl = RcvCashFlowSchedule + 4 * NRcvCF;

	long* Pay_ForwardStartExl = PayCashFlowSchedule;
	long* Pay_ForwardEndExl = PayCashFlowSchedule + NPayCF;
	long* Pay_StartDateExl = PayCashFlowSchedule + 2 * NPayCF;
	long* Pay_EndDateExl = PayCashFlowSchedule + 3 * NPayCF;
	long* Pay_PayDateExl = PayCashFlowSchedule + 4 * NPayCF;

	//////////////////
	// SOFR History //
	long Rcv_NHistory = NHistory[0];
	long Pay_NHistory = NHistory[1];
	long* Rcv_HistoryDateExl = HistoryDateExl;
	long* Pay_HistoryDateExl = HistoryDateExl + Rcv_NHistory;
	double* Rcv_HistoryRate = HistoryRate;
	double* Pay_HistoryRate = HistoryRate + Rcv_NHistory;

	long* Rcv_HistoryRelDate = (long*)malloc(sizeof(long) * Rcv_NHistory);
	long* Pay_HistoryRelDate = (long*)malloc(sizeof(long) * Pay_NHistory);
	for (i = 0; i < Rcv_NHistory; i++) Rcv_HistoryRelDate[i] = Rcv_HistoryDateExl[i] - PriceDate_Exl;
	for (i = 0; i < Pay_NHistory; i++) Pay_HistoryRelDate[i] = Pay_HistoryDateExl[i] - PriceDate_Exl;

	long WeekCheckStart;
	long WeekCheckEnd;

	long NWeekend;

	long* NRcv_Weekend = (long*)malloc(sizeof(long) * NRcvCF);
	long** Rcv_Weekend = (long**)malloc(sizeof(long*) * NRcvCF);
	
	for (i = 0; i < NRcvCF; i++)
	{
		WeekCheckStart = Rcv_ForwardStartExl[i] - 100;
		WeekCheckEnd = Rcv_ForwardEndExl[i];
		NWeekend = 0;
		for (j = WeekCheckStart; j < WeekCheckEnd; j++)
		{
			if (isweekend(j) == 1)
			{
				NWeekend += 1;
			}
		}
		NRcv_Weekend[i] = NWeekend;
		Rcv_Weekend[i] = (long*)malloc(sizeof(long) * max(1,NWeekend) );
		k = 0;
		for (j = WeekCheckStart; j < WeekCheckEnd; j++)
		{
			if (isweekend(j) == 1)
			{
				Rcv_Weekend[i][k] =  j - PriceDate_Exl;
				k += 1;
			}
		}
	}

	long* NPay_Weekend = (long*)malloc(sizeof(long) * NPayCF);
	long** Pay_Weekend = (long**)malloc(sizeof(long*) * NPayCF);

	for (i = 0; i < NPayCF; i++)
	{
		WeekCheckStart = Pay_ForwardStartExl[i] - 100;
		WeekCheckEnd = Pay_ForwardEndExl[i];
		NWeekend = 0;
		for (j = WeekCheckStart; j < WeekCheckEnd; j++)
		{
			if (isweekend(j) == 1)
			{
				NWeekend += 1;
			}
		}
		NPay_Weekend[i] = NWeekend;
		Pay_Weekend[i] = (long*)malloc(sizeof(long) * max(1, NWeekend));
		k = 0;
		for (j = WeekCheckStart; j < WeekCheckEnd; j++)
		{
			if (isweekend(j) == 1)
			{
				Pay_Weekend[i][k] = j - PriceDate_Exl;
				k += 1;
			}
		}
	}

	//////////////////////
	// CRS ���� �� ����
	//////////////////////
	long CalcCRSFlag = CRS_Flag[0];
	long Rcv_NTermFX = CRS_Flag[1];
	long Pay_NTermFX = CRS_Flag[1];

	double* Rcv_TermFX = CRS_Info;
	double* Pay_TermFX = CRS_Info;
	double* Rcv_FX = CRS_Info + CRS_Flag[1];
	double* Pay_FX = CRS_Info + 2 * CRS_Flag[1];

	long NHolidays_Array[2] = { 0,0};
	long MaxDay;
	long MinDay;

	////////////////////////////////////////////
	// Excel Date Type -> C Date Type
	////////////////////////////////////////////

	long PriceDate_C = ExcelDateToCDate(PriceDate_Exl);

	long* Rcv_ForwardStart_C = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Rcv_ForwardEnd_C = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Rcv_StartDate_C = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Rcv_EndDate_C = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Rcv_PayDate_C = (long*)malloc(sizeof(long) * (NRcvCF));
	for (i = 0; i < NRcvCF; i++)
	{
		Rcv_ForwardStart_C[i] = ExcelDateToCDate(Rcv_ForwardStartExl[i]);
		Rcv_ForwardEnd_C[i] = ExcelDateToCDate(Rcv_ForwardEndExl[i]);
		Rcv_StartDate_C[i] = ExcelDateToCDate(Rcv_StartDateExl[i]);
		Rcv_EndDate_C[i] = ExcelDateToCDate(Rcv_EndDateExl[i]);
		Rcv_PayDate_C[i] = ExcelDateToCDate(Rcv_PayDateExl[i]);
	}


	long* Pay_ForwardStart_C = (long*)malloc(sizeof(long) * (NPayCF));
	long* Pay_ForwardEnd_C = (long*)malloc(sizeof(long) * (NPayCF));
	long* Pay_StartDate_C = (long*)malloc(sizeof(long) * (NPayCF));
	long* Pay_EndDate_C = (long*)malloc(sizeof(long) * (NPayCF));
	long* Pay_PayDate_C = (long*)malloc(sizeof(long) * (NPayCF));
	for (i = 0; i < NPayCF; i++)
	{
		Pay_ForwardStart_C[i] = ExcelDateToCDate(Pay_ForwardStartExl[i]);
		Pay_ForwardEnd_C[i] = ExcelDateToCDate(Pay_ForwardEndExl[i]);
		Pay_StartDate_C[i] = ExcelDateToCDate(Pay_StartDateExl[i]);
		Pay_EndDate_C[i] = ExcelDateToCDate(Pay_EndDateExl[i]);
		Pay_PayDate_C[i] = ExcelDateToCDate(Pay_PayDateExl[i]);
	}

	long* Days_Rcv_ForwardStart = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Days_Rcv_ForwardEnd = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Days_Rcv_StartDate = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Days_Rcv_EndDate = (long*)malloc(sizeof(long) * (NRcvCF));
	long* Days_Rcv_PayDate = (long*)malloc(sizeof(long) * (NRcvCF));
	for (i = 0; i < NRcvCF; i++)
	{
		Days_Rcv_ForwardStart[i] = DayCountAtoB(PriceDate_C, Rcv_ForwardStart_C[i]);
		Days_Rcv_ForwardEnd[i] = DayCountAtoB(PriceDate_C, Rcv_ForwardEnd_C[i]);
		Days_Rcv_StartDate[i] = DayCountAtoB(PriceDate_C, Rcv_StartDate_C[i]);
		Days_Rcv_EndDate[i] = DayCountAtoB(PriceDate_C, Rcv_EndDate_C[i]);
		Days_Rcv_PayDate[i] = DayCountAtoB(PriceDate_C, Rcv_PayDate_C[i]);
	}

	long* Days_Pay_ForwardStart = (long*)malloc(sizeof(long) * (NPayCF));
	long* Days_Pay_ForwardEnd = (long*)malloc(sizeof(long) * (NPayCF));
	long* Days_Pay_StartDate = (long*)malloc(sizeof(long) * (NPayCF));
	long* Days_Pay_EndDate = (long*)malloc(sizeof(long) * (NPayCF));
	long* Days_Pay_PayDate = (long*)malloc(sizeof(long) * (NPayCF));
	for (i = 0; i < NPayCF; i++)
	{
		Days_Pay_ForwardStart[i] = DayCountAtoB(PriceDate_C, Pay_ForwardStart_C[i]);
		Days_Pay_ForwardEnd[i] = DayCountAtoB(PriceDate_C, Pay_ForwardEnd_C[i]);
		Days_Pay_StartDate[i] = DayCountAtoB(PriceDate_C, Pay_StartDate_C[i]);
		Days_Pay_EndDate[i] = DayCountAtoB(PriceDate_C, Pay_EndDate_C[i]);
		Days_Pay_PayDate[i] = DayCountAtoB(PriceDate_C, Pay_PayDate_C[i]);
	}

	long RcvRef_HolidayCalc = HolidayCalcFlag[0];
	long PayRef_HolidayCalc = HolidayCalcFlag[1];

	long RcvRef_NHolidays = NHolidays[0];
	long PayRef_NHolidays = NHolidays[1];

	long* RcvRef_Holidays_Exl = Holidays ;
	long* PayRef_Holidays_Exl = Holidays + RcvRef_HolidayCalc ;

	////////////////////////////////////////////
	// Holiday ����, �ʹ� ���ų� �ʹ� �̷��ΰ� ����
	////////////////////////////////////////////

	MaxDay = max(Rcv_PayDateExl[NPayCF - 1], Pay_PayDateExl[NPayCF - 1]);
	MinDay = min(min(PriceDate_Exl, min(Rcv_ForwardStartExl[0], Pay_ForwardStartExl[0])), min(Rcv_StartDateExl[0], Pay_StartDateExl[0]));
	// RcvRef
	NHolidays_Array[0] = 0;
	for (i = 0; i < RcvRef_NHolidays; i++) {
		if (RcvRef_Holidays_Exl[i] <= MaxDay && RcvRef_Holidays_Exl[i] > MinDay - 100) { NHolidays_Array[0] += 1; }
	}
	long* RcvRef_Days_Holiday = (long*)malloc(sizeof(long) * max(1,NHolidays_Array[0]));
	NHolidays_Array[0] = 0;
	for (i = 0; i < RcvRef_NHolidays; i++) {
		if (RcvRef_Holidays_Exl[i] <= MaxDay && RcvRef_Holidays_Exl[i] > MinDay - 100) {
			RcvRef_Days_Holiday[NHolidays_Array[0]] = DayCountAtoB(PriceDate_C, ExcelDateToCDate(RcvRef_Holidays_Exl[i]));
			NHolidays_Array[0] += 1;
		}
	}
	// PayRef
	NHolidays_Array[1] = 0;
	for (i = 0; i < PayRef_NHolidays; i++) {
		if (PayRef_Holidays_Exl[i] <= MaxDay && PayRef_Holidays_Exl[i] > MinDay - 100) { NHolidays_Array[1] += 1; }
	}
	long* PayRef_Days_Holiday = (long*)malloc(sizeof(long) * max(1,NHolidays_Array[1]));
	NHolidays_Array[1] = 0;
	for (i = 0; i < PayRef_NHolidays; i++) {
		if (PayRef_Holidays_Exl[i] <= MaxDay && PayRef_Holidays_Exl[i] > MinDay - 100) {
			PayRef_Days_Holiday[NHolidays_Array[1]] = DayCountAtoB(PriceDate_C, ExcelDateToCDate(PayRef_Holidays_Exl[i]));
			NHolidays_Array[1] += 1;
		}
	}

	////////////////////////////////
	// Receive Schedule Mapping  ///
	////////////////////////////////

	SCHD* Rcv_Schedule = new SCHD;
	Rcv_Schedule->HolidayFlag_Ref = RcvRef_HolidayCalc;
	Rcv_Schedule->NHolidays_Ref = NHolidays_Array[0];
	Rcv_Schedule->Days_Holidays_Ref = RcvRef_Days_Holiday;
	Rcv_Schedule->NCF = NRcvCF;
	Rcv_Schedule->ForwardStart_C = Rcv_ForwardStart_C;
	Rcv_Schedule->ForwardEnd_C = Rcv_ForwardEnd_C;
	Rcv_Schedule->StartDate_C = Rcv_StartDate_C;
	Rcv_Schedule->EndDate_C = Rcv_EndDate_C;
	Rcv_Schedule->PayDate_C = Rcv_PayDate_C;
	Rcv_Schedule->Days_ForwardStart = Days_Rcv_ForwardStart;
	Rcv_Schedule->Days_ForwardEnd = Days_Rcv_ForwardEnd;
	Rcv_Schedule->Days_StartDate = Days_Rcv_StartDate;
	Rcv_Schedule->Days_EndDate = Days_Rcv_EndDate;
	Rcv_Schedule->Days_PayDate = Days_Rcv_PayDate;
	Rcv_Schedule->NotionalPayDate_C = ExcelDateToCDate(Rcv_NotionalPayDate);
	Rcv_Schedule->Days_Notional = DayCountAtoB(PriceDate_C, Rcv_Schedule->NotionalPayDate_C);
	Rcv_Schedule->ReferenceType = Rcv_RefRateType;
	Rcv_Schedule->FixedFlotype = Rcv_FixFloFlag;
	Rcv_Schedule->DayCount = Rcv_DayCount;
	Rcv_Schedule->NotionalAmount = Rcv_NotionalAMT;
	Rcv_Schedule->NWeekendDate = NRcv_Weekend;
	Rcv_Schedule->WeekendList = Rcv_Weekend;
	Rcv_Schedule->NRefHistory = Rcv_NHistory;
	Rcv_Schedule->RefHistoryDate = Rcv_HistoryRelDate;
	Rcv_Schedule->RefHistory = Rcv_HistoryRate;

	if (Rcv_RefRateType != 0 && Rcv_RefRateType != 2) Rcv_Schedule->RefSwapFlag = 1;
	Rcv_Schedule->NSwapPayAnnual = Rcv_SwapYearlyNPayment;
	Rcv_Schedule->RefSwapMaturity = Rcv_SwapMaturity;
	Rcv_Schedule->FixedRefRate = Rcv_FixedRefRate;
	Rcv_Schedule->Slope = Rcv_Slope;
	Rcv_Schedule->CPN = Rcv_CPN;
	Rcv_Schedule->NAFlag = NAFlag;
	Rcv_Schedule->PriceDate_C = PriceDate_C;
	
	Rcv_Schedule->LockOutRef = SOFRConv[0];
	Rcv_Schedule->LookBackRef = SOFRConv[1];
	Rcv_Schedule->ObservationShift = SOFRConv[2];

	////////////////////////////
	// Pay Schedule Mapping  ///
	////////////////////////////

	SCHD* Pay_Schedule = new SCHD;
	Pay_Schedule->HolidayFlag_Ref = RcvRef_HolidayCalc;
	Pay_Schedule->NHolidays_Ref = NHolidays_Array[1];
	Pay_Schedule->Days_Holidays_Ref = PayRef_Days_Holiday;
	Pay_Schedule->NCF = NPayCF;
	Pay_Schedule->ForwardStart_C = Pay_ForwardStart_C;
	Pay_Schedule->ForwardEnd_C = Pay_ForwardEnd_C;
	Pay_Schedule->StartDate_C = Pay_StartDate_C;
	Pay_Schedule->EndDate_C = Pay_EndDate_C;
	Pay_Schedule->PayDate_C = Pay_PayDate_C;
	Pay_Schedule->Days_ForwardStart = Days_Pay_ForwardStart;
	Pay_Schedule->Days_ForwardEnd = Days_Pay_ForwardEnd;
	Pay_Schedule->Days_StartDate = Days_Pay_StartDate;
	Pay_Schedule->Days_EndDate = Days_Pay_EndDate;
	Pay_Schedule->Days_PayDate = Days_Pay_PayDate;
	Pay_Schedule->NotionalPayDate_C = ExcelDateToCDate(Pay_NotionalPayDate);
	Pay_Schedule->Days_Notional = DayCountAtoB(PriceDate_C, Pay_Schedule->NotionalPayDate_C);
	Pay_Schedule->ReferenceType = Pay_RefRateType;
	Pay_Schedule->FixedFlotype = Pay_FixFloFlag;
	Pay_Schedule->DayCount = Pay_DayCount;
	Pay_Schedule->NotionalAmount = Pay_NotionalAMT;
	Pay_Schedule->NWeekendDate = NPay_Weekend;
	Pay_Schedule->WeekendList = Pay_Weekend;
	Pay_Schedule->NRefHistory = Pay_NHistory;
	Pay_Schedule->RefHistoryDate = Pay_HistoryRelDate;
	Pay_Schedule->RefHistory = Pay_HistoryRate;

	if (Pay_RefRateType != 0 && Pay_RefRateType != 2) Pay_Schedule->RefSwapFlag = 1;
	Pay_Schedule->NSwapPayAnnual = Pay_SwapYearlyNPayment;
	Pay_Schedule->RefSwapMaturity = Pay_SwapMaturity;
	Pay_Schedule->FixedRefRate = Pay_FixedRefRate;
	Pay_Schedule->Slope = Pay_Slope;
	Pay_Schedule->CPN = Pay_CPN;
	Pay_Schedule->NAFlag = NAFlag;
	Pay_Schedule->PriceDate_C = PriceDate_C;

	Pay_Schedule->LockOutRef = SOFRConv[3];
	Pay_Schedule->LookBackRef = SOFRConv[4];
	Pay_Schedule->ObservationShift = SOFRConv[5];

	long PricingOnly = 1;
	ResultCode = SwapPricer(PriceDate_C, CalcCRSFlag, Rcv_Schedule, RcvDisc_NTerm, RcvDisc_Term,
							RcvDisc_Rate, RcvRef_NTerm, RcvRef_Term, RcvRef_Rate,
							Pay_Schedule, PayDisc_NTerm, PayDisc_Term, PayDisc_Rate,
							PayRef_NTerm, PayRef_Term, PayRef_Rate,
							Rcv_NTermFX, Rcv_TermFX, Rcv_FX, Pay_NTermFX, Pay_TermFX, Pay_FX,
							GreekFlag, ResultPrice, ResultRefRate, ResultCPN, ResultDF,
							PV01, KeyRateRcvPV01, KeyRatePayPV01, PricingOnly);

	free(Rcv_HistoryRelDate);
	free(Pay_HistoryRelDate);

	free(NRcv_Weekend);
	for (i = 0; i < NRcvCF; i++) free(Rcv_Weekend[i]);
	free(Rcv_Weekend);

	free(NPay_Weekend);
	for (i = 0; i < NPayCF; i++) free(Pay_Weekend[i]);
	free(Pay_Weekend);


	free(Rcv_ForwardStart_C);
	free(Rcv_ForwardEnd_C);
	free(Rcv_StartDate_C);
	free(Rcv_EndDate_C);
	free(Rcv_PayDate_C);

	free(Pay_ForwardStart_C);
	free(Pay_ForwardEnd_C);
	free(Pay_StartDate_C);
	free(Pay_EndDate_C);
	free(Pay_PayDate_C);
	free(RcvRef_Days_Holiday);
	free(PayRef_Days_Holiday);

	delete (Rcv_Schedule);
	delete (Pay_Schedule);

	return ResultCode;
}