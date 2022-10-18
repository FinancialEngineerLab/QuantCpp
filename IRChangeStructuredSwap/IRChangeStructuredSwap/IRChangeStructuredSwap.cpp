#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "HW_Calibration.h"
#include "ErrorCheck.h"

#ifndef UTILITY
#include "Util.h"
#endif

#include "CalcDate.h"
#include "Structure.h"
#include <crtdbg.h>


DLLEXPORT(long) isin(long x, long* array, long narray)
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

long isinFindIndex(long x, long* array, long narray, long& Idx)
{
    long i;
    long s = 0;
    Idx = -1;
    for (i = 0; i < narray; i++)
    {
        if (x == array[i])
        {
            s = 1;
            Idx = i;
            break;
        }
    }
    return s;
}

long isbetweenFindIndex(long x, long* array, long narray, long& idx)
{
    long i;
    long s = 0;
    long temp = -1;

    for (i = idx; i < narray - 1; i++)
    {
        if (x >= array[i] && x < array[i + 1])
        {
            s = 1;
            temp = 1;
            idx = temp;
            break;
        }
    }
    return s;
}

long FindIndex(
    long x,
    long* array,
    long narray
)
{
    long i;
    for (i = 0; i < narray; i++)
    {
        if (x == array[i])
        {
            return i;
        }
    }
    return i;
}

//레퍼런스금리정보
typedef struct ReferenceInfo {
    long CurveNum;              // 커브 번호
    long RefRateType;           // 기초금리 타입 0: CD/Libor (Not Accrual) 1: CD/Libor 또는 이자율스왑 2: SOFR
    long RefSwapNCPN_Ann;       // 기초금리 타입이 스왑일 경우 연 이자지급 수
    double RefSwapMaturity;     // 기초금리 타입이 스왑일 경우 스왑 만기

    long DayCountFlag;          // 0: Day1Y = 365, 1: Day1Y = 360
    double Day1Y;
    double RefRateCondMultiple; // 기준금리 결정조건 Multiple
    long PowerSpreadFlag;       // 파워스프레드 사용여부
    double PowerSpreadMat1;     // 파워스프레드 만기1 
    double PowerSpreadMat2;     // 파워스프레드 만기2
    double RangeUp;             // 쿠폰 상방배리어
    double RangeDn;             // 쿠폰 하방배리어
    double PayoffMultiple;      // 페이오프 결정조건 Multiple

} REFERENCE_INFO;

typedef struct LegInfo {
    long FixFloFlag;                    // 0: 고정금리, 1:변동금리
    long AccrualFlag;                   // Range Accrual 사용여부
    double MaxLossY;                    // 최대손실
    double MaxRetY;                     // 최대이익
    long NReference;                    // 기초금리개수
    REFERENCE_INFO* Reference_Inform;   // 레퍼런스금리정보
    long DiscCurveNum;                  // 할인 커브 번호
    long DiscDayCountFlag;              // 할인 DayCountYear 0:365 1:360

    long NCashFlow;                     // 현금흐름지급개수
    long* ForwardStart_C;               // 금리 추정(또는 리셋) 시작일 Ctype
    long* ForwardEnd_C;                 // 금리 추정(또는 리셋) 종료일 Ctype
    long* DaysForwardStart;             // 평가일 to 금리 추정(또는 리셋) 시작일
    long* DaysForwardEnd;               // 평가일 to 금리 추정(또는 리셋) 종료일
    long* FractionStart_C;              // 기산일 Ctype
    long* FractionEnd_C;                // 기말일 Ctype
    long* PayDate_C;                    // 지급일 Ctype
    long* DaysPayDate_C;                // 평가일 to 지급일
    double* CouponRate;                 // 쿠폰Rate
    double* FixedRate;                  // 과거 확정금리 (PlainVanila, SOFR)
    double* RangeCoupon;                // Range충족시 쿠폰
    long PayoffStructure;               // 페이오프 조건 0: &조건 1: 합조건

    long nUsingCurve;                   // 사용되는커브숫자
    long* HolidayCalcFlag;              // Holiday계산Flag
    long* HolidayCount;                 // Holiday수
    long** HolidayDays;                 // 평가일 to 각 레퍼런스별 Holiday까지 날짜 수

    long* NDayHistory;                  // 과거 레퍼런스 History 개수
    long** RateHistoryDateMatrix;       // 과거 레퍼런스 History Date Matrix
    double** RateHistoryMatrix;         // 과거 레퍼런스 History Rate Matrix

    // 옵션관련 데이터
    long OptionUseFlag;                 // 옵션 사용여부
    long NOption;                       // 옵션 개수
    long CallConditionFlag;             // 행사조건Flag 0:Range1&Range2&Range3 1:Range1||Range2||Range3
    long OptionType;                    // 옵션타입 0:Payer가 Call옵션 1:Receiver가 풋옵션
    long* DaysOptionDate;               // Pricing Date To Option Date
    double* StrikeRate;                 // 옵션행사비율
    double** RangeUp;                   // Ref 123 Range 상한
    double** RangeDn;                   // Ref 123 Range 하한


}LEG_INFO;



typedef struct SimulationInfo {
    long NSimul;
    long NDays;
    long NAsset;
    long DailySimulFlag;
    long* DaysForSimul;
    double* dt_Array;
    double* T_Array;
    double** FixedRandn;
    long* NHWVol;
    double* HWKappa;
    double** HWVolTerm;
    double** HWVol;
    long* NRateTerm;
    double** RateTerm;
    double** Rate;
    long* SimulCurveIdx;

}SIMUL_INFO;

void bubble_sort_long(long* arr, long count, long ascending)
{
    long temp;
    long i, j;
    if (ascending == 1)
    {
        for (i = 0; i < count; i++)
        {
            for (j = 0; j < count - 1; j++)
            {
                if (arr[j] > arr[j + 1])          // 이전 값이 더 크면
                {                                 // 이전 값을 다음 요소로 보내고 다음 요소를 이전 요소 자리로
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < count; i++)
        {
            for (j = 0; j < count - 1; j++)
            {
                if (arr[j] < arr[j + 1])
                {
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
    }
}

long* Make_Unique_Array(long n_array, long* Sorted_Array, long& target_length)
{
    long i;
    long n;
    long unqvalue = Sorted_Array[0];
    n = 1;
    for (i = 1; i < n_array; i++)
    {
        if (Sorted_Array[i] != unqvalue)
        {
            n += 1;
            unqvalue = Sorted_Array[i];
        }
    }

    target_length = n;
    long* ResultArray = (long*)malloc(sizeof(long) * n);
    ResultArray[0] = Sorted_Array[0];
    unqvalue = Sorted_Array[0];
    n = 1;
    for (i = 1; i < n_array; i++)
    {
        if (Sorted_Array[i] != unqvalue)
        {
            ResultArray[n] = Sorted_Array[i];
            unqvalue = Sorted_Array[i];
            n += 1;
        }
    }

    return ResultArray;
}

double FSR(
    double* Term,
    double* Rate,
    long NTerm,
    double T_Option,
    double Tenor,
    double FreqMonth,
    long ndates,
    long* dates
)
{
    long i;
    double Swap_Rate;

    double* P0T = (double*)malloc(sizeof(double) * (ndates + 1));
    P0T[0] = Calc_Discount_Factor(Term, Rate, NTerm, T_Option); // 옵션만기시점
    for (i = 1; i < ndates + 1; i++)
    {
        P0T[i] = Calc_Discount_Factor(Term, Rate, NTerm, (double)dates[i - 1] / 365.0);
    }

    double a, b, dt;
    a = P0T[0] - P0T[ndates];
    b = 0.0;
    for (i = 0; i < ndates; i++)
    {
        if (i == 0)
            dt = (double)dates[0] / 365.0 - T_Option;
        else
            dt = ((double)(dates[i] - dates[i - 1])) / 365.0;
        b += dt * P0T[i + 1];
    }
    Swap_Rate = a / b;



    if (P0T) free(P0T);
    return Swap_Rate;
}

//A_k = exp(-kappa * (t_(k+1) - t_(k))
double XA(
    double kappa,
    double t0,
    double t1
)
{
    double A_k;
    A_k = exp(-kappa * (t1 - t0));
    return A_k;
}

void Set_XA(
    double kappa,
    long Length_XA_Array,
    double* XA_Array,
    double* dt_Array
)
{
    long i;
    double T, t0, t1;
    T = dt_Array[0];
    for (i = 0; i < Length_XA_Array; i++)
    {
        t0 = T;
        t1 = t0 + dt_Array[i];
        XA_Array[i] = XA(kappa, t0, t1);
        T = t1;
    }
}

// XV = sqrt( int_t0^t1 [ e^(-2kappa*(t1-s)) * sig(s)^2 * ds   ])
double XV(
    double kappa,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    double t0,
    double t1
)
{
    long i;
    long nSquare = 4; // 적분구간을 nSquare만큼 나눔
    double ds = (t1 - t0) / (double)nSquare;
    double s = ds;
    double vol = 0.0;
    double B_k, B_k_Square;

    B_k_Square = 0.0;
    for (i = 0; i < nSquare; i++)
    {
        vol = Interpolate_Linear(HWVolTerm, HWVol, NHWVol, s);
        B_k_Square += exp(-2.0 * kappa * (t1 - s)) * vol * vol * ds;
        s += ds;
    }
    B_k = sqrt(B_k_Square);
    return B_k;
}

void Set_XV(
    double kappa,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    long Length_XV_Array,
    double* XV_Array,
    double* dt_Array
)
{
    long i;
    double T, t0, t1;
    T = dt_Array[0];
    for (i = 0; i < Length_XV_Array; i++)
    {
        t0 = T;
        t1 = t0 + dt_Array[i];
        XV_Array[i] = XV(kappa, NHWVol, HWVolTerm, HWVol, t0, t1);
        T = t1;
    }
}

double Simulated_ShortRate(
    double XA,
    double XV,
    double ShortRate_Prev,
    double epsilon
)
{
    double xt = ShortRate_Prev;
    double x_next = XA * xt + XV * epsilon;
    return x_next;
}

double B_s_to_t(
    double kappa,
    double s,
    double t
)
{
    return (1.0 - exp(-kappa * (t - s))) / kappa;
}


double P_t_to_T(
    double t,
    double T,
    double kappa,
    double ShortRate,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    double Pt,
    double PT,
    long nInteg
)
{
    long i;
    double PtT;
    double Bst, BsT, BtT;

    double s;
    double vol;
    double RHS = 0.0;
    //long nInteg = 4;
    double ds;
    double xt = ShortRate;

    ds = t / (double)nInteg;
    s = ds;

    for (i = 0; i < nInteg; i++)
    {
        vol = Interpolate_Linear(HWVolTerm, HWVol, NHWVol, s);
        Bst = B_s_to_t(kappa, s, t);
        BsT = B_s_to_t(kappa, s, T);
        RHS += 0.5 * vol * vol * (Bst * Bst - BsT * BsT) * ds;
        s += ds;
    }

    BtT = B_s_to_t(kappa, t, T);
    PtT = PT / Pt * exp(-xt * BtT + RHS);
    return PtT;
}

double HullWhiteQVTerm(
    double t,
    double T,
    double kappa,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol
)
{
    long i;
    double Bst, BsT, vol;
    double RHS = 0.0;
    long nInteg = 4;
    double s, ds;

    ds = t / (double)nInteg;
    s = ds;

    for (i = 0; i < nInteg; i++)
    {
        vol = Interpolate_Linear(HWVolTerm, HWVol, NHWVol, s);
        Bst = B_s_to_t(kappa, s, t);
        BsT = B_s_to_t(kappa, s, T);
        RHS += 0.5 * vol * vol * (Bst * Bst - BsT * BsT) * ds;
        s += ds;
    }
    return RHS;
}

void Set_HW_Params_RefType0(
    SIMUL_INFO* Simul,
    LEG_INFO* Leg,
    long ReferenceIdx,
    long Curveidx,
    double** PV_t_T,
    double** QVTerm,
    double** B_t_T,
    double** PV_t_T_PowerSpread,
    double** QVTerm_PowerSpread,
    double** B_t_T_PowerSpread
)
{
    long i;
    long Dayidx = -1;
    double Day1Y;
    if (Leg->Reference_Inform[ReferenceIdx].DayCountFlag == 0) Day1Y = 365.0;
    else Day1Y = 360.0;

    double T1, T2, P1, P2;

    if (Leg->Reference_Inform[ReferenceIdx].PowerSpreadFlag == 0)
    {
        if (Simul->DailySimulFlag == 0)
        {
            for (i = 0; i < Simul->NDays; i++)
            {

                T1 = ((double)Leg->DaysForwardStart[i]) / Day1Y;
                T2 = ((double)Leg->DaysForwardEnd[i]) / Day1Y;
                P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                PV_t_T[i][0] = P2 / P1;
                QVTerm[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                B_t_T[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);
            }
        }
        else
        {
            for (i = 0; i < Simul->NDays; i++)
            {
                if (isinFindIndex(i, Leg->DaysForwardStart, Leg->NCashFlow, Dayidx))
                {
                    T1 = ((double)Leg->DaysForwardStart[Dayidx]) / Day1Y;
                    T2 = ((double)Leg->DaysForwardEnd[Dayidx]) / Day1Y;
                    P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                    P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                    PV_t_T[i][0] = P2 / P1;
                    QVTerm[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                    B_t_T[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);
                }
            }
        }
    }
    else
    {
        if (Simul->DailySimulFlag == 0)
        {
            for (i = 0; i < Simul->NDays; i++)
            {
                T1 = ((double)Leg->DaysForwardStart[i]) / Day1Y;
                T2 = T1 + Leg->Reference_Inform[ReferenceIdx].PowerSpreadMat1;
                P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                PV_t_T[i][0] = P2 / P1;
                QVTerm[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                B_t_T[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);

                T1 = ((double)Leg->DaysForwardStart[i]) / Day1Y;
                T2 = T1 + Leg->Reference_Inform[ReferenceIdx].PowerSpreadMat2;
                P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                PV_t_T_PowerSpread[i][0] = P2 / P1;
                QVTerm_PowerSpread[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                B_t_T_PowerSpread[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);
            }
        }
        else
        {
            for (i = 0; i < Simul->NDays; i++)
            {
                if (isinFindIndex(i, Leg->DaysForwardStart, Leg->NCashFlow, Dayidx))
                {
                    T1 = ((double)Leg->DaysForwardStart[Dayidx]) / Day1Y;
                    T2 = T1 + Leg->Reference_Inform[ReferenceIdx].PowerSpreadMat1;
                    P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                    P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                    PV_t_T[i][0] = P2 / P1;
                    QVTerm[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                    B_t_T[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);

                    T1 = ((double)Leg->DaysForwardStart[Dayidx]) / Day1Y;
                    T2 = T1 + Leg->Reference_Inform[ReferenceIdx].PowerSpreadMat2;
                    P1 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T1);
                    P2 = Calc_Discount_Factor(Simul->RateTerm[Curveidx], Simul->Rate[Curveidx], Simul->NRateTerm[Curveidx], T2);
                    PV_t_T_PowerSpread[i][0] = P2 / P1;
                    QVTerm_PowerSpread[i][0] = HullWhiteQVTerm(T1, T2, Simul->HWKappa[Curveidx], Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx]);
                    B_t_T_PowerSpread[i][0] = B_s_to_t(Simul->HWKappa[Curveidx], T1, T2);
                }
            }
        }
    }

}

void Set_HW_Params_RefType1(
    long NDays,
    double* T_Array,
    long NCF_Swap,
    long NCPN_Ann,
    double ZeroMaturity,
    double kappa,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    long NCurve,
    double* TCurve,
    double* Curve,
    double** PV_t_T,
    double** QVTerm,
    double** B_t_T
)
{
    long i, j;
    double term = 1.0 / ((double)NCPN_Ann);
    double T1, T2;
    double P1, P2;
    if (NCPN_Ann > 0)
    {
        for (i = 0; i < NDays; i++)
        {
            T1 = T_Array[i];
            for (j = 0; j < NCF_Swap; j++)
            {
                T2 = T1 + ((double)j + 1.0) * term;
                P1 = Calc_Discount_Factor(TCurve, Curve, NCurve, T1);
                P2 = Calc_Discount_Factor(TCurve, Curve, NCurve, T2);
                PV_t_T[i][j] = P2 / P1;
                QVTerm[i][j] = HullWhiteQVTerm(T1, T2, kappa, NHWVol, HWVolTerm, HWVol);
                B_t_T[i][j] = B_s_to_t(kappa, T1, T2);
            }
        }
    }
    else
    {
        // ZeroRate의 경우
        term = ZeroMaturity;
        for (i = 0; i < NDays; i++)
        {
            T1 = T_Array[i];
            T2 = T1 + term;
            P1 = Calc_Discount_Factor(TCurve, Curve, NCurve, T1);
            P2 = Calc_Discount_Factor(TCurve, Curve, NCurve, T2);
            PV_t_T[i][0] = P2 / P1;
            QVTerm[i][0] = HullWhiteQVTerm(T1, T2, kappa, NHWVol, HWVolTerm, HWVol);
            B_t_T[i][0] = B_s_to_t(kappa, T1, T2);
        }
    }

}

void Set_HW_Params_RefType2(
    long NDays,
    long NDays_t_to_T,
    double* T_Array,
    double kappa,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    long NCurve,
    double* TCurve,
    double* Curve,
    double** PV_t_T,
    double** QVTerm,
    double** B_t_T
)
{
    long i, j;
    double term = 1.0 / 365.0;
    double T;
    double T1, T2;
    double P1, P2;
    for (i = 0; i < NDays; i++)
    {
        T = T_Array[i];
        for (j = 0; j < NDays_t_to_T; j++)
        {
            T1 = T + ((double)j * term);
            T2 = T + ((double)j + 1.0) * term;
            P1 = Calc_Discount_Factor(TCurve, Curve, NCurve, T1);
            P2 = Calc_Discount_Factor(TCurve, Curve, NCurve, T2);
            PV_t_T[i][j] = P2 / P1;
            QVTerm[i][j] = HullWhiteQVTerm(T1, T2, kappa, NHWVol, HWVolTerm, HWVol);
            B_t_T[i][j] = B_s_to_t(kappa, T1, T2);
        }
    }

}

double Calc_Forward_t_to_T_HW(
    double t,
    double T,
    double kappa,
    double ShortRate,
    long NHWVol,
    double* HWVolTerm,
    double* HWVol,
    long NRate,
    double* RateTerm,
    double* Rate
)
{

    double Pt = Calc_Discount_Factor(RateTerm, Rate, NRate, t);
    double PT = Calc_Discount_Factor(RateTerm, Rate, NRate, T);
    double PtT;
    double xt = ShortRate;
    double r;

    if (t != 0.0)
    {
        PtT = P_t_to_T(t, T, kappa, xt, NHWVol, HWVolTerm, HWVol, Pt, PT, 4);
    }
    else
    {
        PtT = Calc_Discount_Factor(RateTerm, Rate, NRate, 1.0 / 365.0);
    }
    r = 365.0 * (1.0 / PtT - 1.0);
    return r;
}

double HW_Rate(
    long ReferenceType,
    double t,
    double T,
    long NRateTerm,
    double* RateTerm,
    double* Rate,
    double ShortRate,
    long NCfSwap,
    long NCPN_Ann,
    double* PV_t_T,
    double* QVTerm,
    double* B_t_T
)
{
    long i;
    double PtT;
    double term = 1.0 / ((double)NCPN_Ann);
    double dt = T - t;
    double ResultRate;
    double A, B;
    if (ReferenceType == 0)
    {
        if (t >= 0.0) PtT = PV_t_T[0] * exp(-ShortRate * B_t_T[0] + QVTerm[0]);
        else PtT = Calc_Discount_Factor(RateTerm, Rate, NRateTerm, T);

        if (dt < 1.0) ResultRate = (1.0 - PtT) / (dt * PtT);
        else ResultRate = -1.0 / dt * log(PtT);
    }
    else
    {
        if (NCPN_Ann > 0)
        {
            B = 0.0;
            for (i = 0; i < NCfSwap; i++)
            {
                PtT = PV_t_T[i] * exp(-ShortRate * B_t_T[i] + QVTerm[i]);
                B += term * PtT;
            }

            if (i != 0)
            {
                A = 1.0 - PtT;
            }
            else
            {
                A = 1.0 - PV_t_T[0];
                B = term * PV_t_T[0];
            }
            ResultRate = A / B;
        }
        else
        {
            if (t >= 0.0) PtT = PV_t_T[0] * exp(-ShortRate * B_t_T[0] + QVTerm[0]);
            else PtT = Calc_Discount_Factor(RateTerm, Rate, NRateTerm, T);

            if (dt < 1.0) ResultRate = (1.0 - PtT) / (dt * PtT);
            else ResultRate = -1.0 / dt * log(PtT);
        }
    }
    return ResultRate;
}

long Count_Accrual_Fast(
    long NDays,
    double Rate1,
    double Rate2,
    double LowerBound,
    double UpperBound
)
{
    long nAccrual;
    double Ratio = 0.0;
    double a, b;
    double r_min = min(Rate1, Rate2);
    double r_max = max(Rate1, Rate2);
    if (LowerBound <= r_min && UpperBound >= r_min && LowerBound <= r_max && UpperBound >= r_max) nAccrual = NDays; // 둘 다 Range 만족한 경우
    else if ((LowerBound > r_max) || (UpperBound < r_min)) nAccrual = 0;   // 둘 다 Upper 위나 Lower 아래에 있는 경우
    else if ((LowerBound < r_min) && (UpperBound < r_max))
    {
        a = UpperBound - r_min;
        b = r_max - r_min;
        Ratio = a / b;
        nAccrual = (long)(Ratio * (double)NDays);
    }
    else if ((LowerBound > r_min) && (UpperBound > r_max))
    {
        a = UpperBound - r_max;
        b = r_max - r_min;
        Ratio = a / b;
        nAccrual = (long)(Ratio * (double)NDays);
    }
    else
    {
        Ratio = ((UpperBound - LowerBound) / (r_max - r_min));
        nAccrual = (long)(Ratio * (double)NDays);
    }
    return nAccrual;
}

double PayoffStructure(
    long DailySimulFlag,
    LEG_INFO* Leg_Inform,
    double** SimulatedRate,  // Reset Start Rate
    double** SimulatedRate2, // Reset End Rate
    double** DailyRate,
    double** SimulOverNightRate,
    long CashFlowIdx,
    long Flag_Required_History,
    double* OutputRate,
    long& NumAccrual
)
{
    long i;
    long j;

    long N;
    long CalcFastFlag = 1;

    double s;
    double dt;
    double Rate, Rate1, Rate2;
    long DayBefore;
    long DayBeforeIdx;
    long ExistHistoryFlag = 0;

    double LowerBound;
    double UpperBound;
    double CondMultiple;
    double PayoffMultiple;
    double Payoff;
    double Alternative_Rate;
    double CPN_Ratio = 1.0;
    long Cond = 1;
    long NDays;
    NumAccrual = 0;
    if (DailySimulFlag == 0 && Leg_Inform->AccrualFlag == 1 && Flag_Required_History == 0)
    {
        NDays = Leg_Inform->DaysForwardEnd[CashFlowIdx] - Leg_Inform->DaysForwardStart[CashFlowIdx];
        N = NDays;
        for (i = 0; i < Leg_Inform->NReference; i++)
        {
            if (Leg_Inform->Reference_Inform[i].RefRateType == 1)
            {
                Rate1 = SimulatedRate[i][CashFlowIdx];
                Rate2 = SimulatedRate2[i][CashFlowIdx];
                LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                N = min(N, Count_Accrual_Fast(NDays, Rate1, Rate2, LowerBound, UpperBound));
            }
        }
        NumAccrual = N;
        CPN_Ratio = ((double)NumAccrual) / ((double)NDays);
    }
    else if (DailySimulFlag == 1 && Leg_Inform->AccrualFlag == 1 && Flag_Required_History == 0)
    {
        NDays = Leg_Inform->DaysForwardEnd[CashFlowIdx] - Leg_Inform->DaysForwardStart[CashFlowIdx];
        NumAccrual = NDays;
        for (i = 0; i < Leg_Inform->NReference; i++)
        {
            if (Leg_Inform->Reference_Inform[i].RefRateType == 1)
            {
                N = 0;
                for (j = Leg_Inform->DaysForwardStart[CashFlowIdx]; j < Leg_Inform->DaysForwardEnd[CashFlowIdx]; j++)
                {
                    Rate = DailyRate[i][j];
                    LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                    UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                    if (Rate >= LowerBound && Rate <= UpperBound) N += 1;
                }
                NumAccrual = min(NumAccrual, N);
            }
        }
        CPN_Ratio = ((double)NumAccrual) / ((double)NDays);
    }

    Payoff = 0.0;
    dt = ((double)DayCountAtoB(Leg_Inform->FractionStart_C[CashFlowIdx], Leg_Inform->FractionEnd_C[CashFlowIdx])) / 365.0;


    if (Leg_Inform->PayoffStructure == 0)
    {
        if (Flag_Required_History == 0)
        {
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                PayoffMultiple = Leg_Inform->Reference_Inform[i].PayoffMultiple;
                OutputRate[i] = SimulatedRate[i][CashFlowIdx];
                if (LowerBound <= CondMultiple * SimulatedRate[i][CashFlowIdx] && UpperBound >= CondMultiple * SimulatedRate[i][CashFlowIdx])
                {
                    Cond = Cond * 1;
                }
                else
                {
                    Cond = 0;
                }
                Payoff += PayoffMultiple * SimulatedRate[i][CashFlowIdx];
            }
            Payoff = Payoff * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;

        }
        else
        {
            Alternative_Rate = 0.0;
            DayBefore = Leg_Inform->DaysForwardStart[CashFlowIdx];

            ExistHistoryFlag = 1;
            Payoff = 0.0;
            NDays = Leg_Inform->DaysForwardEnd[CashFlowIdx] - Leg_Inform->DaysForwardStart[CashFlowIdx];
            N = NDays;
            NumAccrual = NDays;
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                ExistHistoryFlag = isinFindIndex(DayBefore, Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                if (ExistHistoryFlag == 0)
                {
                    Payoff = 0.0;
                    NumAccrual = 0;
                    break;
                }
                else
                {
                    LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                    UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                    CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                    PayoffMultiple = Leg_Inform->Reference_Inform[i].PayoffMultiple;
                    Alternative_Rate = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                    if (Leg_Inform->AccrualFlag == 1)
                    {
                        Rate2 = 0.0;
                        N = 0;
                        for (j = Leg_Inform->DaysForwardStart[CashFlowIdx]; j < Leg_Inform->DaysForwardEnd[CashFlowIdx]; j++)
                        {
                            ExistHistoryFlag = isinFindIndex(j, Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                            if (ExistHistoryFlag != 0 && j < 0)
                            {
                                Rate2 = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                                if (Rate2 >= LowerBound && Rate2 <= UpperBound) N += 1;
                            }
                            else
                            {
                                Rate2 = SimulatedRate2[i][CashFlowIdx];
                                if (Rate2 >= LowerBound && Rate2 <= UpperBound) N += 1;
                            }
                        }
                        NumAccrual = min(NumAccrual, N);

                    }

                    OutputRate[i] = Alternative_Rate;
                    if (LowerBound <= CondMultiple * Alternative_Rate && UpperBound >= CondMultiple * Alternative_Rate)
                    {
                        Cond = Cond * 1;
                    }
                    else
                    {
                        Cond = 0;
                    }

                    Payoff += PayoffMultiple * Alternative_Rate;
                }
            }

            CPN_Ratio = ((double)NumAccrual) / ((double)NDays);
            Payoff = Payoff * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;

        }
    }
    else
    {
        if (Flag_Required_History == 0)
        {
            s = 0.0;
            LowerBound = Leg_Inform->Reference_Inform[0].RangeDn;
            UpperBound = Leg_Inform->Reference_Inform[0].RangeUp;

            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                LowerBound = max(LowerBound, Leg_Inform->Reference_Inform[i].RangeDn);
                UpperBound = min(UpperBound, Leg_Inform->Reference_Inform[i].RangeUp);
                CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                OutputRate[i] = SimulatedRate[i][CashFlowIdx];
                s += CondMultiple * SimulatedRate[i][CashFlowIdx];
            }

            for (i = 0; i < Leg_Inform->NReference; i++)
            {

                if (LowerBound <= s && UpperBound >= s)
                {
                    Cond = Cond * 1;
                }
                else
                {
                    Cond = 0;
                }
            }

            s = 0.0;
            PayoffMultiple = Leg_Inform->Reference_Inform[0].PayoffMultiple;
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                PayoffMultiple = min(PayoffMultiple, Leg_Inform->Reference_Inform[i].PayoffMultiple);
                s += PayoffMultiple * SimulatedRate[i][CashFlowIdx];
            }
            Payoff = s * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;

        }
        else
        {
            Alternative_Rate = 0.0;
            s = 0.0;
            DayBefore = Leg_Inform->DaysForwardStart[CashFlowIdx];

            ExistHistoryFlag = 1;
            Payoff = 0.0;
            LowerBound = Leg_Inform->Reference_Inform[0].RangeDn;
            UpperBound = Leg_Inform->Reference_Inform[0].RangeUp;
            PayoffMultiple = Leg_Inform->Reference_Inform[0].PayoffMultiple;

            NDays = Leg_Inform->DaysForwardEnd[CashFlowIdx] - Leg_Inform->DaysForwardStart[CashFlowIdx];
            N = NDays;
            NumAccrual = NDays;
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                LowerBound = max(LowerBound, Leg_Inform->Reference_Inform[i].RangeDn);
                UpperBound = min(UpperBound, Leg_Inform->Reference_Inform[i].RangeUp);
                PayoffMultiple = min(PayoffMultiple, Leg_Inform->Reference_Inform[i].PayoffMultiple);
                ExistHistoryFlag = isinFindIndex(DayBefore, Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                if (ExistHistoryFlag == 0)
                {
                    Payoff = 0.0;
                    NumAccrual = 0;
                    break;
                }
                else
                {
                    CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                    Alternative_Rate += CondMultiple * Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                    if (Leg_Inform->AccrualFlag == 1)
                    {
                        Rate2 = 0.0;
                        N = 0;
                        for (j = Leg_Inform->DaysForwardStart[CashFlowIdx]; j < Leg_Inform->DaysForwardEnd[CashFlowIdx]; j++)
                        {
                            ExistHistoryFlag = isinFindIndex(j, Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                            if (ExistHistoryFlag != 0 && j < 0)
                            {
                                Rate2 = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                                if (Rate2 >= LowerBound && Rate2 <= UpperBound) N += 1;
                            }
                            else
                            {
                                Rate2 = SimulatedRate2[i][CashFlowIdx];
                                if (Rate2 >= LowerBound && Rate2 <= UpperBound) N += 1;
                            }
                        }
                        NumAccrual = min(NumAccrual, N);
                    }
                    OutputRate[i] = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                    s += PayoffMultiple * Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                }
            }

            for (i = 0; i < Leg_Inform->NReference; i++)
            {

                if (LowerBound <= Alternative_Rate && UpperBound >= Alternative_Rate)
                {
                    Cond = Cond * 1;
                }
                else
                {
                    Cond = 0;
                }
            }

            Payoff = s * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;
        }
    }
    return dt * Payoff;
}

curveinfo* CurveUpDn_Ref(
    curveinfo* rf_curves,
    long n_curve,
    long UpFlag
)
{
    long i;
    long j;
    long n;
    double bp;
    if (UpFlag == 1) bp = 0.0001;
    else bp = -0.0001;

    double* TempRateTerm;
    double* TempRate;
    curveinfo* result_rf_curves = new curveinfo[n_curve];
    for (i = 0; i < n_curve; i++)
    {
        n = (rf_curves + i)->nterm();
        TempRateTerm = (rf_curves + i)->Term;
        TempRate = (double*)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++)
        {
            TempRate[j] = (rf_curves + i)->Rate[j] + bp;
        }
        (result_rf_curves + i)->initialize(n, TempRateTerm, TempRate);
        free(TempRate);
    }

    return result_rf_curves;
}

double SOFR_ForwardRate_Compound_DiscreteSimul(
    long DailySimulFlag,
    long NShortRate,
    double* TShortRate,
    double* ShortRate,
    long Ref_Days_0,
    long Ref_Days_1,
    long FixedRateFlag,
    double FixedRate,
    long HolidayCalcType,
    long NHolidays,
    long* Days_Holidays,
    double* PV_t_T,
    double* QVTerm,
    double* B_t_T,
    double& AnnualizedOISRate
)
{
    long i;
    long j;

    long HolidayFlag = 0;
    long HolidayFlag2 = 0;
    long nHoliday = 0;

    long idx;
    long maxidx;

    double PtT;
    double dt = 1.0 / 365.0;
    double Prev_PI;
    double TodayRate;
    double ForwardRate;
    double PrevRate2;

    double PI_0;
    double t;
    double T;
    double xt;

    if (Ref_Days_0 < 0)
    {
        if (FixedRateFlag == 1)
        {
            if (FixedRate > 0.0 || FixedRate < -0.00002) // 입력된 FixedRate이 존재하는경우
            {
                Prev_PI = pow((1.0 + FixedRate), (-(double)Ref_Days_0) / 365.0);
            }
            else
            {
                t = dt;
                xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                PtT = PV_t_T[0] * exp(-xt * B_t_T[0] + QVTerm[0]);
                ForwardRate = (1.0 - PtT) / (dt * PtT);
                Prev_PI = pow((1.0 + ForwardRate), (-(double)Ref_Days_0) / 365.0);
            }
        }
        else // Flag가 없어도 과거며칠의 금리 대신 가장 유사한 금리 계산
        {
            t = dt;
            xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
            PtT = PV_t_T[0] * exp(-xt * B_t_T[0] + QVTerm[0]);
            ForwardRate = (1.0 - PtT) / (dt * PtT);
            Prev_PI = exp(ForwardRate * (-(double)Ref_Days_0) / 365.0);
        }
        // 여기까지 과거금리 Compound
        xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, 1.0 / 365.0);
        PtT = PV_t_T[0] * exp(-xt * B_t_T[0] + QVTerm[0]);
        TodayRate = (1.0 - PtT) / (dt * PtT);

        // 오늘금리 Compound
        PI_0 = Prev_PI * (1.0 + TodayRate * dt);

        // 내일부터 금리 Compound
        if (HolidayCalcType == 0) // Forward Fill NA
        {
            maxidx = Ref_Days_1 - 1 - 1;
            for (i = 1; i < Ref_Days_1; i++)
            {
                idx = i - 1;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == 1) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                    ForwardRate = (1.0 - PtT) / (dt * PtT);
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
            maxidx = Ref_Days_1 - 1 - 1;
            for (i = 1; i < Ref_Days_1; i++)
            {
                idx = i - 1;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == 1) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    nHoliday = 0;
                    for (j = i + 1; j < Ref_Days_1; j++)
                    {
                        HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
                        if (HolidayFlag2 == 0) break;
                        else nHoliday += 1;
                    }
                    if (nHoliday == 0)
                    {
                        PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                        ForwardRate = (1.0 - PtT) / (dt * PtT);
                    }
                    else
                    {
                        xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                        PtT = PV_t_T[min(maxidx, idx + nHoliday + 1)] * exp(-xt * B_t_T[min(maxidx, idx + nHoliday + 1)] + QVTerm[min(maxidx, idx + nHoliday + 1)]);
                        ForwardRate = (1.0 - PtT) / (dt * PtT);
                    }


                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
                }
            }
        }
        else if (HolidayCalcType == 2)
        {
            maxidx = Ref_Days_1 - 1 - 1;
            for (i = 1; i < Ref_Days_1; i++)
            {
                idx = i - 1;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == 1) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                    ForwardRate = (1.0 - PtT) / (dt * PtT);
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
                    }
                    else
                    {
                        xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                        PtT = PV_t_T[min(maxidx, idx + nHoliday + 1)] * exp(-xt * B_t_T[min(maxidx, idx + nHoliday + 1)] + QVTerm[min(maxidx, idx + nHoliday + 1)]);
                        PrevRate2 = (1.0 - PtT) / (dt * PtT);
                        PI_0 *= (1.0 + (0.5 * ForwardRate + 0.5 * PrevRate2) * (1.0 + (double)nHoliday) * dt);
                    }
                }
            }
        }
        else
        {

            for (i = 1; i < Ref_Days_1; i++)
            {
                idx = i - 1;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                ForwardRate = (1.0 - PtT) / (dt * PtT);
                PI_0 *= (1.0 + ForwardRate * dt);
            }

        }
    }
    else
    {
        PI_0 = 1.0;
        if (HolidayCalcType == 0)
        {
            maxidx = Ref_Days_1 - Ref_Days_0 - 1;
            for (i = Ref_Days_0; i < Ref_Days_1; i++)
            {
                idx = i - Ref_Days_0;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == Ref_Days_0) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                    ForwardRate = (1.0 - PtT) / (dt * PtT);
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
            maxidx = Ref_Days_1 - Ref_Days_0 - 1;
            for (i = Ref_Days_0; i < Ref_Days_1; i++)
            {
                idx = i - Ref_Days_0;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == Ref_Days_0) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    nHoliday = 0;
                    for (j = i + 1; j < Ref_Days_1; j++)
                    {
                        HolidayFlag2 = isin(j, Days_Holidays, NHolidays);
                        if (HolidayFlag2 == 0) break;
                        else nHoliday += 1;
                    }
                    if (nHoliday == 0)
                    {
                        PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                        ForwardRate = (1.0 - PtT) / (dt * PtT);
                    }
                    else
                    {
                        xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                        PtT = PV_t_T[min(maxidx, idx + nHoliday + 1)] * exp(-xt * B_t_T[min(maxidx, idx + nHoliday + 1)] + QVTerm[min(maxidx, idx + nHoliday + 1)]);
                        ForwardRate = (1.0 - PtT) / (dt * PtT);
                    }

                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)nHoliday) * dt);
                }
            }
        }
        else if (HolidayCalcType == 2)
        {
            maxidx = Ref_Days_1 - Ref_Days_0 - 1;
            for (i = Ref_Days_0; i < Ref_Days_1; i++)
            {
                idx = i - Ref_Days_0;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                if (i == Ref_Days_0) HolidayFlag = 0;
                else HolidayFlag = isin(i, Days_Holidays, NHolidays);

                if (HolidayFlag == 0)
                {
                    PtT = PV_t_T[i] * exp(-xt * B_t_T[i] + QVTerm[i]);
                    ForwardRate = (1.0 - PtT) / (dt * PtT);
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
                    }
                    else
                    {
                        xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                        PtT = PV_t_T[min(maxidx, idx + nHoliday + 1)] * exp(-xt * B_t_T[min(maxidx, idx + nHoliday + 1)] + QVTerm[min(maxidx, idx + nHoliday + 1)]);
                        PrevRate2 = (1.0 - PtT) / (dt * PtT);
                        PI_0 *= (1.0 + (0.5 * ForwardRate + 0.5 * PrevRate2) * (1.0 + (double)nHoliday) * dt);
                    }
                }
            }
        }
        else
        {
            for (i = Ref_Days_0; i < Ref_Days_1; i++)
            {
                idx = i - 1;

                t = ((double)i) / 365.0;
                if (DailySimulFlag == 0) xt = Interpolate_Linear(TShortRate, ShortRate, NShortRate, t);
                else xt = ShortRate[i];

                HolidayFlag = isin(i, Days_Holidays, NHolidays);

                PtT = PV_t_T[idx] * exp(-xt * B_t_T[idx] + QVTerm[idx]);
                ForwardRate = (1.0 - PtT) / (dt * PtT);
                PI_0 *= (1.0 + ForwardRate * dt);
            }
        }
    }
    T = (double)(Ref_Days_1 - Ref_Days_0) * dt;
    AnnualizedOISRate = (PI_0 - 1.0) * 1.0 / T;
    return PI_0 - 1.0;
}


long Simulate_HW(
    long PricingOnly,
    long PricingDateC,
    long NAFlag,
    double Notional,
    LEG_INFO* RcvLeg,
    LEG_INFO* PayLeg,
    SIMUL_INFO* Simul,
    long Length_XA_Array,
    double** XA_Array,
    long Length_XV_Array,
    double** XV_Array,
    double** SimulatedRateRcv,
    double** SimulatedRateRcv2,
    double** SimulatedRatePay,
    double** SimulatedRatePay2,
    double** SimulShortRate,
    double** SimulOverNightRate,
    long* CurveIdx_Rcv,
    long* CurveIdx_Pay,
    double* ResultPrice,                // 
    double* ResultRcv,                  // 
    double* ResultPay                   // 
)
{
    long i;
    long j;
    long k;

    long n;
    long nAccrual = 0;
    long ndates;

    double Rcv_DF_Day1Y;
    double Pay_DF_Day1Y;
    if (RcvLeg->DiscDayCountFlag == 0) Rcv_DF_Day1Y = 365.0;
    else Rcv_DF_Day1Y = 360.0;

    if (PayLeg->DiscDayCountFlag == 0) Pay_DF_Day1Y = 365.0;
    else Pay_DF_Day1Y = 360.0;

    long Curveidx;

    long xt_idx;
    long Day1, Day2;
    long Dayidx;
    long RateHistoryUseFlag = 0;
    double PricePath_Rcv, PricePath_Pay;
    double RcvPrice, PayPrice;
    double t, T, T1, T2, Pt, PT, PtT;
    double kappa;
    double xt = 0.0;
    double Rate = 0.0, Rate1 = 0.0, Rate2 = 0.0;


    double** RcvDailyRate = (double**)malloc(sizeof(double*) * RcvLeg->NReference);
    double** PayDailyRate = (double**)malloc(sizeof(double*) * PayLeg->NReference);
    if (Simul->DailySimulFlag == 1)
    {
        for (i = 0; i < RcvLeg->NReference; i++) RcvDailyRate[i] = (double*)malloc(sizeof(double) * (Simul->DaysForSimul[Simul->NDays - 1] + 1));
        for (i = 0; i < PayLeg->NReference; i++) PayDailyRate[i] = (double*)malloc(sizeof(double) * (Simul->DaysForSimul[Simul->NDays - 1] + 1));
    }

    double* RcvPayoff = (double*)calloc((RcvLeg->NCashFlow + 1), sizeof(double));
    double* PayPayoff = (double*)calloc((PayLeg->NCashFlow + 1), sizeof(double));
    long* RcvCFFlag = (long*)calloc(RcvLeg->NCashFlow, sizeof(long));
    long* PayCFFlag = (long*)calloc(RcvLeg->NCashFlow, sizeof(long));
    long CurveIdx_DiscRcv = FindIndex(RcvLeg->DiscCurveNum, Simul->SimulCurveIdx, Simul->NAsset);
    long CurveIdx_DiscPay = FindIndex(PayLeg->DiscCurveNum, Simul->SimulCurveIdx, Simul->NAsset);
    double* Rcv_DF_0_t = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    double* Rcv_DF_0_T = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    for (i = 0; i < RcvLeg->NCashFlow * 6; i++) ResultRcv[i] = 0.0;
    for (i = 0; i < PayLeg->NCashFlow * 6; i++) ResultPay[i] = 0.0;

    for (i = 0; i < RcvLeg->NCashFlow; i++)
    {
        t = ((double)RcvLeg->DaysForwardStart[i]) / Rcv_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)RcvLeg->DaysPayDate_C[i]) / Rcv_DF_Day1Y;
        Rcv_DF_0_t[i] = Calc_Discount_Factor(Simul->RateTerm[CurveIdx_DiscRcv], Simul->Rate[CurveIdx_DiscRcv], Simul->NRateTerm[CurveIdx_DiscRcv], t);
        Rcv_DF_0_T[i] = Calc_Discount_Factor(Simul->RateTerm[CurveIdx_DiscRcv], Simul->Rate[CurveIdx_DiscRcv], Simul->NRateTerm[CurveIdx_DiscRcv], T);
    }
    double* Pay_DF_0_t = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    double* Pay_DF_0_T = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    for (i = 0; i < PayLeg->NCashFlow; i++)
    {
        t = ((double)PayLeg->DaysForwardStart[i]) / Pay_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)PayLeg->DaysPayDate_C[i]) / Pay_DF_Day1Y;
        Pay_DF_0_t[i] = Calc_Discount_Factor(Simul->RateTerm[CurveIdx_DiscPay], Simul->Rate[CurveIdx_DiscPay], Simul->NRateTerm[CurveIdx_DiscPay], t);
        Pay_DF_0_T[i] = Calc_Discount_Factor(Simul->RateTerm[CurveIdx_DiscPay], Simul->Rate[CurveIdx_DiscPay], Simul->NRateTerm[CurveIdx_DiscPay], T);
    }

    double* B_t_T_RcvDisc = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    double* B_t_T_PayDisc = (double*)calloc(PayLeg->NCashFlow, sizeof(double));
    double* QVTerm_RcvDisc = (double*)calloc(RcvLeg->NCashFlow, sizeof(double));
    double* QVTerm_PayDisc = (double*)calloc(PayLeg->NCashFlow, sizeof(double));

    for (i = 0; i < RcvLeg->NCashFlow; i++)
    {
        t = ((double)RcvLeg->DaysForwardStart[i]) / Rcv_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)RcvLeg->DaysPayDate_C[i]) / Rcv_DF_Day1Y;
        B_t_T_RcvDisc[i] = B_s_to_t(Simul->HWKappa[CurveIdx_DiscRcv], t, T);
        QVTerm_RcvDisc[i] = HullWhiteQVTerm(t, T, Simul->HWKappa[CurveIdx_DiscRcv], Simul->NHWVol[CurveIdx_DiscRcv], Simul->HWVolTerm[CurveIdx_DiscRcv], Simul->HWVol[CurveIdx_DiscRcv]);
    }

    for (i = 0; i < PayLeg->NCashFlow; i++)
    {
        t = ((double)PayLeg->DaysForwardStart[i]) / Pay_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)PayLeg->DaysPayDate_C[i]) / Pay_DF_Day1Y;
        B_t_T_PayDisc[i] = B_s_to_t(Simul->HWKappa[CurveIdx_DiscPay], t, T);
        QVTerm_PayDisc[i] = HullWhiteQVTerm(t, T, Simul->HWKappa[CurveIdx_DiscPay], Simul->NHWVol[CurveIdx_DiscPay], Simul->HWVolTerm[CurveIdx_DiscPay], Simul->HWVol[CurveIdx_DiscPay]);
    }

    double* RcvOutputRate = (double*)calloc(RcvLeg->NReference, sizeof(double));
    double* PayOutputRate = (double*)calloc(PayLeg->NReference, sizeof(double));

    double*** PV_t_T_Rcv = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    double*** QVTerm_Rcv = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    double*** B_t_T_Rcv = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    long* ndates_Rcv = (long*)malloc(sizeof(long) * RcvLeg->NReference);

    double*** PV_t_T_Rcv_PowerSpread = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    double*** QVTerm_Rcv_PowerSpread = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    double*** B_t_T_Rcv_PowerSpread = (double***)malloc(sizeof(double**) * RcvLeg->NReference);
    long* ndates_Rcv_PowerSpread = (long*)malloc(sizeof(long) * RcvLeg->NReference);

    for (i = 0; i < RcvLeg->NReference; i++)
    {
        Curveidx = CurveIdx_Rcv[i];
        PV_t_T_Rcv[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        QVTerm_Rcv[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        B_t_T_Rcv[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

        if (RcvLeg->Reference_Inform[i].RefRateType == 0)
        {
            if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                ndates = 2;
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }

                ndates_Rcv[i] = ndates;
                Set_HW_Params_RefType0(Simul, RcvLeg, i, Curveidx, PV_t_T_Rcv[i], QVTerm_Rcv[i],
                    B_t_T_Rcv[i], PV_t_T_Rcv[i], QVTerm_Rcv[i], B_t_T_Rcv[i]);
            }
            else
            {
                ndates = 2;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Rcv[i] = ndates;

                PV_t_T_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                QVTerm_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                B_t_T_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                ndates = 2;
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Rcv_PowerSpread[i] = ndates;

                Set_HW_Params_RefType0(Simul, RcvLeg, i, Curveidx, PV_t_T_Rcv[i], QVTerm_Rcv[i],
                    B_t_T_Rcv[i], PV_t_T_Rcv_PowerSpread[i], QVTerm_Rcv_PowerSpread[i], B_t_T_Rcv_PowerSpread[i]);
            }
        }
        else if (RcvLeg->Reference_Inform[i].RefRateType == 1)
        {
            if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].RefSwapMaturity, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }

                ndates_Rcv[i] = ndates;
                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann, RcvLeg->Reference_Inform[i].RefSwapMaturity, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv[i], QVTerm_Rcv[i], B_t_T_Rcv[i]);
            }
            else
            {
                if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].PowerSpreadMat1, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Rcv[i] = ndates;

                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann, RcvLeg->Reference_Inform[i].PowerSpreadMat1, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv[i], QVTerm_Rcv[i], B_t_T_Rcv[i]);

                PV_t_T_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                QVTerm_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                B_t_T_Rcv_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

                if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].PowerSpreadMat2, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Rcv_PowerSpread[i] = ndates;

                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann, RcvLeg->Reference_Inform[i].PowerSpreadMat2, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv_PowerSpread[i], QVTerm_Rcv_PowerSpread[i], B_t_T_Rcv_PowerSpread[i]);
            }
        }
        else if (RcvLeg->Reference_Inform[i].RefRateType == 2)
        {
            if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                ndates = 0;
                for (j = 0; j < RcvLeg->NCashFlow; j++) ndates = max(ndates, RcvLeg->DaysForwardEnd[j] - RcvLeg->DaysForwardStart[j]);

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv[i], QVTerm_Rcv[i], B_t_T_Rcv[i]);
            }
            else
            {
                ndates = (long)(RcvLeg->Reference_Inform[i].PowerSpreadMat1 * 365.0 + 0.5);
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv[i], QVTerm_Rcv[i], B_t_T_Rcv[i]);

                ndates = (long)(RcvLeg->Reference_Inform[i].PowerSpreadMat2 * 365.0 + 0.5);
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Rcv[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Rcv_PowerSpread[i], QVTerm_Rcv_PowerSpread[i], B_t_T_Rcv_PowerSpread[i]);
            }
        }
    }

    double*** PV_t_T_Pay = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    double*** QVTerm_Pay = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    double*** B_t_T_Pay = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    long* ndates_Pay = (long*)malloc(sizeof(long) * PayLeg->NReference);

    double*** PV_t_T_Pay_PowerSpread = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    double*** QVTerm_Pay_PowerSpread = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    double*** B_t_T_Pay_PowerSpread = (double***)malloc(sizeof(double**) * PayLeg->NReference);
    long* ndates_Pay_PowerSpread = (long*)malloc(sizeof(long) * PayLeg->NReference);
    for (i = 0; i < PayLeg->NReference; i++)
    {
        Curveidx = CurveIdx_Pay[i];
        PV_t_T_Pay[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        QVTerm_Pay[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        B_t_T_Pay[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

        if (PayLeg->Reference_Inform[i].RefRateType == 0)
        {
            if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                ndates = 2;
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }

                ndates_Pay[i] = ndates;
                Set_HW_Params_RefType0(Simul, PayLeg, i, Curveidx, PV_t_T_Pay[i], QVTerm_Pay[i],
                    B_t_T_Pay[i], PV_t_T_Pay[i], QVTerm_Pay[i], B_t_T_Pay[i]);
            }
            else
            {
                ndates = 2;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Pay[i] = ndates;

                PV_t_T_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                QVTerm_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                B_t_T_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                ndates = 2;
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Pay_PowerSpread[i] = ndates;

                Set_HW_Params_RefType0(Simul, PayLeg, i, Curveidx, PV_t_T_Pay[i], QVTerm_Pay[i],
                    B_t_T_Pay[i], PV_t_T_Pay_PowerSpread[i], QVTerm_Pay_PowerSpread[i], B_t_T_Pay_PowerSpread[i]);
            }
        }
        else if (PayLeg->Reference_Inform[i].RefRateType == 1)
        {
            if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].RefSwapMaturity, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }

                ndates_Pay[i] = ndates;
                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, PayLeg->Reference_Inform[i].RefSwapNCPN_Ann, PayLeg->Reference_Inform[i].RefSwapMaturity, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay[i], QVTerm_Pay[i], B_t_T_Pay[i]);
            }
            else
            {
                if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].PowerSpreadMat1, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Pay[i] = ndates;

                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, PayLeg->Reference_Inform[i].RefSwapNCPN_Ann, PayLeg->Reference_Inform[i].PowerSpreadMat1, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay[i], QVTerm_Pay[i], B_t_T_Pay[i]);

                PV_t_T_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                QVTerm_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
                B_t_T_Pay_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

                if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) ndates = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].PowerSpreadMat2, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                else ndates = 1;

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                ndates_Pay_PowerSpread[i] = ndates;

                Set_HW_Params_RefType1(Simul->NDays, Simul->T_Array, ndates, PayLeg->Reference_Inform[i].RefSwapNCPN_Ann, PayLeg->Reference_Inform[i].PowerSpreadMat2, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay_PowerSpread[i], QVTerm_Pay_PowerSpread[i], B_t_T_Pay_PowerSpread[i]);
            }
        }
        else if (PayLeg->Reference_Inform[i].RefRateType == 2)
        {
            if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
            {
                ndates = 0;
                for (j = 0; j < PayLeg->NCashFlow; j++) ndates = max(ndates, PayLeg->DaysForwardEnd[j] - PayLeg->DaysForwardStart[j]);

                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay[i], QVTerm_Pay[i], B_t_T_Pay[i]);
            }
            else
            {
                ndates = (long)(PayLeg->Reference_Inform[i].PowerSpreadMat1 * 365.0 + 0.5);
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay[i], QVTerm_Pay[i], B_t_T_Pay[i]);

                ndates = (long)(PayLeg->Reference_Inform[i].PowerSpreadMat2 * 365.0 + 0.5);
                for (j = 0; j < Simul->NDays; j++)
                {
                    PV_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    QVTerm_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                    B_t_T_Pay[i][j] = (double*)malloc(sizeof(double) * ndates);
                }
                Set_HW_Params_RefType2(Simul->NDays, ndates, Simul->T_Array, Simul->HWKappa[Curveidx],
                    Simul->NHWVol[Curveidx], Simul->HWVolTerm[Curveidx], Simul->HWVol[Curveidx], Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                    Simul->Rate[Curveidx], PV_t_T_Pay_PowerSpread[i], QVTerm_Pay_PowerSpread[i], B_t_T_Pay_PowerSpread[i]);
            }
        }
    }

    RcvPrice = 0.0;
    PayPrice = 0.0;
    double SOFRDailyCompound = 0.0;

    ////////////////
    // LSMC Data  //
    ////////////////

    double*** X;
    long ShapeX[2] = { Simul->NSimul, Simul->NAsset + 1};
    double** Y;
    long** OptRangeFlag;
    long LengthY = Simul->NSimul;
    long idxoption;
    long idxstart;
    double T_Opt;
    double** Y_Hat;
    double* Beta;

    double* Value_By_OptTime;
    double* Estimated_Value_By_OptTime;
    long* OptimalIdx;
    double** InterestRate_Opt;
    double EstOptValue;
    if (RcvLeg->OptionUseFlag == 1)
    {
        Value_By_OptTime = (double*)malloc(sizeof(double) * Simul->NSimul);
        Estimated_Value_By_OptTime = (double*)malloc(sizeof(double*) * Simul->NSimul);
        OptimalIdx = (long*)malloc(sizeof(double) * Simul->NSimul);
        for (i = 0; i < Simul->NSimul; i++) OptimalIdx[i] = 0;

        InterestRate_Opt = (double**)malloc(sizeof(double**) * Simul->NSimul);
        for (i = 0; i < Simul->NSimul; i++) InterestRate_Opt[i] = (double*)malloc(sizeof(double) * (Simul->NAsset + 1) );
        for (i = 0; i < Simul->NSimul; i++) InterestRate_Opt[i][Simul->NAsset] = 1.0;

        X = (double***)malloc(sizeof(double**) * RcvLeg->NOption);
        for (i = 0; i < RcvLeg->NOption; i++) X[i] = make_array(Simul->NSimul, Simul->NAsset + 1);

        Y = (double**)malloc(sizeof(double*) * RcvLeg->NOption);
        for (i = 0; i < RcvLeg->NOption; i++) Y[i] = make_array(Simul->NSimul);

        OptRangeFlag = (long**)malloc(sizeof(long*) * RcvLeg->NOption);
        for (i = 0; i < RcvLeg->NOption; i++) OptRangeFlag[i] = (long*)malloc(sizeof(long) * Simul->NSimul);

        Y_Hat = (double**)malloc(sizeof(double*) * RcvLeg->NOption);
        for (i = 0; i < RcvLeg->NOption; i++) Y_Hat[i] = make_array(Simul->NSimul);

        for (i = 0; i < RcvLeg->NOption; i++)
            for (j = 0; j < Simul->NSimul; j++) X[i][j][Simul->NAsset] = 1.0;


    }
    else
    {
        RcvLeg->NOption = 0;
        X = (double***)malloc(sizeof(double**) * 1);
        Y = (double**)malloc(sizeof(double*) * 1);
        Value_By_OptTime = (double*)malloc(sizeof(double) * 1);
        Estimated_Value_By_OptTime = (double*)malloc(sizeof(double) * 1);
        OptimalIdx = (long*)malloc(sizeof(long) * 1);
        InterestRate_Opt = (double**)malloc(sizeof(double*) * 1);
        OptRangeFlag = (long**)malloc(sizeof(long*) * 1);
        Y_Hat = (double**)malloc(sizeof(double*) * 1);

    }



    long lastresetend = -1;
    for (i = 0; i < Simul->NSimul; i++)
    {
        PricePath_Rcv = 0.0;
        PricePath_Pay = 0.0;

        for (j = 0; j < Simul->NDays; j++)
        {
            for (n = 0; n < Simul->NAsset; n++)
            {
                if (j == 0) SimulShortRate[n][j] = 0.0;
                else  SimulShortRate[n][j] = Simulated_ShortRate(XA_Array[n][j], XV_Array[n][j], SimulShortRate[n][j - 1], Simul->FixedRandn[i * Simul->NDays + j][n]);
            }
        }

        if (RcvLeg->OptionUseFlag == 1)
        {
            idxstart = 0;
            for (idxoption = 0; idxoption < RcvLeg->NOption; idxoption++)
            {
                T_Opt = ((double)RcvLeg->DaysOptionDate[idxoption]) / 365.0;
                for (j = idxstart; j < Simul->NDays; j++)
                {
                    if (Simul->T_Array[j] >= T_Opt)
                    {
                        for (n = 0; n < Simul->NAsset; n++)
                        {
                            X[idxoption][i][n] = SimulShortRate[n][j];
                        }
                        idxstart = j;
                        break;
                    }
                }
            }
        }

        // RcvLeg Reference Rate 결정
        for (n = 0; n < RcvLeg->NReference; n++)
        {
            Curveidx = CurveIdx_Rcv[n];
            lastresetend = -1;
            Rate = 0.0;
            Rate1 = 0.0;
            Rate2 = 0.0;

            if (Simul->DailySimulFlag == 0)
            {
                //Libor, Zero의 경우 
                if (RcvLeg->Reference_Inform[n].RefRateType == 0)
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;

                                Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);

                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);

                                Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv_PowerSpread[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv_PowerSpread[n][j],
                                    QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv[n][Dayidx] = Rate;
                        }
                    }
                }
                else if (RcvLeg->Reference_Inform[n].RefRateType == 1) // Libor Swap
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = t + RcvLeg->Reference_Inform[n].RefSwapMaturity;

                                if (Simul->DaysForSimul[j] != lastresetend) // 이전 Reset End와 Reset Start가 같으면 중복계산하지않는다.(연산속도)
                                {
                                    Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                        QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);
                                }

                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                                if (Simul->DaysForSimul[j] != lastresetend)
                                {
                                    Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                        QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);

                                    Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Rcv_PowerSpread[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv_PowerSpread[n][j],
                                        QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j]);
                                }

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv[n][Dayidx] = Rate;

                        }
                        if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardEnd, RcvLeg->NCashFlow, Dayidx))
                        {
                            lastresetend = Simul->DaysForSimul[j];
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardEnd[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = t + RcvLeg->Reference_Inform[n].RefSwapMaturity;

                                Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);
                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardEnd[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);

                                Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv_PowerSpread[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv_PowerSpread[n][j],
                                    QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv2[n][Dayidx] = Rate;

                        }
                    }
                }
                else if (RcvLeg->Reference_Inform[n].RefRateType == 2) // SOFR
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0.0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                SOFRDailyCompound = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardEnd[Dayidx], RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv[n][j], QVTerm_Rcv[n][j], B_t_T_Rcv[n][j], Rate);
                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                Rate = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardStart[Dayidx] + (long)(RcvLeg->Reference_Inform[n].PowerSpreadMat1 * 365.0), RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv[n][j], QVTerm_Rcv[n][j], B_t_T_Rcv[n][j], Rate1);

                                Rate2 = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardStart[Dayidx] + (long)(RcvLeg->Reference_Inform[n].PowerSpreadMat2 * 365.0), RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv_PowerSpread[n][j], QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j], Rate2);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv[n][Dayidx] = Rate;

                        }
                    }
                }
            }
            else
            {

                //Libor, Zero의 경우 
                if (RcvLeg->Reference_Inform[n].RefRateType == 0)
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(j, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;

                                Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);
                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                    QVTerm_Rcv[n][Dayidx], B_t_T_Rcv[n][j]);

                                Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Rcv_PowerSpread[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv_PowerSpread[n][j],
                                    QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv[n][Dayidx] = Rate;

                        }
                    }
                }
                else if (RcvLeg->Reference_Inform[n].RefRateType == 1) // Libor Swap
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        xt = SimulShortRate[Curveidx][j];
                        t = Simul->T_Array[j];
                        if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            T = t + RcvLeg->Reference_Inform[n].RefSwapMaturity;
                            Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);
                        }
                        else
                        {
                            T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;
                            Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv[n][j],
                                QVTerm_Rcv[n][j], B_t_T_Rcv[n][j]);

                            Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Rcv_PowerSpread[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Rcv_PowerSpread[n][j],
                                QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j]);
                            Rate = Rate1 - Rate2;
                        }
                        Dayidx = 0;
                        RcvDailyRate[n][j] = Rate;
                        if (isinFindIndex(j, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx)) SimulatedRateRcv[n][Dayidx] = Rate;
                        if (isinFindIndex(j, RcvLeg->DaysForwardEnd, RcvLeg->NCashFlow, Dayidx)) SimulatedRateRcv2[n][Dayidx] = Rate;

                    }
                }
                else if (RcvLeg->Reference_Inform[n].RefRateType == 2) // SOFR
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(j, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0.0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                SOFRDailyCompound = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardEnd[Dayidx], RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv[n][j], QVTerm_Rcv[n][j], B_t_T_Rcv[n][j], Rate);
                            }
                            else
                            {
                                t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                Rate = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardStart[Dayidx] + (long)(RcvLeg->Reference_Inform[n].PowerSpreadMat1 * 365.0), RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv[n][j], QVTerm_Rcv[n][j], B_t_T_Rcv[n][j], Rate1);

                                Rate2 = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    RcvLeg->DaysForwardStart[Dayidx], RcvLeg->DaysForwardStart[Dayidx] + (long)(RcvLeg->Reference_Inform[n].PowerSpreadMat2 * 365.0), RateHistoryUseFlag, RcvLeg->FixedRate[Dayidx], RcvLeg->HolidayCalcFlag[n],
                                    RcvLeg->HolidayCount[n], RcvLeg->HolidayDays[n], PV_t_T_Rcv_PowerSpread[n][j], QVTerm_Rcv_PowerSpread[n][j], B_t_T_Rcv_PowerSpread[n][j], Rate2);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRateRcv[n][Dayidx] = Rate;

                        }
                    }
                }
            }
        }
        for (j = 0; j < RcvLeg->NCashFlow; j++) RcvCFFlag[j] = 0;

        for (j = 0; j < RcvLeg->NCashFlow; j++)
        {
            RateHistoryUseFlag = 0;
            Day1 = RcvLeg->DaysForwardStart[j];
            Day2 = RcvLeg->DaysPayDate_C[j];
            t = ((double)Day1) / 365.0;
            T = ((double)RcvLeg->DaysPayDate_C[j]) / 365.0;
            Pt = Rcv_DF_0_t[j];
            PT = Rcv_DF_0_T[j];
            if (Day1 >= 0)
            {
                RcvCFFlag[j] = 1;
                RateHistoryUseFlag = 0;
            }
            else if (Day2 > 0)
            {
                RcvCFFlag[j] = 1;
                RateHistoryUseFlag = 1;
            }
            else
            {
                RcvCFFlag[j] = 0;
                RateHistoryUseFlag = 0;
            }

            if (RcvCFFlag[j] == 1)
            {
                RcvPayoff[j] = Notional * PayoffStructure(Simul->DailySimulFlag, RcvLeg, SimulatedRateRcv, SimulatedRateRcv2, RcvDailyRate, SimulOverNightRate, j, RateHistoryUseFlag, RcvOutputRate, nAccrual);
                xt_idx = FindIndex(Day1, Simul->DaysForSimul, Simul->NDays);
                xt = SimulShortRate[CurveIdx_DiscRcv][xt_idx];
                if (t > 0.0)
                {
                    PtT = PT / Pt * exp(-xt * B_t_T_RcvDisc[j] + QVTerm_RcvDisc[j]);
                }
                else
                {
                    Pt = 1.0;
                    PtT = Rcv_DF_0_T[j];
                }
            }
            else
            {
                RcvPayoff[j] = 0.0;
                for (n = 0; n < RcvLeg->NReference; n++) RcvOutputRate[n] = 0.0;
                nAccrual = 0;
                Pt = 1.0;
                PtT = Rcv_DF_0_T[j];
            }

            PricePath_Rcv += Pt * PtT * RcvPayoff[j];

            for (n = 0; n < RcvLeg->NReference; n++)
            {
                ResultRcv[j + n * RcvLeg->NCashFlow] += RcvOutputRate[n] / (double)Simul->NSimul;
            }
            ResultRcv[3 * RcvLeg->NCashFlow + j] += ((double)nAccrual) / (double)Simul->NSimul;
            ResultRcv[4 * RcvLeg->NCashFlow + j] += RcvPayoff[j] / (double)Simul->NSimul;
            ResultRcv[5 * RcvLeg->NCashFlow + j] += Pt * PtT / (double)Simul->NSimul;

            if (RcvLeg->OptionUseFlag == 1)
            {
                for (n = 0; n < RcvLeg->NOption; n++)
                {
                    if ( RcvLeg->DaysOptionDate[n] < Day2)
                    {
                        Y[n][i] += Pt * PtT * RcvPayoff[j];
                    }
                }
            }
        }

        if (NAFlag == 1) PricePath_Rcv += Notional * Pt * PtT;

        if (RcvLeg->OptionUseFlag == 1)
        {
            for (n = 0; n < RcvLeg->NOption; n++)
            {
                Y[n][i] += Pt * PtT * Notional * (double)(NAFlag == 1) ;
                if (RcvLeg->CallConditionFlag == 0)
                {
                    OptRangeFlag[n][i] = 1;
                    for (k = 0; k < RcvLeg->NReference; k++)
                    {
                        if (RcvOutputRate[k] > RcvLeg->RangeUp[k][n] || RcvOutputRate[k] < RcvLeg->RangeDn[k][n])
                        {
                            OptRangeFlag[n][i] = 0;
                            break;
                        }
                    }
                }
                else
                {
                    OptRangeFlag[n][i] = 0;
                    for (k = 0; k < RcvLeg->NReference; k++)
                    {
                        if (RcvOutputRate[k] <= RcvLeg->RangeUp[k][n] && RcvOutputRate[k] >= RcvLeg->RangeDn[k][n])
                        {
                            OptRangeFlag[n][i] = 1;
                            break;
                        }
                    }
                }
            }
        }


        RcvPrice += PricePath_Rcv / (double)Simul->NSimul;


        // Rcv 끝
        // PayLeg 시작
        // PayLeg Reference Rate 결정
        for (n = 0; n < PayLeg->NReference; n++)
        {
            Curveidx = CurveIdx_Pay[n];
            lastresetend = -1;
            Rate = 0.0;
            Rate1 = 0.0;
            Rate2 = 0.0;
            if (Simul->DailySimulFlag == 0)
            {
                //Libor, Zero의 경우 
                if (PayLeg->Reference_Inform[n].RefRateType == 0)
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;

                                Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][j], B_t_T_Pay[n][j]);

                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][j], B_t_T_Pay[n][j]);

                                Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay_PowerSpread[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay_PowerSpread[n][j],
                                    QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay[n][Dayidx] = Rate;
                        }
                    }
                }
                else if (PayLeg->Reference_Inform[n].RefRateType == 1) // Libor Swap
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = t + PayLeg->Reference_Inform[n].RefSwapMaturity;

                                if (Simul->DaysForSimul[j] != lastresetend) // 이전 Reset End와 Reset Start가 같으면 중복계산하지않는다.(연산속도)
                                {
                                    Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                        QVTerm_Pay[n][j], B_t_T_Pay[n][j]);
                                }

                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                                if (Simul->DaysForSimul[j] != lastresetend)
                                {
                                    Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                        QVTerm_Pay[n][j], B_t_T_Pay[n][j]);

                                    Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                        Simul->Rate[Curveidx], xt, ndates_Pay_PowerSpread[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay_PowerSpread[n][j],
                                        QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j]);
                                }

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay[n][Dayidx] = Rate;
                        }
                        if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardEnd, PayLeg->NCashFlow, Dayidx))
                        {
                            lastresetend = Simul->DaysForSimul[j];
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardEnd[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = t + PayLeg->Reference_Inform[n].RefSwapMaturity;

                                Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][j], B_t_T_Pay[n][j]);
                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardEnd[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][j], B_t_T_Pay[n][j]);

                                Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay_PowerSpread[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay_PowerSpread[n][j],
                                    QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay2[n][Dayidx] = Rate;
                        }
                    }
                }
                else if (PayLeg->Reference_Inform[n].RefRateType == 2) // SOFR
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0.0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                SOFRDailyCompound = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardEnd[Dayidx], RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay[n][j], QVTerm_Pay[n][j], B_t_T_Pay[n][j], Rate);
                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                Rate = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardStart[Dayidx] + (long)(PayLeg->Reference_Inform[n].PowerSpreadMat1 * 365.0), RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay[n][j], QVTerm_Pay[n][j], B_t_T_Pay[n][j], Rate1);

                                Rate2 = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardStart[Dayidx] + (long)(PayLeg->Reference_Inform[n].PowerSpreadMat2 * 365.0), RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay_PowerSpread[n][j], QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j], Rate2);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay[n][Dayidx] = Rate;
                        }
                    }
                }
            }
            else
            {
                //Libor, Zero의 경우 
                if (PayLeg->Reference_Inform[n].RefRateType == 0)
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(j, PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;

                                Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][j], B_t_T_Pay[n][j]);
                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                    QVTerm_Pay[n][Dayidx], B_t_T_Pay[n][j]);

                                Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, ndates_Pay_PowerSpread[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay_PowerSpread[n][j],
                                    QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay[n][Dayidx] = Rate;
                        }
                    }
                }
                else if (PayLeg->Reference_Inform[n].RefRateType == 1) // Libor Swap
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        xt = SimulShortRate[Curveidx][j];
                        t = Simul->T_Array[j];
                        if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            T = t + PayLeg->Reference_Inform[n].RefSwapMaturity;
                            Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                QVTerm_Pay[n][j], B_t_T_Pay[n][j]);
                        }
                        else
                        {
                            T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;
                            Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay[n][j],
                                QVTerm_Pay[n][j], B_t_T_Pay[n][j]);

                            Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, ndates_Pay_PowerSpread[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, PV_t_T_Pay_PowerSpread[n][j],
                                QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j]);
                            Rate = Rate1 - Rate2;
                        }
                        Dayidx = 0;
                        PayDailyRate[n][j] = Rate;
                        if (isinFindIndex(j, PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx)) SimulatedRatePay[n][Dayidx] = Rate;
                        if (isinFindIndex(j, PayLeg->DaysForwardEnd, PayLeg->NCashFlow, Dayidx)) SimulatedRatePay2[n][Dayidx] = Rate;
                    }
                }
                else if (PayLeg->Reference_Inform[n].RefRateType == 2) // SOFR
                {
                    for (j = 0; j < Simul->NDays; j++)
                    {
                        if (isinFindIndex(j, PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx))
                        {
                            xt = SimulShortRate[Curveidx][j];
                            // PowerSpread가 아닌 경우
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0.0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                SOFRDailyCompound = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardEnd[Dayidx], RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay[n][j], QVTerm_Pay[n][j], B_t_T_Pay[n][j], Rate);
                            }
                            else
                            {
                                t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;
                                kappa = Simul->HWKappa[Curveidx];
                                if (t < 0) RateHistoryUseFlag = 1;
                                else RateHistoryUseFlag = 0;

                                Rate = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardStart[Dayidx] + (long)(PayLeg->Reference_Inform[n].PowerSpreadMat1 * 365.0), RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay[n][j], QVTerm_Pay[n][j], B_t_T_Pay[n][j], Rate1);

                                Rate2 = SOFR_ForwardRate_Compound_DiscreteSimul(Simul->DailySimulFlag, Simul->NDays, Simul->T_Array, SimulShortRate[Curveidx],
                                    PayLeg->DaysForwardStart[Dayidx], PayLeg->DaysForwardStart[Dayidx] + (long)(PayLeg->Reference_Inform[n].PowerSpreadMat2 * 365.0), RateHistoryUseFlag, PayLeg->FixedRate[Dayidx], PayLeg->HolidayCalcFlag[n],
                                    PayLeg->HolidayCount[n], PayLeg->HolidayDays[n], PV_t_T_Pay_PowerSpread[n][j], QVTerm_Pay_PowerSpread[n][j], B_t_T_Pay_PowerSpread[n][j], Rate2);

                                Rate = Rate1 - Rate2;
                            }
                            SimulatedRatePay[n][Dayidx] = Rate;
                        }
                    }
                }
            }
        }
        for (j = 0; j < PayLeg->NCashFlow; j++) PayCFFlag[j] = 0;

        for (j = 0; j < PayLeg->NCashFlow; j++)
        {
            RateHistoryUseFlag = 0;
            Day1 = PayLeg->DaysForwardStart[j];
            Day2 = PayLeg->DaysPayDate_C[j];
            t = ((double)Day1) / 365.0;
            T = ((double)PayLeg->DaysPayDate_C[j]) / 365.0;
            Pt = Pay_DF_0_t[j];
            PT = Pay_DF_0_T[j];
            if (Day1 >= 0)
            {
                PayCFFlag[j] = 1;
                RateHistoryUseFlag = 0;
            }
            else if (Day2 > 0)
            {
                PayCFFlag[j] = 1;
                RateHistoryUseFlag = 1;
            }
            else
            {
                PayCFFlag[j] = 0;
                RateHistoryUseFlag = 0;
            }

            if (PayCFFlag[j] == 1)
            {
                PayPayoff[j] = Notional * PayoffStructure(Simul->DailySimulFlag, PayLeg, SimulatedRatePay, SimulatedRatePay2, PayDailyRate, SimulOverNightRate, j, RateHistoryUseFlag, PayOutputRate, nAccrual);
                xt_idx = FindIndex(Day1, Simul->DaysForSimul, Simul->NDays);
                xt = SimulShortRate[CurveIdx_DiscPay][xt_idx];
                if (t > 0.0)
                {
                    //PtT = P_t_to_T(t, T, Simul->HWKappa[CurveIdx_DiscPay], xt, Simul->NHWVol[CurveIdx_DiscPay], Simul->HWVolTerm[CurveIdx_DiscPay], Simul->HWVol[CurveIdx_DiscPay], Pt, PT, 4);
                    PtT = PT / Pt * exp(-xt * B_t_T_PayDisc[j] + QVTerm_PayDisc[j]);
                }
                else
                {
                    Pt = 1.0;
                    PtT = Rcv_DF_0_T[j];
                }

            }
            else
            {
                PayPayoff[j] = 0.0;
                for (n = 0; n < PayLeg->NReference; n++) PayOutputRate[n] = 0.0;
                nAccrual = 0;
                Pt = 1.0;
                PtT = Rcv_DF_0_T[j];
            }

            PricePath_Pay += Pt * PtT * PayPayoff[j];

            for (n = 0; n < PayLeg->NReference; n++)
            {
                ResultPay[j + n * PayLeg->NCashFlow] += PayOutputRate[n] / (double)Simul->NSimul;
            }
            ResultPay[3 * PayLeg->NCashFlow + j] += ((double)nAccrual) / (double)Simul->NSimul;
            ResultPay[4 * PayLeg->NCashFlow + j] += PayPayoff[j] / (double)Simul->NSimul;
            ResultPay[5 * PayLeg->NCashFlow + j] += Pt * PtT / (double)Simul->NSimul;
            if (RcvLeg->OptionUseFlag == 1)
            {
                for (n = 0; n < RcvLeg->NOption; n++)
                {
                    if (RcvLeg->DaysOptionDate[n] < Day2)
                    {
                        Y[n][i] -= Pt * PtT * PayPayoff[j];
                    }
                }
            }
        }

        PayPrice += PricePath_Pay / (double)Simul->NSimul;
        
        if (RcvLeg->OptionUseFlag == 1)
        {
            for (n = 0; n < RcvLeg->NOption; n++)
            {
                Y[n][i] -= Pt * PtT * Notional * (double)(NAFlag == 1);
            }

            if (OptRangeFlag[0][i] > 0)
            {
                if (RcvLeg->OptionType == 1) {
                    if (Y[0][i] < 0.0) Value_By_OptTime[i] = fabs(Y[0][i]);
                    else Value_By_OptTime[i] = 0.0;
                }
                else
                {
                    if (Y[0][i] > 0.0) Value_By_OptTime[i] = Y[0][i];
                    else Value_By_OptTime[i] = 0.0;
                }
            }
            else Value_By_OptTime[i] = 0.0;             
            OptimalIdx[i] = 0;

            for (n = 1; n < RcvLeg->NOption; n++)
            {
                if (RcvLeg->OptionType == 1)
                {
                    if (Y[n][i] < 0.0 && -Y[n][i] > Value_By_OptTime[i] & OptRangeFlag[n][i] > 0)
                    {
                        OptimalIdx[i] = n;
                        Value_By_OptTime[i] = -Y[n][i];
                    }
                }
                else
                {
                    if (Y[n][i] > 0.0 && Y[n][i] > Value_By_OptTime[i] & OptRangeFlag[n][i] > 0)
                    {
                        OptimalIdx[i] = n;
                        Value_By_OptTime[i] = Y[n][i];
                    }
                }
            }

            for (n = 0; n < Simul->NAsset; n++)
            {
                InterestRate_Opt[i][n] = X[OptimalIdx[i]][i][n];
            }

        }
        
        n = 0;

    }

    //////////////
    // LSMC OLS //
    //////////////
    double OptionPrice = 0.0;
    double MaxLegValue;
    double ExerciseValue = 0.0;
    if (RcvLeg->OptionUseFlag == 1)
    {
        Beta = OLSBeta(Value_By_OptTime, LengthY, InterestRate_Opt, ShapeX);
        for (i = 0; i < Simul->NSimul; i++)
        {
            for (idxoption = 0; idxoption < RcvLeg->NOption; idxoption++)
            {
                EstOptValue = Beta[Simul->NAsset];
                for (j = 0; j < Simul->NAsset; j++)
                {
                    EstOptValue += Beta[j] * X[idxoption][i][j];
                }
                
                if (RcvLeg->OptionType == 1)
                {
                    if (Y[idxoption][i] < 0.0 && -Y[idxoption][i] > EstOptValue)
                    {
                        OptionPrice += -Y[idxoption][i] / (double)Simul->NSimul;
                        break;
                    }
                }
                else
                {
                    if (Y[idxoption][i] > 0.0 && Y[idxoption][i] > EstOptValue)
                    {
                        OptionPrice += Y[idxoption][i] / (double)Simul->NSimul;
                        break;
                    }
                }
            }
        }
        
        free(Beta);
    }
    

    if (PricingOnly == 1)
    {
        ResultPrice[0] = RcvPrice;
        ResultPrice[1] = PayPrice;
        ResultPrice[2] = RcvPrice - PayPrice;
        if (RcvLeg->OptionUseFlag == 1)
        {
            ResultPrice[3] = OptionPrice;
            if (RcvLeg->OptionType == 1) ResultPrice[2] += OptionPrice;
            else ResultPrice[2] -= OptionPrice;
        }
    }

    if (Simul->DailySimulFlag == 1)
    {
        for (i = 0; i < RcvLeg->NReference; i++) free(RcvDailyRate[i]);
        for (i = 0; i < PayLeg->NReference; i++) free(PayDailyRate[i]);
    }
    free(RcvDailyRate);
    free(PayDailyRate);
    free(RcvPayoff);
    free(PayPayoff);
    free(RcvCFFlag);
    free(PayCFFlag);
    free(Rcv_DF_0_t);
    free(Pay_DF_0_t);
    free(Rcv_DF_0_T);
    free(Pay_DF_0_T);
    free(B_t_T_RcvDisc);
    free(B_t_T_PayDisc);
    free(QVTerm_RcvDisc);
    free(QVTerm_PayDisc);
    free(RcvOutputRate);
    free(PayOutputRate);

    for (i = 0; i < RcvLeg->NReference; i++)
    {
        if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Rcv[i][j]);
                free(QVTerm_Rcv[i][j]);
                free(B_t_T_Rcv[i][j]);
            }
        }
        else
        {

            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Rcv[i][j]);
                free(QVTerm_Rcv[i][j]);
                free(B_t_T_Rcv[i][j]);
            }

            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Rcv_PowerSpread[i][j]);
                free(QVTerm_Rcv_PowerSpread[i][j]);
                free(B_t_T_Rcv_PowerSpread[i][j]);
            }
            free(PV_t_T_Rcv_PowerSpread[i]);
            free(QVTerm_Rcv_PowerSpread[i]);
            free(B_t_T_Rcv_PowerSpread[i]);
        }

        free(PV_t_T_Rcv[i]);
        free(QVTerm_Rcv[i]);
        free(B_t_T_Rcv[i]);
    }
    free(PV_t_T_Rcv);
    free(QVTerm_Rcv);
    free(B_t_T_Rcv);
    free(PV_t_T_Rcv_PowerSpread);
    free(QVTerm_Rcv_PowerSpread);
    free(B_t_T_Rcv_PowerSpread);
    free(ndates_Rcv);
    free(ndates_Rcv_PowerSpread);

    for (i = 0; i < PayLeg->NReference; i++)
    {
        if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Pay[i][j]);
                free(QVTerm_Pay[i][j]);
                free(B_t_T_Pay[i][j]);
            }
        }
        else
        {

            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Pay[i][j]);
                free(QVTerm_Pay[i][j]);
                free(B_t_T_Pay[i][j]);
            }

            for (j = 0; j < Simul->NDays; j++)
            {
                free(PV_t_T_Pay_PowerSpread[i][j]);
                free(QVTerm_Pay_PowerSpread[i][j]);
                free(B_t_T_Pay_PowerSpread[i][j]);
            }
            free(PV_t_T_Pay_PowerSpread[i]);
            free(QVTerm_Pay_PowerSpread[i]);
            free(B_t_T_Pay_PowerSpread[i]);
        }

        free(PV_t_T_Pay[i]);
        free(QVTerm_Pay[i]);
        free(B_t_T_Pay[i]);
    }
    free(PV_t_T_Pay);
    free(QVTerm_Pay);
    free(B_t_T_Pay);
    free(PV_t_T_Pay_PowerSpread);
    free(QVTerm_Pay_PowerSpread);
    free(B_t_T_Pay_PowerSpread);
    free(ndates_Pay);
    free(ndates_Pay_PowerSpread);


    // 제거중
    if (RcvLeg->OptionUseFlag == 1)
    {
        free(Value_By_OptTime);
        free(Estimated_Value_By_OptTime);
        free(OptimalIdx);

        for (i = 0; i < Simul->NSimul; i++) free(InterestRate_Opt[i]); 
        free(InterestRate_Opt);

        for (i = 0; i < RcvLeg->NOption; i++)
        {
            for (j = 0; j < Simul->NSimul; j++) free(X[i][j]);
            free(X[i]);
        }
        free(X);

        for (i = 0; i < RcvLeg->NOption; i++) free(Y[i]);
        free(Y);

        for (i = 0; i < RcvLeg->NOption; i++) free(OptRangeFlag[i]);
        free(OptRangeFlag);

        for (i = 0; i < RcvLeg->NOption; i++) free(Y_Hat[i]); 

    }
    else
    {
        RcvLeg->NOption = 0;
        free(X); 
        free(Y); 
        free(Value_By_OptTime); 
        free(Estimated_Value_By_OptTime); 
        free(OptimalIdx); 
        free(InterestRate_Opt); 
        free(OptRangeFlag); 
        free(Y_Hat);

    }

    return 1;
}

long IRStructuredSwap(
    long PricingDateC,
    long NAFlag,
    double Notional,
    LEG_INFO* RcvLeg,
    LEG_INFO* PayLeg,
    SIMUL_INFO* Simul,
    long GreekFlag,                     // 
    double* ResultPrice,                // 
    double* ResultRcv,                  // 
    double* ResultPay                   // 
)
{
    long i;

    long ResultCode = 0;

    double xt = 0.0;
    double Rate = 0.0, Rate1 = 0.0, Rate2 = 0.0;

    long Length_XA_Array = Simul->NDays;
    long Length_XV_Array = Simul->NDays;
    double** XA_Array = (double**)malloc(sizeof(double*) * Simul->NAsset);
    double** XV_Array = (double**)malloc(sizeof(double*) * Simul->NAsset);
    for (i = 0; i < Simul->NAsset; i++)
    {
        XA_Array[i] = (double*)malloc(sizeof(double) * Simul->NDays);
        XV_Array[i] = (double*)malloc(sizeof(double) * Simul->NDays);
        Set_XA(Simul->HWKappa[i], Length_XA_Array, XA_Array[i], Simul->dt_Array);
        Set_XV(Simul->HWKappa[i], Simul->NHWVol[i], Simul->HWVolTerm[i], Simul->HWVol[i], Length_XV_Array, XV_Array[i], Simul->dt_Array);
    }

    double** SimulatedRateRcv = make_array(RcvLeg->NReference, RcvLeg->NCashFlow);
    double** SimulatedRateRcv2 = make_array(RcvLeg->NReference, RcvLeg->NCashFlow);

    double** SimulatedRatePay = make_array(PayLeg->NReference, PayLeg->NCashFlow);
    double** SimulatedRatePay2 = make_array(PayLeg->NReference, PayLeg->NCashFlow);

    double** SimulShortRate = make_array(Simul->NAsset, Simul->NDays);
    double** SimulOverNightRate = make_array(Simul->NAsset, Simul->NDays);

    long* CurveIdx_Rcv = (long*)malloc(sizeof(long) * RcvLeg->NReference);
    long* CurveIdx_Pay = (long*)malloc(sizeof(long) * PayLeg->NReference);
    for (i = 0; i < RcvLeg->NReference; i++)
    {
        CurveIdx_Rcv[i] = FindIndex(RcvLeg->Reference_Inform[i].CurveNum, Simul->SimulCurveIdx, Simul->NAsset);
    }

    for (i = 0; i < PayLeg->NReference; i++)
    {
        CurveIdx_Pay[i] = FindIndex(PayLeg->Reference_Inform[i].CurveNum, Simul->SimulCurveIdx, Simul->NAsset);
    }

    ResultCode = Simulate_HW(1, PricingDateC, NAFlag, Notional, RcvLeg, PayLeg,
        Simul, Length_XA_Array, XA_Array, Length_XV_Array, XV_Array,
        SimulatedRateRcv, SimulatedRateRcv2, SimulatedRatePay, SimulatedRatePay2, SimulShortRate, SimulOverNightRate, CurveIdx_Rcv,
        CurveIdx_Pay, ResultPrice, ResultRcv, ResultPay);



    for (i = 0; i < Simul->NAsset; i++)
    {
        free(XA_Array[i]);
        free(XV_Array[i]);
    }
    free(XA_Array);
    free(XV_Array);

    for (i = 0; i < RcvLeg->NReference; i++)
    {
        free(SimulatedRateRcv[i]);
    }
    free(SimulatedRateRcv);

    for (i = 0; i < RcvLeg->NReference; i++)
    {
        free(SimulatedRateRcv2[i]);
    }
    free(SimulatedRateRcv2);

    for (i = 0; i < PayLeg->NReference; i++)
    {
        free(SimulatedRatePay[i]);
    }
    free(SimulatedRatePay);

    for (i = 0; i < PayLeg->NReference; i++)
    {
        free(SimulatedRatePay2[i]);
    }
    free(SimulatedRatePay2);

    for (i = 0; i < Simul->NAsset; i++)
    {
        free(SimulShortRate[i]);
    }
    free(SimulShortRate);

    for (i = 0; i < Simul->NAsset; i++)
    {
        free(SimulOverNightRate[i]);
    }
    free(SimulOverNightRate);

    free(CurveIdx_Rcv);
    free(CurveIdx_Pay);
    return 1;
}

DLLEXPORT(long) Pricing_IRStructuredSwap_Excel(
    long PriceDateExcel,                // 1 가격계산일 Excel long Type
    long NAFlag,                        // 2 Notional Amount 사용 여부 0: No ExChange Notional Amount 1: ExChange Notioanl Amount
    double Notional,                    // 3 Notional Amount
    long* PayoffStructure,              // 4 [0]: RcvLeg 페이오프조건 [1]: PayLeg 페이오프조건 {페이오프조건 0:Range(R1)&Range(R2)&Range(R3) 1:Range(R1+R2+R3)}
    long* HolidayFlagCount,             // 5 [0~2]: R1,R2,R3 Rcv HolidayFlag, [3~5]: R1, R2, R3 Pay HolidayFlag, [6~8]: R1, R2, R3 Rcv HolidayCount, [9~11]: R1, R2, R3 Pay HolidayCount
    long* Holidays_Excel,               // 6 [0~sum(HolidayFlagCount[6~11])]: Holidays Rcv R1R2R3, Pay R1R2R3

    // Receive Leg Information
    long* RcvFlag,                      // 7 [0]: NReference, [1]: Fix/Flo Flag, [2]: Accrual Flag, [3]: DailySimulFlag
    double* RcvMaxLossReturn,           // 8 [0]: Max Loss [1] Max Return
    long* RcvRefCurveNumber,            // 9 Reference Curve 번호 Array [Len = NReference]
    long* RcvRefRateType,               // 10 기초금리 타입 0: CD/Libor (AccrualX) 1: CD/Libor, 이자율Swap 2: SOFR
    double* RcvSwapInfo,                // 11 [짝수] 기초금리 타입이 1일 경우 연 이자 지급횟수 [홀수] 만기
    long* RcvDayCountFlag,              // 12 DayCount 0: 365 ,  1: 360
    double* RcvRefCondMultiple,         // 13 금리결정조건 Multiple
    long* RcvRefPowerSpreadFlag,        // 14 파워스프레드 사용여부
    double* RcvRefPowerSpreadMat,       // 15 [짝수] 파워스프레드 만기1 [홀수] 파워스프레드 만기2
    double* RcvRefRange,                // 16 [짝수] Range 상한 [홀수] Range 하한
    double* RcvPayoffMultiple,          // 17 페이오프 조건 Multiple
    long NRcvCashFlow,                  // 18 Receive Leg CashFlow 개수
    long* RcvCashFlowScheduleExcelDate, // 19 [0~NCF-1]: ResetStart, [NCF~2*NCF-1]: ResetEnd, [2*NCF~3*NCF-1]: 기산일, [3*NCF~4*NCF-1]: 기말일, [4*NCF~5*NCF-1]: 지급일
    double* RcvCouponFixRate,           // 20 [0~NCF-1]: 고정 쿠폰, [NCF~2*NCF-1]: 확정과거금리 [2*NCF~3*NCF-1]: Range충족쿠폰Rate

    // Receive Leg Information
    long* PayFlag,                      // 21 [0]: NReference, [1]: Fix/Flo Flag, [2]: Accrual Flag
    double* PayMaxLossReturn,           // 22 [0]: Max Loss [1] Max Return
    long* PayRefCurveNumber,            // 23 Reference Curve 번호 Array [Len = NReference]
    long* PayRefRateType,               // 24 기초금리 타입 0: CD/Libor (AccrualX) 1: CD/Libor, 이자율Swap 2: SOFR
    double* PaySwapInfo,                // 25 [짝수] 기초금리 타입이 1일 경우 연 이자 지급횟수 [홀수] 만기
    long* PayDayCountFlag,              // 26 DayCount 0: 365 ,  1: 360
    double* PayRefCondMultiple,         // 27 금리결정조건 Multiple
    long* PayRefPowerSpreadFlag,        // 28 파워스프레드 사용여부
    double* PayRefPowerSpreadMat,       // 29 [짝수] 파워스프레드 만기1 [홀수] 파워스프레드 만기2
    double* PayRefRange,                // 30 [짝수] Range 상한 [홀수] Range 하한
    double* PayPayoffMultiple,          // 31 페이오프 조건 Multiple
    long NPayCashFlow,                  // 32 Receive Leg CashFlow 개수
    long* PayCashFlowScheduleExcelDate, // 33 [0~NCF-1]: ResetStart, [NCF~2*NCF-1]: ResetEnd, [2*NCF~3*NCF-1]: 기산일, [3*NCF~4*NCF-1]: 기말일, [4*NCF~5*NCF-1]: 지급일
    double* PayCouponFixRate,           // 34 [0~NCF-1]: 고정 쿠폰, [NCF~2*NCF-1]: 확정과거금리 [2*NCF~3*NCF-1]: Range충족쿠폰Rate

    // 커브 및 HW정보
    long* NHWVolCount,                  // 35 [0~3] Hull White Vol 개수
    double* HWVolTermArray,             // 36 [0~sum(NHWVolCount)-1] Hull White Vol Term Array
    double* HWVolArray,                 // 37 [0~sum(NHWVolCount)-1] Hull White Vol Array
    double* HWKappaArray,               // 38 [0~3] Hull White Kappa
    long* NZeroRate,                    // 39 [0~3] ZeroRate 개수
    double* ZeroTerm,                   // 40 [0~sum(NZeroRate)-1] Zero Rate Term
    double* ZeroRate,                   // 41 [0~sum(NZeroRate)-1] Zero Rate

    // 상관계수 및 히스토리 
    double* Corr,                       // 42 Correlation Matrix.rehaped(-1)
    long* NDayHistRcv,                  // 43 [0~3] Rcv Rate History 개수
    long* RateHistDateRcv,              // 44 [0~sum(NDayHistRcv)-1] RateHistory의 Date
    double* RateHistRcv,                // 45 [0~sum(NDayHistRcv)-1] RateHistory의 Rate
    long* NDayHistPay,                  // 46 [0~3] Pay Rate History 개수
    long* RateHistDatePay,              // 47 [0~sum(NDayHistPay)-1] RateHistory의 Date
    double* RateHistPay,                // 48 [0~sum(NDayHistPay)-1] RateHistory의 Rate

    // 시뮬레이션
    long NSimul,                        // 49 시뮬레이션횟수
    long GreekFlag,                     // 50 Greek산출여부

    // 옵션관련
    long NOption,                       // 51 옵션개수
    long* OptionDateAndFlag,            // 52 [0~NOption-1]옵션행사일, [NOption~2NOption-1]행사조건Flag, [2NOption~3NOption-1]옵션Type
    double* OptionKAndRange,            // 53 [0~NOption-1]옵션행사가격% [NOption~7NOption-1] Ref1~3 Range상한 및 하한 

    double* ResultPrice,                // 54 산출된 가격 [0] Rcv [1] Pay [2] Price
    double* ResultRcv,                  // 55 [0~NCF-1]금리1, [NCF~2NCF-1]금리2, [2NCF~3NCF-1]금리3, [3NCF~4NCF-1] E(Accrual수), [4NCF~5NCF-1] E(CPN), [5NCF~6NCF-1] E(DF)
    double* ResultPay,                  // 56 [0~NCF-1]금리1, [NCF~2NCF-1]금리2, [2NCF~3NCF-1]금리3, [3NCF~4NCF-1] E(Accrual수), [4NCF~5NCF-1] E(CPN), [5NCF~6NCF-1] E(DF)
    char* Error                         // 57 에러메시지
)
{

    long i;
    //_CrtSetBreakAlloc(13172);
    long ResultCode = 0;
    long ErrorCode = 0;
    long j;
    long k;
    long n;

    ErrorCode = ErrorCheck_IRStructuredSwap_Excel(
        PriceDateExcel, NAFlag, Notional, PayoffStructure, HolidayFlagCount,
        Holidays_Excel, RcvFlag, RcvMaxLossReturn, RcvRefCurveNumber, RcvRefRateType,
        RcvSwapInfo, RcvDayCountFlag, RcvRefCondMultiple, RcvRefPowerSpreadFlag, RcvRefPowerSpreadMat,
        RcvRefRange, RcvPayoffMultiple, NRcvCashFlow, RcvCashFlowScheduleExcelDate, RcvCouponFixRate,
        PayFlag, PayMaxLossReturn, PayRefCurveNumber, PayRefRateType, PaySwapInfo,
        PayDayCountFlag, PayRefCondMultiple, PayRefPowerSpreadFlag, PayRefPowerSpreadMat, PayRefRange,
        PayPayoffMultiple, NPayCashFlow, PayCashFlowScheduleExcelDate, PayCouponFixRate, NHWVolCount,
        HWVolTermArray, HWVolArray, HWKappaArray, NZeroRate, ZeroTerm,
        ZeroRate, Corr, NDayHistRcv, RateHistDateRcv, RateHistRcv,
        NDayHistPay, RateHistDatePay, RateHistPay, NSimul, GreekFlag, Error);

    if (ErrorCode < 0) return -1;

    long N_Curve_Max = 4;       // 커브 최대 개수
    long N_Ref_Max = 3;         // 레퍼런스 최대 개수

    //////////////////////////////////////////////////////////////////////////
    // HWVolMatrix, ZeroRateMatrix 매핑 Shape = (N_Curve_Max, 각 Term 개수) //
    //////////////////////////////////////////////////////////////////////////

    double** HWVolTermMatrix = (double**)malloc(sizeof(double*) * N_Curve_Max);
    double** HWVolMatrix = (double**)malloc(sizeof(double*) * N_Curve_Max);

    n = 0;
    for (i = 0; i < N_Curve_Max; i++)
    {
        HWVolTermMatrix[i] = HWVolTermArray + n;
        HWVolMatrix[i] = HWVolArray + n;
        n = n + NHWVolCount[i];
    }

    double** ZeroRateTermMatrix = (double**)malloc(sizeof(double*) * N_Curve_Max);
    double** ZeroRateMatrix = (double**)malloc(sizeof(double*) * N_Curve_Max);
    n = 0;
    for (i = 0; i < N_Curve_Max; i++)
    {
        ZeroRateTermMatrix[i] = ZeroTerm + n;
        ZeroRateMatrix[i] = ZeroRate + n;
        n = n + NZeroRate[i];
    }

    long DailySimulFlag = RcvFlag[3];

    long PriceDateC = ExcelDateToCDate(PriceDateExcel);
    long RcvPayoffStructure = PayoffStructure[0];
    long PayPayoffStructure = PayoffStructure[1];

    ///////////////
    // Holiday Mapping 시작
    ///////////////

    long* RcvHolidayCalcFlag = HolidayFlagCount;
    long* PayHolidayCalcFlag = HolidayFlagCount + N_Ref_Max;
    long* RcvHolidayCount = HolidayFlagCount + 2 * N_Ref_Max;
    long* PayHolidayCount = HolidayFlagCount + 3 * N_Ref_Max;
    long* RcvHolidays_Excel = Holidays_Excel;
    n = 0;
    for (i = 0; i < N_Ref_Max; i++)
    {
        n = n + RcvHolidayCount[i];
    }
    long* PayHolidays_Excel = Holidays_Excel + n;
    long** RcvHolidayDays = (long**)malloc(sizeof(long*) * RcvFlag[0]);
    long** PayHolidayDays = (long**)malloc(sizeof(long*) * PayFlag[0]);
    long* RcvHolidayCountAdj = (long*)malloc(sizeof(long) * N_Ref_Max);
    long* PayHolidayCountAdj = (long*)malloc(sizeof(long) * N_Ref_Max);
    for (i = 0; i < N_Ref_Max; i++)
    {
        RcvHolidayCountAdj[i] = RcvHolidayCount[i];
    }
    for (i = 0; i < N_Ref_Max; i++)
    {
        PayHolidayCountAdj[i] = PayHolidayCount[i];
    }


    k = 0;
    n = 0;
    for (i = 0; i < RcvFlag[0]; i++)
    {
        RcvHolidayDays[i] = (long*)malloc(sizeof(long) * max(1, RcvHolidayCount[i]));
        n = 0;
        for (j = 0; j < RcvHolidayCount[i]; j++)
        {
            if (RcvHolidays_Excel[k] >= RcvCashFlowScheduleExcelDate[0] && RcvHolidays_Excel[k] <= RcvCashFlowScheduleExcelDate[NRcvCashFlow * 2 - 1] + 100)
            {
                RcvHolidayDays[i][n] = (RcvHolidays_Excel[k] - PriceDateExcel);
                n += 1;
            }
            else
            {
                RcvHolidayCountAdj[i] -= 1;
            }
            k += 1;
        }
    }
    k = 0;
    n = 0;
    for (i = 0; i < PayFlag[0]; i++)
    {
        PayHolidayDays[i] = (long*)malloc(sizeof(long) * max(1, PayHolidayCount[i]));
        n = 0;
        for (j = 0; j < PayHolidayCount[i]; j++)
        {
            if (PayHolidays_Excel[k] >= PayCashFlowScheduleExcelDate[0] && PayHolidays_Excel[k] <= PayCashFlowScheduleExcelDate[NPayCashFlow * 2 - 1] + 100)
            {
                PayHolidayDays[i][n] = (PayHolidays_Excel[k] - PriceDateExcel);
                n += 1;
            }
            else
            {
                PayHolidayCountAdj[i] -= 1;
            }

            k += 1;
        }
    }

    //////////////////////////////
    // 옵션 관련 Date
    //////////////////////////////

    long* OptionDateCtype;
    long* DaysOptionDate;
    OptionDateCtype = (long*)malloc(sizeof(long) * max(1, NOption));
    DaysOptionDate = (long*)malloc(sizeof(long) * max(1, NOption));
    for (i = 0; i < NOption; i++)
    {
        OptionDateCtype[i] = ExcelDateToCDate(OptionDateAndFlag[i]);
        DaysOptionDate[i] = DayCountAtoB(PriceDateC, OptionDateCtype[i]);
    }
    long OptionUseFlag;
    if (NOption > 0) OptionUseFlag = 1;
    else OptionUseFlag = 0;

    long CallConditionFlag = OptionDateAndFlag[NOption];
    long OptionType = OptionDateAndFlag[NOption + 1];

    double* OptionStrikeRate = (double*)malloc(sizeof(double) * max(1, NOption));
    for (i = 0; i < NOption; i++) OptionStrikeRate[i] = OptionKAndRange[i];

    double** RangeUp;
    double** RangeDn;

    RangeUp = (double**)malloc(sizeof(double*) * max(1, RcvFlag[0]));
    RangeDn = (double**)malloc(sizeof(double*) * max(1, RcvFlag[0]));
    for (i = 0; i < RcvFlag[0]; i++)
    {
        RangeUp[i] = OptionKAndRange + (2 * i + 1) * NOption;
        RangeDn[i] = OptionKAndRange + (2 * i + 2) * NOption;
    }


    //////////////////////////////
    // 시뮬레이션에 사용되는 커브 번호 Mapping 시작
    //////////////////////////////

    long nCurveUsing = 0;
    long nUsingCurveRcv = RcvFlag[0];
    long nUsingCurvePay = PayFlag[0];

    long nTotalCurve = nUsingCurveRcv + nUsingCurvePay + 2;
    long* RcvPayCurveNumber = (long*)malloc(sizeof(long) * nTotalCurve);
    for (i = 0; i < nUsingCurveRcv; i++)
    {
        RcvPayCurveNumber[i] = RcvRefCurveNumber[i];
    }
    for (i = 0; i < nUsingCurvePay; i++)
    {
        RcvPayCurveNumber[i + nUsingCurveRcv] = PayRefCurveNumber[i];
    }
    RcvPayCurveNumber[nUsingCurveRcv + nUsingCurvePay] = RcvRefCurveNumber[3];
    RcvPayCurveNumber[nUsingCurveRcv + nUsingCurvePay + 1] = PayRefCurveNumber[3];

    bubble_sort_long(RcvPayCurveNumber, nTotalCurve, 1);

    long nSimulateCurve = 0;
    long* SimulateCurveIdx = Make_Unique_Array(nTotalCurve, RcvPayCurveNumber, nSimulateCurve);

    //////////////////////////////
    // 각 Leg의 정보 Mapping
    //////////////////////////////

    LEG_INFO* RcvLeg = new LEG_INFO;
    LEG_INFO* PayLeg = new LEG_INFO;
    RcvLeg->NReference = nUsingCurveRcv;
    PayLeg->NReference = nUsingCurvePay;
    RcvLeg->FixFloFlag = RcvFlag[1];
    PayLeg->FixFloFlag = PayFlag[1];
    RcvLeg->AccrualFlag = RcvFlag[2];
    PayLeg->AccrualFlag = PayFlag[2];
    RcvLeg->MaxLossY = RcvMaxLossReturn[0];
    PayLeg->MaxLossY = PayMaxLossReturn[0];
    RcvLeg->MaxRetY = RcvMaxLossReturn[1];
    PayLeg->MaxRetY = PayMaxLossReturn[1];
    RcvLeg->DiscCurveNum = RcvRefCurveNumber[3];
    PayLeg->DiscCurveNum = PayRefCurveNumber[3];
    RcvLeg->Reference_Inform = new REFERENCE_INFO[nUsingCurveRcv];
    PayLeg->Reference_Inform = new REFERENCE_INFO[nUsingCurvePay];
    RcvLeg->PayoffStructure = RcvPayoffStructure;
    PayLeg->PayoffStructure = PayPayoffStructure;
    for (i = 0; i < nUsingCurveRcv; i++)
    {
        RcvLeg->Reference_Inform[i].CurveNum = RcvRefCurveNumber[i];
        RcvLeg->Reference_Inform[i].RefRateType = RcvRefRateType[i];
        RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann = (long)RcvSwapInfo[2 * i];
        RcvLeg->Reference_Inform[i].RefSwapMaturity = RcvSwapInfo[2 * i + 1];
        RcvLeg->Reference_Inform[i].DayCountFlag = RcvDayCountFlag[i];
        if (RcvLeg->Reference_Inform[i].DayCountFlag == 0) RcvLeg->Reference_Inform[i].Day1Y = 365.0;
        else RcvLeg->Reference_Inform[i].Day1Y = 360.0;
        RcvLeg->Reference_Inform[i].RefRateCondMultiple = RcvRefCondMultiple[i];
        RcvLeg->Reference_Inform[i].PowerSpreadFlag = RcvRefPowerSpreadFlag[i];
        RcvLeg->Reference_Inform[i].PowerSpreadMat1 = RcvRefPowerSpreadMat[2 * i];
        RcvLeg->Reference_Inform[i].PowerSpreadMat2 = RcvRefPowerSpreadMat[2 * i + 1];
        RcvLeg->Reference_Inform[i].RangeUp = RcvRefRange[2 * i];
        RcvLeg->Reference_Inform[i].RangeDn = RcvRefRange[2 * i + 1];
        RcvLeg->Reference_Inform[i].PayoffMultiple = RcvPayoffMultiple[i];
    }

    RcvLeg->NCashFlow = NRcvCashFlow;
    RcvLeg->ForwardStart_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->ForwardEnd_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->DaysForwardStart = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->DaysForwardEnd = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->FractionStart_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->FractionEnd_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->PayDate_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->DaysPayDate_C = (long*)malloc(sizeof(long) * NRcvCashFlow);
    RcvLeg->CouponRate = (double*)malloc(sizeof(double) * NRcvCashFlow);
    RcvLeg->FixedRate = (double*)malloc(sizeof(double) * NRcvCashFlow);
    RcvLeg->RangeCoupon = (double*)malloc(sizeof(double) * NRcvCashFlow);

    RcvLeg->OptionUseFlag = OptionUseFlag;
    RcvLeg->NOption = NOption;
    RcvLeg->DaysOptionDate = DaysOptionDate;
    RcvLeg->CallConditionFlag = CallConditionFlag;
    RcvLeg->OptionType = OptionType;
    RcvLeg->StrikeRate = OptionStrikeRate;
    RcvLeg->RangeUp = RangeUp;
    RcvLeg->RangeDn = RangeDn;

    for (i = 0; i < NRcvCashFlow; i++)
    {
        RcvLeg->ForwardStart_C[i] = ExcelDateToCDate(RcvCashFlowScheduleExcelDate[i]);
        RcvLeg->ForwardEnd_C[i] = ExcelDateToCDate(RcvCashFlowScheduleExcelDate[i + NRcvCashFlow]);
        RcvLeg->DaysForwardStart[i] = DayCountAtoB(PriceDateC, RcvLeg->ForwardStart_C[i]);
        RcvLeg->DaysForwardEnd[i] = DayCountAtoB(PriceDateC, RcvLeg->ForwardEnd_C[i]);
        if (RcvLeg->DaysForwardStart[i] == RcvLeg->DaysForwardEnd[i]) RcvLeg->DaysForwardEnd[i] = DayPlus(RcvLeg->DaysForwardEnd[i], 1);
        RcvLeg->FractionStart_C[i] = ExcelDateToCDate(RcvCashFlowScheduleExcelDate[i + 2 * NRcvCashFlow]);
        RcvLeg->FractionEnd_C[i] = ExcelDateToCDate(RcvCashFlowScheduleExcelDate[i + 3 * NRcvCashFlow]);
        RcvLeg->PayDate_C[i] = ExcelDateToCDate(RcvCashFlowScheduleExcelDate[i + 4 * NRcvCashFlow]);
        RcvLeg->DaysPayDate_C[i] = DayCountAtoB(PriceDateC, RcvLeg->PayDate_C[i]);

        RcvLeg->CouponRate[i] = RcvCouponFixRate[i];
        RcvLeg->FixedRate[i] = RcvCouponFixRate[i + NRcvCashFlow];
        RcvLeg->RangeCoupon[i] = RcvCouponFixRate[i + 2 * NRcvCashFlow];
    }

    for (i = 0; i < nUsingCurvePay; i++)
    {
        PayLeg->Reference_Inform[i].CurveNum = PayRefCurveNumber[i];
        PayLeg->Reference_Inform[i].RefRateType = PayRefRateType[i];
        PayLeg->Reference_Inform[i].RefSwapNCPN_Ann = (long)PaySwapInfo[2 * i];
        PayLeg->Reference_Inform[i].RefSwapMaturity = PaySwapInfo[2 * i + 1];
        PayLeg->Reference_Inform[i].DayCountFlag = PayDayCountFlag[i];
        if (PayLeg->Reference_Inform[i].DayCountFlag == 0) PayLeg->Reference_Inform[i].Day1Y = 365.0;
        else PayLeg->Reference_Inform[i].Day1Y = 360.0;
        PayLeg->Reference_Inform[i].RefRateCondMultiple = PayRefCondMultiple[i];
        PayLeg->Reference_Inform[i].PowerSpreadFlag = PayRefPowerSpreadFlag[i];
        PayLeg->Reference_Inform[i].PowerSpreadMat1 = PayRefPowerSpreadMat[2 * i];
        PayLeg->Reference_Inform[i].PowerSpreadMat2 = PayRefPowerSpreadMat[2 * i + 1];
        PayLeg->Reference_Inform[i].RangeUp = PayRefRange[2 * i];
        PayLeg->Reference_Inform[i].RangeDn = PayRefRange[2 * i + 1];
        PayLeg->Reference_Inform[i].PayoffMultiple = PayPayoffMultiple[i];
    }

    PayLeg->NCashFlow = NPayCashFlow;
    PayLeg->ForwardStart_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->ForwardEnd_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->DaysForwardStart = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->DaysForwardEnd = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->FractionStart_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->FractionEnd_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->PayDate_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->DaysPayDate_C = (long*)malloc(sizeof(long) * NPayCashFlow);
    PayLeg->CouponRate = (double*)malloc(sizeof(double) * NPayCashFlow);
    PayLeg->FixedRate = (double*)malloc(sizeof(double) * NPayCashFlow);
    PayLeg->RangeCoupon = (double*)malloc(sizeof(double) * NPayCashFlow);
    for (i = 0; i < NPayCashFlow; i++)
    {
        PayLeg->ForwardStart_C[i] = ExcelDateToCDate(PayCashFlowScheduleExcelDate[i]);
        PayLeg->ForwardEnd_C[i] = ExcelDateToCDate(PayCashFlowScheduleExcelDate[i + NPayCashFlow]);
        PayLeg->DaysForwardStart[i] = DayCountAtoB(PriceDateC, PayLeg->ForwardStart_C[i]);
        PayLeg->DaysForwardEnd[i] = DayCountAtoB(PriceDateC, PayLeg->ForwardEnd_C[i]);
        if (PayLeg->DaysForwardStart[i] == PayLeg->DaysForwardEnd[i]) PayLeg->DaysForwardEnd[i] = DayPlus(PayLeg->DaysForwardEnd[i], 1);
        PayLeg->FractionStart_C[i] = ExcelDateToCDate(PayCashFlowScheduleExcelDate[i + 2 * NPayCashFlow]);
        PayLeg->FractionEnd_C[i] = ExcelDateToCDate(PayCashFlowScheduleExcelDate[i + 3 * NPayCashFlow]);
        PayLeg->PayDate_C[i] = ExcelDateToCDate(PayCashFlowScheduleExcelDate[i + 4 * NPayCashFlow]);
        PayLeg->DaysPayDate_C[i] = DayCountAtoB(PriceDateC, PayLeg->PayDate_C[i]);
        PayLeg->CouponRate[i] = PayCouponFixRate[i];
        PayLeg->FixedRate[i] = PayCouponFixRate[i + NPayCashFlow];
        PayLeg->RangeCoupon[i] = PayCouponFixRate[i + 2 * NPayCashFlow];
    }

    RcvLeg->nUsingCurve = nUsingCurveRcv;
    PayLeg->nUsingCurve = nUsingCurvePay;
    RcvLeg->HolidayCalcFlag = RcvHolidayCalcFlag;
    PayLeg->HolidayCalcFlag = PayHolidayCalcFlag;
    RcvLeg->HolidayCount = RcvHolidayCountAdj;
    PayLeg->HolidayCount = PayHolidayCountAdj;
    RcvLeg->HolidayDays = RcvHolidayDays;
    PayLeg->HolidayDays = PayHolidayDays;

    PayLeg->OptionUseFlag = 0;
    PayLeg->NOption = 0;

    double* CorrMatReshaped = (double*)malloc(sizeof(double) * (nSimulateCurve * nSimulateCurve));
    k = 0;
    n = 0;
    for (i = 0; i < N_Curve_Max; i++)
    {
        for (j = 0; j < N_Curve_Max; j++)
        {
            if (isin(i + 1, SimulateCurveIdx, N_Curve_Max) && isin(j + 1, SimulateCurveIdx, N_Curve_Max))
            {
                CorrMatReshaped[n] = Corr[k];
                n = n + 1;
            }
            k += 1;
        }
    }

    double** CorrMatAdj = (double**)malloc(sizeof(double*) * nSimulateCurve);
    for (i = 0; i < nSimulateCurve; i++)
    {
        CorrMatAdj[i] = CorrMatReshaped + i * nSimulateCurve;
    }
    ////////////////////////////////////////////////////////////
    // 시뮬레이션에 필요한 날짜 Node 결정 : 두 Leg의 ForwardStartDate, ForwardEndDate를 Unique Sort하자.
    ////////////////////////////////////////////////////////////

    long NDays = 0;
    long* AllDays = (long*)malloc(sizeof(long) * (RcvLeg->NCashFlow * 2 + PayLeg->NCashFlow * 2));
    for (i = 0; i < RcvLeg->NCashFlow; i++)
    {
        if (PriceDateC <= RcvLeg->ForwardStart_C[i]) AllDays[2 * i] = DayCountAtoB(PriceDateC, RcvLeg->ForwardStart_C[i]);
        else AllDays[2 * i] = 0;

        if (PriceDateC <= RcvLeg->ForwardEnd_C[i]) AllDays[2 * i + 1] = DayCountAtoB(PriceDateC, RcvLeg->ForwardEnd_C[i]);
        else AllDays[2 * i + 1] = 0;
    }
    for (i = 0; i < PayLeg->NCashFlow; i++)
    {
        if (PriceDateC <= PayLeg->ForwardStart_C[i]) AllDays[RcvLeg->NCashFlow * 2 + i] = DayCountAtoB(PriceDateC, PayLeg->ForwardStart_C[i]);
        else AllDays[RcvLeg->NCashFlow * 2 + i] = 0;

        if (PriceDateC <= PayLeg->ForwardEnd_C[i]) AllDays[RcvLeg->NCashFlow * 2 + PayLeg->NCashFlow + i] = DayCountAtoB(PriceDateC, PayLeg->ForwardEnd_C[i]);
        else AllDays[RcvLeg->NCashFlow * 2 + PayLeg->NCashFlow + i] = 0;
    }
    bubble_sort_long(AllDays, (RcvLeg->NCashFlow * 2 + PayLeg->NCashFlow * 2), 1);
    long* DaysForSimul;
    DaysForSimul = Make_Unique_Array((RcvLeg->NCashFlow * 2 + PayLeg->NCashFlow * 2), AllDays, NDays);

    long MaxDaysSimul;
    if (DailySimulFlag == 0) MaxDaysSimul = NDays;
    else MaxDaysSimul = DaysForSimul[NDays - 1];

    if (DailySimulFlag != 0)
    {
        free(DaysForSimul);
        DaysForSimul = (long*)malloc(sizeof(long) * MaxDaysSimul);
        for (i = 0; i < MaxDaysSimul; i++)
        {
            DaysForSimul[i] = i;
        }
    }

    //////////////////////////////////
    // Random Number Generate
    //////////////////////////////////

    long seed = 0;
    randnorm(seed);
    double** FixedRandn = random_mvrn(NSimul * MaxDaysSimul, nSimulateCurve, CorrMatAdj);

    /////////////////
    // Rate History
    /////////////////

    long** RateHistoryDateMatrixRcv = (long**)malloc(sizeof(long*) * nUsingCurveRcv);
    double** RateHistoryMatrixRcv = (double**)malloc(sizeof(double*) * nUsingCurveRcv);
    k = 0;
    n = 0;
    for (i = 0; i < nUsingCurveRcv; i++)
    {
        RateHistoryDateMatrixRcv[i] = (long*)malloc(sizeof(long) * NDayHistRcv[i]);
        RateHistoryMatrixRcv[i] = RateHistRcv + n;
        for (j = 0; j < NDayHistRcv[i]; j++)
        {
            RateHistoryDateMatrixRcv[i][j] = RateHistDateRcv[k] - PriceDateExcel;
            k += 1;
        }
        n = n + NDayHistRcv[i];
    }

    long** RateHistoryDateMatrixPay = (long**)malloc(sizeof(long*) * nUsingCurvePay);
    double** RateHistoryMatrixPay = (double**)malloc(sizeof(double*) * nUsingCurvePay);
    k = 0;
    n = 0;
    for (i = 0; i < nUsingCurvePay; i++)
    {
        RateHistoryDateMatrixPay[i] = (long*)malloc(sizeof(long) * NDayHistPay[i]);
        RateHistoryMatrixPay[i] = RateHistPay + n;
        for (j = 0; j < NDayHistPay[i]; j++)
        {
            RateHistoryDateMatrixPay[i][j] = RateHistDatePay[k] - PriceDateExcel;
            k += 1;
        }
        n = n + NDayHistPay[i];
    }
    RcvLeg->NDayHistory = NDayHistRcv;
    PayLeg->NDayHistory = NDayHistPay;
    RcvLeg->RateHistoryDateMatrix = RateHistoryDateMatrixRcv;
    PayLeg->RateHistoryDateMatrix = RateHistoryDateMatrixPay;
    RcvLeg->RateHistoryMatrix = RateHistoryMatrixRcv;
    PayLeg->RateHistoryMatrix = RateHistoryMatrixPay;

    //////////////////////
    // 시뮬레이션 정보  //
    //////////////////////

    SIMUL_INFO* Simul = new SIMUL_INFO;
    Simul->NSimul = NSimul;
    Simul->DailySimulFlag = DailySimulFlag;
    Simul->NDays = MaxDaysSimul;
    Simul->NAsset = nSimulateCurve;
    Simul->DaysForSimul = DaysForSimul;

    Simul->dt_Array = (double*)malloc(sizeof(double) * Simul->NDays);
    Simul->T_Array = (double*)malloc(sizeof(double) * Simul->NDays);
    for (i = 0; i < Simul->NDays; i++)
    {
        if (i == 0)
        {
            Simul->dt_Array[0] = ((double)DaysForSimul[0]) / 365.0;
        }
        else
        {
            Simul->dt_Array[i] = ((double)(DaysForSimul[i] - DaysForSimul[i - 1])) / 365.0;
        }
        Simul->T_Array[i] = ((double)DaysForSimul[i]) / 365.0;
    }

    Simul->FixedRandn = FixedRandn;
    Simul->NHWVol = (long*)malloc(sizeof(long) * Simul->NAsset);
    Simul->HWKappa = (double*)malloc(sizeof(double) * Simul->NAsset);
    Simul->HWVolTerm = (double**)malloc(sizeof(double*) * Simul->NAsset);
    Simul->HWVol = (double**)malloc(sizeof(double*) * Simul->NAsset);
    Simul->NRateTerm = (long*)malloc(sizeof(long) * Simul->NAsset);
    Simul->RateTerm = (double**)malloc(sizeof(double*) * Simul->NAsset);
    Simul->Rate = (double**)malloc(sizeof(double*) * Simul->NAsset);
    Simul->SimulCurveIdx = SimulateCurveIdx;
    k = 0;
    for (i = 0; i < N_Curve_Max; i++)
    {
        if (isin(i + 1, SimulateCurveIdx, nSimulateCurve))
        {
            Simul->NHWVol[k] = NHWVolCount[i];
            Simul->HWVolTerm[k] = HWVolTermMatrix[i];
            Simul->HWVol[k] = HWVolMatrix[i];
            Simul->HWKappa[k] = HWKappaArray[i];
            Simul->NRateTerm[k] = NZeroRate[i];
            Simul->RateTerm[k] = ZeroRateTermMatrix[i];
            Simul->Rate[k] = ZeroRateMatrix[i];
            k = k + 1;
        }
    }

    ResultCode = IRStructuredSwap(PriceDateC, NAFlag, Notional, RcvLeg, PayLeg, Simul, GreekFlag, ResultPrice, ResultRcv, ResultPay);


    free(HWVolTermMatrix);
    free(HWVolMatrix);
    free(ZeroRateTermMatrix);
    free(ZeroRateMatrix);

    for (i = 0; i < RcvFlag[0]; i++)
    {
        free(RcvHolidayDays[i]);
    }
    free(RcvHolidayDays);

    for (i = 0; i < PayFlag[0]; i++)
    {
        free(PayHolidayDays[i]);
    }
    free(PayHolidayDays);

    free(OptionDateCtype);
    free(OptionStrikeRate);
    free(RangeUp);
    free(RangeDn);
    free(DaysOptionDate);
    free(RcvPayCurveNumber);
    free(SimulateCurveIdx);

    free(RcvLeg->ForwardStart_C);
    free(RcvLeg->ForwardEnd_C);
    free(RcvLeg->DaysForwardStart);
    free(RcvLeg->DaysForwardEnd);
    free(RcvLeg->FractionStart_C);
    free(RcvLeg->FractionEnd_C);
    free(RcvLeg->PayDate_C);
    free(RcvLeg->DaysPayDate_C);
    free(RcvLeg->CouponRate);
    free(RcvLeg->FixedRate);
    free(RcvLeg->RangeCoupon);

    free(PayLeg->ForwardStart_C);
    free(PayLeg->ForwardEnd_C);
    free(PayLeg->DaysForwardStart);
    free(PayLeg->DaysForwardEnd);
    free(PayLeg->FractionStart_C);
    free(PayLeg->FractionEnd_C);
    free(PayLeg->PayDate_C);
    free(PayLeg->DaysPayDate_C);
    free(PayLeg->CouponRate);
    free(PayLeg->FixedRate);
    free(PayLeg->RangeCoupon);
    free(RcvHolidayCountAdj);
    free(PayHolidayCountAdj);
    free(CorrMatReshaped);
    free(CorrMatAdj);
    free(AllDays);
    free(DaysForSimul);

    for (i = 0; i < NSimul * MaxDaysSimul; i++)
    {
        free(FixedRandn[i]);
    }
    free(FixedRandn);
    delete (RcvLeg->Reference_Inform);
    delete (PayLeg->Reference_Inform);
    delete (RcvLeg);
    delete (PayLeg);

    for (i = 0; i < nUsingCurveRcv; i++)
    {
        free(RateHistoryDateMatrixRcv[i]);
    }
    free(RateHistoryDateMatrixRcv);
    free(RateHistoryMatrixRcv);

    for (i = 0; i < nUsingCurvePay; i++)
    {
        free(RateHistoryDateMatrixPay[i]);
    }
    free(RateHistoryDateMatrixPay);
    free(RateHistoryMatrixPay);

    free(Simul->dt_Array);
    free(Simul->T_Array);
    free(Simul->NHWVol);
    free(Simul->HWKappa);
    free(Simul->HWVolTerm);
    free(Simul->HWVol);
    free(Simul->NRateTerm);
    free(Simul->RateTerm);
    free(Simul->Rate);
    delete(Simul);
    _CrtDumpMemoryLeaks();
    return 1;
}