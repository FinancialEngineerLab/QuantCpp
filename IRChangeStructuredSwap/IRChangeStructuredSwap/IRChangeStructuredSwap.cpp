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

long isin(long x, long* array, long narray, long& startidx)
{
    long i;
    long s = 0;
    for (i = startidx; i < narray; i++)
    {
        if (x == array[i])
        {
            s = 1;
            startidx = i;
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

long isinFindIndex(long x, long* array, long narray, long& Idx, long* startidx)
{
    long i;
    long ii = *startidx;
    long s = 0;
    Idx = -1;
    for (i = ii; i < narray; i++)
    {
        if (x == array[i])
        {
            *startidx = i;
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

long isweekend(long ExlDate)
{
    // 나머지 1이면 일요일, 2이면 월요일, 3이면 화요일, 4이면 수요일, 5이면 목요일, 6이면 금요일, 0이면 토요일
    long MOD7;
    if (ExlDate > 0)
    {
        MOD7 = ExlDate % 7;
        if (MOD7 == 1 || MOD7 == 0) return 1;
    }
    return 0;
}

void LockOutCheck(long& LockOutFlag, long LockOutDay, long Today, double& LockOutDayRate, double& ForwardRate)
{
    if (LockOutFlag == 1) ForwardRate = LockOutDayRate + 0.0;
    else
    {
        if (LockOutDay <= Today)
        {
            LockOutDayRate = ForwardRate + 0.0;
            LockOutFlag = 1;
        }
    }
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

    if (i == narray) i = -1;

    return i;
}

long FindIndex(
    long x,
    long* array,
    long narray,
    long* loopstart
)
{
    long i;
    long idxstart = *loopstart;
    for (i = max(0, idxstart); i < narray; i++)
    {
        if (x == array[i])
        {
            *loopstart = i;
            return i;
        }
    }

    if (i == narray) i = -1;

    return i;
}

double sumation_array(long n, double* myarray)
{
    double res = 0;
    for (long i = 0; i < n; i++) res += myarray[i];
    return res;
}

typedef struct HW_Parameter {
    long NAsset;                            // 시뮬레이션 자산 개수
    long NDays;                             // 시뮬레이션 Day 개수

    double** XA;                            // Hull White XA Term = exp(-kappa * (T2-T1)) --> Shape : (NAsset, NDays)
    double** XV;                            // Hull White XV Term = sqrt( integral_(t1)^(t2) exp(-2kappa(t2-t1)) sigma_t^2 dt ) --> Shape : (NAsset, NDays)

    long CurveIdx_DiscRcv;                  // 시뮬레이션 기준 Rcv Discount 커브 index
    long CurveIdx_DiscPay;                  // 시뮬레이션 기준 Pay Discount 커브 index
    long NRcv_Cf;                           // Rcv Leg Cashflow 개수
    long NPay_Cf;                           // Pay Leg Cashflow 개수
    double* Rcv_DF_0_t;                     // Rcv Disc P(0,t)
    double* Rcv_DF_0_T;                     // Rcv Disc P(0,T)
    double* Rcv_DF_t_T;                     // Rcv Disc P(t,T)
    double* B_t_T_RcvDisc;                  // Rcv disc Hull White B(t,T) = (1-exp(-kappa(T-t))/kappa
    double* QVTerm_RcvDisc;                 // Rcv disc Hull White QVTerm(t,T) = integ_0^t sigma^2 * [B(s,t)^2 - B(s,T)^2] ds

    double* Pay_DF_0_t;                     // Pay Disc P(0,t)
    double* Pay_DF_0_T;                     // Pay Disc P(0,T)
    double* Pay_DF_t_T;                     // Pay Disc P(t,T)
    double* B_t_T_PayDisc;                  // Pay Disc Hull White B(t,T) = (1-exp(-kappa(T-t))/kappa
    double* QVTerm_PayDisc;                 // Pay Disc Hull White QVTerm(t,T) = integ_0^t sigma^2 * [B(s,t)^2 - B(s,T)^2] ds

    long NRcvCurve;                         // Rcv Leg Reference Curve 개수
    long NPayCurve;                         // Pay Leg Reference Curve 개수

    long* ndates_cpn_rcv;                   // Rcv Reference 금리별 쿠폰개수
    double*** RcvRef_DF_t_T;                // Rcv Ref P(t,T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_B_t_T;                 // Rcv Ref B(t,T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_QVTerm;                // Rcv Ref QVTerm(t, T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_dt;                    // Rcv Ref DeltaT(t, T)--> Shape = NReference, NDays, NCpn

    long* ndates_cpn_powerspread_rcv;       // Rcv Reference 반대쪽 파워스프레드 금리별 쿠폰개수
    double*** RcvRef_DF_t_T_PowerSpread;    // Rcv Ref 반대쪽 파워스프레드 P(t,T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_B_t_T_PowerSpread;     // Rcv Ref 반대쪽 파워스프레드 B(t,T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_QVTerm_PowerSpread;    // Rcv Ref 반대쪽 파워스프레드 QVTerm(t, T) --> Shape = NReference, NDays, NCpn
    double*** RcvRef_dt_PowerSpread;        // Rcv Ref 반대쪽 파워스프레드 DeltaT(t, T)--> Shape = NReference, NDays, NCpn

    long* ndates_cpn_pay;                   // Pay Reference 금리별 쿠폰개수
    double*** PayRef_DF_t_T;                // Pay Ref P(t,T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_B_t_T;                 // Pay Ref B(t,T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_QVTerm;                // Pay Ref QVTerm(t, T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_dt;                    // Pay Ref DeltaT(t, T)--> Shape = NReference, NDays, NCpn

    long* ndates_cpn_powerspread_pay;       // Pay Reference 반대쪽 파워스프레드 금리별 쿠폰개수
    double*** PayRef_DF_t_T_PowerSpread;    // Pay Ref 반대쪽 파워스프레드 P(t,T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_B_t_T_PowerSpread;     // Pay Ref 반대쪽 파워스프레드 B(t,T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_QVTerm_PowerSpread;    // Pay Ref 반대쪽 파워스프레드 QVTerm(t, T) --> Shape = NReference, NDays, NCpn
    double*** PayRef_dt_PowerSpread;        // Pay Ref 반대쪽 파워스프레드 DeltaT(t, T)--> Shape = NReference, NDays, NCpn

} HW_INFO;

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
    long CurrentIdx;

    long nUsingCurve;                   // 사용되는커브숫자
    long* HolidayCalcFlag;              // Holiday계산Flag
    long* HolidayCount;                 // Holiday수
    long** HolidayDays;                 // 평가일 to 각 레퍼런스별 Holiday까지 날짜 수

    long* NDayHistory;                  // 과거 레퍼런스 History 개수
    long** RateHistoryDateMatrix;       // 과거 레퍼런스 History Date Matrix
    double** RateHistoryMatrix;         // 과거 레퍼런스 History Rate Matrix

    long* NWeekend;                     // CashFlow 시점별 주말개수
    long** Weekend;                     // 토, 일 Array shape = NCashFlow * NWeekend[i]
    double* SOFR_Annualized_R0;         // 첫 Cpn의 Annualized Rate
    double* SOFR_Compound0;             // 첫 Cpn의 Compound Value

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
    long NSimul;                    // 시뮬레이션 개수
    long NDays;                     // 시뮬레이션 날짜 개수
    long NAsset;                    // 기초자산 중복 제외하고 개수
    long DailySimulFlag;            // DailySimul여부
    long* DaysForSimul;             // 시뮬 날짜 Array
    double* dt_Array;               // 시뮬레이션 Day들의 간격
    double* T_Array;                // 시뮬레이션 Day를 T로 변환
    double** FixedRandn;            // 고정된랜덤넘버
    long* NHWVol;                   // HWVol 개수를 NAsset만큼 배열
    double* HWKappa;                // HW Kappa값을 NAsset만큼 배열
    double** HWVolTerm;             // HW Vol Term을 --> NAsset, NHWVol 만큼 배열
    double** HWVol;                 // HW Vol을 --> NAsset, NHWVol 만큼 배열
    long* NRateTerm;                // 기초금리 Term 개수 Array
    double** RateTerm;              // 기초금리 Term Array --> Shape = NAsset, NRateTerm
    double** Rate;                  // 기초금리 Array --> Shape = NAsset, NRateTerm
    long* SimulCurveIdx;            // 시뮬레이션되는 커브 번호
    long SOFRBusinessDaySimul;      // SOFR로 인해 BusinessDay Simul되는지 여부
}SIMUL_INFO;

// Linear Interpolation (X변수, Y변수, X길이, 타겟X)
double Interpolate_Linear_Pointer(double* x, double* fx, long nx, double targetx, long* idx)
{
    long i;
    long startidx = *idx;
    double result = 0.0;

    if (nx == 1 || targetx == x[0])
        return fx[0];
    else if (targetx == x[nx - 1])
        return fx[nx - 1];


    if (targetx < x[0]) return fx[0];
    else if (targetx > x[nx - 1]) return fx[nx - 1];
    else
    {
        for (i = max(1, startidx); i < nx; i++)
        {
            if (targetx < x[i])
            {
                *idx = i - 1;
                result = (fx[i] - fx[i - 1]) / (x[i] - x[i - 1]) * (targetx - x[i - 1]) + fx[i - 1];
                break;
            }
        }
        return result;
    }
}

double SOFR_ForwardRate_Compound(
    long NRefCrvTerm,
    double* RefCrvTerm,
    double* RefCrvRate,
    long ForwardStartIdx,
    long ForwardEndIdx,
    long LockOutDays,
    long LookBackDays,
    long ObservShift,
    long HolidayFlag,
    long NHoliday,
    long* Holiday,
    long NSaturSunDay,
    long* SaturSunDay,
    long UseHistorySOFR,
    long NRefHistory,
    long* RefHistoryDate,
    double* RefHistory,
    double denominator,
    double& AnnualizedOISRate,
    long SOFRUseFlag
)
{
    long i;
    long j;
    long k;
    long n;

    long nday = ForwardEndIdx - ForwardStartIdx;
    long isinflag = 0;
    long HistDayIdx = 0;
    long AverageFlag = 0;
    if (SOFRUseFlag == 3) AverageFlag = 1;

    // 빠른 Holiday찾기용 TimePos Pointer
    long HoliStartIdx = 0;
    long HoliStartIdx2 = 0;
    long SatSunIdx = 0;
    long SatSunIdx2 = 0;

    long isHolidayFlag = 0;
    long isHolidayFlag2 = 0;

    long CountHoliday = 0;
    long CountHoliday2 = 0;

    long TimePos = 0;
    long TimePos2 = 0;

    double dt = 1.0 / 365.0;
    double Prev_PI;
    double ForwardRate;

    double t = 0.0;
    double PI_0;
    double T;
    T = (double)(ForwardEndIdx - ForwardStartIdx) / denominator;

    double CurrentRate;
    CurrentRate = Interpolate_Linear(RefCrvTerm, RefCrvRate, NRefCrvTerm, dt);

    long ObservShiftFlag = 0;
    if (LookBackDays > 0 && ObservShift > 0)     ObservShiftFlag = 1;

    ///////////////////////////
    // Average 키움증권 추가 //
    ///////////////////////////
    long NCumpound = 0;
    double AverageRate = 0.0;

    ///////////////////////////
    // N영업일 전 LockOutDay 계산
    ///////////////////////////
    long LockOutDay = ForwardEndIdx;
    long LockOutFlag = 0;

    if (LockOutDays > 0)
    {
        k = 0;
        for (i = 1; i < 30; i++)
        {
            LockOutDay -= 1;
            // N영업일까지 날짜 뒤로가기
            if (max(isin(LockOutDay, Holiday, NHoliday), isin(LockOutDay, SaturSunDay, NSaturSunDay)) == 0) k += 1;
            if (k == LockOutDays) break;
        }
    }
    double LockOutDayRate = 0.0;


    ///////////////////////////
    // Holiday Rate가 Interpolate일 때 사용할 변수
    ///////////////////////////
    double TargetRate[2] = { 0.0,0.0 };
    double TargetT[2] = { 0.0,0.0 };

    long lookbackday;
    Prev_PI = 1.0;
    PI_0 = 1.0;
    if (UseHistorySOFR == 1 && ObservShiftFlag == 0)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Forward Start 날짜가 가격 계산일보다 앞에 있을 경우 History Rate 뒤져봐야함 + Observe Shift 적용 안하는 경우 //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        HistDayIdx = 0;
        if (LookBackDays < 1)
        {
            ////////////////////////////
            // LookBack 적용안할 경우 //
            ////////////////////////////
            if (ForwardEndIdx >= 0)
            {
                if (HolidayFlag == 3)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        isinflag= isinFindIndex(i, RefHistoryDate, NRefHistory, HistDayIdx);

                        if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                        else ForwardRate = CurrentRate;

                        Prev_PI *= 1.0 + ForwardRate * 1.0 / denominator;
                        AverageRate += ForwardRate;

                        NCumpound += 1;
                    }
                }
                else if (HolidayFlag == 0)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        if (i == ForwardStartIdx) isHolidayFlag = 0;
                        else isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                        if (isHolidayFlag == 0)
                        {
                            isinflag= isinFindIndex(i, RefHistoryDate, NRefHistory, HistDayIdx);

                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = i + 1; j < 0; j++)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }
                            }
                            Prev_PI *= 1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator;
                            AverageRate += ForwardRate;

                            NCumpound += 1;
                        }

                    }
                }
                else if (HolidayFlag == 1)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        if (i == ForwardStartIdx) isHolidayFlag = 0;
                        else isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                        if (isHolidayFlag == 0)
                        {
                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = i + 1; j < ForwardEndIdx; j++)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }
                            }
                            isinflag = isinFindIndex(i + CountHoliday + 1, RefHistoryDate, NRefHistory, HistDayIdx);

                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            Prev_PI *= 1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator;
                            AverageRate += ForwardRate;

                            NCumpound += 1;
                        }
                    }
                }
                if (HolidayFlag == 2)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        if (i == ForwardStartIdx) isHolidayFlag = 0;
                        else isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                        if (isHolidayFlag == 0)
                        {
                            ///////////////////////////////
                            // 영업일인 경우 당일 이자율 //
                            ///////////////////////////////
                            isinflag= isinFindIndex(i, RefHistoryDate, NRefHistory, HistDayIdx);
                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;
                            CountHoliday = 0;
                        }
                        else if (isHolidayFlag == 1 && NRefHistory > 0)
                        {
                            ///////////////////////////////////////////////////////////
                            // 영업일이 아닌 경우 직전, 직후 영업일 ForwardRate 찾기 //
                            ///////////////////////////////////////////////////////////

                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = j - 1; j < ForwardStartIdx; j--)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }
                            }
                            isinflag= isinFindIndex(i - CountHoliday - 1, RefHistoryDate, NRefHistory, HistDayIdx);
                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            TargetRate[0] = ForwardRate;
                            TargetT[0] = (-(double)(CountHoliday)-1.0);

                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = i + 1; j < ForwardEndIdx; j++)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }

                            }
                            isinflag = isinFindIndex(i + CountHoliday + 1, RefHistoryDate, NRefHistory, HistDayIdx);
                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            TargetRate[1] = ForwardRate;
                            TargetT[1] = ((double)(CountHoliday)+1.0);
                            ForwardRate = Interpolate_Linear(TargetT, TargetRate, 2, 0.0);
                        }
                        else ForwardRate = CurrentRate;

                        Prev_PI *= 1.0 + ForwardRate * 1.0 / denominator;
                        AverageRate += ForwardRate;

                        NCumpound += 1;
                    }
                }
            }
            ////////////////////////////////////////////////////////
            // 과거 반영 다 끝났으면 ForwardStartIdx 0으로 바꾸기 //
            ////////////////////////////////////////////////////////
            ForwardStartIdx = 0;
        }
        else if (LookBackDays > 0)
        {
            //////////////////////////
            // LookBack 적용할 경우 //
            //////////////////////////
            if (ForwardEndIdx >= 0)
            {
                if (HolidayFlag == 3)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        lookbackday = i;

                        k = 0;

                        for (n = 1; n < 120; n++)
                        {
                            ///////////////////////////
                            // N영업일 전 LookBackDay 계산
                            ///////////////////////////
                            lookbackday -= 1;
                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        isinflag = isinFindIndex(lookbackday, RefHistoryDate, NRefHistory, HistDayIdx);
                        if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                        else ForwardRate = CurrentRate;
                        Prev_PI *= 1.0 + ForwardRate * 1.0 / denominator;
                        AverageRate += ForwardRate;

                        NCumpound += 1;
                    }
                }
                else if (HolidayFlag == 0)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        if (i == ForwardStartIdx) isHolidayFlag = 0;
                        else isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                        if (isHolidayFlag == 0)
                        {
                            lookbackday = i;

                            k = 0;

                            for (n = 1; n < 120; n++)
                            {
                                ///////////////////////////
                                // N영업일 전 LookBackDay 계산
                                ///////////////////////////
                                lookbackday -= 1;

                                if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                                if (k == LookBackDays) break;
                            }
                            isinflag= isinFindIndex(lookbackday, RefHistoryDate, NRefHistory, HistDayIdx);
                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = i + 1; j < ForwardEndIdx; j++)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }

                            }
                            Prev_PI *= 1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator;
                            AverageRate += ForwardRate;

                            NCumpound += 1;
                        }
                    }
                }
                else if (HolidayFlag == 1)
                {
                    for (i = ForwardStartIdx; i < 0; i++)
                    {
                        if (i == ForwardStartIdx) isHolidayFlag = 0;
                        else isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                        if (isHolidayFlag == 0)
                        {
                            CountHoliday = 0;
                            HoliStartIdx = 0;
                            SatSunIdx = 0;

                            for (j = i + 1; j < ForwardEndIdx; j++)
                            {
                                isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                                if (isHolidayFlag2 == 0) break;
                                else
                                {
                                    CountHoliday += 1;
                                }
                            }

                            if (CountHoliday == 0)
                            {
                                lookbackday = i;
                            }
                            else
                            {
                                lookbackday = i + CountHoliday + 1;
                            }

                            k = 0;
                            for (n = 1; n < 120; n++)
                            {
                                ///////////////////////////
                                // N영업일 전 LookBackDay 계산
                                ///////////////////////////
                                lookbackday -= 1;
                                if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                                if (k == LookBackDays) break;
                            }
                            isinflag = isinFindIndex(lookbackday, RefHistoryDate, NRefHistory, HistDayIdx);
                            if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                            else ForwardRate = CurrentRate;

                            Prev_PI *= 1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator;
                            AverageRate += ForwardRate;

                            NCumpound += 1;
                        }
                    }
                }
            }
            ForwardStartIdx = 0;
        }
    }
    else if (UseHistorySOFR == 1 && ObservShiftFlag == 1)
    {
        ///////////////////////////
        // Forward Start 날짜가 가격 계산일보다 앞에 있을 경우
        // 또한 ObservShift가 있는 경우
        ///////////////////////////
        HistDayIdx = 0;
        if (ForwardEndIdx >= 0)
        {
            for (i = ForwardStartIdx; i < 0; i++)
            {
                lookbackday = i;
                isHolidayFlag = max(isin(i, Holiday, NHoliday), isin(i, SaturSunDay, NSaturSunDay));

                if (isHolidayFlag == 0)
                {
                    //////////////////
                    // Holiday가 아니라면 LookBack영업일 Back
                    //////////////////

                    k = 0;
                    if (LookBackDays > 0)
                    {
                        for (n = 1; n < 120; n++)
                        {
                            lookbackday -= 1;
                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                    }
                    isinflag= isinFindIndex(lookbackday, RefHistoryDate, NRefHistory, HistDayIdx);
                    if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                    else ForwardRate = Interpolate_Linear(RefCrvTerm, RefCrvRate, NRefCrvTerm, dt);

                    ///////////////////////////
                    // Observation Dt 계산
                    ///////////////////////////
                    CountHoliday = 0;
                    HoliStartIdx = 0;
                    SatSunIdx = 0;

                    for (j = lookbackday + 1; j < 0; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx));
                        if (isHolidayFlag2 == 0) break;
                        else CountHoliday += 1;
                    }

                    if (CountHoliday == 0) Prev_PI *= 1.0 + ForwardRate * 1.0 / denominator;
                    else Prev_PI *= 1.0 + ForwardRate * ((double)(1 + CountHoliday) / denominator);
                    AverageRate += ForwardRate;

                    NCumpound += 1;
                }
            }
        }
        ForwardStartIdx = 0;
    }

    PI_0 = Prev_PI;
    HoliStartIdx = 0;
    HoliStartIdx2 = 0;
    TimePos = 0;
    TimePos2 = 0;
    LockOutFlag = 0;
    NCumpound = 0;
    SatSunIdx = 0;
    SatSunIdx2 = 0;

    if (ObservShiftFlag == 0)
    {
        if (HolidayFlag == 0)
        {
            //////////////////
            // HolidayRate Forward Fill
            //////////////////
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                //////////////////////////
                // 첫 날은 영업일로 처리
                //////////////////////////
                if (i == ForwardStartIdx)
                {
                    if (UseHistorySOFR == 0) isHolidayFlag = 0;
                    else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));
                }
                else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    lookbackday = i + 0;

                    ////////////////////
                    // Forward Rate Time 추정 LookBack 반영
                    ////////////////////
                    if (LookBackDays < 1) t = ((double)i) / denominator;
                    else
                    {
                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;

                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }

                    ////////////////////
                    // 해당일 이후 Holiday개수 체크하여 dt설정
                    ////////////////////
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;
                    for (j = i + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday += 1;
                        }
                    }

                    if (CountHoliday == 0)
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);
                    }
                    else
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos, CountHoliday);
                    }

                    LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator);
                    AverageRate += ForwardRate;

                    NCumpound += 1;
                }
            }
        }
        else if (HolidayFlag == 1)
        {
            //////////////////
            // HolidayRate Backward Fill
            //////////////////
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                //////////////////////////
                // 첫 날은 영업일로 처리
                //////////////////////////
                if (i == ForwardStartIdx)
                {
                    if (UseHistorySOFR == 0) isHolidayFlag = 0;
                    else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));
                }
                else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    ////////////////////
                    // 해당일 이후 Holiday개수 체크하여 dt설정
                    ////////////////////
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    for (j = i + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday += 1;
                        }
                    }

                    if (CountHoliday == 0)
                    {
                        lookbackday = i + 0;
                    }
                    else
                    {
                        lookbackday = i + CountHoliday + 1;
                    }

                    ////////////////////
                    // Forward Rate Time 추정 LookBack 반영
                    ////////////////////
                    if (LookBackDays < 1) t = ((double)lookbackday) / denominator;
                    else
                    {
                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;

                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }

                    CountHoliday2 = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    for (j = lookbackday + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday2 += 1;
                        }
                    }

                    if (CountHoliday2 == 0)
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);
                    }
                    else
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos, CountHoliday2);
                    }

                    LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator);
                    AverageRate += ForwardRate;

                    NCumpound += 1;
                }
            }
        }
        else if (HolidayFlag == 2)
        {
            //////////////////
            // HolidayRate Interpolated Fill
            //////////////////
            TimePos2 = 0;
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                if (i == ForwardStartIdx)
                {
                    if (UseHistorySOFR == 0) isHolidayFlag = 0;
                    else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));
                }
                else isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    lookbackday = i + 0;

                    if (LookBackDays < 1) t = ((double)i) / denominator;
                    else
                    {
                        ///////////////////////////
                        // N영업일 전 LookBackDay 계산
                        ///////////////////////////

                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;

                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }
                    ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                    TargetRate[0] = ForwardRate;
                    TargetT[0] = 0.0;

                    ////////////////////
                    // 해당일 이후 Holiday개수 체크하여 dt설정
                    ////////////////////
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    for (j = i + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));

                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday += 1;
                        }
                    }

                    if (CountHoliday == 0)
                    {
                        LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                        PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                        AverageRate += ForwardRate;
                        NCumpound += 1;
                    }
                    else
                    {
                        ////////////////////
                        // Holiday 이후 이자율추정
                        ////////////////////
                        lookbackday = i + CountHoliday + 1;

                        if (LookBackDays < 1) t = ((double)lookbackday) / denominator;
                        else
                        {
                            k = 0;
                            for (n = 1; n < 120; n++)
                            {
                                /////////////////////
                                // N영업일 전까지 빼기
                                /////////////////////
                                lookbackday -= 1;

                                if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                                if (k == LookBackDays) break;
                            }
                            t = ((double)lookbackday) / denominator;
                        }
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos2);

                        LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                        TargetRate[1] = ForwardRate;
                        TargetT[1] = ((double)(CountHoliday + 1));

                        PI_0 *= (1.0 + TargetRate[0] * 1.0 / denominator);
                        AverageRate += TargetRate[0];
                        NCumpound += 1;
                        for (j = 0; j < CountHoliday; j++)
                        {
                            //////////////////////////////////////
                            // Holiday동안 Linterp Rate로 Compound
                            //////////////////////////////////////
                            if (LockOutFlag == 0) ForwardRate = Interpolate_Linear(TargetT, TargetRate, 2, (double)(j + 1));
                            PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                            AverageRate += ForwardRate;
                            NCumpound += 1;
                        }
                    }
                }
            }
        }
        else
        {
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                lookbackday = i + 0;

                ////////////////////
                // Forward Rate Time 추정 LookBack 반영
                ////////////////////
                if (LookBackDays < 1) t = ((double)i) / denominator;
                else
                {
                    k = 0;
                    for (n = 1; n < 120; n++)
                    {
                        /////////////////////
                        // N영업일 전까지 빼기
                        /////////////////////
                        lookbackday -= 1;
                        if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                        if (k == LookBackDays) break;
                    }
                    t = ((double)lookbackday) / denominator;
                }
                ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                AverageRate += ForwardRate;
                NCumpound += 1;
            }
        }
    }
    else if (ObservShiftFlag == 1)
    {
        if (HolidayFlag == 0)
        {
            //////////////////
            // HolidayRate Forward Fill
            //////////////////
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    lookbackday = i + 0;

                    ////////////////////
                    // Forward Rate Time 추정 LookBack 반영
                    ////////////////////
                    if (LookBackDays < 1) t = ((double)i) / denominator;
                    else
                    {
                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;

                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }

                    ////////////////////
                    // Observe Shift Dt 추정 LookBack 반영
                    ////////////////////
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;
                    for (j = lookbackday + 1; j < ForwardEndIdx; j++)
                    {

                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else CountHoliday += 1;
                    }

                    if (CountHoliday == 0)
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);
                    }
                    else
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos, CountHoliday);
                    }

                    if (lookbackday < 0)
                    {
                        isinflag = isinFindIndex(lookbackday, RefHistoryDate, NRefHistory, HistDayIdx);
                        if (HistDayIdx >= 0) ForwardRate = RefHistory[HistDayIdx];
                        else ForwardRate = CurrentRate;
                    }

                    LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator);
                    AverageRate += ForwardRate;

                    NCumpound += 1;
                }
            }

        }
        else if (HolidayFlag == 1)
        {
            //////////////////
            // HolidayRate Backward Fill
            //////////////////
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    for (j = i + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday += 1;
                        }
                    }

                    if (CountHoliday == 0)
                    {
                        lookbackday = i + 0;
                    }
                    else
                    {
                        lookbackday = i + CountHoliday + 1;
                    }

                    ////////////////////
                    // Forward Rate Time 추정 LookBack 반영
                    ////////////////////
                    if (LookBackDays < 1) t = ((double)lookbackday) / denominator;
                    else
                    {
                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;
                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }

                    CountHoliday2 = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    for (j = lookbackday + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else
                        {
                            CountHoliday2 += 1;
                        }
                    }

                    if (CountHoliday2 == 0)
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);
                    }
                    else
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos, CountHoliday2);
                    }
                    LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                    ////////////////////
                    // Observe Shift Dt 추정 LookBack 반영
                    ////////////////////                    
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;

                    if (lookbackday != i)
                    {
                        for (j = lookbackday + 1; j < ForwardEndIdx; j++)
                        {
                            isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                            if (isHolidayFlag2 == 0) break;
                            else CountHoliday += 1;
                        }
                    }

                    PI_0 *= (1.0 + ForwardRate * (1.0 + (double)CountHoliday) / denominator);
                    AverageRate += ForwardRate;

                    NCumpound += 1;
                }
            }

        }
        else if (HolidayFlag == 2)
        {
            //////////////////
            // HolidayRate Interpolated Fill
            //////////////////
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                isHolidayFlag = max(isin(i, Holiday, NHoliday, HoliStartIdx), isin(i, SaturSunDay, NSaturSunDay, SatSunIdx));

                if (isHolidayFlag == 0)
                {
                    lookbackday = i + 0;

                    ////////////////////
                    // Forward Rate Time1 추정 LookBack 반영
                    ////////////////////
                    if (LookBackDays < 1) t = ((double)i) / denominator;
                    else
                    {
                        k = 0;
                        for (n = 1; n < 120; n++)
                        {
                            /////////////////////
                            // N영업일 전까지 빼기
                            /////////////////////
                            lookbackday -= 1;
                            if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                            if (k == LookBackDays) break;
                        }
                        t = ((double)lookbackday) / denominator;
                    }

                    ////////////////////
                    // Observ Shift
                    ////////////////////
                    CountHoliday = 0;
                    HoliStartIdx2 = 0;
                    SatSunIdx2 = 0;
                    for (j = lookbackday + 1; j < ForwardEndIdx; j++)
                    {
                        isHolidayFlag2 = max(isin(j, Holiday, NHoliday, HoliStartIdx2), isin(j, SaturSunDay, NSaturSunDay, SatSunIdx2));
                        if (isHolidayFlag2 == 0) break;
                        else CountHoliday += 1;
                    }

                    if (CountHoliday == 0)
                    {
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                        LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                        PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                        AverageRate += ForwardRate;

                    }
                    else
                    {
                        TargetT[0] = 0.0;
                        ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                        LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                        TargetRate[0] = ForwardRate;
                        PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                        AverageRate += ForwardRate;
                        NCumpound += 1;


                        ////////////////////
                        // Forward Rate Time2 추정 LookBack 반영
                        ////////////////////
                        t = ((double)(lookbackday + CountHoliday + 1)) / denominator;
                        TargetT[1] = (double)(CountHoliday + 1);
                        TargetRate[1] = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                        LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                        for (j = 0; j < CountHoliday; j++)
                        {
                            if (LockOutFlag == 0) ForwardRate = Interpolate_Linear(TargetT, TargetRate, 2, ((double)(j + 1)));
                            PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                            AverageRate += ForwardRate;
                            NCumpound += 1;
                        }
                    }
                }
            }
        }
        else
        {
            for (i = ForwardStartIdx; i < ForwardEndIdx; i++)
            {
                lookbackday = i + 0;

                ////////////////////
                // Forward Rate Time 추정 LookBack 반영
                ////////////////////
                if (LookBackDays < 1) t = ((double)i) / denominator;
                else
                {
                    k = 0;
                    for (n = 1; n < 120; n++)
                    {
                        /////////////////////
                        // N영업일 전까지 빼기
                        /////////////////////
                        lookbackday -= 1;
                        if (max(isin(lookbackday, Holiday, NHoliday), isin(lookbackday, SaturSunDay, NSaturSunDay)) == 0) k += 1;
                        if (k == LookBackDays) break;
                    }
                    t = ((double)lookbackday) / denominator;
                }
                ForwardRate = Calc_Forward_Rate_Daily(RefCrvTerm, RefCrvRate, NRefCrvTerm, t, &TimePos);

                LockOutCheck(LockOutFlag, LockOutDay, i, LockOutDayRate, ForwardRate);

                PI_0 *= (1.0 + ForwardRate * 1.0 / denominator);
                AverageRate += ForwardRate;
                NCumpound += 1;
            }
        }
    }

    if (AverageFlag == 1)
    {
        AverageRate = AverageRate / (double)NCumpound;
        PI_0 = pow(1.0 + AverageRate / denominator, (int)nday);
    }

    AnnualizedOISRate = (PI_0 - 1.0) * 1.0 / T;

    return PI_0 - 1.0;
}

double Calc_Discount_Factor_Pointer(
    double* TermArray,
    double* RateArray,
    long LengthArray,
    double T,
    long* idx
)
{
    double DF;
    DF = exp(-Interpolate_Linear_Pointer(TermArray, RateArray, LengthArray, T, idx) * T);

    return DF;
}

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
    T = 0.0;
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
    T = 0.0;
    for (i = 0; i < Length_XV_Array; i++)
    {
        t0 = T;
        t1 = t0 + dt_Array[i];
        XV_Array[i] = XV(kappa, NHWVol, HWVolTerm, HWVol, t0, t1);
        T = t1;
    }
}

double B_s_to_t(
    double kappa,
    double s,
    double t
)
{
    return (1.0 - exp(-kappa * (t - s))) / kappa;
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
    double* B_t_T,
    double* dt
)
{
    long i;
    double PtT;
    double term = 1.0 / ((double)NCPN_Ann);
    double ResultRate;
    double A, B;
    if (ReferenceType == 0) NCfSwap = 1;
    if (NCPN_Ann > 0)
    {
        B = 0.0;
        for (i = 0; i < NCfSwap; i++)
        {
            term = dt[i];
            PtT = PV_t_T[i] * exp(-ShortRate * B_t_T[i] + QVTerm[i]);
            B += term * PtT;
        }

        ResultRate = (1. - PtT) / B;
    }
    else
    {
        if (t >= 0.0) PtT = PV_t_T[0] * exp(-ShortRate * B_t_T[0] + QVTerm[0]);
        else PtT = Calc_Discount_Factor(RateTerm, Rate, NRateTerm, T);

        if (dt[0] < 1.0) ResultRate = (1.0 - PtT) / (dt[0] * PtT);
        else ResultRate = -1.0 / dt[0] * log(PtT);
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
    long SimulatedRateShape0,
    long SimulatedRateShape1,
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
    double RateArrayTemp[2];
    double RateTermTemp[2];
    long DayBefore;
    long DayBeforeIdx;
    long ExistHistoryFlag = 0;

    double LowerBound;
    double UpperBound;
    double CondMultiple;
    double PayoffMultiple;
    double Payoff = 0.0;

    double CPN_Ratio = 1.0;
    long Cond = 1;
    long NDays;
    NumAccrual = 0;
    NDays = Leg_Inform->DaysForwardEnd[CashFlowIdx] - Leg_Inform->DaysForwardStart[CashFlowIdx];
    dt = ((double)NDays) / 365.0;
    N = NDays; // 일단 100% Accrual부터 시작

    if (DailySimulFlag == 0 && Leg_Inform->AccrualFlag == 1 && Flag_Required_History == 0)
    {
        //////////////////////////////////////////////////
        // dt간격 시뮬레이션 and 히스토리 필요없을경우  //
        //////////////////////////////////////////////////

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
        for (i = 0; i < Leg_Inform->NReference; i++)
        {
            if (Leg_Inform->Reference_Inform[i].RefRateType == 1)
            {
                N = 0;
                LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                for (j = Leg_Inform->DaysForwardStart[CashFlowIdx]; j < Leg_Inform->DaysForwardEnd[CashFlowIdx]; j++)
                {
                    Rate = DailyRate[i][j];
                    if (Rate >= LowerBound && Rate <= UpperBound) N += 1;
                }
                NumAccrual = min(NumAccrual, N);
            }
        }
        CPN_Ratio = ((double)NumAccrual) / ((double)NDays);
    }

    if (Leg_Inform->PayoffStructure == 0)
    {
        ////////////////////////////
        // 기초금리 AND 조건일 때 //
        ////////////////////////////
        if (Flag_Required_History == 0)
        {
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                PayoffMultiple = Leg_Inform->Reference_Inform[i].PayoffMultiple;
                OutputRate[i] = SimulatedRate[i][CashFlowIdx];
                if (LowerBound <= CondMultiple * SimulatedRate[i][CashFlowIdx] && UpperBound >= CondMultiple * SimulatedRate[i][CashFlowIdx]) Cond = Cond * 1;
                else Cond = 0;
                Payoff += PayoffMultiple * SimulatedRate[i][CashFlowIdx];
            }
            Payoff = Payoff * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;

        }
        else
        {

            DayBefore = -1;
            ExistHistoryFlag = 1;

            RateTermTemp[0] = (double)Leg_Inform->DaysForwardStart[CashFlowIdx];
            RateTermTemp[1] = (double)Leg_Inform->DaysForwardStart[min(CashFlowIdx + 1, SimulatedRateShape1 - 1)];

            Payoff = 0.0;
            NumAccrual = NDays;     // 100% Accrual부터 시작
            for (i = 0; i < Leg_Inform->NReference; i++)
            {
                RateArrayTemp[0] = SimulatedRate[i][CashFlowIdx];
                RateArrayTemp[1] = SimulatedRate[i][min(CashFlowIdx + 1, SimulatedRateShape1 - 1)];
                if (RateArrayTemp[0] == 0.0)
                {
                    ExistHistoryFlag = isinFindIndex(Leg_Inform->DaysForwardStart[CashFlowIdx], Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                    RateArrayTemp[0] = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                    Rate2 = RateArrayTemp[0];
                }

                LowerBound = Leg_Inform->Reference_Inform[i].RangeDn;
                UpperBound = Leg_Inform->Reference_Inform[i].RangeUp;
                CondMultiple = Leg_Inform->Reference_Inform[i].RefRateCondMultiple;
                PayoffMultiple = Leg_Inform->Reference_Inform[i].PayoffMultiple;
                N = NDays;

                if (Leg_Inform->AccrualFlag == 1)
                {
                    for (j = Leg_Inform->DaysForwardStart[CashFlowIdx]; j < Leg_Inform->DaysForwardEnd[CashFlowIdx]; j++)
                    {
                        if (j <= 0)
                        {
                            ExistHistoryFlag = isinFindIndex(j, Leg_Inform->RateHistoryDateMatrix[i], Leg_Inform->NDayHistory[i], DayBeforeIdx);
                            if (ExistHistoryFlag != 0) Rate2 = Leg_Inform->RateHistoryMatrix[i][DayBeforeIdx];
                            else Rate2 = SimulatedRate[i][CashFlowIdx];
                        }
                        else
                        {
                            if (DailySimulFlag == 0) Rate2 = Interpolate_Linear(RateTermTemp, RateArrayTemp, 2, (double)j);
                            else Rate2 = DailyRate[i][j];
                        }
                        if (Rate2 >= LowerBound && Rate2 <= UpperBound) N += 1;
                    }
                }
                NumAccrual = min(NumAccrual, N);
                if (Leg_Inform->Reference_Inform[i].RefRateType == 2)
                {
                    if (Leg_Inform->SOFR_Annualized_R0[i] != 0.0)
                    {
                        OutputRate[i] = Leg_Inform->SOFR_Annualized_R0[i];
                        Payoff += PayoffMultiple * Leg_Inform->SOFR_Compound0[i]/dt;
                    }
                    else
                    {
                        OutputRate[i] = Rate2;
                        Payoff += PayoffMultiple * Rate2;
                    }
                }
                else
                {
                    OutputRate[i] = Rate2;
                    Payoff += PayoffMultiple * Rate2;
                }
                if (LowerBound <= CondMultiple * Rate2 && UpperBound >= CondMultiple * Rate2) Cond = Cond * 1;
                else Cond = 0;

            }

            CPN_Ratio = ((double)NumAccrual) / ((double)NDays);
            Payoff = Payoff * (double)Cond + Leg_Inform->CouponRate[CashFlowIdx] + (double)Cond * Leg_Inform->RangeCoupon[CashFlowIdx] * CPN_Ratio;

        }
    }

    return Payoff * dt;
}

long Simulate_HW(
    long PricingOnly,
    long PricingDateC,
    long NAFlag,
    double Notional,
    LEG_INFO* RcvLeg,
    LEG_INFO* PayLeg,
    SIMUL_INFO* Simul,
    HW_INFO* HW_Information,
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
    long idx0, idx1, idx2, idx3, idx4;
    long idxrcv, idxrcv2, idxrcv3, idxpay, idxpay2, idxpay3;

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


    double** RcvDailyRate = (double**)malloc(sizeof(double*) * RcvLeg->NReference);             // 1
    double** PayDailyRate = (double**)malloc(sizeof(double*) * PayLeg->NReference);             // 2
    if (Simul->DailySimulFlag == 1)
    {
        for (i = 0; i < RcvLeg->NReference; i++) RcvDailyRate[i] = (double*)malloc(sizeof(double) * (Simul->DaysForSimul[Simul->NDays - 1] + 1));
        for (i = 0; i < PayLeg->NReference; i++) PayDailyRate[i] = (double*)malloc(sizeof(double) * (Simul->DaysForSimul[Simul->NDays - 1] + 1));
    }

    double* RcvPayoff = (double*)malloc((RcvLeg->NCashFlow + 1) * sizeof(double));             // 3
    double* PayPayoff = (double*)malloc((PayLeg->NCashFlow + 1) * sizeof(double));             // 4
    long* RcvCFFlag = (long*)malloc(RcvLeg->NCashFlow * sizeof(long));             // 5
    long* PayCFFlag = (long*)malloc(RcvLeg->NCashFlow * sizeof(long));             // 6

    RcvPrice = 0.0;
    PayPrice = 0.0;
    double SOFRDailyCompound = 0.0;

    double* RcvOutputRate = (double*)malloc(RcvLeg->NReference * sizeof(double));
    double* PayOutputRate = (double*)malloc(PayLeg->NReference * sizeof(double));


    ////////////////
    // LSMC Data  //
    ////////////////

    double*** X;
    long ShapeX[2] = { Simul->NSimul, Simul->NAsset + 1 + Simul->NAsset };
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
        Value_By_OptTime = (double*)malloc(sizeof(double) * Simul->NSimul);             // 7
        Estimated_Value_By_OptTime = (double*)malloc(sizeof(double*) * Simul->NSimul);             // 8
        OptimalIdx = (long*)malloc(sizeof(double) * Simul->NSimul);             // 9
        for (i = 0; i < Simul->NSimul; i++) OptimalIdx[i] = 0;

        InterestRate_Opt = (double**)malloc(sizeof(double**) * Simul->NSimul);             // 10
        for (i = 0; i < Simul->NSimul; i++) InterestRate_Opt[i] = (double*)malloc(sizeof(double) * (Simul->NAsset + 1 + Simul->NAsset));
        for (i = 0; i < Simul->NSimul; i++) InterestRate_Opt[i][Simul->NAsset] = 1.0;

        X = (double***)malloc(sizeof(double**) * RcvLeg->NOption);             // 11
        for (i = 0; i < RcvLeg->NOption; i++) X[i] = make_array(Simul->NSimul, Simul->NAsset + 1);

        Y = (double**)malloc(sizeof(double*) * RcvLeg->NOption);             // 12
        for (i = 0; i < RcvLeg->NOption; i++) Y[i] = make_array(Simul->NSimul);

        OptRangeFlag = (long**)malloc(sizeof(long*) * RcvLeg->NOption);             // 13
        for (i = 0; i < RcvLeg->NOption; i++) OptRangeFlag[i] = (long*)malloc(sizeof(long) * Simul->NSimul);

        Y_Hat = (double**)malloc(sizeof(double*) * RcvLeg->NOption);             // 14
        for (i = 0; i < RcvLeg->NOption; i++) Y_Hat[i] = (double*)malloc(sizeof(double) * (Simul->NSimul));

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
    double xt_termlsmc[2] = { 0.0,0.0 };
    double xt_lsmc[2] = { 0.0,0.0 };
    double xt_interp;
    for (i = 0; i < Simul->NSimul; i++)
    {
        PricePath_Rcv = 0.0;
        PricePath_Pay = 0.0;

        for (j = 0; j < Simul->NDays; j++)
        {
            for (n = 0; n < Simul->NAsset; n++)
            {
                if (j == 0) SimulShortRate[n][j] = 0.0;
                else  SimulShortRate[n][j] = Simulated_ShortRate(HW_Information->XA[n][j], HW_Information->XV[n][j], SimulShortRate[n][j - 1], Simul->FixedRandn[i * Simul->NDays + j][n]);
            }
        }

        for (idxoption = 0; idxoption < RcvLeg->NOption; idxoption++)
        {
            T_Opt = ((double)RcvLeg->DaysOptionDate[idxoption]) / 365.0;
            for (j = 0; j < Simul->NDays - 1; j++)
            {
                if (T_Opt < Simul->T_Array[j + 1] && T_Opt >= Simul->T_Array[j])
                {
                    xt_termlsmc[0] = Simul->T_Array[j];
                    xt_termlsmc[1] = Simul->T_Array[j + 1];
                    for (n = 0; n < Simul->NAsset; n++)
                    {
                        xt_lsmc[0] = SimulShortRate[n][j];
                        xt_lsmc[1] = SimulShortRate[n][j + 1];
                        xt_interp = Interpolate_Linear(xt_termlsmc, xt_lsmc, 2, T_Opt);
                        X[idxoption][i][n] = xt_interp;
                    }
                    break;
                }
                else if (T_Opt >= Simul->T_Array[Simul->NDays - 1])
                {
                    for (n = 0; n < Simul->NAsset; n++) X[idxoption][i][n] = SimulShortRate[n][Simul->NDays - 1];
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
            idxrcv = 0;
            idxrcv2 = 0;
            idxrcv3 = 0;
            if (Simul->DailySimulFlag == 0)
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx, &idxrcv))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            t = ((double)RcvLeg->DaysForwardStart[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                            T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                            Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);
                        }
                        else
                        {
                            t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                            T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                            Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);

                            Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T_PowerSpread[n][j],
                                HW_Information->RcvRef_QVTerm_PowerSpread[n][j], HW_Information->RcvRef_B_t_T_PowerSpread[n][j], HW_Information->RcvRef_dt_PowerSpread[n][j]);

                            Rate = Rate1 - Rate2;
                        }
                        SimulatedRateRcv[n][Dayidx] = Rate;
                    }

                }
            }
            else if (Simul->DailySimulFlag == 1 && Simul->SOFRBusinessDaySimul == 0)
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(j, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx, &idxrcv))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        t = Simul->T_Array[j];
                        if (RcvLeg->Reference_Inform[n].RefRateType <= 2)
                        {
                            if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                                Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                    HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);
                            }
                            else
                            {
                                T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                    HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);

                                Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T_PowerSpread[n][j],
                                    HW_Information->RcvRef_QVTerm_PowerSpread[n][j], HW_Information->RcvRef_B_t_T_PowerSpread[n][j], HW_Information->RcvRef_dt_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            Dayidx = 0;
                            RcvDailyRate[n][j] = Rate;
                            if (isinFindIndex(j, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx, &idxrcv2)) SimulatedRateRcv[n][Dayidx] = Rate;
                            if (isinFindIndex(j, RcvLeg->DaysForwardEnd, RcvLeg->NCashFlow, Dayidx, &idxrcv3)) SimulatedRateRcv2[n][Dayidx] = Rate;
                        }
                    }
                }
            }
            else
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, Dayidx, &idxrcv))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        if (RcvLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            t = ((double)RcvLeg->DaysForwardStart[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                            T = ((double)RcvLeg->DaysForwardEnd[Dayidx]) / RcvLeg->Reference_Inform[n].Day1Y;
                            Rate = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);
                        }
                        else
                        {
                            t = ((double)RcvLeg->DaysForwardStart[Dayidx] / RcvLeg->Reference_Inform[n].Day1Y);
                            T1 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + RcvLeg->Reference_Inform[n].PowerSpreadMat2;

                            Rate1 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T[n][j],
                                HW_Information->RcvRef_QVTerm[n][j], HW_Information->RcvRef_B_t_T[n][j], HW_Information->RcvRef_dt[n][j]);

                            Rate2 = HW_Rate(RcvLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_rcv[n], RcvLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->RcvRef_DF_t_T_PowerSpread[n][j],
                                HW_Information->RcvRef_QVTerm_PowerSpread[n][j], HW_Information->RcvRef_B_t_T_PowerSpread[n][j], HW_Information->RcvRef_dt_PowerSpread[n][j]);

                            Rate = Rate1 - Rate2;
                        }
                        SimulatedRateRcv[n][Dayidx] = Rate;
                    }

                }
            }
        }

        for (j = 0; j < RcvLeg->NCashFlow; j++) RcvCFFlag[j] = 0;

        idx1 = 0;
        for (j = 0; j < RcvLeg->NCashFlow; j++)
        {
            RateHistoryUseFlag = 0;
            Day1 = RcvLeg->DaysForwardStart[j];
            Day2 = RcvLeg->DaysPayDate_C[j];
            t = ((double)Day1) / 365.0;
            T = ((double)RcvLeg->DaysPayDate_C[j]) / 365.0;
            Pt = HW_Information->Rcv_DF_0_t[j];
            PT = HW_Information->Rcv_DF_0_T[j];

            if (Day1 >= 0)
            {
                ////////////////////////////////
                // Forward Start 아직  안지남 //
                ////////////////////////////////
                RcvCFFlag[j] = 1;
                RateHistoryUseFlag = 0;
            }
            else if (Day2 > 0)
            {
                ////////////////////////////////////
                // Forward Start지나고 Pay 안지남 //
                ////////////////////////////////////
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
                RcvPayoff[j] = Notional * PayoffStructure(Simul->DailySimulFlag, RcvLeg, RcvLeg->NReference, Simul->NDays, SimulatedRateRcv, SimulatedRateRcv2, RcvDailyRate, SimulOverNightRate, j, RateHistoryUseFlag, RcvOutputRate, nAccrual);
                xt_idx = FindIndex(Day1, Simul->DaysForSimul, Simul->NDays, &idx1);
                xt = SimulShortRate[HW_Information->CurveIdx_DiscRcv][xt_idx];
                if (t > 0.0) PtT = HW_Information->Rcv_DF_t_T[j] * exp(-xt * HW_Information->B_t_T_RcvDisc[j] + HW_Information->QVTerm_RcvDisc[j]);
                else PtT = HW_Information->Rcv_DF_0_T[j];

            }
            else
            {
                RcvPayoff[j] = 0.0;
                for (n = 0; n < RcvLeg->NReference; n++) RcvOutputRate[n] = 0.0;
                nAccrual = 0;
                Pt = 1.0;
                PtT = HW_Information->Rcv_DF_0_T[j];
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
                ///////////////////////////////////////////////////////////////////
                // 해당 Payoff를 받을 수 있는 옵션날짜에 다 Payoff 할인하여 담음 //
                ///////////////////////////////////////////////////////////////////
                for (n = 0; n < RcvLeg->NOption; n++)
                {
                    if (RcvLeg->DaysOptionDate[n] < Day2)
                    {
                        Y[n][i] += Pt * PtT * RcvPayoff[j];
                    }
                }

                if (j == RcvLeg->NCashFlow - 1 && NAFlag == 1)
                {
                    for (n = 0; n < RcvLeg->NOption; n++) Y[n][i] += Pt * PtT * Notional;
                }
            }
        }

        RcvPrice += PricePath_Rcv / (double)Simul->NSimul;

        if (RcvLeg->OptionUseFlag == 1)
        {
            for (n = 0; n < RcvLeg->NOption; n++)
            {
                if (RcvLeg->CallConditionFlag == 0)
                {
                    OptRangeFlag[n][i] = 1;
                    for (k = 0; k < RcvLeg->NReference; k++)
                    {
                        //////////////////////////////////
                        // 하나라도 벗어나면 행사 못함  //
                        //////////////////////////////////
                        if (RcvOutputRate[k] > RcvLeg->RangeUp[k][n] || RcvOutputRate[k] < RcvLeg->RangeDn[k][n])
                        {
                            OptRangeFlag[n][i] = 0;
                            break;
                        }
                    }
                }
                else if (RcvLeg->CallConditionFlag == 1)
                {
                    OptRangeFlag[n][i] = 0;
                    for (k = 0; k < RcvLeg->NReference; k++)
                    {
                        /////////////////////////////////////
                        // Range 하나라도 만족하면 행사함  //
                        /////////////////////////////////////
                        if (RcvOutputRate[k] <= RcvLeg->RangeUp[k][n] && RcvOutputRate[k] >= RcvLeg->RangeDn[k][n])
                        {
                            OptRangeFlag[n][i] = 1;
                            break;
                        }
                    }
                }
                else
                {
                    OptRangeFlag[n][i] = 1;
                }
            }
        }

        // PayLeg Reference Rate 결정
        for (n = 0; n < PayLeg->NReference; n++)
        {
            Curveidx = CurveIdx_Pay[n];
            lastresetend = -1;
            Rate = 0.0;
            Rate1 = 0.0;
            Rate2 = 0.0;
            idxpay = 0;
            idxpay2 = 0;
            idxpay3 = 0;
            if (Simul->DailySimulFlag == 0)
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx, &idxpay))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            t = ((double)PayLeg->DaysForwardStart[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                            T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                            Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);
                        }
                        else
                        {
                            t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                            T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                            Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);

                            Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T_PowerSpread[n][j],
                                HW_Information->PayRef_QVTerm_PowerSpread[n][j], HW_Information->PayRef_B_t_T_PowerSpread[n][j], HW_Information->PayRef_dt_PowerSpread[n][j]);

                            Rate = Rate1 - Rate2;
                        }
                        SimulatedRatePay[n][Dayidx] = Rate;
                    }

                }
            }
            else if (Simul->DailySimulFlag == 1 && Simul->SOFRBusinessDaySimul == 0)
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(j, PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx, &idxpay))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        t = Simul->T_Array[j];
                        if (PayLeg->Reference_Inform[n].RefRateType <= 2)
                        {
                            if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                            {
                                T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                                Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                    HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);
                            }
                            else
                            {
                                T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                                T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                                Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                    HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);

                                Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                    Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T_PowerSpread[n][j],
                                    HW_Information->PayRef_QVTerm_PowerSpread[n][j], HW_Information->PayRef_B_t_T_PowerSpread[n][j], HW_Information->PayRef_dt_PowerSpread[n][j]);

                                Rate = Rate1 - Rate2;
                            }
                            Dayidx = 0;
                            PayDailyRate[n][j] = Rate;
                            if (isinFindIndex(j, PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx, &idxpay2)) SimulatedRatePay[n][Dayidx] = Rate;
                            if (isinFindIndex(j, PayLeg->DaysForwardEnd, PayLeg->NCashFlow, Dayidx, &idxpay3)) SimulatedRatePay2[n][Dayidx] = Rate;
                        }
                    }
                }
            }
            else
            {
                for (j = 0; j < Simul->NDays; j++)
                {
                    if (isinFindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, Dayidx, &idxpay))
                    {
                        xt = SimulShortRate[Curveidx][j];
                        if (PayLeg->Reference_Inform[n].PowerSpreadFlag == 0)
                        {
                            t = ((double)PayLeg->DaysForwardStart[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                            T = ((double)PayLeg->DaysForwardEnd[Dayidx]) / PayLeg->Reference_Inform[n].Day1Y;
                            Rate = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);
                        }
                        else
                        {
                            t = ((double)PayLeg->DaysForwardStart[Dayidx] / PayLeg->Reference_Inform[n].Day1Y);
                            T1 = t + PayLeg->Reference_Inform[n].PowerSpreadMat1;
                            T2 = t + PayLeg->Reference_Inform[n].PowerSpreadMat2;

                            Rate1 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T1, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T[n][j],
                                HW_Information->PayRef_QVTerm[n][j], HW_Information->PayRef_B_t_T[n][j], HW_Information->PayRef_dt[n][j]);

                            Rate2 = HW_Rate(PayLeg->Reference_Inform[n].RefRateType, t, T2, Simul->NRateTerm[Curveidx], Simul->RateTerm[Curveidx],
                                Simul->Rate[Curveidx], xt, HW_Information->ndates_cpn_powerspread_pay[n], PayLeg->Reference_Inform[n].RefSwapNCPN_Ann, HW_Information->PayRef_DF_t_T_PowerSpread[n][j],
                                HW_Information->PayRef_QVTerm_PowerSpread[n][j], HW_Information->PayRef_B_t_T_PowerSpread[n][j], HW_Information->PayRef_dt_PowerSpread[n][j]);

                            Rate = Rate1 - Rate2;
                        }
                        SimulatedRatePay[n][Dayidx] = Rate;
                    }

                }
            }
        }

        for (j = 0; j < PayLeg->NCashFlow; j++) PayCFFlag[j] = 0;

        idx2 = 0;
        for (j = 0; j < PayLeg->NCashFlow; j++)
        {
            RateHistoryUseFlag = 0;
            Day1 = PayLeg->DaysForwardStart[j];
            Day2 = PayLeg->DaysPayDate_C[j];
            t = ((double)Day1) / 365.0;
            T = ((double)PayLeg->DaysPayDate_C[j]) / 365.0;
            Pt = HW_Information->Pay_DF_0_t[j];
            PT = HW_Information->Pay_DF_0_T[j];

            if (Day1 >= 0)
            {
                ////////////////////////////////
                // Forward Start 아직  안지남 //
                ////////////////////////////////
                PayCFFlag[j] = 1;
                RateHistoryUseFlag = 0;
            }
            else if (Day2 > 0)
            {
                ////////////////////////////////////
                // Forward Start지나고 Pay 안지남 //
                ////////////////////////////////////
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
                PayPayoff[j] = Notional * PayoffStructure(Simul->DailySimulFlag, PayLeg, PayLeg->NReference, Simul->NDays, SimulatedRatePay, SimulatedRatePay2, PayDailyRate, SimulOverNightRate, j, RateHistoryUseFlag, PayOutputRate, nAccrual);
                xt_idx = FindIndex(Day1, Simul->DaysForSimul, Simul->NDays, &idx2);
                xt = SimulShortRate[HW_Information->CurveIdx_DiscPay][xt_idx];
                if (t > 0.0) PtT = HW_Information->Pay_DF_t_T[j] * exp(-xt * HW_Information->B_t_T_PayDisc[j] + HW_Information->QVTerm_PayDisc[j]);
                else PtT = HW_Information->Pay_DF_0_T[j];

            }
            else
            {
                PayPayoff[j] = 0.0;
                for (n = 0; n < PayLeg->NReference; n++) PayOutputRate[n] = 0.0;
                nAccrual = 0;
                Pt = 1.0;
                PtT = HW_Information->Pay_DF_0_T[j];
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
                ///////////////////////////////////////////////////////////////////
                // 해당 Payoff를 받을 수 있는 옵션날짜에 다 Payoff 할인하여 담음 //
                ///////////////////////////////////////////////////////////////////
                for (n = 0; n < RcvLeg->NOption; n++)
                {
                    if (RcvLeg->DaysOptionDate[n] < Day2)
                    {
                        Y[n][i] -= Pt * PtT * PayPayoff[j];
                    }
                }

                if (j == PayLeg->NCashFlow - 1 && NAFlag == 1)
                {
                    for (n = 0; n < PayLeg->NOption; n++) Y[n][i] -= Pt * PtT * Notional;
                }
            }

        }

        PayPrice += PricePath_Pay / (double)Simul->NSimul;

        if (RcvLeg->OptionUseFlag == 1)
        {
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
                    if (Y[n][i] < 0.0 && -Y[n][i] > Value_By_OptTime[i] && OptRangeFlag[n][i] > 0)
                    {
                        OptimalIdx[i] = n;
                        Value_By_OptTime[i] = -Y[n][i];
                    }
                }
                else
                {
                    if (Y[n][i] > 0.0 && Y[n][i] > Value_By_OptTime[i] && OptRangeFlag[n][i] > 0)
                    {
                        OptimalIdx[i] = n;
                        Value_By_OptTime[i] = Y[n][i];
                    }
                }
            }

            for (n = 0; n < Simul->NAsset; n++)
            {
                InterestRate_Opt[i][n] = X[OptimalIdx[i]][i][n];
                InterestRate_Opt[i][Simul->NAsset + 1 + n] = X[OptimalIdx[i]][i][n] * X[OptimalIdx[i]][i][n];
            }

        }
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
                    EstOptValue += Beta[Simul->NAsset + 1 + j] * (X[idxoption][i][j] * X[idxoption][i][j]);
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
    free(RcvDailyRate);                                  // 1
    free(PayDailyRate);             // 2
    free(RcvPayoff);             // 3
    free(PayPayoff);             // 4
    free(RcvCFFlag);             // 5
    free(PayCFFlag);             // 6
    free(RcvOutputRate);
    free(PayOutputRate);
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
        free(Y_Hat);

    }
    else
    {
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
    long i, j, k;
    long idx0, idx1, idx2, idx3;
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
    for (i = 0; i < RcvLeg->NReference; i++) CurveIdx_Rcv[i] = FindIndex(RcvLeg->Reference_Inform[i].CurveNum, Simul->SimulCurveIdx, Simul->NAsset);

    for (i = 0; i < PayLeg->NReference; i++) CurveIdx_Pay[i] = FindIndex(PayLeg->Reference_Inform[i].CurveNum, Simul->SimulCurveIdx, Simul->NAsset);

    double Rcv_DF_Day1Y;
    double Pay_DF_Day1Y;
    if (RcvLeg->DiscDayCountFlag == 0) Rcv_DF_Day1Y = 365.0;
    else Rcv_DF_Day1Y = 360.0;

    if (PayLeg->DiscDayCountFlag == 0) Pay_DF_Day1Y = 365.0;
    else Pay_DF_Day1Y = 360.0;

    long CurveIdx_DiscRcv = FindIndex(RcvLeg->DiscCurveNum, Simul->SimulCurveIdx, Simul->NAsset);
    long CurveIdx_DiscPay = FindIndex(PayLeg->DiscCurveNum, Simul->SimulCurveIdx, Simul->NAsset);

    double* Rcv_DF_0_t = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* Rcv_DF_0_T = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* Rcv_DF_t_T = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double t, T;
    idx1 = 0;
    idx2 = 0;
    for (i = 0; i < RcvLeg->NCashFlow; i++)
    {
        t = ((double)RcvLeg->DaysForwardStart[i]) / Rcv_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)RcvLeg->DaysPayDate_C[i]) / Rcv_DF_Day1Y;
        if (T < 0.0) T = 0.0;
        Rcv_DF_0_t[i] = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdx_DiscRcv], Simul->Rate[CurveIdx_DiscRcv], Simul->NRateTerm[CurveIdx_DiscRcv], t, &idx1);
        Rcv_DF_0_T[i] = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdx_DiscRcv], Simul->Rate[CurveIdx_DiscRcv], Simul->NRateTerm[CurveIdx_DiscRcv], T, &idx2);
        Rcv_DF_t_T[i] = Rcv_DF_0_T[i] / Rcv_DF_0_t[i];
    }

    double* Pay_DF_0_t = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* Pay_DF_0_T = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* Pay_DF_t_T = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));

    idx1 = 0;
    idx2 = 0;
    for (i = 0; i < PayLeg->NCashFlow; i++)
    {
        t = ((double)PayLeg->DaysForwardStart[i]) / Pay_DF_Day1Y;
        if (t < 0.0) t = 0.0;
        T = ((double)PayLeg->DaysPayDate_C[i]) / Pay_DF_Day1Y;
        if (T < 0.0) T = 0.0;
        Pay_DF_0_t[i] = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdx_DiscPay], Simul->Rate[CurveIdx_DiscPay], Simul->NRateTerm[CurveIdx_DiscPay], t, &idx1);
        Pay_DF_0_T[i] = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdx_DiscPay], Simul->Rate[CurveIdx_DiscPay], Simul->NRateTerm[CurveIdx_DiscPay], T, &idx2);
        Pay_DF_t_T[i] = Pay_DF_0_T[i] / Pay_DF_0_t[i];
    }

    double* B_t_T_RcvDisc = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* B_t_T_PayDisc = (double*)malloc(PayLeg->NCashFlow * sizeof(double));
    double* QVTerm_RcvDisc = (double*)malloc(RcvLeg->NCashFlow * sizeof(double));
    double* QVTerm_PayDisc = (double*)malloc(PayLeg->NCashFlow * sizeof(double));

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

    long ndates, ndates2 = 1;
    long CurveIdxRef;
    long idxT1, DayT1;
    long ncpn;
    double T1, T_PowerSpread;
    double DF_0_t, DF_0_T, term;

    long* ndates_Rcv = (long*)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_DF_t_T = (double***)malloc(sizeof(double) * RcvLeg->NReference);        // Shape = NReference, NDays, NCpn
    double*** RcvRef_B_t_T = (double***)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_QVTerm = (double***)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_dt = (double***)malloc(sizeof(double**) * RcvLeg->NReference);

    long* ndates_Rcv_PowerSpread = (long*)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_DF_t_T_PowerSpread = (double***)malloc(sizeof(double) * RcvLeg->NReference);        // Shape = NReference, NDays, NCpn
    double*** RcvRef_B_t_T_PowerSpread = (double***)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_QVTerm_PowerSpread = (double***)malloc(sizeof(double) * RcvLeg->NReference);
    double*** RcvRef_dt_PowerSpread = (double***)malloc(sizeof(double**) * RcvLeg->NReference);

    for (i = 0; i < RcvLeg->NReference; i++)
    {
        if (RcvLeg->Reference_Inform[i].RefRateType == 2)
        {
            // SOFR면 대충 아무거나 오류 안나게 할당하게 세팅
            RcvLeg->Reference_Inform[i].RefSwapMaturity = 0.25;
            RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann = 4;
        }

        if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann >= 1 && RcvLeg->Reference_Inform[i].RefSwapMaturity >= 1.0)
            {
                ndates = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].RefSwapMaturity, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
            }
            else if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann == 1 && RcvLeg->Reference_Inform[i].RefSwapMaturity < 1.0)
            {
                ndates = 1;
                RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann = (long)(1.0 / RcvLeg->Reference_Inform[i].RefSwapMaturity + 0.05);
            }
            else ndates = 1;
            ndates_Rcv[i] = ndates;
        }
        else
        {
            if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0)
            {
                ndates = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].PowerSpreadMat1, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                ndates2 = Number_Of_Payment(0.0, RcvLeg->Reference_Inform[i].PowerSpreadMat2, 12.0 / (double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                RcvLeg->Reference_Inform[i].RefRateType = 1;
            }
            else
            {
                ndates = 1;
                ndates2 = 1;
            }
            ndates_Rcv[i] = ndates;
            ndates_Rcv_PowerSpread[i] = ndates2;
        }
    }

    long idxfrdstart;
    for (i = 0; i < RcvLeg->NReference; i++)
    {
        CurveIdxRef = CurveIdx_Rcv[i];
        RcvRef_DF_t_T[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_B_t_T[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_QVTerm[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_dt[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

        RcvRef_DF_t_T_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_B_t_T_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_QVTerm_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        RcvRef_dt_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

        if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            if (RcvLeg->Reference_Inform[i].RefRateType == 0 || RcvLeg->Reference_Inform[i].RefRateType == 2)
            {
                ncpn = 1;
                for (j = 0; j < Simul->NDays; j++)
                {
                    RcvRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    RcvRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    RcvRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    RcvRef_dt[i][j] = (double*)malloc(sizeof(double) * ncpn);
                }

                idxT1 = RcvLeg->CurrentIdx;
                DayT1 = RcvLeg->DaysForwardEnd[RcvLeg->CurrentIdx];
                T1 = ((double)DayT1) / Rcv_DF_Day1Y;
                idx1 = 0;
                idx2 = 0;
                idx3 = 0;
                for (j = 0; j < Simul->NDays; j++)
                {
                    t = ((double)Simul->DaysForSimul[j]) / Rcv_DF_Day1Y;
                    idxfrdstart = FindIndex(Simul->DaysForSimul[j], RcvLeg->DaysForwardStart, RcvLeg->NCashFlow, &idx3);
                    if (idxfrdstart >= 0) T1 = ((double)RcvLeg->DaysForwardEnd[idxfrdstart]) / Rcv_DF_Day1Y;;
                    DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx1);
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx2);
                    RcvRef_dt[i][j][0] = T1 - t;
                    RcvRef_DF_t_T[i][j][0] = DF_0_T / DF_0_t;
                    RcvRef_B_t_T[i][j][0] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                    RcvRef_QVTerm[i][j][0] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                }

            }
            else if (RcvLeg->Reference_Inform[i].RefRateType == 1)
            {
                idx0 = 0;
                for (j = 0; j < Simul->NDays; j++)
                {
                    t = ((double)(Simul->DaysForSimul[j])) / Rcv_DF_Day1Y;

                    ndates = ndates_Rcv[i];

                    if (RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) term = 1.0 / ((double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                    else term = RcvLeg->Reference_Inform[i].RefSwapMaturity;

                    RcvRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                    RcvRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                    RcvRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ndates);
                    RcvRef_dt[i][j] = (double*)malloc(sizeof(double) * ndates);
                    DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx0);
                    idx1 = 0;
                    for (k = 0; k < ndates; k++)
                    {
                        T1 = t + term * ((double)k + 1.0);
                        RcvRef_dt[i][j][k] = term;
                        DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx1);
                        RcvRef_DF_t_T[i][j][k] = DF_0_T / DF_0_t;
                        RcvRef_QVTerm[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                        RcvRef_B_t_T[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                    }
                }
            }
        }
        else
        {
            idx0 = 0;
            for (j = 0; j < Simul->NDays; j++)
            {
                t = ((double)(Simul->DaysForSimul[j])) / Rcv_DF_Day1Y;

                ndates = ndates_Rcv[i];
                term = 1.0 / ((double)RcvLeg->Reference_Inform[i].RefSwapNCPN_Ann);

                RcvRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                RcvRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                RcvRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ndates);
                RcvRef_dt[i][j] = (double*)malloc(sizeof(double) * ndates);
                DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx0);
                idx1 = 0;
                for (k = 0; k < ndates; k++)
                {
                    T1 = t + term * ((double)k + 1.0);
                    RcvRef_dt[i][j][k] = term;
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx1);
                    RcvRef_DF_t_T[i][j][k] = DF_0_T / DF_0_t;
                    RcvRef_QVTerm[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                    RcvRef_B_t_T[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                }

                ndates2 = ndates_Rcv_PowerSpread[i];

                RcvRef_DF_t_T_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                RcvRef_B_t_T_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                RcvRef_QVTerm_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                RcvRef_dt_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                idx2 = 0;
                for (k = 0; k < ndates2; k++)
                {
                    T1 = t + term * ((double)k + 1.0);
                    RcvRef_dt_PowerSpread[i][j][k] = term;
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx2);
                    RcvRef_DF_t_T_PowerSpread[i][j][k] = DF_0_T / DF_0_t;
                    RcvRef_QVTerm_PowerSpread[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                    RcvRef_B_t_T_PowerSpread[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                }
            }
        }
    }

    long* ndates_Pay = (long*)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_DF_t_T = (double***)malloc(sizeof(double) * PayLeg->NReference);        // Shape = NReference, NDays, NCpn
    double*** PayRef_B_t_T = (double***)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_QVTerm = (double***)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_dt = (double***)malloc(sizeof(double**) * PayLeg->NReference);

    long* ndates_Pay_PowerSpread = (long*)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_DF_t_T_PowerSpread = (double***)malloc(sizeof(double) * PayLeg->NReference);        // Shape = NReference, NDays, NCpn
    double*** PayRef_B_t_T_PowerSpread = (double***)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_QVTerm_PowerSpread = (double***)malloc(sizeof(double) * PayLeg->NReference);
    double*** PayRef_dt_PowerSpread = (double***)malloc(sizeof(double**) * PayLeg->NReference);

    for (i = 0; i < PayLeg->NReference; i++)
    {
        if (PayLeg->Reference_Inform[i].RefRateType == 2)
        {
            // SOFR면 대충 아무거나 오류 안나게 할당하게 세팅
            PayLeg->Reference_Inform[i].RefSwapMaturity = 0.25;
            PayLeg->Reference_Inform[i].RefSwapNCPN_Ann = 4;
        }

        if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann >= 1 && PayLeg->Reference_Inform[i].RefSwapMaturity >= 1.0)
            {
                ndates = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].RefSwapMaturity, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
            }
            else if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann == 1 && PayLeg->Reference_Inform[i].RefSwapMaturity < 1.0)
            {
                ndates = 1;
                PayLeg->Reference_Inform[i].RefSwapNCPN_Ann = (long)(1.0 / PayLeg->Reference_Inform[i].RefSwapMaturity + 0.05);
            }
            else ndates = 1;
            ndates_Pay[i] = ndates;
        }
        else
        {
            if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0)
            {
                ndates = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].PowerSpreadMat1, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                ndates2 = Number_Of_Payment(0.0, PayLeg->Reference_Inform[i].PowerSpreadMat2, 12.0 / (double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                PayLeg->Reference_Inform[i].RefRateType = 1;
            }
            else
            {
                ndates = 1;
                ndates2 = 1;
            }
            ndates_Pay[i] = ndates;
            ndates_Pay_PowerSpread[i] = ndates2;
        }
    }

    for (i = 0; i < PayLeg->NReference; i++)
    {
        CurveIdxRef = CurveIdx_Pay[i];
        PayRef_DF_t_T[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_B_t_T[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_QVTerm[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_dt[i] = (double**)malloc(sizeof(double*) * Simul->NDays);

        PayRef_DF_t_T_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_B_t_T_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_QVTerm_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        PayRef_dt_PowerSpread[i] = (double**)malloc(sizeof(double*) * Simul->NDays);
        if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            if (PayLeg->Reference_Inform[i].RefRateType == 0 || PayLeg->Reference_Inform[i].RefRateType == 2)
            {
                ncpn = 1;
                for (j = 0; j < Simul->NDays; j++)
                {
                    PayRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    PayRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    PayRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ncpn);
                    PayRef_dt[i][j] = (double*)malloc(sizeof(double) * ncpn);
                }

                idxT1 = PayLeg->CurrentIdx;
                DayT1 = PayLeg->DaysForwardEnd[PayLeg->CurrentIdx];
                T1 = ((double)DayT1) / Pay_DF_Day1Y;
                idx1 = 0;
                idx2 = 0;
                idx3 = 0;
                for (j = 0; j < Simul->NDays; j++)
                {
                    t = ((double)Simul->DaysForSimul[j]) / Pay_DF_Day1Y;
                    idxfrdstart = FindIndex(Simul->DaysForSimul[j], PayLeg->DaysForwardStart, PayLeg->NCashFlow, &idx3);
                    if (idxfrdstart >= 0) T1 = ((double)PayLeg->DaysForwardEnd[idxfrdstart]) / Pay_DF_Day1Y;;
                    DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx1);
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx2);
                    PayRef_dt[i][j][0] = T1 - t;
                    PayRef_DF_t_T[i][j][0] = DF_0_T / DF_0_t;
                    PayRef_B_t_T[i][j][0] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                    PayRef_QVTerm[i][j][0] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                }

            }
            else if (PayLeg->Reference_Inform[i].RefRateType == 1)
            {
                idx0 = 0;
                for (j = 0; j < Simul->NDays; j++)
                {
                    t = ((double)(Simul->DaysForSimul[j])) / Pay_DF_Day1Y;

                    ndates = ndates_Pay[i];
                    if (PayLeg->Reference_Inform[i].RefSwapNCPN_Ann > 0) term = 1.0 / ((double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);
                    else term = PayLeg->Reference_Inform[i].RefSwapMaturity;

                    PayRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                    PayRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                    PayRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ndates);
                    PayRef_dt[i][j] = (double*)malloc(sizeof(double) * ndates);

                    DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx0);
                    idx1 = 0;
                    for (k = 0; k < ndates; k++)
                    {
                        T1 = t + term * ((double)k + 1.0);
                        PayRef_dt[i][j][k] = term;
                        DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx1);
                        PayRef_DF_t_T[i][j][k] = DF_0_T / DF_0_t;
                        PayRef_QVTerm[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                        PayRef_B_t_T[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                    }
                }
            }
        }
        else
        {
            idx0 = 0;
            for (j = 0; j < Simul->NDays; j++)
            {
                t = ((double)(Simul->DaysForSimul[j])) / Pay_DF_Day1Y;

                ndates = ndates_Pay[i];
                term = 1.0 / ((double)PayLeg->Reference_Inform[i].RefSwapNCPN_Ann);

                PayRef_DF_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                PayRef_B_t_T[i][j] = (double*)malloc(sizeof(double) * ndates);
                PayRef_QVTerm[i][j] = (double*)malloc(sizeof(double) * ndates);
                PayRef_dt[i][j] = (double*)malloc(sizeof(double) * ndates);
                DF_0_t = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], t, &idx0);
                idx1 = 0;
                for (k = 0; k < ndates; k++)
                {
                    T1 = t + term * ((double)k + 1.0);
                    PayRef_dt[i][j][k] = term;
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx1);
                    PayRef_DF_t_T[i][j][k] = DF_0_T / DF_0_t;
                    PayRef_QVTerm[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                    PayRef_B_t_T[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                }

                ndates2 = ndates_Pay_PowerSpread[i];

                PayRef_DF_t_T_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                PayRef_B_t_T_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                PayRef_QVTerm_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                PayRef_dt_PowerSpread[i][j] = (double*)malloc(sizeof(double) * ndates2);
                idx2 = 0;
                for (k = 0; k < ndates2; k++)
                {
                    T1 = t + term * ((double)k + 1.0);
                    PayRef_dt_PowerSpread[i][j][k] = term;
                    DF_0_T = Calc_Discount_Factor_Pointer(Simul->RateTerm[CurveIdxRef], Simul->Rate[CurveIdxRef], Simul->NRateTerm[CurveIdxRef], T1, &idx2);
                    PayRef_DF_t_T_PowerSpread[i][j][k] = DF_0_T / DF_0_t;
                    PayRef_QVTerm_PowerSpread[i][j][k] = HullWhiteQVTerm(t, T1, Simul->HWKappa[CurveIdxRef], Simul->NHWVol[CurveIdxRef], Simul->HWVolTerm[CurveIdxRef], Simul->HWVol[CurveIdxRef]);
                    PayRef_B_t_T_PowerSpread[i][j][k] = B_s_to_t(Simul->HWKappa[CurveIdxRef], t, T1);
                }
            }
        }
    }

    HW_INFO* HW_information = new HW_INFO;
    HW_information->NAsset = Simul->NAsset;
    HW_information->NDays = Simul->NDays;

    HW_information->XA = XA_Array;
    HW_information->XV = XV_Array;

    HW_information->CurveIdx_DiscRcv = CurveIdx_DiscRcv;
    HW_information->CurveIdx_DiscPay = CurveIdx_DiscPay;
    HW_information->NRcv_Cf = RcvLeg->NCashFlow;
    HW_information->NPay_Cf = PayLeg->NCashFlow;

    HW_information->Rcv_DF_0_t = Rcv_DF_0_t;
    HW_information->Pay_DF_0_t = Pay_DF_0_t;
    HW_information->Rcv_DF_0_T = Rcv_DF_0_T;
    HW_information->Pay_DF_0_T = Pay_DF_0_T;
    HW_information->Rcv_DF_t_T = Rcv_DF_t_T;
    HW_information->Pay_DF_t_T = Pay_DF_t_T;

    HW_information->B_t_T_RcvDisc = B_t_T_RcvDisc;
    HW_information->B_t_T_PayDisc = B_t_T_PayDisc;
    HW_information->QVTerm_RcvDisc = QVTerm_RcvDisc;
    HW_information->QVTerm_PayDisc = QVTerm_PayDisc;
    HW_information->NRcvCurve = RcvLeg->NReference;
    HW_information->NPayCurve = RcvLeg->NReference;

    HW_information->ndates_cpn_rcv = ndates_Rcv;
    HW_information->RcvRef_DF_t_T = RcvRef_DF_t_T;
    HW_information->RcvRef_B_t_T = RcvRef_B_t_T;
    HW_information->RcvRef_QVTerm = RcvRef_QVTerm;
    HW_information->RcvRef_dt = RcvRef_dt;

    HW_information->ndates_cpn_powerspread_rcv = ndates_Rcv_PowerSpread;
    HW_information->RcvRef_DF_t_T_PowerSpread = RcvRef_DF_t_T_PowerSpread;
    HW_information->RcvRef_B_t_T_PowerSpread = RcvRef_B_t_T_PowerSpread;
    HW_information->RcvRef_QVTerm_PowerSpread = RcvRef_QVTerm_PowerSpread;
    HW_information->RcvRef_dt_PowerSpread = RcvRef_dt_PowerSpread;

    HW_information->ndates_cpn_pay = ndates_Pay;
    HW_information->PayRef_DF_t_T = PayRef_DF_t_T;
    HW_information->PayRef_B_t_T = PayRef_B_t_T;
    HW_information->PayRef_QVTerm = PayRef_QVTerm;
    HW_information->PayRef_dt = PayRef_dt;

    HW_information->ndates_cpn_powerspread_pay = ndates_Pay_PowerSpread;
    HW_information->PayRef_DF_t_T_PowerSpread = PayRef_DF_t_T_PowerSpread;
    HW_information->PayRef_B_t_T_PowerSpread = PayRef_B_t_T_PowerSpread;
    HW_information->PayRef_QVTerm_PowerSpread = PayRef_QVTerm_PowerSpread;
    HW_information->PayRef_dt_PowerSpread = PayRef_dt_PowerSpread;

    ////////////////////////////
    // SOFR 첫 번째 Cpn계산 //
    ////////////////////////////
    long idx;
    double SOFRComp0_Rcv[3] = { 0.0,0.0,0.0 };
    double SOFR_AnnOISRate0_Rcv[3] = { 0.0,0.0,0.0 };
    double SOFRComp0_Pay[3] = { 0.0,0.0,0.0 };
    double SOFR_AnnOISRate0_Pay[3] = { 0.0,0.0,0.0 };
    for (i = 0; i < RcvLeg->NReference; i++)
    {
        idx = CurveIdx_Rcv[i];
        if (RcvLeg->DaysForwardStart[RcvLeg->CurrentIdx] <= 0 && RcvLeg->DaysForwardEnd[RcvLeg->CurrentIdx] > 0)
        {
            SOFRComp0_Rcv[i] = SOFR_ForwardRate_Compound(Simul->NRateTerm[idx], Simul->RateTerm[idx], Simul->Rate[idx], RcvLeg->DaysForwardStart[RcvLeg->CurrentIdx], RcvLeg->DaysForwardEnd[RcvLeg->CurrentIdx],
                0, 0, 0, RcvLeg->HolidayCalcFlag[i], RcvLeg->HolidayCount[i], RcvLeg->HolidayDays[i], RcvLeg->NWeekend[i], RcvLeg->Weekend[i], 1, RcvLeg->NDayHistory[i],
                RcvLeg->RateHistoryDateMatrix[i], RcvLeg->RateHistoryMatrix[i], Rcv_DF_Day1Y, SOFR_AnnOISRate0_Rcv[i], 1);
        }
    }
    RcvLeg->SOFR_Annualized_R0 = SOFR_AnnOISRate0_Rcv;
    RcvLeg->SOFR_Compound0 = SOFRComp0_Rcv;
    PayLeg->SOFR_Annualized_R0 = SOFR_AnnOISRate0_Pay;
    PayLeg->SOFR_Compound0 = SOFRComp0_Pay;



    ResultCode = Simulate_HW(1, PricingDateC, NAFlag, Notional, RcvLeg, PayLeg,
        Simul, HW_information,
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
    free(Rcv_DF_0_t);
    free(Rcv_DF_0_T);
    free(Pay_DF_0_t);
    free(Pay_DF_0_T);
    free(Rcv_DF_t_T);
    free(Pay_DF_t_T);
    free(B_t_T_RcvDisc);
    free(B_t_T_PayDisc);
    free(QVTerm_RcvDisc);
    free(QVTerm_PayDisc);
    for (i = 0; i < RcvLeg->NReference; i++)
    {
        if (RcvLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(RcvRef_DF_t_T[i][j]);
                free(RcvRef_B_t_T[i][j]);
                free(RcvRef_QVTerm[i][j]);
                free(RcvRef_dt[i][j]);
            }
        }
        else
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(RcvRef_DF_t_T[i][j]);
                free(RcvRef_B_t_T[i][j]);
                free(RcvRef_QVTerm[i][j]);
                free(RcvRef_dt[i][j]);

                free(RcvRef_DF_t_T_PowerSpread[i][j]);
                free(RcvRef_B_t_T_PowerSpread[i][j]);
                free(RcvRef_QVTerm_PowerSpread[i][j]);
                free(RcvRef_dt_PowerSpread[i][j]);
            }

        }
        free(RcvRef_DF_t_T[i]);
        free(RcvRef_B_t_T[i]);
        free(RcvRef_QVTerm[i]);
        free(RcvRef_dt[i]);

        free(RcvRef_DF_t_T_PowerSpread[i]);
        free(RcvRef_B_t_T_PowerSpread[i]);
        free(RcvRef_QVTerm_PowerSpread[i]);
        free(RcvRef_dt_PowerSpread[i]);
    }

    free(RcvRef_DF_t_T);
    free(RcvRef_B_t_T);
    free(RcvRef_QVTerm);
    free(RcvRef_dt);

    free(RcvRef_DF_t_T_PowerSpread);
    free(RcvRef_B_t_T_PowerSpread);
    free(RcvRef_QVTerm_PowerSpread);
    free(RcvRef_dt_PowerSpread);
    //
    for (i = 0; i < PayLeg->NReference; i++)
    {
        if (PayLeg->Reference_Inform[i].PowerSpreadFlag == 0)
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(PayRef_DF_t_T[i][j]);
                free(PayRef_B_t_T[i][j]);
                free(PayRef_QVTerm[i][j]);
                free(PayRef_dt[i][j]);
            }
        }
        else
        {
            for (j = 0; j < Simul->NDays; j++)
            {
                free(PayRef_DF_t_T[i][j]);
                free(PayRef_B_t_T[i][j]);
                free(PayRef_QVTerm[i][j]);
                free(PayRef_dt[i][j]);

                free(PayRef_DF_t_T_PowerSpread[i][j]);
                free(PayRef_B_t_T_PowerSpread[i][j]);
                free(PayRef_QVTerm_PowerSpread[i][j]);
                free(PayRef_dt_PowerSpread[i][j]);
            }

        }
        free(PayRef_DF_t_T[i]);
        free(PayRef_B_t_T[i]);
        free(PayRef_QVTerm[i]);
        free(PayRef_dt[i]);

        free(PayRef_DF_t_T_PowerSpread[i]);
        free(PayRef_B_t_T_PowerSpread[i]);
        free(PayRef_QVTerm_PowerSpread[i]);
        free(PayRef_dt_PowerSpread[i]);
    }

    free(PayRef_DF_t_T);
    free(PayRef_B_t_T);
    free(PayRef_QVTerm);
    free(PayRef_dt);

    free(PayRef_DF_t_T_PowerSpread);
    free(PayRef_B_t_T_PowerSpread);
    free(PayRef_QVTerm_PowerSpread);
    free(PayRef_dt_PowerSpread);

    free(ndates_Rcv);
    free(ndates_Pay);
    free(ndates_Rcv_PowerSpread);
    free(ndates_Pay_PowerSpread);
    delete(HW_information);

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

    double* ResultPrice,                // 54 산출된 가격 [0] Rcv [1] Pay [2] Price [3] OptionPrice [4~
    double* ResultRcv,                  // 55 [0~NCF-1]금리1, [NCF~2NCF-1]금리2, [2NCF~3NCF-1]금리3, [3NCF~4NCF-1] E(Accrual수), [4NCF~5NCF-1] E(CPN), [5NCF~6NCF-1] E(DF)
    double* ResultPay,                  // 56 [0~NCF-1]금리1, [NCF~2NCF-1]금리2, [2NCF~3NCF-1]금리3, [3NCF~4NCF-1] E(Accrual수), [4NCF~5NCF-1] E(CPN), [5NCF~6NCF-1] E(DF)
    char* Error                         // 57 에러메시지
)
{

    long i;
    //_CrtSetBreakAlloc(44399);
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
    
    double* ResultPV01 = ResultPrice + 4;
    double* ResultKeyPV01 = ResultPV01 + 8;
    

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
    // 주말 처리 관련
    //////////////////////////////
    long SOFRUsingFlag = 0;
    for (i = 0; i < 3; i++)
    {
        if (RcvRefRateType[i] >= 2 || PayRefRateType[i] >= 2)
        {
            SOFRUsingFlag = 1;
            break;
        }
    }

    long WeekCheckStart;
    long WeekCheckEnd;
    long NWeekend;

    long* NRcv_Weekend = (long*)malloc(sizeof(long) * NRcvCashFlow);
    long** Rcv_Weekend = (long**)malloc(sizeof(long*) * NRcvCashFlow);
    long* NPay_Weekend = (long*)malloc(sizeof(long) * NPayCashFlow);
    long** Pay_Weekend = (long**)malloc(sizeof(long*) * NPayCashFlow);

    if (SOFRUsingFlag == 1)
    {
        for (i = 0; i < NRcvCashFlow; i++)
        {
            WeekCheckStart = RcvCashFlowScheduleExcelDate[i] - 100;
            WeekCheckEnd = RcvCashFlowScheduleExcelDate[i + NRcvCashFlow] + 365;

            NWeekend = 0;
            for (j = WeekCheckStart; j < WeekCheckEnd; j++)
            {
                if (isweekend(j) == 1)
                {
                    NWeekend += 1;
                }
            }
            NRcv_Weekend[i] = NWeekend;
            Rcv_Weekend[i] = (long*)malloc(sizeof(long) * max(1, NWeekend));
            k = 0;
            for (j = WeekCheckStart; j < WeekCheckEnd; j++)
            {
                if (isweekend(j) == 1)
                {
                    Rcv_Weekend[i][k] = j - PriceDateExcel;
                    k += 1;
                }
            }
        }

        for (i = 0; i < NPayCashFlow; i++)
        {
            WeekCheckStart = PayCashFlowScheduleExcelDate[i] - 100;
            WeekCheckEnd = PayCashFlowScheduleExcelDate[i + NPayCashFlow] + 365;

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
                    Pay_Weekend[i][k] = j - PriceDateExcel;
                    k += 1;
                }
            }
        }
    }

    /////////////////////////////////////////////////
    // SOFR 하나라도 있으면 
    // 그리고 Daily시뮬레이션이라면
    // 주말 빼고 시뮬레이션
    ///////////////////////////////////////////
    long nSOFR = 0;
    if (DailySimulFlag == 1 && SOFRUsingFlag == 1) nSOFR = 1;


    //////////////////////////////
    // 옵션 관련 Date
    //////////////////////////////
    long NOptionOrg = NOption;
    long NPrevOpt = 0;
    long* OptionDateCtype;
    long* DaysOptionDate;
    for (i = 0; i < NOption; i++)
    {
        if (OptionDateAndFlag[i] <= PriceDateExcel) NPrevOpt += 1;
    }

    OptionDateCtype = (long*)malloc(sizeof(long) * max(1, NOption));
    DaysOptionDate = (long*)malloc(sizeof(long) * max(1, NOption));
    k = 0;
    for (i = 0; i < NOption; i++)
    {
        if (OptionDateAndFlag[i] > PriceDateExcel)
        {
            OptionDateCtype[k] = ExcelDateToCDate(OptionDateAndFlag[i]);
            DaysOptionDate[k] = DayCountAtoB(PriceDateC, OptionDateCtype[k]);
            k += 1;
        }
    }
    NOption = NOption - NPrevOpt;

    long OptionUseFlag;
    if (NOption > 0) OptionUseFlag = 1;
    else OptionUseFlag = 0;

    long CallConditionFlag = OptionDateAndFlag[NOptionOrg];
    long OptionType = OptionDateAndFlag[NOptionOrg + 1];

    double* OptionStrikeRate = (double*)malloc(sizeof(double) * max(1, NOption));
    for (i = 0; i < NOption; i++) OptionStrikeRate[i] = OptionKAndRange[i + NPrevOpt];

    double** RangeUp;
    double** RangeDn;

    RangeUp = (double**)malloc(sizeof(double*) * max(1, RcvFlag[0]));
    RangeDn = (double**)malloc(sizeof(double*) * max(1, RcvFlag[0]));
    for (i = 0; i < RcvFlag[0]; i++)
    {
        RangeUp[i] = OptionKAndRange + (2 * i + 1) * NOptionOrg + NPrevOpt;
        RangeDn[i] = OptionKAndRange + (2 * i + 2) * NOptionOrg + NPrevOpt;
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
    RcvLeg->NWeekend = NRcv_Weekend;
    RcvLeg->Weekend = Rcv_Weekend;

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
    long IdxCurrentRcvLeg = 0;
    for (i = 0; i < NRcvCashFlow; i++)
    {
        if (PriceDateC < RcvLeg->ForwardEnd_C[i] && PriceDateC >= RcvLeg->ForwardStart_C[i])
        {
            IdxCurrentRcvLeg = i;

        }
    }
    RcvLeg->CurrentIdx = IdxCurrentRcvLeg;

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

    PayLeg->NWeekend = NPay_Weekend;
    PayLeg->Weekend = Pay_Weekend;

    long IdxCurrentPayLeg = 0;
    for (i = 0; i < NPayCashFlow; i++)
    {
        if (PriceDateC < PayLeg->ForwardEnd_C[i] && PriceDateC >= PayLeg->ForwardStart_C[i])
        {
            IdxCurrentPayLeg = i;

        }
    }
    PayLeg->CurrentIdx = IdxCurrentPayLeg;


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
    for (i = 0; i < nSimulateCurve; i++) CorrMatAdj[i] = CorrMatReshaped + i * nSimulateCurve;

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

    long MaxDaysSimul = 1000;

    if (DailySimulFlag != 0)
    {
        if (nSOFR == 0)
        {
            MaxDaysSimul = DaysForSimul[NDays - 1];
            free(DaysForSimul);
            DaysForSimul = (long*)malloc(sizeof(long) * MaxDaysSimul);
            for (i = 0; i < MaxDaysSimul; i++) DaysForSimul[i] = i;
        }
        else
        {
            n = DaysForSimul[NDays - 1];
            k = 0;
            for (i = 0; i < n; i++)
            {
                if (isweekend(PriceDateExcel + i) == 0 || isin(i, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow) || isin(i, RcvLeg->DaysForwardEnd, RcvLeg->NCashFlow) || isin(i, PayLeg->DaysForwardStart, PayLeg->NCashFlow) || isin(i, PayLeg->DaysForwardEnd, PayLeg->NCashFlow))
                {
                    k = k + 1;
                }
            }

            MaxDaysSimul = k;
            free(DaysForSimul);
            
            DaysForSimul = (long*)malloc(sizeof(long) * MaxDaysSimul);
            k = 0;
            for (i = 0; i < n; i++)
            {
                if (isweekend(PriceDateExcel + i) == 0 || isin(i, RcvLeg->DaysForwardStart, RcvLeg->NCashFlow) || isin(i, RcvLeg->DaysForwardEnd, RcvLeg->NCashFlow) || isin(i, PayLeg->DaysForwardStart, PayLeg->NCashFlow) || isin(i, PayLeg->DaysForwardEnd, PayLeg->NCashFlow))
                {
                    DaysForSimul[k] = i;
                    k = k + 1;
                }
            }
        }
    }
    else
    {
        MaxDaysSimul = NDays;
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
    Simul->SOFRBusinessDaySimul = nSOFR;

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

    if (GreekFlag == 1)
    {
        long idxtemp = 0;
        long idxtermup = 0;
        long kk = 0;
        long h = 0;
        long NTotalZeroRate = 0;
        long PointerZero[5] = { 0,0,0,0,0 };
        long idxstart, idxend;
        long UselessTermFlag = 0;
        double RcvValue, PayValue;
        double PV01Rcv, PV01Pay;
        double TempResultPrice[5] = { 0.0, 0.0,0.0,0.0 };
        double* TempResultRcv = (double*)malloc(sizeof(double) * (NRcvCashFlow * 7));
        double* TempResultPay = (double*)malloc(sizeof(double) * (NPayCashFlow * 7));
        double* KeyPV01Rcv;
        double* KeyPV01Pay;
        RcvValue = ResultPrice[0];
        PayValue = ResultPrice[1];
        for (i = 0; i < N_Curve_Max; i++)
        {
            NTotalZeroRate += NZeroRate[i];
            PointerZero[i + 1] = NTotalZeroRate;
        }

        double* ZeroRateParallelUp = (double*)malloc(sizeof(double) * NTotalZeroRate);
        double* ZeroRateKeyUp = (double*)malloc(sizeof(double) * NTotalZeroRate);
        double** ZeroRateMatrixParralelUp = (double**)malloc(sizeof(double*) * N_Curve_Max);
        n = 0;
        for (i = 0; i < N_Curve_Max; i++)
        {
            ZeroRateMatrixParralelUp[i] = ZeroRateParallelUp + n;
            n = n + NZeroRate[i];
        }

        SIMUL_INFO* Simul_ForGreek = new SIMUL_INFO;
        Simul_ForGreek->NSimul = Simul->NSimul;
        Simul_ForGreek->DailySimulFlag= Simul->DailySimulFlag;
        Simul_ForGreek->NDays = Simul->NDays;
        Simul_ForGreek->NAsset = Simul->NAsset;
        Simul_ForGreek->DaysForSimul = Simul->DaysForSimul;
        Simul_ForGreek->SOFRBusinessDaySimul = nSOFR;

        Simul_ForGreek->dt_Array = Simul->dt_Array;
        Simul_ForGreek->T_Array = Simul->T_Array;
        Simul_ForGreek->FixedRandn = FixedRandn;
        Simul_ForGreek->NHWVol = Simul->NHWVol ;
        Simul_ForGreek->HWKappa = Simul->HWKappa;
        Simul_ForGreek->HWVolTerm  = Simul->HWVolTerm ;
        Simul_ForGreek->HWVol = Simul->HWVol  ;
        Simul_ForGreek->NRateTerm  = Simul->NRateTerm ;
        Simul_ForGreek->RateTerm = Simul->RateTerm ;
        Simul_ForGreek->Rate = (double**)malloc(sizeof(double*) * Simul->NAsset);
        Simul_ForGreek->SimulCurveIdx = SimulateCurveIdx;
        k = 0;
        for (i = 1; i < N_Curve_Max + 1; i++)
        {
            if (isin(i, Simul->SimulCurveIdx, Simul->NAsset))
            {
                idxstart = PointerZero[i - 1];
                idxend = PointerZero[i];

                for (j = 0; j < NTotalZeroRate; j++)
                {
                    if (j >= idxstart && j < idxend) ZeroRateParallelUp[j] = ZeroRate[j] + 0.0001;
                    else ZeroRateParallelUp[j] = ZeroRate[j];                  
                }


                kk = 0;
                for (idxtemp = 0; idxtemp < N_Curve_Max; idxtemp++)
                {
                    if (isin(idxtemp + 1, SimulateCurveIdx, nSimulateCurve))
                    {
                        Simul_ForGreek->Rate[kk] = ZeroRateMatrixParralelUp[idxtemp];
                        kk = kk + 1;
                    }
                }

                ResultCode = IRStructuredSwap(PriceDateC, NAFlag, Notional, RcvLeg, PayLeg, Simul_ForGreek, GreekFlag, TempResultPrice, TempResultRcv, TempResultPay);

                PV01Rcv = TempResultPrice[0] - RcvValue;
                PV01Pay = TempResultPrice[1] - PayValue;
                ResultPV01[i - 1] = PV01Rcv;
                ResultPV01[i - 1 + 4] = PV01Pay;
                k = k + 1;
            }
        }

        k = 0;
        for (i = 1; i < N_Curve_Max + 1; i++)
        {
            KeyPV01Rcv = ResultKeyPV01 + PointerZero[i - 1];
            KeyPV01Pay = ResultKeyPV01 + NTotalZeroRate + PointerZero[i - 1];
            if (isin(i, Simul->SimulCurveIdx, Simul->NAsset))
            {
                idxstart = PointerZero[i - 1];
                idxend = PointerZero[i];
                for (idxtermup = idxstart; idxtermup < idxend; idxtermup++)
                {
                    for (j = 0; j < NTotalZeroRate; j++)
                    {
                        if (j == idxtermup) ZeroRateParallelUp[j] = ZeroRate[j] + 0.0001;
                        else ZeroRateParallelUp[j] = ZeroRate[j];
                    }

                    kk = 0;
                    for (idxtemp = 0; idxtemp < N_Curve_Max; idxtemp++)
                    {
                        if (isin(idxtemp + 1, SimulateCurveIdx, nSimulateCurve))
                        {
                            Simul_ForGreek->Rate[kk] = ZeroRateMatrixParralelUp[idxtemp];
                            kk = kk + 1;
                        }
                    }

                    UselessTermFlag = 0;
                    if (NZeroRate[i] > 2)
                    {
                        ////////////////////////////////////////////////////
                        // 스왑 만기보다 너무 Term이 길면 노트 가격 영향 없다
                        ////////////////////////////////////////////////////
                        if (ZeroTerm[idxtermup - 1] > Simul_ForGreek->T_Array[Simul_ForGreek->NDays - 1]) UselessTermFlag = 1;
                    }

                    if (UselessTermFlag == 1)
                    {
                        KeyPV01Rcv[idxtermup - idxstart] = 0.0;
                        KeyPV01Pay[idxtermup - idxstart] = 0.0;
                    }
                    else
                    {
                        ResultCode = IRStructuredSwap(PriceDateC, NAFlag, Notional, RcvLeg, PayLeg, Simul_ForGreek, GreekFlag, TempResultPrice, TempResultRcv, TempResultPay);

                        PV01Rcv = TempResultPrice[0] - RcvValue;
                        PV01Pay = TempResultPrice[1] - PayValue;
                        KeyPV01Rcv[idxtermup - idxstart] = PV01Rcv;
                        KeyPV01Pay[idxtermup - idxstart] = PV01Pay;
                    }
                    k = k + 1;
                }
            }
        }

        free(TempResultRcv);
        free(TempResultPay);
        free(ZeroRateParallelUp);
        free(ZeroRateKeyUp);
        free(ZeroRateMatrixParralelUp);
        free(Simul_ForGreek->Rate);
        delete(Simul_ForGreek);
    }

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
    
    if (SOFRUsingFlag == 1)
    {
        for (i = 0; i < NRcvCashFlow; i++) free(Rcv_Weekend[i]);
        for (i = 0; i < NPayCashFlow; i++) free(Pay_Weekend[i]);
    }
    free(NRcv_Weekend);
    free(Rcv_Weekend);
    free(NPay_Weekend);
    free(Pay_Weekend);

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