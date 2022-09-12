#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "Util.h"

#include <crtdbg.h>


#ifndef DLLEXPORT(A)
#ifdef WIN32
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A __stdcall 
#elif _WIN64
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A __stdcall 
#elif __linux__
#define DLLEXPORT(A) extern "C" __declspec
#elif __hpux
#define DLLEXPORT(A) extern "C" __declspec
#elif __unix__
#define DLLEXPORT(A) extern "C" __declspec
#else 
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A __stdcall 
#endif
#endif

#define NTAB                32
#define Notional_Amount     10000.0
#define Param1_Min          0.002         //Kappa최소값
#define Param1_Max          0.1          //Kappa최대값
#define Interval1           0.002         //Kappa간격
#define Param2_Min          0.001         //Vol최소값
#define Param2_Max          0.04          //Vol최대값
#define Interval2           0.0002          //Vol간격

#define Tiny_Value                1.0e-7        // 변동성 최소값
#define Max_Value                1.0            // 변동성 최대값

long Number_Of_Payment(double T_OptionMaturity, double T_SwapMaturity, double PayFreqOfMonth)
{
    double dT = (T_SwapMaturity - T_OptionMaturity + 0.00001);
    double Num_AnnPayment = 12.0 / PayFreqOfMonth;
    long N;
    N = (long)(dT * Num_AnnPayment + 0.01);
    return N;
}

long PaymentDatesMapping(double T_SwapMaturity, double FreqMonth, long* TempDatesArray, long NDates)
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

double FSR(
    double* Term,
    double* Rate,
    long NTerm,
    double T_Option,
    double Tenor,
    double FreqMonth
)
{
    long i;
    double Swap_Rate;

    long ndates = Number_Of_Payment(T_Option, T_Option + Tenor, FreqMonth);
    long* dates = (long*)malloc(sizeof(long) * (ndates));
    ndates = PaymentDatesMapping(T_Option + Tenor, FreqMonth, dates, ndates);

    double* PT = (double*)malloc(sizeof(double) * (ndates + 1));
    PT[0] = Calc_Discount_Factor(Term, Rate, NTerm, T_Option); // 옵션만기시점
    for (i = 1; i < ndates + 1; i++)
    {
        PT[i] = Calc_Discount_Factor(Term, Rate, NTerm, (double)dates[i - 1] / 365.0);
    }

    double a, b, dt;
    a = PT[0] - PT[ndates];
    b = 0.0;
    for (i = 0; i < ndates; i++)
    {
        if (i == 0)
            dt = (double)dates[0] / 365.0 - T_Option;
        else
            dt = ((double)(dates[i] - dates[i - 1])) / 365.0;
        b += dt * PT[i + 1];
    }
    Swap_Rate = a / b;


    if (dates) free(dates);
    if (PT) free(PT);
    return Swap_Rate;
}


// 적분 계산 공통함수
// I(t) = Int_0^t sigma(s)^2 A exp(Bs) ds
double Integ(
    double t,
    double A,
    double kappa,
    double* tVol,
    double* Vol,
    long nVol
)
{
    long i;
    long NodeNum = 3;
    double ds = t / (double)NodeNum;
    double s;
    double value = 0.0;
    double sigma;

    for (i = 0; i < NodeNum; i++)
    {
        s = (double)(i + 1) * ds;
        sigma = Interpolate_Linear(tVol, Vol, nVol, s);
        value += sigma * sigma * A * exp(kappa * s) * ds;
    }
    return value;
}

double B(double s, double t, double kappa)
{
    return (1.0 - exp(-kappa * (t - s))) / kappa;
}

// 1-factor 모형의 Fixed Payer Swaption 가격 계산
double HW_Swaption(
    double NA,            // 액면금액
    double kappa,        // 회귀속도 
    double* tVol,        // 변동성 구간 종점
    double* Vol,            // 구간 변동성
    long nVol,            // 변동성 구간 개수
    double* t,            // 할인채 만기
    double* r,            // 할인채 가격
    long nr,                // 할인채 개수
    double StrikeRate,    // 고정금리(지급부분)
    long MaturityDate,    // 옵션 만기일까지 일수
    long* Dates,            // 지급일: 계산일로부터 각 중간지급일까지의 일수
    long nDates            // 지급 회수(계산일 이후 남은 회수)
)
{
    long i;
    double T0, PrevT, T, deltaT;
    double* PT = (double*)malloc(sizeof(double) * (nDates + 1));
    double VT0, G, H;
    double d1, d2;
    double value;

    if (kappa < Tiny_Value) kappa = Tiny_Value;

    for (i = 0; i < nVol; i++) {
        if (Vol[i] < 0.0) Vol[i] = -Vol[i];
        if (Vol[i] < Tiny_Value) Vol[i] = Tiny_Value;
    }

    T0 = (double)MaturityDate / 365.0;
    PT[0] = Calc_Discount_Factor(t, r, nr, T0);

    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        PT[i + 1] = Calc_Discount_Factor(t, r, nr, T);
    }

    VT0 = exp(-2.0 * kappa * T0) * Integ(T0, 1.0, 2.0 * kappa, tVol, Vol, nVol);

    G = PT[nDates];
    PrevT = T0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;
        G += StrikeRate * deltaT * PT[i + 1];
        PrevT = T;
    }
    G /= PT[0];

    H = PT[nDates] * B(T0, (double)Dates[nDates - 1] / 365.0, kappa);
    PrevT = T0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;
        H += StrikeRate * deltaT * PT[i + 1] * B(T0, T, kappa);
        PrevT = T;
    }
    H /= G * PT[0];

    d1 = -log(G) / (H * sqrt(VT0)) + 0.5 * H * sqrt(VT0);
    d2 = -log(G) / (H * sqrt(VT0)) - 0.5 * H * sqrt(VT0);

    value = PT[0] * (CDF_N(d1) - G * CDF_N(d2));

    if (PT) free(PT);

    return NA * value;
}

// 1-factor 모형의 Fixed Payer Swaption 가격 빠른 계산
double HW_Swaption(
    double NA,            // 액면금액
    double kappa,        // 회귀속도 
    double* tVol,        // 변동성 구간 종점
    double* Vol,            // 구간 변동성
    long nVol,            // 변동성 구간 개수
    double* t,            // 할인채 만기
    double* r,            // 할인채 가격
    long nr,                // 할인채 개수
    double StrikeRate,    // 고정금리(지급부분)
    long MaturityDate,    // 옵션 만기일까지 일수
    long* Dates,            // 지급일: 계산일로부터 각 중간지급일까지의 일수
    long nDates,            // 지급 회수(계산일 이후 남은 회수
    double* PT              // DiscountFactor Shape = nDates + 1
)
{
    long i;
    double T0, PrevT, T, deltaT;

    double VT0, G, H;
    double d1, d2;
    double value;

    if (kappa < Tiny_Value) kappa = Tiny_Value;

    for (i = 0; i < nVol; i++) {
        if (Vol[i] < 0.0) Vol[i] = -Vol[i];
        if (Vol[i] < Tiny_Value) Vol[i] = Tiny_Value;
    }

    T0 = (double)MaturityDate / 365.0;

    VT0 = exp(-2.0 * kappa * T0) * Integ(T0, 1.0, 2.0 * kappa, tVol, Vol, nVol);

    G = PT[nDates];
    PrevT = T0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;
        G += StrikeRate * deltaT * PT[i + 1];
        PrevT = T;
    }
    G /= PT[0];

    H = PT[nDates] * B(T0, (double)Dates[nDates - 1] / 365.0, kappa);
    PrevT = T0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;
        H += StrikeRate * deltaT * PT[i + 1] * B(T0, T, kappa);
        PrevT = T;
    }
    H /= G * PT[0];

    d1 = -log(G) / (H * sqrt(VT0)) + 0.5 * H * sqrt(VT0);
    d2 = -log(G) / (H * sqrt(VT0)) - 0.5 * H * sqrt(VT0);

    value = PT[0] * (CDF_N(d1) - G * CDF_N(d2));


    return NA * value;
}

// 1-factor 모형의 Cap 가격 계산
double HW_Cap(
    double NA,            // 액면금액
    double kappa,        // 회귀속도 
    double* tVol,        // 변동성 구간 종점
    double* Vol,            // 구간 변동성
    long nVol,            // 변동성 구간 개수
    double* t,            // 할인채 만기
    double* r,            // 할인채 가격
    long nr,                // 할인채 개수
    double StrikeRate,    // 행사금리
    long* Dates,            // 지급일: 계산일로부터 각 지급일까지의 일수
    long nDates            // 지급 회수(계산일 이후 남은 회수)
)
{
    long i;
    double d1, d2, PrevT, T, deltaT;
    double PrevDisc, Disc;
    double BP, u;
    double value;

    if (kappa < Tiny_Value) kappa = Tiny_Value;

    for (i = 0; i < nVol; i++) {
        if (Vol[i] < 0.0) Vol[i] = -Vol[i];
        if (Vol[i] < Tiny_Value) Vol[i] = Tiny_Value;
    }

    PrevT = 0.0;
    value = 0.0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;

        PrevDisc = Calc_Discount_Factor(t, r, nr, PrevT);
        Disc = Calc_Discount_Factor(t, r, nr, T);

        BP = B(PrevT, T, kappa);
        u = BP * BP * exp(-2.0 * kappa * PrevT) * Integ(PrevT, 1.0, 2.0 * kappa, tVol, Vol, nVol);

        if (i == 0) {
            if (PrevDisc > (1.0 + StrikeRate * deltaT) * Disc) value = PrevDisc - (1.0 + StrikeRate * deltaT) * Disc;
        }
        else {
            d1 = -log((1.0 + StrikeRate * deltaT) * Disc / PrevDisc) / sqrt(u) + 0.5 * sqrt(u);
            d2 = d1 - sqrt(u);

            value += PrevDisc * CDF_N(d1) - (1.0 + StrikeRate * deltaT) * Disc * CDF_N(d2);
        }

        PrevT = T;
    }

    return NA * value;
}

// 1-factor 모형의 Cap 가격 계산
double HW_Cap(
    double NA,            // 액면금액
    double kappa,        // 회귀속도 
    double* tVol,        // 변동성 구간 종점
    double* Vol,            // 구간 변동성
    long nVol,            // 변동성 구간 개수
    double* t,            // 할인채 만기
    double* r,            // 할인채 가격
    long nr,                // 할인채 개수
    double StrikeRate,    // 행사금리
    long* Dates,            // 지급일: 계산일로부터 각 지급일까지의 일수
    long nDates,            // 지급 회수(계산일 이후 남은 회수)
    double* PT
)
{
    long i;
    double d1, d2, PrevT, T, deltaT;
    double PrevDisc, Disc;
    double BP, u;
    double value;

    if (kappa < Tiny_Value) kappa = Tiny_Value;

    for (i = 0; i < nVol; i++) {
        if (Vol[i] < 0.0) Vol[i] = -Vol[i];
        if (Vol[i] < Tiny_Value) Vol[i] = Tiny_Value;
    }

    PrevT = 0.0;
    value = 0.0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        deltaT = T - PrevT;

        PrevDisc = PT[i];
        Disc = PT[i + 1];

        BP = B(PrevT, T, kappa);
        u = BP * BP * exp(-2.0 * kappa * PrevT) * Integ(PrevT, 1.0, 2.0 * kappa, tVol, Vol, nVol);

        if (i == 0) {
            if (PrevDisc > (1.0 + StrikeRate * deltaT) * Disc) value = PrevDisc - (1.0 + StrikeRate * deltaT) * Disc;
        }
        else {
            d1 = -log((1.0 + StrikeRate * deltaT) * Disc / PrevDisc) / sqrt(u) + 0.5 * sqrt(u);
            d2 = d1 - sqrt(u);

            value += PrevDisc * CDF_N(d1) - (1.0 + StrikeRate * deltaT) * Disc * CDF_N(d2);
        }

        PrevT = T;
    }

    return NA * value;
}


double BS_Swaption(
    double NA,
    double Vol,
    double StrikePrice,
    double* Term,
    double* Rate,
    long NTerm,
    double T_Option,
    double Tenor,
    double FreqMonth
)
{
    long i;
    long NDates;
    double dT = FreqMonth / 12.0;
    double Annuity;
    double ForwardSwapRate;
    double d1, d2;
    double Ti;

    ForwardSwapRate = FSR(Term, Rate, NTerm, T_Option, Tenor, FreqMonth);

    NDates = Number_Of_Payment(T_Option, T_Option + Tenor, FreqMonth);
    Annuity = 0.0;
    for (i = 0; i < NDates; i++)
    {
        Ti = T_Option + Tenor - (double)i * dT;
        Annuity += dT * Calc_Discount_Factor(Term, Rate, NTerm, Ti);
    }

    d1 = (log(ForwardSwapRate / StrikePrice) + 0.5 * Vol * Vol * T_Option) / (Vol * sqrt(T_Option));
    d2 = d1 - Vol * sqrt(T_Option);

    double value = Annuity * (ForwardSwapRate * CDF_N(d1) - StrikePrice * CDF_N(d2));
    return NA * value;
}

double BS_Cap(
    double NA,            // 액면금액
    double Vol,            // 변동성
    double* t,            // 할인채 만기
    double* r,            // 할인채 가격
    long nr,                // 할인채 개수
    double StrikeRate,    // 행사금리
    double Swaption_Mat,        // 입력: 각 스왑션(캡)의 만기(연환산) 
    double Swap_Mat,            // 입력: 각 스왑션(캡)의 기초자산(스왑) 만기(연환산)
    double Swap_Period        // 입력: 각 스왑션의 기초자산인 스왑(캡) 지급 주기(월환산)
)
{
    long i;
    double d1, d2, PrevT, T, delta;
    double PrevDisc, Disc;
    double F, value;
    long nDates;
    long* Dates;

    nDates = Number_Of_Payment(Swaption_Mat, Swaption_Mat + Swap_Mat, Swap_Period);
    Dates = (long*)malloc(sizeof(long) * (nDates));
    PaymentDatesMapping(Swaption_Mat + Swap_Mat, Swap_Period, Dates, nDates);

    PrevT = 0.0;
    value = 0.0;
    for (i = 0; i < nDates; i++) {
        T = (double)Dates[i] / 365.0;
        delta = T - PrevT;

        PrevDisc = Calc_Discount_Factor(t, r, nr, PrevT);
        Disc = Calc_Discount_Factor(t, r, nr, T);

        F = (PrevDisc / Disc - 1.0) / delta;        // 선도금리

        if (i == 0) {    // 첫번째 지급일 
            if (F > StrikeRate) value = delta * Disc * (F - StrikeRate);
        }
        else {
            d1 = (log(F / StrikeRate) + 0.5 * Vol * Vol * PrevT) / (Vol * sqrt(PrevT));
            d2 = d1 - Vol * sqrt(PrevT);

            value += delta * Disc * (F * CDF_N(d1) - StrikeRate * CDF_N(d2));
        }

        PrevT = T;
    }
    if (Dates) free(Dates);

    return NA * value;
}

void Calc_Error_Array(
    long nSwapMaturity,
    double* HW_Vol,
    double* tHW_Vol,
    long nHW_Vol,
    double* temp_vol,

    double* t,
    double* r,
    long nr,

    double* Swaption_Mat,
    double* Swap_Mat,
    double* Swap_Period,
    double* Strike_Rate,
    double* Swaption_Price,

    long Number_KappaRange,
    double* KappaRange,
    long Number_VolRange,
    double* VolRange,
    double** ErrorArray
)
{
    long i, j, k;
    long p, q;
    long nodet;
    long nDates;
    long* dates;
    double kappa;
    double Price;
    double* PT;
    long MaturityDate;

    k = 0;
    for (j = 0; j < nSwapMaturity; j++)
    {
        // Swaption의 경우 NodeT계산
        if (Swaption_Mat[k] > 0.0)
        {
            for (nodet = 0; nodet < nHW_Vol - 1; nodet++) {
                if (Swaption_Mat[k] >= tHW_Vol[nodet] && Swaption_Mat[k] < tHW_Vol[nodet + 1]) {
                    break;
                }
            }
        }
        else // Cap의 경우 NodeT계산
        {
            for (nodet = 0; nodet < nHW_Vol - 1; nodet++) {
                if (Swap_Mat[k] >= tHW_Vol[nodet] && Swap_Mat[k] < tHW_Vol[nodet + 1]) {
                    break;
                }
            }
        }

        nDates = Number_Of_Payment(Swaption_Mat[k], Swaption_Mat[k] + Swap_Mat[k], Swap_Period[k]);
        dates = (long*)malloc(nDates * sizeof(long));
        PaymentDatesMapping(Swaption_Mat[k] + Swap_Mat[k], Swap_Period[k], dates, nDates);
        MaturityDate = (long)(Swaption_Mat[k] * 365.0);

        PT = (double*)malloc(sizeof(double) * (nDates + 1));
        PT[0] = Calc_Discount_Factor(t, r, nr, (double)MaturityDate / 365.0);

        for (i = 0; i < nDates; i++) {
            PT[i + 1] = Calc_Discount_Factor(t, r, nr, (double)dates[i] / 365.0);
        }

        for (p = 0; p < nodet; p++) {
            temp_vol[p] = HW_Vol[p];
        }

        for (p = 0; p < Number_KappaRange; p++) {
            for (q = 0; q < Number_VolRange; q++) {
                kappa = KappaRange[p];
                temp_vol[nodet] = VolRange[q];
                if (Swaption_Mat[k] > 0.0) {
                    Price = HW_Swaption(Notional_Amount, kappa, tHW_Vol, temp_vol, nodet + 1, t, r, nr, Strike_Rate[k], MaturityDate, dates, nDates, PT);
                }
                else {
                    Price = HW_Cap(Notional_Amount, kappa, tHW_Vol, temp_vol, nodet + 1, t, r, nr, Strike_Rate[k], dates, nDates);
                }
                ErrorArray[p][q] = ErrorArray[p][q] + fabs(Price - Swaption_Price[k]);
            }
        }


        k += 1;
        free(dates);
        free(PT);
    }
}

double OptimizeHWCap(
    long nCap,
    double* t,
    double* r,
    long nr,
    long CaliFlag,
    double* Cap_Mat,
    double* Cap_Strike,
    double* Cap_Tenor,
    double* Cap_Vol,
    double* Cap_Price,
    long na,
    long* ia,
    double* a,
    double* tHW_Vol,
    double* HW_Vol,
    long nHW_Vol,
    long* Cap_ErrorShape,
    double*** Distance_Array,
    long HW_Kappa_Flag,
    double Fixed_Kappa
)
{
    long i, j, p, q;
    long N, M;
    long N2, M2;

    double* Params1_Array;
    double* Params2_Array;

    double* Params1_Array2;
    double* Params2_Array2;
    double** Distance_Array2;

    double Param1_Min2;
    double Param1_Max2;
    double Param2_Min2;
    double Param2_Max2;

    double minvalue;
    double kappa;
    long break_flag;


    double* Swaption_Mat = (double*)malloc(sizeof(double) * nCap);          // 할당 1
    for (i = 0; i < nCap; i++)
        Swaption_Mat[i] = 0.0;

    for (i = 0; i < nCap; i++)
        tHW_Vol[i] = Cap_Mat[i];

    double* ResultKappa;
    double* temp_vol;
    long nSwapMaturity = 1;

    N = (long)((Param1_Max - Param1_Min) / Interval1) + 1;
    M = (long)((Param2_Max - Param2_Min) / Interval2) + 1;

    Params1_Array = (double*)malloc(sizeof(double) * N);                    // 할당 2
    Params2_Array = (double*)malloc(sizeof(double) * M);                    // 할당 3

    Cap_ErrorShape[0] = nCap;
    Cap_ErrorShape[1] = N;
    Cap_ErrorShape[2] = M;

    for (i = 0; i < nCap; i++) {
        Distance_Array[i] = (double**)malloc(sizeof(double*) * N);          // Result 할당부분
        for (j = 0; j < N; j++)
        {
            Distance_Array[i][j] = (double*)malloc(sizeof(double) * M);
        }
    }

    if (HW_Kappa_Flag != 1) {
        for (i = 0; i < N; i++)
            Params1_Array[i] = Param1_Min + (double)i * Interval1;
    }
    else {
        for (i = 0; i < N; i++)
            Params1_Array[i] = Fixed_Kappa;
    }

    for (i = 0; i < M; i++)
        Params2_Array[i] = Param2_Min + (double)i * Interval2;

    ResultKappa = (double*)calloc(nCap, sizeof(double));                    // 할당 4
    temp_vol = (double*)calloc(nCap, sizeof(double));                       // 할당 5

    for (i = 0; i < nCap; i++)
    {
        for (p = 0; p < N; p++)
            for (q = 0; q < M; q++)
                Distance_Array[i][p][q] = 0.0;

        Calc_Error_Array(nSwapMaturity, HW_Vol, tHW_Vol, nCap, temp_vol,
            t, r, nr, Swaption_Mat + i * nSwapMaturity, Cap_Mat + i * nSwapMaturity,
            Cap_Tenor + i * nSwapMaturity, Cap_Strike + i * nSwapMaturity, Cap_Price + i * nSwapMaturity, N, Params1_Array,
            M, Params2_Array, Distance_Array[i]);

        minvalue = Distance_Array[i][0][0];
        for (p = 1; p < N; p++)
            for (q = 1; q < M; q++)
                minvalue = min(minvalue, Distance_Array[i][p][q]);

        for (p = 0; p < N; p++) {
            break_flag = 0;
            for (q = 0; q < M; q++) {
                if (Distance_Array[i][p][q] == minvalue) {
                    break_flag = 1;
                    break;
                }
            }
            if (break_flag == 1)
                break;
        }

        Param1_Min2 = Params1_Array[max(0, p - 1)];
        Param1_Max2 = Params1_Array[min(N - 1, p + 1)];

        Param2_Min2 = Params2_Array[max(0, q - 1)];
        Param2_Max2 = Params2_Array[min(M - 1, q + 1)];

        N2 = (long)((Param1_Max2 - Param1_Min2) / (Interval1 / 10.0)) + 1;
        M2 = (long)((Param2_Max2 - Param2_Min2) / (Interval2 / 10.0)) + 1;

        Params1_Array2 = (double*)malloc(sizeof(double) * N2);                   // 반복문 내부 할당 1
        Params2_Array2 = (double*)malloc(sizeof(double) * M2);                   // 반복문 내부 할당 2

        Distance_Array2 = (double**)calloc(N2, sizeof(double*));                 // 반복문 내부 할당 3
        for (p = 0; p < N2; p++)
            Distance_Array2[p] = (double*)calloc(M2, sizeof(double));

        if (HW_Kappa_Flag != 1)
        {
            for (p = 0; p < N2; p++) {
                Params1_Array2[p] = Param1_Min2 + (double)p * Interval1 / 10.0;
            }
        }
        else
        {
            for (p = 0; p < N2; p++) {
                Params1_Array2[p] = Fixed_Kappa;
            }
        }

        for (p = 0; p < M2; p++) {
            Params2_Array2[p] = Param2_Min2 + (double)p * Interval2 / 10.0;
        }

        Calc_Error_Array(nSwapMaturity, HW_Vol, tHW_Vol, nCap, temp_vol,
            t, r, nr, Swaption_Mat + i * nSwapMaturity, Cap_Mat + i * nSwapMaturity,
            Cap_Tenor + i * nSwapMaturity, Cap_Strike + i * nSwapMaturity, Cap_Price + i * nSwapMaturity, N2, Params1_Array2,
            M2, Params2_Array2, Distance_Array2);

        minvalue = Distance_Array2[0][0];
        for (p = 1; p < N2; p++)
            for (q = 1; q < M2; q++)
                minvalue = min(minvalue, Distance_Array2[p][q]);

        for (p = 0; p < N2; p++)
        {
            break_flag = 0;
            for (q = 0; q < M2; q++)
            {
                if (Distance_Array2[p][q] == minvalue)
                {
                    break_flag = 1;
                    break;
                }
            }
            if (break_flag == 1)
                break;
        }

        ResultKappa[i] = Params1_Array2[p];
        HW_Vol[i] = Params2_Array2[q];

        free(Params1_Array2);
        free(Params2_Array2);
        for (p = 0; p < N2; p++)
        {
            free(Distance_Array2[p]);
        }
        free(Distance_Array2);

    }


    kappa = 0.0;
    for (i = 0; i < nCap; i++)
    {
        kappa += ResultKappa[i] / (double)nCap;
    }
    a[0] = kappa;

    free(Swaption_Mat);
    free(Params1_Array);
    free(Params2_Array);

    free(ResultKappa);
    free(temp_vol);
    return 1.0;
}

double OptimizeHWSwaption(
    long nSwaption,
    long nCap,
    long* Cap_ErrorShape,
    double*** Cap_ErrorArray,
    double* Cap_Mat,
    double* t,                    //
    double* r,                    // 입력: 무위험 금리 기간 구조
    long nr,                    //
    long CaliFlag,               // 0 Swaption Cap 둘다 1 Cap 2 Swaption
    double* Swaption_Mat,        // 입력: 각 스왑션(캡)의 만기(연환산) 
    double* Swap_Mat,            // 입력: 각 스왑션(캡)의 기초자산(스왑) 만기(연환산)
    double* Strike_Rate,        // 입력: 각 스왑션(캡)의 행사가격
    double* Swap_Period,        // 입력: 각 스왑션의 기초자산인 스왑(캡) 지급 주기(월환산)
    double* Swaption_Vol,        // 입력: 각 스왑션(캡)의 변동성
    double* Swaption_Price,
    long na,
    long* ia,
    double* a,
    double* tHW_Vol,            // 입력: HW 변동성 기간 구조 시점
    double* HW_Vol,            // 입력/결과: HW 변동성
    long nHW_Vol,                // 입력: HW 변동성 기간 구조 개수
    long HW_Kappa_Flag,
    double Fixed_Kappa
)
{
    long i, j, k, k2;
    long p, q;
    long N, M;
    long N2, M2;

    long nMaturity;
    long nOption;

    double* Params1_Array;
    double* Params2_Array;
    double** Distance_Array;

    double* Params1_Array2;
    double* Params2_Array2;
    double** Distance_Array2;

    double Param1_Min2;
    double Param1_Max2;
    double Param2_Min2;
    double Param2_Max2;

    double minvalue;
    double kappa;
    long break_flag;

    double* ResultKappa;
    double* temp_vol;

    if (CaliFlag == 1)
        nMaturity = nSwaption / nHW_Vol;
    else if (CaliFlag == 0)
        nMaturity = nSwaption / (nHW_Vol - nCap);

    nOption = nSwaption / nMaturity;

    N = (long)((Param1_Max - Param1_Min) / Interval1) + 1;
    M = (long)((Param2_Max - Param2_Min) / Interval2) + 1;

    Params1_Array = (double*)malloc(sizeof(double) * N);          // 할당 1
    Params2_Array = (double*)malloc(sizeof(double) * M);          // 할당 2

    Distance_Array = (double**)malloc(sizeof(double*) * N);       // 2차원 할당3
    for (i = 0; i < N; i++)
        Distance_Array[i] = (double*)malloc(sizeof(double) * M);
    if (HW_Kappa_Flag != 1)
    {
        for (i = 0; i < N; i++)
        {
            Params1_Array[i] = Param1_Min + (double)i * Interval1;
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            Params1_Array[i] = Fixed_Kappa;
        }
    }

    for (i = 0; i < M; i++)
        Params2_Array[i] = Param2_Min + (double)i * Interval2;

    ResultKappa = (double*)calloc(nHW_Vol, sizeof(double));     // 할당4
    temp_vol = (double*)calloc(nHW_Vol, sizeof(double));       // 할당5

    k = 0;
    k2 = 0;
    for (i = 0; i < nOption; i++)
    {

        for (p = 0; p < N; p++) {
            for (q = 0; q < M; q++) {
                Distance_Array[p][q] = 0.0;
            }
        }

        if (CaliFlag == 0)
        {
            for (j = 0; j < nCap; j++)
            {
                if (tHW_Vol[i] == Cap_Mat[j])
                {
                    for (p = 0; p < N; p++) {
                        for (q = 0; q < M; q++)
                            Distance_Array[p][q] = Cap_ErrorArray[j][p][q];
                    }
                    break;
                }
            }
        }

        Calc_Error_Array(nMaturity, HW_Vol, tHW_Vol, nOption, temp_vol,
            t, r, nr, Swaption_Mat + i * nMaturity, Swap_Mat + i * nMaturity,
            Swap_Period + i * nMaturity, Strike_Rate + i * nMaturity, Swaption_Price + i * nMaturity, N, Params1_Array,
            M, Params2_Array, Distance_Array);

        minvalue = Distance_Array[0][0];
        for (p = 1; p < N; p++)
            for (q = 1; q < M; q++)
                minvalue = min(minvalue, Distance_Array[p][q]);

        for (p = 0; p < N; p++) {
            break_flag = 0;
            for (q = 0; q < M; q++) {
                if (Distance_Array[p][q] == minvalue) {
                    break_flag = 1;
                    break;
                }
            }
            if (break_flag == 1)
                break;
        }

        Param1_Min2 = Params1_Array[max(0, p - 1)];
        Param1_Max2 = Params1_Array[min(N - 1, p + 1)];

        Param2_Min2 = Params2_Array[max(0, q - 1)];
        Param2_Max2 = Params2_Array[min(M - 1, q + 1)];

        N2 = (long)((Param1_Max2 - Param1_Min2) / (Interval1 / 10.0)) + 1;
        M2 = (long)((Param2_Max2 - Param2_Min2) / (Interval2 / 10.0)) + 1;

        Params1_Array2 = (double*)malloc(sizeof(double) * N2);                   // 반복문 내부 할당 1
        Params2_Array2 = (double*)malloc(sizeof(double) * M2);                   // 반복문 내부 할당 2

        Distance_Array2 = (double**)calloc(N2, sizeof(double*));                 // 반복문 내부 할당 3 (2차원)
        for (p = 0; p < N2; p++)
            Distance_Array2[p] = (double*)calloc(M2, sizeof(double));

        if (HW_Kappa_Flag != 1)
        {
            for (p = 0; p < N2; p++) {
                Params1_Array2[p] = Param1_Min2 + (double)p * Interval1 / 10.0;
            }
        }
        else
        {
            for (p = 0; p < N2; p++) {
                Params1_Array2[p] = Fixed_Kappa;
            }
        }

        for (p = 0; p < M2; p++) {
            Params2_Array2[p] = Param2_Min2 + (double)p * Interval2 / 10.0;
        }

        Calc_Error_Array(nMaturity, HW_Vol, tHW_Vol, nOption, temp_vol,
            t, r, nr, Swaption_Mat + i * nMaturity, Swap_Mat + i * nMaturity,
            Swap_Period + i * nMaturity, Strike_Rate + i * nMaturity, Swaption_Price + i * nMaturity, N2, Params1_Array2,
            M2, Params2_Array2, Distance_Array2);

        minvalue = Distance_Array2[0][0];
        for (p = 1; p < N2; p++)
            for (q = 1; q < M2; q++)
                minvalue = min(minvalue, Distance_Array2[p][q]);

        for (p = 0; p < N2; p++)
        {
            break_flag = 0;
            for (q = 0; q < M2; q++)
            {
                if (Distance_Array2[p][q] == minvalue)
                {
                    break_flag = 1;
                    break;
                }
            }
            if (break_flag == 1)
                break;
        }

        ResultKappa[i] = Params1_Array2[p];
        HW_Vol[i] = Params2_Array2[q];
        free(Params1_Array2);               // 반복문 내부 할당 해제 1
        free(Params2_Array2);               // 반복문 내부 할당 해제 2
        for (p = 0; p < N2; p++)
        {
            free(Distance_Array2[p]);
        }
        free(Distance_Array2);              // 반복문 내부 할당 해제 3 (2차원)

    }

    kappa = 0.0;
    for (i = 0; i < nOption; i++)
    {
        kappa += ResultKappa[i] / (double)nOption;
    }
    a[0] = kappa;

    free(Params1_Array);
    free(Params2_Array);
    for (i = 0; i < N; i++)
        free(Distance_Array[i]);
    free(Distance_Array);


    free(temp_vol);
    free(ResultKappa);

    return 1.0;
}

long Sort_VolData(
    long nHWVol_Swaption,
    double* tHWVol_Swaption,
    double* HWVol_Swaption,
    long nHWVol_Cap,
    double* tHWVol_Cap,
    double* HWVol_Cap,
    double* tHWVol,
    double* HWVol,
    long nHWVol
)
{
    long i;
    long j;
    long k;
    long nSame;

    // 중복되는 Term 개수 Count
    long N_Time_Equal;
    N_Time_Equal = 0;
    for (i = 0; i < nHWVol_Swaption; i++)
        for (j = 0; j < nHWVol_Cap; j++)
        {
            if (tHWVol_Swaption[i] == tHWVol_Cap[j])
            {
                N_Time_Equal += 1;
                break;
            }
        }

    // Term 겹치는 중복 개수 제거
    long NVol = nHWVol_Swaption + nHWVol_Cap - N_Time_Equal;
    double* Term_Vol = (double*)malloc(sizeof(double) * NVol);
    double* Vol = (double*)malloc(sizeof(double) * NVol);

    for (i = 0; i < nHWVol_Swaption; i++)
    {
        Term_Vol[i] = tHWVol_Swaption[i];
        Vol[i] = HWVol_Swaption[i];
    }

    k = nHWVol_Swaption;
    for (i = 0; i < nHWVol_Cap; i++)
    {
        nSame = 0;
        for (j = 0; j < nHWVol_Swaption; j++)
        {
            if (tHWVol_Cap[i] == tHWVol_Swaption[j])
            {
                nSame += 1;
                break;
            }
        }

        if (nSame == 0)
        {
            Term_Vol[k] = tHWVol_Cap[i];
            Vol[k] = HWVol_Cap[i];
            k += 1;
        }
    }

    double temp_idx;
    double temp_value;
    for (i = 0; i < NVol; i++)
    {
        for (j = 0; j < NVol - 1; j++)
        {
            if (Term_Vol[j] > Term_Vol[j + 1])
            {
                temp_idx = Term_Vol[j];
                temp_value = Vol[i];
                Term_Vol[j] = Term_Vol[j + 1];
                Vol[j] = Vol[j + 1];
                Term_Vol[j + 1] = temp_idx;
                Vol[j + 1] = temp_value;
            }
        }
    }

    for (i = 0; i < NVol; i++)
    {
        tHWVol[i] = Term_Vol[i];
        HWVol[i] = Vol[i];
    }

    for (i = NVol; i < nHWVol; i++)
    {
        tHWVol[i] = tHWVol[i - 1] + 1.0e-7;
        HWVol[i] = HWVol[i - 1];
    }

    free(Term_Vol);
    free(Vol);
    return NVol;
}

long Calibration(
    double* t,                    //
    double* r,                    // 입력: 무위험 금리 기간 구조
    long nr,                    //
    long CaliFlag,               // 0 Swaption Cap 둘다 1 Cap 2 Swaption
    double* Swaption_Mat,        // 입력: 각 스왑션(캡)의 만기(연환산) 
    double* Swap_Mat,            // 입력: 각 스왑션(캡)의 기초자산(스왑) 만기(연환산)
    double* Strike_Rate,        // 입력: 각 스왑션(캡)의 행사가격
    double* Swap_Period,        // 입력: 각 스왑션의 기초자산인 스왑(캡) 지급 주기(월환산)
    double* Swaption_Vol,        // 입력: 각 스왑션(캡)의 변동성
    long nSwaption,                // 입력: 스왑션(캡) 개수
    double* Cap_Mat,            // 입력: 각 캡의 만기(연환산)
    double* Cap_Tenor,            // 입력: 각 캡의 지급 주기(월환산)
    double* Cap_Strike,
    double* Cap_Vol,            // 입력: 각 캡의 변동성
    long nCap,                    // 입력: 캡 개수
    double* HW_Kappa,            // 결과: HW 모형 모수 Kappa
    long HW_Kappa_Flag,            // 입력: HW_Kappa 고정 여부(0:고정, 1:변동)
    double Fixed_Kappa,
    double* tHW_Vol,            // 입력: HW 변동성 기간 구조 시점
    double* HW_Vol,            // 입력/결과: HW 변동성
    long nHW_Vol                // 입력: HW 변동성 기간 구조 개수

)
{
    long i;
    long j;

    long ma;

    double* tVol_g = (double*)malloc(sizeof(double) * (nHW_Vol));

    for (i = 0; i < nHW_Vol; i++) tVol_g[i] = tHW_Vol[i];
    long nVol_g = nHW_Vol;

    // BS 모형으로 스왑션(캡) 가격을 계산
    double* Swaption_Price = (double*)malloc(sizeof(double) * (nSwaption));
    double* Cap_Price = (double*)malloc(sizeof(double) * (nCap));

    for (i = 0; i < nSwaption; i++) {

        if (Swaption_Mat[i] > 0.0) {
            Swaption_Price[i] = BS_Swaption(Notional_Amount, Swaption_Vol[i], Strike_Rate[i], t, r, nr, Swaption_Mat[i], Swap_Mat[i], Swap_Period[i]);

        }
        else {
            Swaption_Price[i] = BS_Cap(Notional_Amount, Swaption_Vol[i], t, r, nr,
                Strike_Rate[i], Swaption_Mat[i], Swap_Mat[i], Swap_Period[i]);
        }

    }


    for (i = 0; i < nCap; i++) {
        Cap_Price[i] = BS_Cap(Notional_Amount, Cap_Vol[i], t, r, nr, Cap_Strike[i], 0.0, Cap_Mat[i], Cap_Tenor[i]);

    }

    ma = nHW_Vol + 1;
    double* a = (double*)malloc(sizeof(double) * (ma));
    long* ia = (long*)malloc(sizeof(long) * (ma));

    ia[0] = HW_Kappa_Flag;
    a[0] = *HW_Kappa;
    for (i = 1; i < ma; i++) {
        ia[i] = 1;
        a[i] = HW_Vol[i - 1];
    }

    double Temp;
    double*** Cap_ErrorArray = (double***)malloc(sizeof(double**) * max(nCap, 1));
    long Cap_ErrorShape[3] = { 0,0,0 };
    long nSwaptionOption = 0;

    if (CaliFlag == 0)
    {
        nSwaptionOption = (nHW_Vol - nCap);
    }


    for (i = 0; i < nCap; i++) {
        tHW_Vol[i + nSwaptionOption] = Cap_Mat[i];
    }

    if (CaliFlag == 0 || CaliFlag == 2)
    {
        Temp = OptimizeHWCap(nCap, t, r, nr, CaliFlag,
            Cap_Mat, Cap_Strike, Cap_Tenor, Cap_Vol, Cap_Price,
            ma, ia, a, tHW_Vol + nSwaptionOption, HW_Vol + nSwaptionOption,
            nCap, Cap_ErrorShape, Cap_ErrorArray, HW_Kappa_Flag, Fixed_Kappa);

        if (CaliFlag == 0)
        {
            Temp = OptimizeHWSwaption(nSwaption, nCap, Cap_ErrorShape, Cap_ErrorArray, Cap_Mat, t, r, nr, CaliFlag,
                Swaption_Mat, Swap_Mat, Strike_Rate, Swap_Period, Swaption_Vol,
                Swaption_Price, ma, ia, a, tHW_Vol, HW_Vol, nHW_Vol, HW_Kappa_Flag, Fixed_Kappa);
        }

        for (i = 0; i < Cap_ErrorShape[0]; i++)
        {
            for (j = 0; j < Cap_ErrorShape[1]; j++)
            {
                if (Cap_ErrorArray[i][j]) free(Cap_ErrorArray[i][j]);
            }
            if (Cap_ErrorArray[i]) free(Cap_ErrorArray[i]);
        }

    }

    if (CaliFlag == 1)
    {
        Temp = OptimizeHWSwaption(nSwaption, nCap, Cap_ErrorShape, Cap_ErrorArray, Cap_Mat, t, r, nr, CaliFlag,
            Swaption_Mat, Swap_Mat, Strike_Rate, Swap_Period, Swaption_Vol,
            Swaption_Price, ma, ia, a, tHW_Vol, HW_Vol, nHW_Vol, HW_Kappa_Flag, Fixed_Kappa);

    }

    if (CaliFlag == 0)
        nHW_Vol = Sort_VolData(nSwaptionOption, tHW_Vol, HW_Vol, nCap, tHW_Vol + nSwaptionOption, HW_Vol + nSwaptionOption, tHW_Vol, HW_Vol, nHW_Vol);

    HW_Kappa[0] = a[0];


    free(tVol_g);
    free(Swaption_Price);
    free(Cap_Price);
    free(a);
    free(ia);
    if (Cap_ErrorArray) free(Cap_ErrorArray);
    return 1;
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

long ErrorCheck_HW_Calibration(
    double* t,                    //
    double* r,                    // 입력: 무위험 금리 기간 구조
    long nr,                        //
    long CaliFlag,               // 0 Swaption Cap 둘다 1 Swaption 2 Cap
    double* Cap_Mat,            // 입력: 각 캡의 만기(연환산)
    double Caplet_Freq,            // 입력: 각 캡의 지급 주기(월환산)
    double* Cap_Vol,            // 입력: 각 캡의 변동성
    long nCap,                    // 입력: 캡 개수
    double* Swaption_Mat,        // 입력: 각 스왑션의 만기(연환산) 
    double* Swap_Mat,            // 입력: 각 스왑 만기(연환산)
    double* Swaption_Vol,        // 입력: 각 스왑션(캡)의 변동성
    double SwapPay_Freq,
    long nSwaption,                // 입력: 스왑션(캡) 개수
    long HW_Kappa_Flag,            // 입력: HW_Kappa 고정 여부(0:고정, 1:변동)
    double Fixed_Kappa,
    long nHW_Vol,                // 입력: HW 변동성 기간 구조 개수
    long ErrorCheckFlag,
    char* Error
)
{
    long i;

    char ErrorName[100];

    if (nCap + nSwaption < 0 || nHW_Vol < 1)
    {
        strcpy_s(ErrorName, "Check nCap nSwaption nHWVol");
        return SaveErrorName(Error, ErrorName);
    }

    if (ErrorCheckFlag > 0)
    {
        for (i = 1; i < nr; i++)
        {
            if (t[i] < t[i - 1])
            {
                strcpy_s(ErrorName, "Curve t1 is bigger than t2");
                return SaveErrorName(Error, ErrorName);
            }
        }

        if (CaliFlag < 0 || CaliFlag > 2)
        {
            strcpy_s(ErrorName, "Check CaliFlag");
            return SaveErrorName(Error, ErrorName);
        }

        if (nCap < 0)
        {
            strcpy_s(ErrorName, "nCap is zero");
            return SaveErrorName(Error, ErrorName);
        }

        for (i = 1; i < nCap; i++)
        {
            if (Cap_Mat[i] < Cap_Mat[i - 1]) {
                strcpy_s(ErrorName, "Cap Mat t1 is bigger than t2");
                return SaveErrorName(Error, ErrorName);
            }
            if (Cap_Vol[i] < 0.0)
            {
                strcpy_s(ErrorName, "CapVol must be Positive");
                return SaveErrorName(Error, ErrorName);
            }
        }
        if (Cap_Vol[0] < 0.0)
        {
            strcpy_s(ErrorName, "CapVol must be Positive");
            return SaveErrorName(Error, ErrorName);
        }

        if (Cap_Mat[0] < (double)Caplet_Freq / 12.0)
        {
            strcpy_s(ErrorName, "First Cap term must be bigger than Cap Frequency");
            return SaveErrorName(Error, ErrorName);
        }

        if (nSwaption < 0)
        {
            strcpy_s(ErrorName, "nSwaption must be Positive");
            return SaveErrorName(Error, ErrorName);
        }

        for (i = 1; i < nSwaption; i++)
        {
            if (Swaption_Vol[i] < 0.0)
            {
                strcpy_s(ErrorName, "SwaptionVol must be Positive");
                return SaveErrorName(Error, ErrorName);
            }
        }
        if (Swaption_Vol[0] < 0.0)
        {
            strcpy_s(ErrorName, "SwaptionVol must be Positive");
            return SaveErrorName(Error, ErrorName);
        }

        if (HW_Kappa_Flag < 0 || HW_Kappa_Flag > 2)
        {
            strcpy_s(ErrorName, "Check the Kappa Flag");
            return SaveErrorName(Error, ErrorName);
        }

    }


    return 1;
}

// 1-factor 모형의 Calibration 함수
// ATM 변동성만 사용
DLLEXPORT(long) HW_Calibration_ATM(
    double* t,                    //
    double* r,                    // 입력: 무위험 금리 기간 구조
    long nr,                        //
    long CaliFlag,               // 0 Swaption Cap 둘다 1 Swaption 2 Cap
    double* Cap_Mat,            // 입력: 각 캡의 만기(연환산)
    double Caplet_Freq,            // 입력: 각 캡의 지급 주기(월환산)
    double* Cap_Vol,            // 입력: 각 캡의 변동성
    long nCap,                    // 입력: 캡 개수
    double* Swaption_Mat,        // 입력: 각 스왑션의 만기(연환산) 
    double* Swap_Mat,            // 입력: 각 스왑 만기(연환산)
    double* Swaption_Vol,        // 입력: 각 스왑션(캡)의 변동성
    double SwapPay_Freq,
    long nSwaption,                // 입력: 스왑션(캡) 개수
    double* HW_Kappa,            // 결과: HW 모형 모수 Kappa
    long HW_Kappa_Flag,            // 입력: HW_Kappa 고정 여부(0:고정, 1:변동)
    double Fixed_Kappa,
    double* tHW_Vol,            // 입력: HW 변동성 기간 구조 시점
    double* HW_Vol,                // 입력/결과: HW 변동성
    long nHW_Vol,                // 입력: HW 변동성 기간 구조 개수

    char* Error
)
{

    long i;
    long Result = 0;
    long ErrorCheckFlag = 1;

    Result = ErrorCheck_HW_Calibration(t, r, nr, CaliFlag, Cap_Mat, Caplet_Freq, Cap_Vol, nCap, Swaption_Mat, Swap_Mat, Swaption_Vol, SwapPay_Freq, nSwaption, HW_Kappa_Flag, Fixed_Kappa, nHW_Vol, ErrorCheckFlag, Error);
    if (CaliFlag == 1) {
        nHW_Vol = nHW_Vol - nCap;
        nCap = 0;
    }

    if (CaliFlag == 2) {
        nHW_Vol = nCap;
        nSwaption = 0;
    }

    double* Cap_Tenor = (double*)malloc(sizeof(double) * nCap);             //할당1
    double* Swap_Period = (double*)malloc(nSwaption * sizeof(double));      //할당2
    for (i = 0; i < nCap; i++)
        Cap_Tenor[i] = Caplet_Freq;

    for (i = 0; i < nSwaption; i++)
        Swap_Period[i] = SwapPay_Freq;

    double* Cap_Strike = (double*)malloc(max(nCap, 1) * sizeof(double));     //할당3
    if (CaliFlag != 1)
    {
        for (i = 0; i < nCap; i++) {
            Cap_Strike[i] = FSR(t, r, nr, 0.0, Cap_Mat[i], Cap_Tenor[i]);
        }
    }

    double* Swap_Rate = (double*)malloc(max(nSwaption, 1) * sizeof(double)); //할당4
    if (CaliFlag != 2)
    {
        for (i = 0; i < nSwaption; i++) {
            Swap_Rate[i] = FSR(t, r, nr, Swaption_Mat[i], Swap_Mat[i], Swap_Period[i]);
        }
    }

    if (Result > 0)
    {
        Result = Calibration(t, r, nr, CaliFlag, Swaption_Mat, Swap_Mat, Swap_Rate, Swap_Period, Swaption_Vol, nSwaption, Cap_Mat, Cap_Tenor, Cap_Strike, Cap_Vol, nCap,
            HW_Kappa, HW_Kappa_Flag, Fixed_Kappa, tHW_Vol, HW_Vol, nHW_Vol);
    }

    if (Cap_Tenor) free(Cap_Tenor);
    if (Swap_Period) free(Swap_Period);
    if (Cap_Strike) free(Cap_Strike);
    if (Swap_Rate) free(Swap_Rate);


    //_CrtDumpMemoryLeaks();

    return Result;
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

void Levenberg_Marquardt(long NParams, long NResidual, double* NextParams, double* CurrentParams, double* lambda, double** Jacov, double* Residual, double& ParamSum, double** JT_J, double** Inverse_JT_J, double** JT_Res, double** ResultMatrix)
{
    long i;
    double mu = *lambda;
    double MIN_kappa = 0.001;
    double MIN_Vol = 0.00001;
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
        if (i == 0)
        {
            if (NextParams[0] < MIN_kappa) NextParams[0] = MIN_kappa;
        }
        else
        {
            if (NextParams[i] < MIN_Vol) NextParams[i] = MIN_Vol;
        }
    }

    double s = 0.0;
    for (i = 0; i < NParams; i++) s += fabs(ResultMatrix[i][0]);
    ParamSum = s;
}

void Gradient_Desent(long NParams, long NResidual, double* NextParams, double* CurrentParams, double* lambda, double** Jacov, double* Residual, double& ParamSum, double** JT_Res)
{
    long i;
    double mu = *lambda;
    double MIN_kappa = 0.001;
    double MIN_Vol = 0.00001;
    long n = NResidual, m = NParams;
    long Shape_J[2] = { n,m };


    // J' dot Res               Shape = m * n dot n * 1 = m * 1
    long Shape_JT_Res[2] = { m,1 };
    XprimeY(Jacov, Shape_J, Residual, n, JT_Res);

    for (i = 0; i < NParams; i++)
    {
        NextParams[i] = CurrentParams[i] - 2.0 * mu * JT_Res[i][0];
        if (i == 0)
        {
            if (NextParams[0] < MIN_kappa) NextParams[0] = MIN_kappa;
        }
        else
        {
            if (NextParams[i] < MIN_Vol) NextParams[i] = MIN_Vol;
        }
    }

    double s = 0.0;
    for (i = 0; i < NParams; i++) s += 2.0 * mu * JT_Res[i][0];
    ParamSum = fabs(s);
}

void make_Residual_HWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double& absErrorSum,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    double kappa;
    double* HWVol;
    if (FixedKappaFlag != 1)
    {
        kappa = Params[0];
        HWVol = Params + 1;
    }
    else
    {
        kappa = FixedKappa;
        HWVol = Params;
    }
    double s = 0.0;
    for (i = 0; i < NResidual; i++) TempHWSwaptionPrice[i] = HW_Swaption(1.0, kappa, HWVolTerm, HWVol, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
    for (i = 0; i < NResidual; i++) ResidualArray[i] = (BSSwaptionPrice[i] - TempHWSwaptionPrice[i]);
    for (i = 0; i < NResidual; i++) s += fabs(ResidualArray[i]);
    absErrorSum = s;
}

void make_Residual_HWCapHWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    long* SwaptionFlag,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double& absErrorSum,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    double kappa;
    double* HWVol;
    if (FixedKappaFlag != 1)
    {
        kappa = Params[0];
        HWVol = Params + 1;
    }
    else
    {
        kappa = FixedKappa;
        HWVol = Params;
    }

    double s = 0.0;
    for (i = 0; i < NResidual; i++)
    {
        if (SwaptionFlag[i] == 0) TempHWSwaptionPrice[i] = HW_Cap(1.0, kappa, HWVolTerm, HWVol, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], dates[i], nDates[i], PT[i]);
        else TempHWSwaptionPrice[i] = HW_Swaption(1.0, kappa, HWVolTerm, HWVol, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
    }
    for (i = 0; i < NResidual; i++) ResidualArray[i] = (BSSwaptionPrice[i] - TempHWSwaptionPrice[i]);
    for (i = 0; i < NResidual; i++) s += fabs(ResidualArray[i]);
    absErrorSum = s;
}

void make_Jacov_HWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double* ParamsUp,
    double* ParamsDn,
    double** TempJacovMatrix,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    long j;
    long n;

    double dhwvol_up = 0.0001;
    double dkappa_up = 0.0001;

    double T;
    double HWT;

    double kappa_up;
    double kappa_dn;

    double* HWVol_up;
    double* HWVol_dn;

    double Pup = 0.0;
    double Pdn = 0.0;
    double ErrorUp = 0.0;
    double ErrorDn = 0.0;
    double dParams;
    double temp = 0.0;

    for (i = 0; i < NResidual; i++)
    {
        for (j = 0; j < NParams; j++)
        {
            // 파라미터 초기화
            for (n = 0; n < NParams; n++) {
                ParamsUp[n] = Params[n];
                ParamsDn[n] = Params[n];
            }

            // 파라미터 Up And Dn
            if (FixedKappaFlag != 1)
            {
                if (j == 0)
                {
                    ParamsUp[j] = Params[j] + dkappa_up;
                    ParamsDn[j] = Params[j] - dkappa_up;
                    dParams = dkappa_up;
                }
                else
                {
                    ParamsUp[j] = Params[j] + dhwvol_up;
                    ParamsDn[j] = max(0.0001, Params[j] - dhwvol_up);
                    dParams = dhwvol_up;
                }

                kappa_up = ParamsUp[0];
                kappa_dn = ParamsDn[0];
                HWVol_up = ParamsUp + 1;
                HWVol_dn = ParamsDn + 1;
            }
            else
            {
                ParamsUp[j] = Params[j] + dhwvol_up;
                ParamsDn[j] = max(0.0001, Params[j] - dhwvol_up);
                dParams = dhwvol_up;

                kappa_up = FixedKappa;
                kappa_dn = FixedKappa;
                HWVol_up = ParamsUp;
                HWVol_dn = ParamsDn;
            }

            T = (double)dates[i][nDates[i] - 1] / 365.0;
            HWT = HWVolTerm[j];
            if (HWT > T)
            {
                TempJacovMatrix[i][j] = 0.0;
            }
            else
            {
                Pup = HW_Swaption(1.0, kappa_up, HWVolTerm, HWVol_up, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
                ErrorUp = BSSwaptionPrice[i] - Pup;
                Pdn = HW_Swaption(1.0, kappa_dn, HWVolTerm, HWVol_dn, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
                ErrorDn = BSSwaptionPrice[i] - Pdn;
                TempJacovMatrix[i][j] = (ErrorUp - ErrorDn) / (2.0 * dParams);
            }
        }
    }
}

void make_Jacov_HWCapHWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    long* SwaptionFlag,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double* ParamsUp,
    double* ParamsDn,
    double** TempJacovMatrix,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    long j;
    long n;

    double dhwvol_up = 0.0001;
    double dkappa_up = 0.0001;

    double kappa_up;
    double kappa_dn;

    double* HWVol_up;
    double* HWVol_dn;

    double Pup = 0.0;
    double Pdn = 0.0;
    double ErrorUp = 0.0;
    double ErrorDn = 0.0;
    double dParams;
    double temp = 0.0;

    for (i = 0; i < NResidual; i++)
    {
        for (j = 0; j < NParams; j++)
        {
            // 파라미터 초기화
            for (n = 0; n < NParams; n++) {
                ParamsUp[n] = Params[n];
                ParamsDn[n] = Params[n];
            }

            // 파라미터 Up And Dn
            if (FixedKappaFlag != 1)
            {
                if (j == 0)
                {
                    ParamsUp[j] = Params[j] + dkappa_up;
                    ParamsDn[j] = Params[j] - dkappa_up;
                    dParams = dkappa_up;
                }
                else
                {
                    ParamsUp[j] = Params[j] + dhwvol_up;
                    ParamsDn[j] = max(0.0001, Params[j] - dhwvol_up);
                    dParams = dhwvol_up;
                }

                kappa_up = ParamsUp[0];
                kappa_dn = ParamsDn[0];
                HWVol_up = ParamsUp + 1;
                HWVol_dn = ParamsDn + 1;
            }
            else
            {
                ParamsUp[j] = Params[j] + dhwvol_up;
                ParamsDn[j] = max(0.0001, Params[j] - dhwvol_up);
                dParams = dhwvol_up;
                kappa_up = FixedKappa;
                kappa_dn = FixedKappa;
                HWVol_up = ParamsUp;
                HWVol_dn = ParamsDn;
            }

            if (SwaptionFlag[i] == 1) Pup = HW_Swaption(1.0, kappa_up, HWVolTerm, HWVol_up, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
            else Pup = HW_Cap(1.0, kappa_up, HWVolTerm, HWVol_up, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], dates[i], nDates[i], PT[i]);

            ErrorUp = BSSwaptionPrice[i] - Pup;

            if (SwaptionFlag[i] == 1) Pdn = HW_Swaption(1.0, kappa_dn, HWVolTerm, HWVol_dn, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], OptMaturityDates[i], dates[i], nDates[i], PT[i]);
            else Pdn = HW_Cap(1.0, kappa_dn, HWVolTerm, HWVol_dn, NHWVol, ZeroTerm, ZeroRate, NZero, StrikePrice[i], dates[i], nDates[i], PT[i]);

            ErrorDn = BSSwaptionPrice[i] - Pdn;

            TempJacovMatrix[i][j] = (ErrorUp - ErrorDn) / (2.0 * dParams);
        }
    }
}

void Gradient_Desent_HWCapHWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    long* SwaptionFlag,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double FreqMonth,
    double FreqMonthCap,
    double* ParamsUp,
    double* ParamsDn,
    double** TempJacovMatrix,
    long FixedKappaFlag,
    double FixedKappa
)
{

    long i;
    long n;
    long Shape0 = NResidual;
    long Shape1 = NParams;
    long BreakFlag = 0;


    double StopCondition = 0.00001;

    double absErrorSum = 100000.0;
    double PrevAbsErrorSum = 0.0;
    double ParamSum = 10000.0;
    double lambda[1] = { 0.1 };
    double* NextParams = make_array(NParams);
    double** JT_Res = make_array(NParams, 1);
    for (n = 0; n < 100; n++)
    {
        make_Jacov_HWCapHWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
            NHWVol, HWVolTerm, NResidual, SwaptionFlag, BSSwaptionPrice, StrikePrice,
            TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
            nDates, dates, PT, ParamsUp, ParamsDn, TempJacovMatrix, FixedKappaFlag, FixedKappa);

        lambda[0] *= 0.995;
        make_Residual_HWCapHWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
            NHWVol, HWVolTerm, NResidual, SwaptionFlag, BSSwaptionPrice, StrikePrice,
            TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
            nDates, dates, PT, absErrorSum, FixedKappaFlag, FixedKappa);

        Gradient_Desent(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_Res);
        for (i = 0; i < NParams; i++) Params[i] = NextParams[i];

        if (ParamSum < StopCondition && n > 10)
        {
            break;
        }

        PrevAbsErrorSum = absErrorSum;
    }

    for (i = 0; i < NParams; i++) free(JT_Res[i]);
    free(JT_Res);
    free(NextParams);
}

void Levenberg_Marquardt_HWSwaption(
    long NParams,
    double* Params,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,
    long NHWVol,
    double* HWVolTerm,
    long NResidual,
    double* BSSwaptionPrice,
    double* StrikePrice,
    double* TempHWSwaptionPrice,
    double* ResidualArray,
    double* TermSwapNew,
    double* TermOptNew,
    long* OptMaturityDates,
    long* nDates,
    long** dates,
    double** PT,
    double* ParamsUp,
    double* ParamsDn,
    double** TempJacovMatrix,
    long FixedKappaFlag,
    double FixedKappa
)
{

    long i;
    long n;
    long Shape0 = NResidual;
    long Shape1 = NParams;
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
    if (Levenberg == 1)
    {
        for (n = 0; n < 50; n++)
        {

            make_Jacov_HWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
                NHWVol, HWVolTerm, NResidual, BSSwaptionPrice, StrikePrice,
                TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
                nDates, dates, PT, ParamsUp, ParamsDn, TempJacovMatrix, FixedKappaFlag, FixedKappa);

            make_Residual_HWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
                NHWVol, HWVolTerm, NResidual, BSSwaptionPrice, StrikePrice,
                TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
                nDates, dates, PT, absErrorSum, FixedKappaFlag, FixedKappa);

            if (n >= 1) NextLambda(absErrorSum, PrevAbsErrorSum, lambda, BreakFlag);

            Levenberg_Marquardt(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_J, Inverse_JT_J, JT_Res, ResultMatrix);
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
            make_Jacov_HWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
                NHWVol, HWVolTerm, NResidual, BSSwaptionPrice, StrikePrice,
                TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
                nDates, dates, PT, ParamsUp, ParamsDn, TempJacovMatrix, FixedKappaFlag, FixedKappa);


            make_Residual_HWSwaption(NParams, Params, NZero, ZeroTerm, ZeroRate,
                NHWVol, HWVolTerm, NResidual, BSSwaptionPrice, StrikePrice,
                TempHWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
                nDates, dates, PT, absErrorSum, FixedKappaFlag, FixedKappa);

            if (n >= 1) NextLambda(absErrorSum, PrevAbsErrorSum, lambda, BreakFlag);

            Levenberg_Marquardt(NParams, NResidual, NextParams, Params, lambda, TempJacovMatrix, ResidualArray, ParamSum, JT_J, Inverse_JT_J, JT_Res, ResultMatrix);
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



long HW_GradDecent_Calibration_SwaptionCap(
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,

    long NOpt,
    double* TermOpt,
    long NSwap,
    double* TermSwap,
    double* SwaptionVol,
    double* Strike,

    long NCap,
    double* TermCapVol,
    double* CapVol,
    double FreqMonth,
    double FreqMonthCap,

    long NHW,
    double* HWTerm,
    double* HWVol,
    double& kappa,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    long j;

    long k;
    long n;
    n = 0;
    for (i = 0; i < NHW; i++)
    {
        if (TermCapVol[NCap - 1] < HWTerm[i] || TermOpt[NOpt - 1] < HWTerm[i]) n += 1;
    }
    NHW = NHW - n;

    long* SwaptionFlag = (long*)malloc(sizeof(long) * (NSwap * NOpt + NCap));
    for (i = 0; i < (NSwap * NOpt); i++) SwaptionFlag[i] = 1;
    for (i = NSwap * NOpt; i < (NSwap * NOpt + NCap); i++) SwaptionFlag[i] = 0;

    double* TermSwapNew = (double*)malloc(sizeof(double) * (NSwap * NOpt + NCap));
    for (i = 0; i < NSwap * NOpt; i++) TermSwapNew[i] = TermSwap[i % NSwap];
    for (i = NSwap * NOpt; i < NSwap * NOpt + NCap; i++) TermSwapNew[i] = TermCapVol[i - NSwap * NOpt];

    double* TermOptNew = (double*)malloc(sizeof(double) * (NSwap * NOpt + NCap));
    for (i = 0; i < NSwap * NOpt; i++) TermOptNew[i] = TermOpt[i / NSwap];
    for (i = NSwap * NOpt; i < NSwap * NOpt + NCap; i++) TermOptNew[i] = 0.0;

    double* BSSwaptionPrice = make_array(NSwap * NOpt + NCap);
    double* HWSwaptionPrice = make_array(NSwap * NOpt + NCap);

    for (k = 0; k < NSwap * NOpt; k++)
    {
        Strike[k] = FSR(ZeroTerm, ZeroRate, NZero, TermOptNew[k], TermSwapNew[k], FreqMonth);
        BSSwaptionPrice[k] = BS_Swaption(1.0, SwaptionVol[k], Strike[k], ZeroTerm, ZeroRate, NZero, TermOptNew[k], TermSwapNew[k], FreqMonth);
    }
    for (k = NSwap * NOpt; k < NSwap * NOpt + NCap; k++)
    {
        Strike[k] = FSR(ZeroTerm, ZeroRate, NZero, TermOptNew[k], TermSwapNew[k], FreqMonthCap);
        BSSwaptionPrice[k] = BS_Cap(1.0, CapVol[k - NSwap * NOpt], ZeroTerm, ZeroRate, NZero, Strike[k], TermOptNew[k], TermSwapNew[k], FreqMonthCap);
    }

    long* OptMaturityDates = (long*)malloc(sizeof(long) * (NSwap * NOpt + NCap));
    for (i = 0; i < (NSwap * NOpt + NCap); i++) OptMaturityDates[i] = (long)(TermOptNew[i] * 365.0);

    long* nDates = (long*)malloc(sizeof(long) * (NSwap * NOpt + NCap));
    for (i = 0; i < (NSwap * NOpt + NCap); i++)
    {
        if (i < NSwap * NOpt) nDates[i] = Number_Of_Payment(TermOptNew[i], TermOptNew[i] + TermSwapNew[i], FreqMonth);
        else nDates[i] = Number_Of_Payment(TermOptNew[i], TermOptNew[i] + TermSwapNew[i], FreqMonthCap);
    }

    long** dates = (long**)malloc(sizeof(long*) * (NSwap * NOpt + NCap));
    for (i = 0; i < (NSwap * NOpt + NCap); i++)
    {
        dates[i] = (long*)malloc(sizeof(long) * nDates[i]);
        if (i < NSwap * NOpt) PaymentDatesMapping(TermOptNew[i] + TermSwapNew[i], FreqMonth, dates[i], nDates[i]);
        else PaymentDatesMapping(TermOptNew[i] + TermSwapNew[i], FreqMonthCap, dates[i], nDates[i]);
    }

    double** PT = (double**)malloc(sizeof(double*) * (NSwap * NOpt + NCap));
    for (i = 0; i < (NSwap * NOpt + NCap); i++)
    {
        PT[i] = (double*)malloc(sizeof(double) * (nDates[i] + 1));
        PT[i][0] = Calc_Discount_Factor(ZeroTerm, ZeroRate, NZero, TermOptNew[i]);
        for (n = 0; n < nDates[i]; n++)  PT[i][n + 1] = Calc_Discount_Factor(ZeroTerm, ZeroRate, NZero, (double)dates[i][n] / 365.0);
    }

    long NResidual = NOpt * NSwap + NCap;
    double* ResidualArray = make_array(NResidual);

    long nparams;
    if (FixedKappaFlag != 1) nparams = 1 + NHW;
    else nparams = NHW;

    double* params = make_array(nparams);
    double* paramsup = make_array(nparams);
    double* paramsdn = make_array(nparams);
    double** tempjacov = make_array(NResidual, nparams);

    if (FixedKappaFlag != 1)
    {
        params[0] = kappa;
        for (i = 0; i < NHW; i++) params[i + 1] = HWVol[i];
    }
    else
    {
        for (i = 0; i < NHW; i++) params[i] = HWVol[i];
    }

    Gradient_Desent_HWCapHWSwaption(nparams, params, NZero, ZeroTerm, ZeroRate,
        NHW, HWTerm, NResidual, SwaptionFlag, BSSwaptionPrice, Strike,
        HWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
        nDates, dates, PT, FreqMonth, FreqMonthCap, paramsup, paramsdn,
        tempjacov, FixedKappaFlag, FixedKappa);

    if (FixedKappaFlag != 1)
    {
        kappa = params[0];
        for (i = 0; i < NHW; i++) HWVol[i] = params[i+1];
    }
    else
    {
        kappa = FixedKappa;
        for (i = 0; i < NHW; i++) HWVol[i] = params[i];
    }

    free(TermSwapNew);
    free(TermOptNew);
    free(BSSwaptionPrice);
    free(HWSwaptionPrice);
    free(OptMaturityDates);
    free(nDates);
    for (i = 0; i < NOpt * NSwap + NCap; i++) free(dates[i]);
    free(dates);
    for (i = 0; i < NOpt * NSwap + NCap; i++) free(PT[i]);
    free(PT);
    free(params);
    free(paramsup);
    free(paramsdn);
    for (i = 0; i < NResidual; i++) free(tempjacov[i]);
    free(tempjacov);
    free(ResidualArray);
    free(SwaptionFlag);

    return 1;
}

long HW_LevMarq_Calibration_Swaption(
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,

    long NOpt,
    double* TermOpt,
    long NSwap,
    double* TermSwap,
    double* SwaptionVol,
    double* Strike,
    double FreqMonth,

    long NHW,
    double* HWTerm,
    double* HWVol,
    double& kappa,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    long j;

    long k;
    long n;

    double* TermSwapNew = (double*)malloc(sizeof(double) * NSwap * NOpt);
    for (i = 0; i < NSwap * NOpt; i++) TermSwapNew[i] = TermSwap[i % NSwap];
    double* TermOptNew = (double*)malloc(sizeof(double) * NSwap * NOpt);
    for (i = 0; i < NSwap * NOpt; i++)
        TermOptNew[i] = TermOpt[i / NSwap];

    double* BSSwaptionPrice = make_array(NSwap * NOpt);
    double* HWSwaptionPrice = make_array(NSwap * NOpt);

    for (k = 0; k < NSwap * NOpt; k++)
    {
        Strike[k] = FSR(ZeroTerm, ZeroRate, NZero, TermOptNew[k], TermSwapNew[k], FreqMonth);
        BSSwaptionPrice[k] = BS_Swaption(1.0, SwaptionVol[k], Strike[k], ZeroTerm, ZeroRate, NZero, TermOptNew[k], TermSwapNew[k], FreqMonth);
    }

    long* OptMaturityDates = (long*)malloc(sizeof(long) * NOpt * NSwap);
    for (i = 0; i < NOpt * NSwap; i++) OptMaturityDates[i] = (long)(TermOptNew[i] * 365.0);

    long* nDates = (long*)malloc(sizeof(long) * NOpt * NSwap);
    for (i = 0; i < NOpt * NSwap; i++)
        nDates[i] = Number_Of_Payment(TermOptNew[i], TermOptNew[i] + TermSwapNew[i], FreqMonth);

    long** dates = (long**)malloc(sizeof(long*) * NOpt * NSwap);
    for (i = 0; i < NOpt * NSwap; i++)
    {
        dates[i] = (long*)malloc(sizeof(long) * nDates[i]);
        PaymentDatesMapping(TermOptNew[i] + TermSwapNew[i], FreqMonth, dates[i], nDates[i]);
    }

    double** PT = (double**)malloc(sizeof(double*) * NOpt * NSwap);
    for (i = 0; i < NOpt * NSwap; i++)
    {
        PT[i] = (double*)malloc(sizeof(double) * (nDates[i] + 1));
        PT[i][0] = Calc_Discount_Factor(ZeroTerm, ZeroRate, NZero, TermOptNew[i]);
        for (n = 0; n < nDates[i]; n++)  PT[i][n + 1] = Calc_Discount_Factor(ZeroTerm, ZeroRate, NZero, (double)dates[i][n] / 365.0);
    }

    long NResidual = NOpt * NSwap;
    double* ResidualArray = make_array(NResidual);

    long nparams;
    if (FixedKappaFlag != 1) nparams = 1 + NHW;
    else nparams = NHW;

    double* params = make_array(nparams);
    double* paramsup = make_array(nparams);
    double* paramsdn = make_array(nparams);
    double** tempjacov = make_array(NResidual, nparams);

    if (FixedKappaFlag != 1)
    {
        params[0] = kappa;
        for (i = 0; i < NHW; i++) params[i + 1] = HWVol[i];
    }
    else
    {
        for (i = 0; i < NHW; i++) params[i] = HWVol[i];
    }

    Levenberg_Marquardt_HWSwaption(nparams, params, NZero, ZeroTerm, ZeroRate,
        NHW, HWTerm, NResidual, BSSwaptionPrice, Strike,
        HWSwaptionPrice, ResidualArray, TermSwapNew, TermOptNew, OptMaturityDates,
        nDates, dates, PT, paramsup, paramsdn,
        tempjacov, FixedKappaFlag, FixedKappa);

    if (FixedKappaFlag != 1)
    {
        kappa = params[0];
        for (i = 0; i < NHW; i++) HWVol[i] = params[i + 1];
    }
    else
    {
        kappa = FixedKappa;
        for (i = 0; i < NHW; i++) HWVol[i] = params[i];
    }

    free(TermSwapNew);
    free(TermOptNew);
    free(BSSwaptionPrice);
    free(HWSwaptionPrice);
    free(OptMaturityDates);
    free(nDates);
    for (i = 0; i < NOpt * NSwap; i++) free(dates[i]);
    free(dates);
    for (i = 0; i < NOpt * NSwap; i++) free(PT[i]);
    free(PT);
    free(params);
    free(paramsup);
    free(paramsdn);
    for (i = 0; i < NResidual; i++) free(tempjacov[i]);
    free(tempjacov);
    free(ResidualArray);

    return 1;
}

DLLEXPORT(long) HWCapHWSwaptionCalibration(
    long CaliFlag,
    long NZero,
    double* ZeroTerm,
    double* ZeroRate,

    long NOpt,
    double* TermOpt,
    long NSwap,
    double* TermSwap,
    double* SwaptionVol,
    double FreqMonthSwaption,

    long NCap,
    double* TermCap,
    double* CapVol,
    double FreqMonthCap,

    long NHW_initial,
    double* HWTerm_initial,
    double* HWVol_initial,
    double* initial_kappa,
    long FixedKappaFlag,
    double FixedKappa
)
{
    long i;
    long ResultCode = 0;
    //_CrtSetBreakAlloc(161);

    long NHW = 4;
    double* HWTerm;
    double* HWVol;

    double* Strike = (double*)malloc(sizeof(double) * (NSwap * NOpt + NCap));
    double kappa = initial_kappa[0];
    if (NCap == 0 && NSwap * NOpt > 0 || CaliFlag == 1 )
    {
        HWTerm = TermOpt;
        NHW = min(NOpt, NHW_initial);
        HWVol = (double*)malloc(sizeof(double) * NHW);
        for (i = 0; i < NHW; i++) HWVol[i] = max(0.001, HWVol_initial[i]);
        ResultCode = HW_LevMarq_Calibration_Swaption(NZero, ZeroTerm, ZeroRate,
            NOpt, TermOpt, NSwap, TermSwap, SwaptionVol, Strike, FreqMonthSwaption,
            NHW, HWTerm, HWVol, kappa, FixedKappaFlag, FixedKappa);
        for (i = 0; i < NHW; i++)
        {
            HWTerm_initial[i] = HWTerm[i];
            HWVol_initial[i] = HWVol[i];
            *initial_kappa = kappa;
        }
        free(HWVol);
    }
    else if (NCap > 0)
    {
        if (CaliFlag == 2) NOpt = 0;
        if (NOpt > NCap)
        {
            HWTerm = TermOpt;
            NHW = min(NHW_initial, NOpt);
        }
        else {
            HWTerm = TermCap;
            NHW = min(NHW_initial, NCap);
        }
        HWVol = (double*)malloc(sizeof(double) * NHW);
        for (i = 0; i < NHW; i++) HWVol[i] = max(0.001, HWVol_initial[i]);

        ResultCode = HW_GradDecent_Calibration_SwaptionCap(NZero, ZeroTerm, ZeroRate,
            NOpt, TermOpt, NSwap, TermSwap, SwaptionVol, Strike,
            NCap, TermCap, CapVol, FreqMonthSwaption, FreqMonthCap,
            NHW, HWTerm, HWVol, kappa, FixedKappaFlag, FixedKappa);
        for (i = 0; i < NHW; i++)
        {
            HWTerm_initial[i] = HWTerm[i];
            HWVol_initial[i] = HWVol[i];
            *initial_kappa = kappa;
        }
        free(HWVol);
    }
    else
    {
        return -1;
    }

    free(Strike);
    _CrtDumpMemoryLeaks();
}