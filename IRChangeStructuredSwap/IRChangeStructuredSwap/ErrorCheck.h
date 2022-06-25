
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

long ErrorCheck_IRStructuredSwap_Excel(
    long PriceDateExcel,                // 1    가격계산일 Excel long Type
    long NAFlag,                        // 2    Notional Amount 사용 여부 0: No ExChange Notional Amount 1: ExChange Notioanl Amount
    double Notional,                    // 3    Notional Amount
    long* PayoffStructure,              // 4    [0]: RcvLeg 페이오프조건 [1]: PayLeg 페이오프조건 {페이오프조건 0:Range(R1)&Range(R2)&Range(R3) 1:Range(R1+R2+R3)}
    long* HolidayFlagCount,             // 5    [0~2]: R1,R2,R3 Receive HolidayFlag, [3~5]: R1, R2, R3 Pay HolidayFlag, [6~8]: R1, R2, R3 Rcv HolidayCount, [9~11]: R1, R2, R3 Pay HolidayCount
    long* Holidays_Excel,               // 6    [0~sum(HolidayFlagCount[6~11])]: Holidays Rcv R1R2R3, Pay R1R2R3
    // Receive Leg Information
    long* RcvFlag,                      // 7    [0]: NReference, [1]: Fix/Flo Flag, [2]: Accrual Flag, [3]: DailySimulFlag
    double* RcvMaxLossReturn,           // 8
    long* RcvRefCurveNumber,            // 9
    long* RcvRefRateType,               // 10
    double* RcvSwapInfo,                // 11
    long* RcvDayCountFlag,              // 12
    double* RcvRefCondMultiple,         // 13
    long* RcvRefPowerSpreadFlag,        // 14
    double* RcvRefPowerSpreadMat,       // 15
    double* RcvRefRange,                // 16
    double* RcvPayoffMultiple,          // 17
    long NRcvCashFlow,                  // 18
    long* RcvCashFlowScheduleExcelDate, // 19   [0~NCF-1]: ResetStart, [NCF~2*NCF-1]: ResetEnd, [2*NCF~3*NCF-1]: 기산일, [3*NCF~4*NCF-1]: 기말일, [4*NCF~5*NCF-1]: 지급일
    double* RcvCouponFixRate,           // 20
    // Receive Leg Information
    long* PayFlag,                      // 21    [0]: NReference, [1]: Fix/Flo Flag, [2]: Accrual Flag
    double* PayMaxLossReturn,           // 22
    long* PayRefCurveNumber,            // 23
    long* PayRefRateType,               // 24
    double* PaySwapInfo,                 // 25
    long* PayDayCountFlag,              // 26
    double* PayRefCondMultiple,         // 27
    long* PayRefPowerSpreadFlag,        // 28
    double* PayRefPowerSpreadMat,       // 29
    double* PayRefRange,                // 30
    double* PayPayoffMultiple,          // 31
    long NPayCashFlow,                  // 32
    long* PayCashFlowScheduleExcelDate, // 33   [0~NCF-1]: ResetStart, [NCF~2*NCF-1]: ResetEnd, [2*NCF~3*NCF-1]: 기산일, [3*NCF~4*NCF-1]: 기말일, [4*NCF~5*NCF-1]: 지급일
    double* PayCouponFixRate,            // 34
    // 커브 및 HW정보

    long* NHWVolCount,                  // 35
    double* HWVolTermArray,             // 36
    double* HWVolArray,                 // 37
    double* HWKappaArray,               // 38
    long* NZeroRate,                    // 39
    double* ZeroTerm,                   // 40
    double* ZeroRate,                   // 41
    // 상관계수 및 히스토리 
    double* Corr,                       // 42
    long* NDayHistRcv,                  // 43
    long* RateHistDateRcv,              // 44
    double* RateHistRcv,                // 45
    long* NDayHistPay,                  // 46
    long* RateHistDatePay,              // 47
    double* RateHistPay,                // 48

    // 시뮬레이션
    long NSimul,                        // 49
    long GreekFlag,                     // 50
    char* Error
)
{
    long i;
    long j;
    char ErrorName[100];
    long* Rcv_ResetStart = RcvCashFlowScheduleExcelDate;
    long* Rcv_ResetEnd = RcvCashFlowScheduleExcelDate + NRcvCashFlow;
    long* Rcv_StartDate = RcvCashFlowScheduleExcelDate + 2 * NRcvCashFlow;
    long* Rcv_EndDate = RcvCashFlowScheduleExcelDate + 3 * NRcvCashFlow;
    long* Rcv_PayDate = RcvCashFlowScheduleExcelDate + 4 * NRcvCashFlow;

    long* Pay_ResetStart = PayCashFlowScheduleExcelDate;
    long* Pay_ResetEnd = PayCashFlowScheduleExcelDate + NPayCashFlow;
    long* Pay_StartDate = PayCashFlowScheduleExcelDate + 2 * NPayCashFlow;
    long* Pay_EndDate = PayCashFlowScheduleExcelDate + 3 * NPayCashFlow;
    long* Pay_PayDate = PayCashFlowScheduleExcelDate + 4 * NPayCashFlow;

    long Rcv_nUnderlying = RcvFlag[0];
    long Pay_nUnderlying = PayFlag[0];

    if (PriceDateExcel < 0 || PriceDateExcel >= RcvCashFlowScheduleExcelDate[5 * NRcvCashFlow - 1] || PriceDateExcel >= PayCashFlowScheduleExcelDate[5 * NPayCashFlow - 1])
    {
        strcpy_s(ErrorName, "Check PriceDate or Rcv Pay ScheduleDate");
        return SaveErrorName2(Error, ErrorName);
    }

    if (NAFlag != 0 && NAFlag != 1)
    {
        strcpy_s(ErrorName, "Check NAFlag");
        return SaveErrorName2(Error, ErrorName);
    }

    if (Notional < 0.0)
    {
        strcpy_s(ErrorName, "Check Notional");
        return SaveErrorName2(Error, ErrorName);
    }

    if ((PayoffStructure[0] != 0 && PayoffStructure[0] != 1) || (PayoffStructure[1] != 0 && PayoffStructure[1] != 1))
    {
        strcpy_s(ErrorName, "Check PayoffStructureFlag");
        return SaveErrorName2(Error, ErrorName);
    }

    for (i = 0; i < 6; i++)
    {
        if (HolidayFlagCount[i] != 0 && HolidayFlagCount[i] != 1 && HolidayFlagCount[i] != 2 && HolidayFlagCount[i] != 3)
        {
            strcpy_s(ErrorName, "Check HolidayFlag");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    for (i = 6; i < 12; i++)
    {
        if (HolidayFlagCount[i] < 0)
        {
            strcpy_s(ErrorName, "Check HolidayCount");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    if (RcvFlag[0] <= 0 || RcvFlag[0] > 3)
    {
        strcpy_s(ErrorName, "Rcv UnderlyingNum is Over Zero or Less Than Three");
        return SaveErrorName2(Error, ErrorName);
    }

    if (PayFlag[0] <= 0 || PayFlag[0] > 3)
    {
        strcpy_s(ErrorName, "Pay UnderlyingNum is Over Zero or Less Than Three");
        return SaveErrorName2(Error, ErrorName);
    }

    if (RcvFlag[1] != 0 && RcvFlag[1] != 1)
    {
        strcpy_s(ErrorName, "Rcv FixeFlotingFlag Error");
        return SaveErrorName2(Error, ErrorName);
    }

    if (PayFlag[1] != 0 && PayFlag[1] != 1)
    {
        strcpy_s(ErrorName, "Pay FixeFlotingFlag Error");
        return SaveErrorName2(Error, ErrorName);
    }

    if (RcvFlag[2] != 0 && RcvFlag[2] != 1)
    {
        strcpy_s(ErrorName, "Rcv AccrualFlag Error");
        return SaveErrorName2(Error, ErrorName);
    }

    if (PayFlag[2] != 0 && PayFlag[2] != 1)
    {
        strcpy_s(ErrorName, "Pay AccrualFlag Error");
        return SaveErrorName2(Error, ErrorName);
    }

    for (i = 0; i < Rcv_nUnderlying; i++)
    {
        if (RcvRefCurveNumber[i] < 1 || RcvRefCurveNumber[i] > 4)
        {
            strcpy_s(ErrorName, "Check Rcv Reference CurveNumber");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    for (i = 0; i < Pay_nUnderlying; i++)
    {
        if (PayRefCurveNumber[i] < 1 || PayRefCurveNumber[i] > 4)
        {
            strcpy_s(ErrorName, "Check Pay Reference CurveNumber");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    if (Rcv_ResetStart[0] > Rcv_ResetEnd[0] || Rcv_PayDate[0] < Rcv_ResetEnd[0])
    {
        strcpy_s(ErrorName, "Sort Rcv Schedule");
        return SaveErrorName2(Error, ErrorName);
    }

    for (i = 1; i < NRcvCashFlow; i++)
    {
        if (Rcv_ResetStart[i] > Rcv_ResetEnd[i] || Rcv_ResetStart[i] < Rcv_ResetStart[i - 1] || Rcv_ResetEnd[i] < Rcv_ResetEnd[i-1] || Rcv_PayDate[i] < Rcv_PayDate[i-1] || Rcv_PayDate[i] < Rcv_ResetEnd[i])
        {
            strcpy_s(ErrorName, "Sort Rcv Schedule");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    if (Pay_ResetStart[0] > Pay_ResetEnd[0] || Pay_PayDate[0] < Pay_ResetEnd[0])
    {
        strcpy_s(ErrorName, "Sort Pay Schedule");
        return SaveErrorName2(Error, ErrorName);
    }

    for (i = 1; i < NPayCashFlow; i++)
    {
        if (Pay_ResetStart[i] > Pay_ResetEnd[i] || Pay_ResetStart[i] < Pay_ResetStart[i - 1] || Pay_ResetEnd[i] < Pay_ResetEnd[i - 1] || Pay_PayDate[i] < Pay_PayDate[i - 1] || Pay_PayDate[i] < Pay_ResetEnd[i])
        {
            strcpy_s(ErrorName, "Sort Pay Schedule");
            return SaveErrorName2(Error, ErrorName);
        }
    }

    return 1;
}