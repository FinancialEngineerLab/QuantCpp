#include <stdio.h>
#include <math.h>

#ifndef DLLEXPORT(A)
#ifdef WIN32
#include <crtdbg.h>
#define DLLEXPORT(A) extern "C" __declspec(dllexport) A _stdcall
#elif _WIN64
#include <crtdbg.h>
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

DLLEXPORT(long) CDateToExcelDate(long Cdate);
DLLEXPORT(long) ExcelDateToCDate(long ExlDate);
DLLEXPORT(long) DayCountAtoB(long Day1, long Day2); 
DLLEXPORT(long) DayPlus(long Cdate, long NDays);

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

// 윤년Check
long LeapCheck(long Year)
{
    if ((Year % 4 == 0 && Year % 100 != 0) || Year % 400 == 0)
        return 1;
    else
        return 0;
}

// 20210627을 엑셀타입 2021-06-27(44374)로 바꿔줌
DLLEXPORT(long) CDateToExcelDate(long Cdate)
{
    long Days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    long CummulativeDays[13] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
    long Leap = 0;

    long Year = (long)Cdate / 10000;
    long Month = (long)(Cdate - Year * 10000) / 100;
    long Day = (long)(Cdate - Year * 10000 - Month * 100);

    Leap = LeapCheck(Year);
    if (Year == 1900) Leap = 1;

    long N4 = (long)Year / 4;
    long N100 = (long)Year / 100;
    long N400 = (long)Year / 400;

    long N4_To_1900 = 1900 / 4;
    long N100_To_1900 = 1900 / 100;
    long N400_To_1900 = 1900 / 400;

    long Result;

    if (Leap == 1)
    {
        if (Month > 2)
            Result = (Year - 1900) * 365 + (N4 - N4_To_1900) - (N100 - N100_To_1900) + (N400 - N400_To_1900) + CummulativeDays[Month - 1] + 1 + Day;
        else
            Result = (Year - 1900) * 365 + (N4 - N4_To_1900) - (N100 - N100_To_1900) + (N400 - N400_To_1900) + CummulativeDays[Month - 1] + Day;
    }

    else
    {
        Result = (Year - 1900) * 365 + (N4 - N4_To_1900) - (N100 - N100_To_1900) + (N400 - N400_To_1900) + CummulativeDays[Month - 1] + 1 + Day;
    }

    return Result;
}


// 엑셀타입 2021-06-27(44374)을 20210627으로 바꿔줌
DLLEXPORT(long) ExcelDateToCDate(long ExlDate)
{
    long i, m;

    long CummulativeDays[13] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
    long CummulativeDays_Leap[13] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };

    // xx00년도부터 xx99년도까지 Day수(xx는 4의 배수)
    long LeapDaysFor4Year[101] = { 0, 366, 731, 1096, 1461, 1827, 2192, 2557, 2922, 3288,
                                   3653, 4018, 4383, 4749, 5114, 5479, 5844, 6210, 6575, 6940,
                                   7305, 7671, 8036, 8401, 8766, 9132, 9497, 9862, 10227, 10593,
                                   10958, 11323, 11688, 12054, 12419, 12784, 13149, 13515, 13880, 14245,
                                   14610, 14976, 15341, 15706, 16071, 16437, 16802, 17167, 17532, 17898,
                                   18263, 18628, 18993, 19359, 19724, 20089, 20454, 20820, 21185, 21550,
                                   21915, 22281, 22646, 23011, 23376, 23742, 24107, 24472, 24837, 25203,
                                   25568, 25933, 26298, 26664, 27029, 27394, 27759, 28125, 28490, 28855,
                                   29220, 29586, 29951, 30316, 30681, 31047, 31412, 31777, 32142, 32508,
                                   32873, 33238, 33603, 33969, 34334, 34699, 35064, 35430, 35795, 36160, 36525 };

    // xx00년도부터 xx99년도까지 Day수(xx는 4의 배수가 아님)
    long LeapDaysFor4YearNot400[101] = { 0,365, 730, 1095, 1460, 1826, 2191, 2556, 2921, 3287,
                                         3652, 4017, 4382, 4748, 5113, 5478, 5843, 6209, 6574, 6939,
                                         7304, 7670, 8035, 8400, 8765, 9131, 9496, 9861, 10226, 10592,
                                         10957, 11322, 11687, 12053, 12418, 12783, 13148, 13514, 13879, 14244,
                                         14609, 14975, 15340, 15705, 16070, 16436, 16801, 17166, 17531, 17897,
                                         18262, 18627, 18992, 19358, 19723, 20088, 20453, 20819, 21184, 21549,
                                         21914, 22280, 22645, 23010, 23375, 23741, 24106, 24471, 24836, 25202,
                                         25567, 25932, 26297, 26663, 27028, 27393, 27758, 28124, 28489, 28854,
                                         29219, 29585, 29950, 30315, 30680, 31046, 31411, 31776, 32141, 32507,
                                         32872, 33237, 33602, 33968, 34333, 34698, 35063, 35429, 35794, 36159, 36524 };

    // 1900~1999.12.31까지의 Day 개수 36525(엑셀기준)
    long ExcelDate;
    long Year, Month, Day;

    if (ExlDate - 36525 > 0)
    {
        ExcelDate = ExlDate - 36525;
        Year = 2000;
    }

    else
    {
        ExcelDate = ExlDate;
        Year = 1900;
    }

    //400년동안 윤년 개수
    //long nLeap_for_400 = 400 / 4 - 400 / 100 + 1;
    //long nDay_for_400 = 365 * 400 + nLeap_for_400;
    long nDay_for_400 = 146097;

    for (i = 0; i < 1000; i++)
    {
        if (ExcelDate - nDay_for_400 <= 0) break;
        else
        {
            ExcelDate -= nDay_for_400;
            Year += 400;
        }
    }
    //400년의 배수만큼 날짜 차감 완료
    
    //100년동안 윤년 개수
    long nDay_for_100_400x = 36525;       //400의 배수 년도에는36525
    long nDay_for_100 = 36524;            //나머지 3년동안은 36524

    for (i = 0; i < 5; i++)
    {
        if (i % 4 == 0)
        {
            if (ExcelDate - nDay_for_100_400x <= 0) break;
            else
            {
                ExcelDate -= nDay_for_100_400x;
                Year += 100;
            }
        }
        else
        {
            if (ExcelDate - nDay_for_100 <= 0) break;
            else
            {
                ExcelDate -= nDay_for_100;
                Year += 100;
            }
        }
    }
    // 100년의 배수만큼 날짜 차감 완료

    long Flag_400;
    long Leap;

    // 이쯤에서 남은 날짜가 0인지 체크
    if (ExcelDate == 0)
    {
        Year -= 1;
        return Year * 10000 + 1231;
    }
    else
    {
        //위에서 년도 100의자리 체크할 때 사용했던 i를 사용
        //년도 백의자리가 4의 배수인지 아닌지 체크 (16xx년 일 경우 해당 100년동안은 36525일)
        Flag_400 = 0;
        if (i % 4 == 0)
        {
            Flag_400 = 1;
        }

        if (Flag_400 == 1)
        {
            for (i = 0; i < 101; i++)
            {
                if (ExcelDate - LeapDaysFor4Year[i + 1] <= 0)
                {
                    ExcelDate = ExcelDate - LeapDaysFor4Year[i];
                    Year += i;
                    break;
                }
            }
        }
        else
        {
            for (i = 0; i < 101; i++)
            {
                if (ExcelDate - LeapDaysFor4YearNot400[i + 1] <= 0)
                {
                    ExcelDate = ExcelDate - LeapDaysFor4YearNot400[i];
                    Year += i;
                    break;
                }
            }
        }

        // 이쯤에서 남은 날짜가 0인지 한번 더 체크
        if (ExcelDate == 0)
        {
            Year -= 1;
            return Year * 10000 + 1231;
        }

        Leap = 0;
        Leap = LeapCheck(Year);
        if (Year == 1900) Leap = 1;

        if (Leap == 1)
        {
            for (m = 0; m < 12; m++)
            {
                if (CummulativeDays_Leap[m + 1] >= ExcelDate && CummulativeDays_Leap[m] < ExcelDate) break;
            }
            if (m < 12)
                Month = m + 1;
            else
                Month = 12;
            Day = ExcelDate - CummulativeDays_Leap[m];
        }
        else
        {
            for (m = 0; m < 12; m++)
            {
                if (CummulativeDays[m + 1] >= ExcelDate && CummulativeDays[m] < ExcelDate) break;
            }

            if (m < 12)
                Month = m + 1;
            else
                Month = 12;

            Day = ExcelDate - CummulativeDays[m];
        }
        return   Year * 10000 + Month * 100 + Day;
    }
}

DLLEXPORT(long) DayCountAtoB(long Day1, long Day2)
{
    return CDateToExcelDate(Day2) - CDateToExcelDate(Day1);
}

DLLEXPORT(long) DayPlus(long Cdate, long NDays)
{
    long i;
    long Days[13] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ,31 };
    long Days_Leap[13] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ,31 };
    long CummulativeDays[14] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ,396 };
    long CummulativeDays_Leap[14] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366,397 };
    long Leap = 0;

    long Year = (long)Cdate / 10000;
    long Month = (long)(Cdate - Year * 10000) / 100;
    long Day = (long)(Cdate - Year * 10000 - Month * 100);
    long EndFlag = 1;
    long NewMonth;
    long NewDay;


    //0보다 작을 때는 귀찮으니 기존 Function 사용하자
    if (NDays < 0) 
        return ExcelDateToCDate(CDateToExcelDate(Cdate) + NDays);
    else if (NDays == 0)
        return Cdate;

    Leap = LeapCheck(Year);
    if (Year == 1900) Leap = 1;

    long LeftDaysofMonth;
    if (Leap == 1)
        LeftDaysofMonth = Days_Leap[Month - 1] - Day;
    else
        LeftDaysofMonth = Days[Month - 1] - Day;

    long NDatesOfYearToday;
    long NDayCountFromJan01;
    if (NDays <= LeftDaysofMonth)
    {
        return Year * 10000 + Month * 100 + (Day + NDays);
    }
    else
    {
        //당해 중 지금까지 날짜 개수
        if (Leap == 1) // 만약 해당 날짜가 윤년이라면
        {
            NDatesOfYearToday = CummulativeDays_Leap[Month - 1] + Day;
            NDayCountFromJan01 = NDatesOfYearToday + NDays;
            if (NDayCountFromJan01 <= 366)
            {
                for (i = 0; i < 13; i++)
                {
                    if (NDayCountFromJan01 > CummulativeDays_Leap[i] && NDayCountFromJan01 < CummulativeDays_Leap[i + 1])
                    {
                        NewMonth = i + 1;
                        NewDay = NDayCountFromJan01 - CummulativeDays_Leap[i];
                        break;
                    }
                    else if (NDayCountFromJan01 == CummulativeDays_Leap[i])
                    {
                        NewMonth = i;
                        NewDay = Days_Leap[NewMonth - 1];
                        break;
                    }
                }
                return Year * 10000 + NewMonth * 100 + NewDay;
            }
            else if (NDayCountFromJan01 == 366)
            {
                return Year * 10000 + 1231;
            }
            else
            {
                Year += 1;
                DayPlus(Year * 10000 + 1 * 100 + 1, NDayCountFromJan01 - 366 - 1);
            }
        }
        else
        {
            NDatesOfYearToday = CummulativeDays[Month - 1] + Day;
            NDayCountFromJan01 = NDatesOfYearToday + NDays;
            if (NDayCountFromJan01 <= 365)
            {
                for (i = 0; i < 13; i++)
                {
                    if (NDayCountFromJan01 > CummulativeDays[i] && NDayCountFromJan01 < CummulativeDays[i + 1])
                    {
                        NewMonth = i + 1;
                        NewDay = NDayCountFromJan01 - CummulativeDays[i];
                        break;
                    }
                    else if (NDayCountFromJan01 == CummulativeDays[i])
                    {
                        NewMonth = i;
                        NewDay = Days[NewMonth - 1];
                        break;
                    }
                }
                return Year * 10000 + NewMonth * 100 + NewDay;
            }
            else if (NDayCountFromJan01 == 365)
            {
                return Year * 10000 + 1231;
            }
            else
            {
                Year += 1;
                DayPlus(Year * 10000 + 1 * 100 + 1, NDayCountFromJan01 - 365 - 1);
            }
        }
    }
}

DLLEXPORT(long) EDate_Cpp(long Cdate, long NMonths)
{
    long Days[13] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ,31 };
    long Days_Leap[13] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ,31 };
    long CummulativeDays[14] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ,396 };
    long CummulativeDays_Leap[14] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366,397 };
    long Leap = 0;

    long Year = (long)Cdate / 10000;
    long Month = (long)(Cdate - Year * 10000) / 100;
    long Day = (long)(Cdate - Year * 10000 - Month * 100);

    long TargetMonth = Month + NMonths;
    long PlusYear;
    long ResultMonth;
    long ResultDay;

    if (TargetMonth > 12)
    {
        PlusYear = (long)(TargetMonth / 12);
        ResultMonth = TargetMonth % 12;
        Leap = LeapCheck(Year + PlusYear);
        if (Leap == 1)
        {
            ResultDay = min(Day, Days_Leap[ResultMonth - 1]);

        }
        else
        {
            ResultDay = min(Day, Days[ResultMonth - 1]);
        }
        if (TargetMonth % 12 != 0)
            return (Year + PlusYear) * 10000 + ResultMonth * 100 + ResultDay;
        else
            return (Year + PlusYear - 1) * 10000 + 12 * 100 + ResultDay;
    }
    else if (TargetMonth > 0)
    {
        ResultMonth = TargetMonth;
        Leap = LeapCheck(Year);
        if (Leap == 1)
        {
            ResultDay = min(Day, Days_Leap[ResultMonth - 1]);

        }
        else
        {
            ResultDay = min(Day, Days[ResultMonth - 1]);
        }
        return Year * 10000 + ResultMonth * 100 + ResultDay;
    }
    else
    {
        PlusYear = (long)(TargetMonth / 12 - 1);
        ResultMonth = (12 + TargetMonth % 12);
        Leap = LeapCheck(Year + PlusYear);
        if (Leap == 1)
        {
            ResultDay = min(Day, Days_Leap[ResultMonth - 1]);

        }
        else
        {
            ResultDay = min(Day, Days[ResultMonth - 1]);
        }
        return (Year + PlusYear) * 10000 + ResultMonth * 100 + ResultDay;
    }
}
