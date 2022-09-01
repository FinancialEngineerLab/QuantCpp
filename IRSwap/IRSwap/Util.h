#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define PI		3.141592653589793238462643383279
#define UTILITY 1

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


long GetMinor(double** src, double** dest, long row, long col, long order);
double CalcDeterminant(double** mat, long order, double*** MinorMatrixList);
long idum = -1;

/////////////////////////////////////////////////
// Numerical Recipy Random Number Generator
// Page 305
/////////////////////////////////////////////////

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
static long iy = 0;
static long iv[NTAB];


double ran1(long* idum, long& iy, long* iv)
{
	long j;
	long k;
	//static long iy = 0;
	//static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0) *idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX) return RNMX;
	else return temp;
}

double randnorm()
{
	static long iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	if (idum < 0) iset = 0;
	if (iset == 0) {
		do {
			v1 = 2.0 * ran1(&idum, iy, iv) - 1.0;
			v2 = 2.0 * ran1(&idum, iy, iv) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

// mu = 0, sigma = 1 RandomNumber Generate (시드값 주면 리턴하지않고 시드셋팅)
void randnorm(long seed)
{
	idum = seed;
	iy = 0;
}

// Make 1d Array
double* make_array(long n)
{
	double* Array = (double*)calloc(n, sizeof(double));
	return Array;
}

// Make 1d Array Random Number
double* make_array_randn(long n)
{
	long i;
	double* Array = (double*)calloc(n, sizeof(double));
	for (i = 0; i < n; i++)
	{
		Array[i] = randnorm();
	}
	return Array;
}

// Make 2d Array
double** make_array(long n, long m)
{
	long i;
	double** Array_2d = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++)
	{
		Array_2d[i] = (double*)calloc(m, sizeof(double));
	}
	return Array_2d;
}

// Make 2d Array Random Number
double** make_array_randn(long n, long m)
{
	long i, j;
	double** Array_2d = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++)
	{
		Array_2d[i] = (double*)calloc(m, sizeof(double));
		for (j = 0; j < m; j++)
			Array_2d[i][j] = randnorm();
	}
	return Array_2d;
}


// Print Array
void Print_Array(double* Array, long n)
{
	long i;
	printf("[");
	for (i = 0; i < n - 1; i++)
		printf("%0.5lf,  ", Array[i]);
	printf("%0.5lf]", Array[n - 1]);
}

// Print Array
void Print_Array(double** Array, long n, long m)
{
	long i;
	printf("[");
	for (i = 0; i < n - 1; i++)
	{
		Print_Array(Array[i], m);
		printf("\n");
	}
	Print_Array(Array[n - 1], m);
	printf("]");
}



//표준정규분포(Standard Normal Distribution) PDF
DLLEXPORT(double) PDF_N(double x)
{
	return exp(-x * x / 2.0) / sqrt(2.0 * PI);
}

//표준정규분포(Standard Normal Distribution) CDF
DLLEXPORT(double) CDF_N(double x)
{
	double y, Exponential, S, S2, value;

	y = fabs(x);
	if (y > 37) value = 0;
	else {
		Exponential = exp(-y * y / 2);
		if (y < 7.07106781186547) {
			S = 3.52624965998911E-02 * y + 0.700383064443688;
			S = S * y + 6.37396220353165;
			S = S * y + 33.912866078383;
			S = S * y + 112.079291497871;
			S = S * y + 221.213596169931;
			S = S * y + 220.206867912376;
			S2 = 8.83883476483184E-02 * y + 1.75566716318264;
			S2 = S2 * y + 16.064177579207;
			S2 = S2 * y + 86.7807322029461;
			S2 = S2 * y + 296.564248779674;
			S2 = S2 * y + 637.333633378831;
			S2 = S2 * y + 793.826512519948;
			S2 = S2 * y + 440.413735824752;
			value = Exponential * S / S2;
		}
		else {
			S = y + 0.65;
			S = y + 4 / S;
			S = y + 3 / S;
			S = y + 2 / S;
			S = y + 1 / S;
			value = Exponential / (S * 2.506628274631);
		}
	}

	if (x > 0) value = 1 - value;

	return value;
}

DLLEXPORT(double) INV_CDF_N(double p)
{

	double a1 = -39.69683028665376;
	double a2 = 220.9460984245205;
	double a3 = -275.9285104469687;
	double a4 = 138.3577518672690;
	double a5 = -30.66479806614716;
	double a6 = 2.506628277459239;

	double b1 = -54.47609879822406;
	double b2 = 161.5858368580409;
	double b3 = -155.6989798598866;
	double b4 = 66.80131188771972;
	double b5 = -13.28068155288572;

	double c1 = -0.007784894002430293;
	double c2 = -0.3223964580411365;
	double c3 = -2.400758277161838;
	double c4 = -2.549732539343734;
	double c5 = 4.374664141464968;
	double c6 = 2.938163982698783;

	double d1 = 0.007784695709041462;
	double d2 = 0.3224671290700398;
	double d3 = 2.445134137142996;
	double d4 = 3.754408661907416;

	//Define break-points.

	double p_low = 0.02425;
	double p_high = 1 - p_low;
	double  q, r, e, u;
	double x = 0.0;


	//Rational approximation for lower region.

	if (0 < p && p < p_low) {
		q = sqrt(-2 * log(p));
		x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	}

	//Rational approximation for central region.

	if (p_low <= p && p <= p_high) {
		q = p - 0.5;
		r = q * q;
		x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
	}

	//Rational approximation for upper region.

	if (p_high < p && p < 1) {
		q = sqrt(-2 * log(1 - p));
		x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	}

	//Pseudo-code algorithm for refinement

	if ((0 < p) && (p < 1)) {
		e = 0.5 * erfc(-x / sqrt(2)) - p;
		u = e * sqrt(2 * PI) * exp(x * x / 2);
		x = x - u / (1 + x * u / 2);
	}

	return x;
}

// Linear Interpolation (X변수, Y변수, X길이, 타겟X)
DLLEXPORT(double) Interpolate_Linear(double* x, double* fx, long nx, double targetx, long extrapolateflag)
{
	long i;
	double result = 0.0;

	if (nx == 1 || targetx == x[0])
		return fx[0];
	else if (targetx == x[nx - 1])
		return fx[nx - 1];


	if (targetx < x[0])
	{
		if (extrapolateflag == 0) return fx[0];
		else return (fx[1] - fx[0]) / (x[1] - x[0]) * (targetx - x[0]) + fx[0];
	}
	else if (targetx > x[nx - 1])
	{
		if (extrapolateflag == 0) return fx[nx - 1];
		else return (fx[nx - 1] - fx[nx - 2]) / (x[nx - 1] - x[nx - 2]) * (targetx - x[nx - 1]) + fx[nx - 1];
	}
	else
	{
		for (i = 1; i < nx; i++)
		{
			if (targetx < x[i])
			{
				result = (fx[i] - fx[i - 1]) / (x[i] - x[i - 1]) * (targetx - x[i - 1]) + fx[i - 1];
				break;
			}
		}
		return result;
	}
}

// Linear Interpolation (X변수, Y변수, X길이, 타겟X)
double Interpolate_Linear(double* x, double* fx, long nx, double targetx)
{
	return Interpolate_Linear(x, fx, nx, targetx, 0);
}

// LogLinear Interpolation (X변수, Y변수, X길이, 타겟X)
DLLEXPORT(double) Interpolate_LogLinear(double* x, double* fx, long nx, double targetx, long extrapolateflag)
{
	long i;
	double result = 0.0;

	if (nx == 1 || targetx == x[0])
		return fx[0];
	else if (targetx == x[nx - 1])
		return fx[nx - 1];


	if (targetx < x[0])
	{
		if (extrapolateflag == 0) return fx[0];
		else return exp((log(fx[1]) - log(fx[0])) / (x[1] - x[0]) * (targetx - x[0]) + log(fx[0]));
	}
	else if (targetx > x[nx - 1])
	{
		if (extrapolateflag == 0) return fx[nx - 1];
		else return exp((log(fx[nx - 1]) - log(fx[nx - 2])) / (x[nx - 1] - x[nx - 2]) * (targetx - x[nx - 1]) + log(fx[nx - 1]));
	}
	else
	{
		for (i = 1; i < nx; i++)
		{
			if (targetx < x[i])
			{
				result = exp((log(fx[i]) - log(fx[i - 1])) / (x[i] - x[i - 1]) * (targetx - x[i - 1]) + log(fx[i - 1]));
				break;
			}
		}
		return result;
	}
}

// LogLinear Interpolation (X변수, Y변수, X길이, 타겟X)
double Interpolate_LogLinear(double* x, double* fx, long nx, double targetx)
{
	return Interpolate_LogLinear(x, fx, nx, targetx, 0);
}

DLLEXPORT(double) Interpolate_ExpLinear(double* x, double* fx, long nx, double targetx, long extrapolateflag)
{
	long i;
	double result = 0.0;

	if (nx == 1 || targetx == x[0])
		return fx[0];
	else if (targetx == x[nx - 1])
		return fx[nx - 1];


	if (targetx < x[0])
	{
		if (extrapolateflag == 0) return fx[0];
		else return log((exp(fx[1]) - exp(fx[0])) / (x[1] - x[0]) * (targetx - x[0]) + exp(fx[0]));
	}
	else if (targetx > x[nx - 1])
	{
		if (extrapolateflag == 0) return fx[nx - 1];
		else return log((exp(fx[nx - 1]) - exp(fx[nx - 2])) / (x[nx - 1] - x[nx - 2]) * (targetx - x[nx - 1]) + exp(fx[nx - 1]));
	}
	else
	{
		for (i = 1; i < nx; i++)
		{
			if (targetx < x[i])
			{
				result = log((exp(fx[i]) - exp(fx[i - 1])) / (x[i] - x[i - 1]) * (targetx - x[i - 1]) + exp(fx[i - 1]));
				break;
			}
		}
		return result;
	}
}

// LogLinear Interpolation (X변수, Y변수, X길이, 타겟X)
double Interpolate_ExpLinear(double* x, double* fx, long nx, double targetx)
{
	return Interpolate_ExpLinear(x, fx, nx, targetx, 0);
}

// Cholesky 분해된 Array를 리턴 (함수내부 동적할당)
double** Cholesky_Decomposition(double** matrix, long n)
{
	long i, j, k;
	double sum;

	double** lower = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++)
	{
		lower[i] = (double*)calloc(n, sizeof(double));
	}

	// Decomposing a matrix into Lower Triangular
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			sum = 0.0;

			if (j == i) // summation for diagonals
			{
				for (k = 0; k < j; k++)
					sum += pow(lower[j][k], 2);
				lower[j][j] = sqrt(matrix[j][j] -
					sum);
			}
			else
			{
				// Evaluating L(i, j) using L(j, j)
				for (k = 0; k < j; k++)
					sum += (lower[i][k] * lower[j][k]);
				lower[i][j] = (matrix[i][j] - sum) /
					lower[j][j];
			}
		}
	}
	return lower;
}

// Cholesky 분해된 Array를 리턴 (CorrMatrix.reshape(-1) , EmptyMatrix.reshape(-1), Length)
DLLEXPORT(long) Cholesky_Decomposition(double* matrix_reshaped, double* temp_array, long n)
{
	long i, j, k;
	double sum;

	double** matrix = (double**)malloc(n * sizeof(double*));
	k = 0;
	for (i = 0; i < n; i++)
	{
		matrix[i] = (double*)malloc(sizeof(double) * n);
		for (j = 0; j < n; j++)
		{
			matrix[i][j] = matrix_reshaped[k];
			k++;
		}
	}

	double** lower = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++)
	{
		lower[i] = (double*)calloc(n, sizeof(double));
	}

	// Decomposing a matrix into Lower Triangular
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			sum = 0.0;

			if (j == i) // summation for diagonals
			{
				for (k = 0; k < j; k++)
					sum += pow(lower[j][k], 2);
				lower[j][j] = sqrt(matrix[j][j] -
					sum);
			}
			else
			{
				// Evaluating L(i, j) using L(j, j)
				for (k = 0; k < j; k++)
					sum += (lower[i][k] * lower[j][k]);
				lower[i][j] = (matrix[i][j] - sum) /
					lower[j][j];
			}
		}
	}


	k = 0;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			temp_array[k] = lower[i][j];
			k++;
		}

	for (i = 0; i < n; i++)
		if (matrix[i]) free(matrix[i]);
	if (matrix) free(matrix);

	for (i = 0; i < n; i++)
		if (lower[i]) free(lower[i]);
	if (lower) free(lower);


	return 1;
}

// Generate Multivariate Random Number
double** random_mvrn(long number_random, long number_variable, double** Corr)
{
	long i, j, n, k;
	double mvrn;

	double** Chol_Array = Cholesky_Decomposition(Corr, number_variable);
	double* RANDN = (double*)malloc(sizeof(double) * number_variable);

	double** MVRN = (double**)malloc(sizeof(double*) * number_random);
	for (i = 0; i < number_random; i++)
		MVRN[i] = (double*)malloc(sizeof(double) * number_variable);

	for (i = 0; i < number_random; i++)
	{
		for (n = 0; n < number_variable; n++)
		{
			RANDN[n] = randnorm();
		}

		for (j = 0; j < number_variable; j++)
		{
			mvrn = 0.0;
			for (k = 0; k <= j; k++)
			{
				mvrn += Chol_Array[j][k] * RANDN[k];
			}
			MVRN[i][j] = mvrn;
		}
	}

	for (i = 0; i < number_variable; i++)
		if (Chol_Array[i]) free(Chol_Array[i]);
	if (Chol_Array) free(Chol_Array);

	if (RANDN) free(RANDN);
	return MVRN;
}

//	tri-diagonal 행렬 방정식을 빠르게 푼다.(임시 Array : TempAlpha, TempBeta, TempGamma를 준비해서 할당을 하는 시간 및 할당 해제하는 시간을 없앤다.)
//  V_T0결과는 V_T1[1], V_T1[2], ... , V_T1[N]에 덮어씌워서 저장된다.
// [[Beta[0], Gamma[0], 0.0,      0.0,       0.0,        ..........      , 0.0]     [V_T0[0]  ]     [V_T1[0]  ]    
//  [Alpha[1],Beta[1],  Gamma[1], 0.0,       0.0,        ..........      , 0.0]     [V_T0[1]  ]     [V_T1[1]  ]
//  [0.0,     Alpha[2], Beta[2],  Gamma[2],  0.0,        ..........      , 0.0]     [V_T0[2]  ]     [V_T1[2]  ]
//  [0.0,     0.0,      Alpha[3], Beta[3],   Gamma[3],   ........        , 0.0]     [V_T0[3]  ]     [V_T1[3]  ]
//  [~~~,     ~~,       ~~~~~,    ~~~,       ~~~~~~,   ........        , ...  ]   x [     ..  ]  =  [   ..    ]
//  [0.0,     0.0,      .....     Alpha[N-3],Beta[N-3],  Gamma[N-3],       0.0]     [V_T0[N-3]]     [V_T1[N-3]]
//  [0.0,     0.0,      .....,    0.0,       Alpha[N-2], Beta[N-2], Gamma[N-2]]     [V_T0[N-2]]     [V_T1[N-2]]
//  [0.0,     0.0,      .....,    0.0,       0.0,        Alpha[N-1], Beta[N-1]]]    [V_T0[N-1]]     [V_T1[N-1]]
DLLEXPORT(long) Tri_diagonal_Fast(
	double* Alpha,    // Tridiagonal Matrix에서 왼쪽 대각선 행렬
	double* Beta,     // Tridiagonal Matrix에서 중앙 대각선 행렬
	double* Gamma,    // Tridiagonal Matrix에서 오른쪽 대각선 행렬
	double* V_T1,     // V(T1)
	long N,           // 대각선 행렬들의 길이
	double* TempAlpha,// 임시 행렬1 
	double* TempBeta, // 임시 행렬2
	double* TempGamma // 임시 행렬3
)
{
	long i;

	for (i = 0; i <= N - 1; i++)
	{
		TempAlpha[i] = Alpha[i];
		TempBeta[i] = Beta[i];
		TempGamma[i] = Gamma[i];
	}

	TempAlpha[0] = 0.0;
	TempGamma[N - 1] = 0.0;

	if (N <= 1) V_T1[0] /= TempBeta[0];
	else {
		for (i = 1; i <= N - 1; i++) {
			TempAlpha[i] /= TempBeta[i - 1];
			TempBeta[i] -= TempAlpha[i] * TempGamma[i - 1];
			V_T1[i] -= TempAlpha[i] * V_T1[i - 1];
		}
		V_T1[N - 1] /= TempBeta[N - 1];

		for (i = N - 2; i >= 0; i--) V_T1[i] = (V_T1[i] - TempGamma[i] * V_T1[i + 1]) / TempBeta[i];
	}
	return 1;
}

// LocalVol Data와 ,S, T를 넣으면 Interpolated Local Vol을 계산해줌. (느린계산)
double Calc_Volatility(
	// LocalVol Data
	long nt,           // Term 개수        
	long nparity,      // Parity 개수
	double* t,         // Term Array     ex) [0.25,  0.5,   1, ...] 
	double* parity,    // Parity Array   exp)[0.5,   0.7,  1.0, ...]
	double** vol_array,// Vol Array
	// Simulation Data
	double T,          // 시점
	double S           // 현재주가(%)
)
{

	long node_T;
	long node_S;
	double slope_T;
	double slope_S;
	double result_vol;


	if (T <= t[0])
	{
		node_T = 0;
		if (S <= parity[0])                      // T 가 minT 보다 작고, S가 minS 보다 작음
		{
			node_S = 0;
			result_vol = vol_array[node_S][node_T];
		}
		else if (S >= parity[nparity - 1])          // T 가 minT 보다 작고, S가 maxS 보다 큼
		{
			node_S = nparity - 1;
			result_vol = vol_array[node_S][node_T];
		}
		else                                    //  T 가 minT 보다 작고, S가 minS ~ maxS 사이에 있음
		{
			for (node_S = 0; node_S < nparity - 1; node_S++)
				if (S >= parity[node_S] && S < parity[node_S + 1]) break;
			slope_S = (vol_array[node_S + 1][node_T] * vol_array[node_S + 1][node_T] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (parity[node_S + 1] - parity[node_S]);
			// V_s * dS + V(S,t)
			result_vol = sqrt(slope_S * (S - parity[node_S]) + vol_array[node_S][node_T] * vol_array[node_S][node_T]);
		}
	}
	else if (T >= t[nt - 1])
	{
		node_T = nt - 1;
		if (S <= parity[0])                     // T 가 maxT 보다 크고, S가 minS 보다 작음
		{
			node_S = 0;
			result_vol = vol_array[node_S][node_T];
		}
		else if (S >= parity[nparity - 1])        // T 가 maxT 보다 크고, S가 maxS 보다 큼
		{
			node_S = nparity - 1;
			result_vol = vol_array[node_S][node_T];
		}
		else                                   // T 가 maxT 보다 크고 S가 minS ~ maxS 사이에 있음
		{
			for (node_S = 0; node_S < nparity - 1; node_S++)
				if (S >= parity[node_S] && S < parity[node_S + 1]) break;
			slope_S = (vol_array[node_S + 1][node_T] * vol_array[node_S + 1][node_T] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (parity[node_S + 1] - parity[node_S]);
			// V_s * dS + V(S,t)
			result_vol = sqrt(slope_S * (S - parity[node_S]) + vol_array[node_S][node_T] * vol_array[node_S][node_T]);
		}
	}
	else {
		for (node_T = 0; node_T < nt - 1; node_T++)
			if (T >= t[node_T] && T < t[node_T + 1]) break;
		if (S <= parity[0])                   // T가 minT~ maxT 사이에 있고, S가 minS 보다 작음
		{
			node_S = 0;
			slope_T = (vol_array[node_S][node_T + 1] * vol_array[node_S][node_T + 1] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (t[node_T + 1] - t[node_T]);
			// V_t * dt + V(S,t)
			result_vol = sqrt(slope_T * (T - t[node_T]) + vol_array[node_S][node_T] * vol_array[node_S][node_T]);
		}
		else if (S >= parity[nparity - 1])       // T가 minT~ maxT 사이에 있고, S가 maxS 보다 큼
		{
			node_S = nparity - 1;
			slope_T = (vol_array[node_S][node_T + 1] * vol_array[node_S][node_T + 1] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (t[node_T + 1] - t[node_T]);
			// V_t * dt + V(S,t)
			result_vol = sqrt(slope_T * (T - t[node_T]) + vol_array[node_S][node_T] * vol_array[node_S][node_T]);
		}
		else                                  // T가 minT~ maxT 사이에 있고, S가 minS~maxS 사이에 있음
		{
			for (node_S = 0; node_S < nparity - 1; node_S++)
				if (S >= parity[node_S] && S < parity[node_S + 1]) break;
			slope_S = (vol_array[node_S + 1][node_T] * vol_array[node_S + 1][node_T] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (parity[node_S + 1] - parity[node_S]);
			slope_T = (vol_array[node_S][node_T + 1] * vol_array[node_S][node_T + 1] - vol_array[node_S][node_T] * vol_array[node_S][node_T]) / (t[node_T + 1] - t[node_T]);
			// V_t * dt + V_s * dS + V(S,t)
			result_vol = sqrt(slope_S * (S - parity[node_S]) + slope_T * (T - t[node_T]) + vol_array[node_S][node_T] * vol_array[node_S][node_T]);
		}
	}

	return result_vol;
}

// LocalVol Data와 ,S, T를 넣으면 Interpolated Local Vol을 계산해줌. 
// (그때그때 가까운 node 위치를 기억하여 빠른계산)
double Calc_Volatility(long N_Term, long N_Parity, double* Term, double* Parity, double** LocalVolMat, double T, double S, long* Prev_NodeT, long* Prev_NodeS)
{

	long Current_Node_T = 0;
	long Current_Node_S = 0;
	double Slope_T;
	double Slope_S;
	double Result_Vol;

	if (T <= Term[0])
	{
		Current_Node_T = 0;
		Prev_NodeT[0] = 0;
		if (S <= Parity[0])
		{
			Current_Node_S = 0;
			Prev_NodeS[0] = 0;
			Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
			return Result_Vol;
		}
		else if (S >= Parity[N_Parity - 1])
		{
			Current_Node_S = N_Parity - 1;
			Prev_NodeS[0] = N_Parity - 1;
			Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
			return Result_Vol;
		}
		else
		{
			if (Parity[Prev_NodeS[0]] <= S)
			{
				for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			else
			{
				for (Current_Node_S = max(Prev_NodeS[0] - 1, 0); Current_Node_S >= 0; Current_Node_S--)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity[Prev_NodeS[0] + 1] - Parity[Prev_NodeS[0]]);
			// V_s * dS + V(S,t)
			Result_Vol = (Slope_S * (S - Parity[Prev_NodeS[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
			return Result_Vol;
		}
	}
	else if (T >= Term[N_Term - 1])
	{
		Current_Node_T = N_Term - 1;
		Prev_NodeT[0] = N_Term - 1;
		if (S <= Parity[0])
		{
			Current_Node_S = 0;
			Prev_NodeS[0] = 0;
			Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
			return Result_Vol;
		}
		else if (S >= Parity[N_Parity - 1])
		{
			Current_Node_S = N_Parity - 1;
			Prev_NodeS[0] = N_Parity - 1;
			Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
			return Result_Vol;
		}
		else
		{
			if (Parity[Prev_NodeS[0]] <= S)
			{
				for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			else
			{
				for (Current_Node_S = max(Prev_NodeS[0] - 1, 0); Current_Node_S >= 0; Current_Node_S--)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity[Prev_NodeS[0] + 1] - Parity[Prev_NodeS[0]]);
			// V_s * dS + V(S,t)
			Result_Vol = (Slope_S * (S - Parity[Prev_NodeS[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
			return Result_Vol;
		}
	}
	else
	{
		for (Current_Node_T = Prev_NodeT[0]; Current_Node_T < N_Term - 1; Current_Node_T++)
			if (T >= Term[Current_Node_T] && T < Term[Current_Node_T + 1])
			{
				Prev_NodeT[0] = Current_Node_T;
				break;
			}

		if (S <= Parity[0])                   // T가 minT~ maxT 사이에 있고, S가 minS 보다 작음
		{
			Current_Node_S = 0;
			Prev_NodeS[0] = 0;
			Result_Vol = LocalVolMat[Current_Node_S][Prev_NodeT[0]];
			return Result_Vol;
		}
		else if (S >= Parity[N_Parity - 1])
		{
			Current_Node_S = N_Parity - 1;
			Prev_NodeS[0] = N_Parity - 1;
			Result_Vol = LocalVolMat[Current_Node_S][Prev_NodeT[0]];
			return Result_Vol;
		}
		else
		{
			if (Parity[Prev_NodeS[0]] <= S)
			{
				for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			else
			{
				for (Current_Node_S = max(Prev_NodeS[0] - 1, 0); Current_Node_S >= 0; Current_Node_S--)
				{
					if (S >= Parity[Current_Node_S] && S < Parity[Current_Node_S + 1])
					{
						Prev_NodeS[0] = Current_Node_S;
						break;
					}
				}
			}
			Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity[Prev_NodeS[0] + 1] - Parity[Prev_NodeS[0]]);
			Slope_T = (LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0] + 1] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Term[Prev_NodeT[0] + 1] - Term[Prev_NodeT[0]]);
			// V_t * dt + V_s * dS + V(S,t)
			Result_Vol = (Slope_S * (S - Parity[Prev_NodeS[0]]) + Slope_T * (T - Term[Prev_NodeT[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
			return Result_Vol;
		}
	}
}

//Discount Factor 계산
DLLEXPORT(double) Calc_Discount_Factor(
	double* TermArray,
	double* RateArray,
	long LengthArray,
	double T
)
{
	double DF;
	DF = exp(-Interpolate_Linear(TermArray, RateArray, LengthArray, T) * T);

	return DF;
}

DLLEXPORT(double) Calc_Zero_Rate(
	double* TermArray,
	double* RateArray,
	long LengthArray,
	double T
)
{
	return Interpolate_ExpLinear(TermArray, RateArray, LengthArray, T);
}

double Calc_Forward_Rate_Fast(
	double* TermArray, // 기간구조의 기간 Array [0.25,  0.5,   1.0, ....]
	double* RateArray, // 기간구조의 Rate Array [0.008, 0.012, 0.014, ...]
	long LengthArray,  // 기간구조 개수
	double T1,         // Forward Start 시점
	double T2,		   // Forward End 시점
	long* TimePos
)
{
	long i;
	long startidx = *TimePos + 0;
	double dt = T2 - T1;
	double r1, r2;
	double DF1, DF2, FRate;

	if (T1 <= TermArray[0])
	{
		r1 = RateArray[0];
	}
	else if (T1 > TermArray[LengthArray - 1])
	{
		r1 = RateArray[LengthArray - 1];
	}
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T1 < TermArray[i])
			{
				*TimePos = i - 1;
				r1 = (RateArray[i] - RateArray[i - 1]) / (TermArray[i] - TermArray[i - 1]) * (T1 - TermArray[i - 1]) + RateArray[i - 1];
				break;
			}
		}
	}

	if (T2 > T1)
	{
		if (T2 <= TermArray[0])
		{
			r2 = RateArray[0];
		}
		else if (T2 > TermArray[LengthArray - 1])
		{
			r2 = RateArray[LengthArray - 1];
		}
		else
		{
			for (i = max(1, startidx); i < LengthArray; i++)
			{
				if (T2 < TermArray[i])
				{
					r2 = (RateArray[i] - RateArray[i - 1]) / (TermArray[i] - TermArray[i - 1]) * (T2 - TermArray[i - 1]) + RateArray[i - 1];
					break;
				}
			}
		}
	}
	else
	{
		r2 = r1;
	}


	if (T2 == T1) return r1;
	else
	{
		DF1 = exp(-r1 * T1);
		DF2 = exp(-r2 * T2);
		FRate = 1.0 / dt * (DF1 / DF2 - 1.0);
		return FRate;
	}
}

double Calc_Forward_Rate_Daily(
	double* Term,
	double* Rate,
	long LengthArray,
	double T1,
	long* TimePos
)
{
	long i;
	long startidx = *TimePos + 0;
	double dt = 0.00273972602739726;
	double T2 = T1 + dt;
	double r1, r2;
	double DF1, DF2, FRate;

	if (T1 <= Term[0]) r1 = Rate[0];
	else if (T1 >= Term[LengthArray - 1]) r1 = Rate[LengthArray - 1];
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T1 < Term[i])
			{
				*TimePos = i - 1;
				r1 = (Rate[i] - Rate[i - 1]) / (Term[i] - Term[i - 1]) * (T1 - Term[i - 1]) + Rate[i - 1];
				break;
			}
		}
	}

	if (T2 <= Term[0]) r2 = Rate[0];
	else if (T2 >= Term[LengthArray - 1]) r2 = Rate[LengthArray - 1];
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T2 < Term[i])
			{
				r2 = (Rate[i] - Rate[i - 1]) / (Term[i] - Term[i - 1]) * (T2 - Term[i - 1]) + Rate[i - 1];
				break;
			}
		}
	}

	DF1 = exp(-r1 * T1);
	DF2 = exp(-r2 * T2);
	FRate = 1.0 / dt * (DF1 / DF2 - 1.0);
	return FRate;
}

double Calc_Forward_Rate_Daily(
	double* Term,
	double* Rate,
	long LengthArray,
	double T1,
	long* TimePos,
	long NHoliday
)
{
	long i;
	long startidx = *TimePos + 0;
	double dt = 0.00273972602739726;
	double T2 = T1 + ((double)(NHoliday + 1)) * dt;
	double r1, r2;
	double DeltaT = T2 - T1;
	double DF1, DF2, FRate;

	if (T1 <= Term[0]) r1 = Rate[0];
	else if (T1 >= Term[LengthArray - 1]) r1 = Rate[LengthArray - 1];
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T1 < Term[i])
			{
				*TimePos = i - 1;
				r1 = (Rate[i] - Rate[i - 1]) / (Term[i] - Term[i - 1]) * (T1 - Term[i - 1]) + Rate[i - 1];
				break;
			}
		}
	}

	if (T2 <= Term[0]) r2 = Rate[0];
	else if (T2 >= Term[LengthArray - 1]) r2 = Rate[LengthArray - 1];
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T2 < Term[i])
			{
				r2 = (Rate[i] - Rate[i - 1]) / (Term[i] - Term[i - 1]) * (T2 - Term[i - 1]) + Rate[i - 1];
				break;
			}
		}
	}

	DF1 = exp(-r1 * T1);
	DF2 = exp(-r2 * T2);
	FRate = 1.0 / DeltaT * (DF1 / DF2 - 1.0);
	return FRate;
}

double Calc_Forward_FXVol_Daily(
	double* TermArray, // 기간구조의 기간 Array [0.25,  0.5,   1.0, ....]
	double* FXVolArray, // 기간구조의 Rate Array [0.008, 0.012, 0.014, ...]
	long LengthArray,  // 기간구조 개수
	double T1,         // Forward Start 시점
	long* TimePos
)
{
	long i;
	long startidx = *TimePos + 0;
	double dt = 0.00273972602739726;
	double T2 = T1 + dt;
	double V1, V2;
	double DF1, DF2, FVar, FVol;

	if (T1 <= TermArray[0])
	{
		V1 = FXVolArray[0];
	}
	else if (T1 > TermArray[LengthArray - 1])
	{
		V1 = FXVolArray[LengthArray - 1];
	}
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T1 < TermArray[i])
			{
				*TimePos = i - 1;
				V1 = (FXVolArray[i] - FXVolArray[i - 1]) / (TermArray[i] - TermArray[i - 1]) * (T1 - TermArray[i - 1]) + FXVolArray[i - 1];
				break;
			}
		}
	}

	if (T2 <= TermArray[0])
	{
		V2 = FXVolArray[0];
	}
	else if (T2 > TermArray[LengthArray - 1])
	{
		V2 = FXVolArray[LengthArray - 1];
	}
	else
	{
		for (i = max(1, startidx); i < LengthArray; i++)
		{
			if (T2 < TermArray[i])
			{
				V2 = (FXVolArray[i] - FXVolArray[i - 1]) / (TermArray[i] - TermArray[i - 1]) * (T2 - TermArray[i - 1]) + FXVolArray[i - 1];
				break;
			}
		}
	}

	FVar = (T2 * V2 * V2 - T1 * V1 * V1) / dt;
	if (FVar > 0.0) FVol = sqrt(FVar);
	else FVol = V1;

	return FVol;
}

//Forward Rate 계산
DLLEXPORT(double) Calc_Forward_Rate(
	double* TermArray, // 기간구조의 기간 Array [0.25,  0.5,   1.0, ....]
	double* RateArray, // 기간구조의 Rate Array [0.008, 0.012, 0.014, ...]
	long LengthArray,  // 기간구조 개수
	double T1,         // Forward Start 시점
	double T2          // Forward End 시점
)
{
	double DFT1;
	double DFT2;
	double FRate;
	//Linear Interpolation
	if (T1 < TermArray[LengthArray - 1] && T1> TermArray[0])
	{
		DFT1 = exp(-Interpolate_Linear(TermArray, RateArray, LengthArray, T1) * T1);
		DFT2 = exp(-Interpolate_Linear(TermArray, RateArray, LengthArray, T2) * T2);
		FRate = 1.0 / (T2 - T1) * (DFT1 / DFT2 - 1.0);
	}
	else if (T1 <= TermArray[0])
	{
		DFT1 = exp(-RateArray[0] * T1);
		DFT2 = exp(-Interpolate_Linear(TermArray, RateArray, LengthArray, T2) * T2);
		FRate = 1.0 / (T2 - T1) * (DFT1 / DFT2 - 1.0);
	}
	else
	{
		DFT1 = exp(-RateArray[LengthArray - 1] * T1);
		DFT2 = exp(-RateArray[LengthArray - 1] * T2);
		FRate = 1.0 / (T2 - T1) * (DFT1 / DFT2 - 1.0);
	}

	return FRate;
}

// matrix inversioon
// the result is put in Y
void MatrixInversion(double** A, long order, double** Y)
{
	long i, j, k, n;
	// MinorMatrixList 추가 2022.08.30 임대선 미리 만들어두고 반복사용하기
	double*** MinorMatrixList = (double***)malloc(sizeof(double**) * order);
	for (i = 0; i < order; i++)
	{
		n = max(1, i);
		MinorMatrixList[i] = (double**)malloc(sizeof(double*) * n);
		for (j = 0; j < n; j++)
		{
			MinorMatrixList[i][j] = (double*)malloc(sizeof(double) * n);
		}
	}
	// get the determinant of a
	double det = 1.0 / CalcDeterminant(A, order, MinorMatrixList);

	// memory allocation
	double* temp = new double[(order - 1) * (order - 1)];
	double** minor = new double* [order - 1];
	for (i = 0; i < order - 1; i++)
		minor[i] = temp + (i * (order - 1));

	for (j = 0; j < order; j++)
	{
		for (i = 0; i < order; i++)
		{
			// get the co-factor (matrix) of A(j,i)
			GetMinor(A, minor, j, i, order);
			Y[i][j] = det * CalcDeterminant(minor, order - 1, MinorMatrixList);
			if ((i + j) % 2 == 1)
				Y[i][j] = -Y[i][j];
		}
	}

	// release memory
	delete[] temp;
	delete[] minor;

	for (i = 0; i < order; i++)
	{
		n = max(1, i);
		for (j = 0; j < n; j++)
		{
			free(MinorMatrixList[i][j]);
		}
		free(MinorMatrixList[i]);
	}
	free(MinorMatrixList);
}

// calculate the cofactor of element (row,col)
long GetMinor(double** src, double** dest, long row, long col, long order)
{
	// indicate which col and row is being copied to dest
	long colCount = 0, rowCount = 0;

	for (long i = 0; i < order; i++)
	{
		if (i != row)
		{
			colCount = 0;
			for (long j = 0; j < order; j++)
			{
				// when j is not the element
				if (j != col)
				{
					dest[rowCount][colCount] = src[i][j];
					colCount++;
				}
			}
			rowCount++;
		}
	}

	return 1;
}

// Calculate the determinant recursively.
double CalcDeterminant(double** mat, long order, double*** MinorMatrixList)
{

	if (order == 1) return mat[0][0];
	else if (order == 2) return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	else if (order == 3) return mat[0][0] * mat[1][1] * mat[2][2] + mat[0][1] * mat[1][2] * mat[2][0] + mat[0][2] * mat[1][0] * mat[2][1] - (mat[0][2] * mat[1][1] * mat[2][0] + mat[0][1] * mat[1][0] * mat[2][2] + mat[0][0] * mat[1][2] * mat[2][1]);
	else if (order == 4)
	{
		double p1, p2, p3, p4, detA;
		p1 = mat[0][0] * (mat[1][1] * mat[2][2] * mat[3][3] + mat[1][2] * mat[2][3] * mat[3][1] + mat[1][3] * mat[2][1] * mat[3][2] -
			mat[1][3] * mat[2][2] * mat[3][1] - mat[1][2] * mat[2][1] * mat[3][3] - mat[1][1] * mat[2][3] * mat[3][2]);

		p2 = -mat[1][0] * (mat[0][1] * mat[2][2] * mat[3][3] + mat[0][2] * mat[2][3] * mat[3][1] + mat[0][3] * mat[2][1] * mat[3][2] -
			mat[0][3] * mat[2][2] * mat[3][1] - mat[0][2] * mat[2][1] * mat[3][3] - mat[0][1] * mat[2][3] * mat[3][2]);

		p3 = mat[2][0] * (mat[0][1] * mat[1][2] * mat[3][3] + mat[0][2] * mat[1][3] * mat[3][1] + mat[0][3] * mat[1][1] * mat[3][2] -
			mat[0][3] * mat[1][2] * mat[3][1] - mat[0][2] * mat[1][1] * mat[3][3] - mat[0][1] * mat[1][3] * mat[3][2]);

		p4 = -mat[3][0] * (mat[0][1] * mat[1][2] * mat[2][3] + mat[0][2] * mat[1][3] * mat[2][1] + mat[0][3] * mat[1][1] * mat[2][2] -
			mat[0][3] * mat[1][2] * mat[2][1] - mat[0][2] * mat[1][1] * mat[2][3] - mat[0][1] * mat[1][3] * mat[2][2]);
		detA = p1 + p2 + p3 + p4;
		return detA;
	}
	double det = 0;

	double** minor;
	minor = MinorMatrixList[max(0, order - 1)];
	//minor = new double* [order - 1];
	//for (long i = 0; i < order - 1; i++)
	//	minor[i] = new double[order - 1];

	for (long i = 0; i < order; i++)
	{
		// get minor of element (0,i)
		GetMinor(mat, minor, 0, i, order);
		// the recusion is here!
		det += (i % 2 == 1 ? -1.0 : 1.0) * mat[0][i] * CalcDeterminant(minor, order - 1, MinorMatrixList);

	}

	// release memory
	//for (long i = 0; i < order - 1; i++)
	//	delete[] minor[i];
	//delete[] minor;

	return det;
}

void bubble_sort(double* arr, long count, long ascending)
{
	double temp;
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

void bubble_sort(double* arr, long count)
{
	bubble_sort(arr, count, 1);
}

void bubble_sort(double* arr, long* Index, long count, long ascending)
{
	double temp;
	long temp_idx;
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
					temp_idx = Index[j];
					arr[j] = arr[j + 1];
					Index[j] = Index[j + 1];
					arr[j + 1] = temp;
					Index[j + 1] = temp_idx;
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
					temp_idx = Index[j];
					arr[j] = arr[j + 1];
					Index[j] = Index[j + 1];
					arr[j + 1] = temp;
					Index[j + 1] = temp_idx;
				}
			}
		}
	}
}

void bubble_sort(double* arr, long* Index, long count)
{
	bubble_sort(arr, Index, count, 1);
}

void bubble_sort_index(long* index_arr, double* value_arr, long count, long ascending)
{
	long temp;
	double temp_value;
	long i, j;

	if (ascending == 1)
	{
		for (i = 0; i < count; i++)
		{
			for (j = 0; j < count - 1; j++)
			{
				if (index_arr[j] > index_arr[j + 1])
				{
					temp = index_arr[j];
					temp_value = value_arr[i];
					index_arr[j] = index_arr[j + 1];
					value_arr[j] = value_arr[j + 1];
					index_arr[j + 1] = temp;
					value_arr[j + 1] = temp_value;
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
				if (index_arr[j] < index_arr[j + 1])
				{
					temp = index_arr[j];
					temp_value = value_arr[i];
					index_arr[j] = index_arr[j + 1];
					value_arr[j] = value_arr[j + 1];
					index_arr[j + 1] = temp;
					value_arr[j + 1] = temp_value;
				}
			}
		}
	}
}

void bubble_sort_index(long* index_arr, double* value_arr, long count)
{
	bubble_sort_index(index_arr, value_arr, count, 1);
}

long* RankArray(double* list, long list_length, long ascending)
{
	long i, j, max_rank = -1;
	long* rank = (long*)malloc(list_length * sizeof(long));
	for (i = 0; i < list_length; i++) rank[i] = 1;

	for (i = 0; i < list_length; i++)
	{
		for (j = i + 1; j < list_length; j++)
		{
			if (list[i] > list[j]) rank[j]++;
			else if (list[i] < list[j]) rank[i]++;
		}
		max_rank = max(max_rank, rank[i]);
	}

	if (ascending != 0)
		for (i = 0; i < list_length; i++)
			rank[i] = max_rank - rank[i] + 1;
	return rank;
}

double** Dot2dArray(double** A, long shape_A[2], double** B, long shape_B[2])
{
	long p, q, i, j, k, maxk;
	double s;
	p = shape_A[0];
	q = shape_B[1];
	double** ResultArray = (double**)malloc(sizeof(double*) * p);
	for (i = 0; i < p; i++)
		ResultArray[i] = (double*)calloc(q, sizeof(double));

	if (shape_A[1] != shape_B[0])
	{
		return ResultArray;
	}
	else
	{
		maxk = shape_A[1];
		for (i = 0; i < p; i++)
		{
			for (j = 0; j < q; j++)
			{
				s = 0.0;
				for (k = 0; k < maxk; k++)
				{
					s += A[i][k] * B[k][j];;
				}
				ResultArray[i][j] = s;
			}
		}
		return ResultArray;
	}
}

void Dot2dArray(double** A, long shape_A[2], double** B, long shape_B[2], double** ResultArray)
{
	long p, q, i, j, k, maxk;
	double s;
	p = shape_A[0];
	q = shape_B[1];

	maxk = shape_A[1];
	for (i = 0; i < p; i++)
	{
		for (j = 0; j < q; j++)
		{
			s = 0.0;
			for (k = 0; k < maxk; k++)
			{
				s += A[i][k] * B[k][j];;
			}
			ResultArray[i][j] = s;
		}
	}
}

double** XprimeDotX(double** X, long ShapeX[2])
{
	long i, j, k;
	long ShapeXpX[2] = { ShapeX[1], ShapeX[1] };
	long n = ShapeX[0];
	long idx_row;
	long idx_col;

	double s = 0.0;
	double a, b;

	double** XpDotX = make_array(ShapeXpX[0], ShapeXpX[1]);
	for (i = 0; i < ShapeXpX[0]; i++)
		for (j = 0; j < ShapeXpX[1]; j++)
		{
			s = 0.0;
			for (k = 0; k < n; k++)
			{
				s += X[k][i] * X[k][j];
			}
			XpDotX[i][j] = s;
		}
	return XpDotX;
}

void XprimeDotX(double** X, long ShapeX[2], double** XpDotX)
{
	long i, j, k;
	long ShapeXpX[2] = { ShapeX[1], ShapeX[1] };
	long n = ShapeX[0];
	long idx_row;
	long idx_col;

	double s = 0.0;
	double a, b;
	for (i = 0; i < ShapeXpX[0]; i++)
		for (j = 0; j < ShapeXpX[1]; j++)
		{
			s = 0.0;
			for (k = 0; k < n; k++)
			{
				s += X[k][i] * X[k][j];
			}
			XpDotX[i][j] = s;
		}
}

long XprimeY(double** X, long shape_X[2], double* Y, long LengthY, double** XprimeYMatrix)
{
	long i, j, k;
	long p = shape_X[1];
	long n = shape_X[0];
	double s = 0.0;
	for (i = 0; i < p; i++)
	{
		s = 0.0;
		for (j = 0; j < n; j++)
		{
			s += X[j][i] * Y[j];
		}
		XprimeYMatrix[i][0] = s;
	}
	return 1;
}

double* OLSBeta(double* Y, long LengthY, double** X, long ShapeX[2])
{
	long i;

	double* Result = make_array(ShapeX[1]);

	double** XpX = XprimeDotX(X, ShapeX);

	long Shape_XpX_Inv[2] = { ShapeX[1], ShapeX[1] };
	double** XpX_Inv = make_array(ShapeX[1], ShapeX[1]);
	MatrixInversion(XpX, ShapeX[1], XpX_Inv);

	long Shape_XpdotY[2] = { ShapeX[1], 1 };
	double** XpdotY = make_array(ShapeX[1], 1);
	i = XprimeY(X, ShapeX, Y, LengthY, XpdotY);

	double** XpXInvXpY = Dot2dArray(XpX_Inv, Shape_XpX_Inv, XpdotY, Shape_XpdotY);

	for (i = 0; i < ShapeX[1]; i++) Result[i] = XpXInvXpY[i][0];

	for (i = 0; i < ShapeX[1]; i++) free(XpX[i]);
	free(XpX);
	for (i = 0; i < ShapeX[1]; i++) free(XpX_Inv[i]);
	free(XpX_Inv);
	for (i = 0; i < ShapeX[1]; i++) free(XpdotY[i]);
	free(XpdotY);
	for (i = 0; i < ShapeX[1]; i++) free(XpXInvXpY[i]);
	free(XpXInvXpY);

	return Result;
}


void Calc_C(long nRate, double* Term, double* Rate, double* CArray)
{
	long i, n;
	double h1, h0;
	long Length = nRate - 2;

	CArray[0] = 0.0;
	CArray[nRate - 1] = 0.0;

	double* RHS = CArray + 1;
	double* Alpha = (double*)malloc(sizeof(double) * (nRate - 2));
	double* Beta = (double*)malloc(sizeof(double) * (nRate - 2));
	double* Gamma = (double*)malloc(sizeof(double) * (nRate - 2));

	n = 0;
	for (i = 1; i < nRate - 1; i++)
	{
		h1 = Term[i + 1] - Term[i];
		h0 = Term[i] - Term[i - 1];
		if (n == 0)
		{
			Alpha[n] = 0.0;
			Beta[n] = 2.0 * (h0 + h1);
			Gamma[n] = h1;
		}
		else if (n == nRate - 2)
		{
			Gamma[n] = 0.0;
			Beta[n] = 2.0 * (h0 + h1);
			Gamma[n] = h1;
		}
		else {
			Alpha[n] = h0;
			Beta[n] = 2.0 * (h0 + h1);
			Gamma[n] = h1;
		}
		RHS[i - 1] = 3.0 * ((Rate[i + 1] - Rate[i]) / h1 - (Rate[i] - Rate[i - 1]) / h0);
		n++;
	}

	for (i = 1; i <= Length - 1; i++) {
		Alpha[i] /= Beta[i - 1];
		Beta[i] -= Alpha[i] * Gamma[i - 1];
		RHS[i] -= Alpha[i] * RHS[i - 1];
	}
	RHS[Length - 1] /= Beta[Length - 1];

	for (i = Length - 2; i >= 0; i--) RHS[i] = (RHS[i] - Gamma[i] * RHS[i + 1]) / Beta[i];

	free(Alpha);
	free(Beta);
	free(Gamma);
}

double CubicSpline(long nRate, double* Term, double* Rate, double TargetTerm)
{
	long i;
	double a, b, d, y, hi, xp;

	if (nRate < 4) return Interpolate_Linear(Term, Rate, nRate, TargetTerm);

	double* CArray = (double*)malloc(sizeof(double) * nRate);
	Calc_C(nRate, Term, Rate, CArray);

	if (Term[0] >= TargetTerm)
	{
		hi = Term[1] - Term[0];
		xp = TargetTerm - Term[0];
		a = Rate[0];
		b = (Rate[1] - Rate[0]) / hi - hi * (2.0 * CArray[0] + CArray[1]) / 3.0;
		d = (CArray[1] - CArray[0]) / (3.0 * hi);
		y = a + b * xp + CArray[0] * xp * xp + d * xp * xp * xp;
	}
	else if (Term[nRate - 1] <= TargetTerm)
	{
		hi = Term[nRate - 1] - Term[nRate - 2];
		xp = TargetTerm - Term[nRate - 2];
		a = Rate[nRate - 2];
		b = (Rate[nRate - 1] - Rate[nRate - 2]) / hi - hi * (2.0 * CArray[nRate - 2] + CArray[nRate - 1]) / 3.0;
		d = (CArray[nRate - 1] - CArray[nRate - 2]) / (3.0 * hi);
		y = a + b * xp + CArray[nRate - 2] * xp * xp + d * xp * xp * xp;
	}
	else
	{
		for (i = 1; i < nRate; i++)
		{
			if (Term[i] > TargetTerm) {
				hi = Term[i] - Term[i - 1];
				xp = TargetTerm - Term[i - 1];
				a = Rate[i - 1];
				b = (Rate[i] - Rate[i - 1]) / hi - hi * (2.0 * CArray[i - 1] + CArray[i]) / 3.0;
				d = (CArray[i] - CArray[i - 1]) / (3.0 * hi);
				y = a + b * xp + CArray[i - 1] * xp * xp + d * xp * xp * xp;
				break;
			}
		}
	}

	free(CArray);
	return y;
}

long Calibration_CubicSpline_Params(
	long nRate,
	double* Term,
	double* Rate,
	double* C_Array	// Out: 2차 계수 Param
)
{
	long ResultCode = 1;

	if (nRate < 4)
	{
		ResultCode = -1;
		return ResultCode;
	}

	Calc_C(nRate, Term, Rate, C_Array);
	return ResultCode;
}

double CubicInterpolation(long nRate, double* Term, double* Rate, double* C_Array, double TargetTerm)
{
	long i;
	double a, b, d, y, hi, xp;

	if (nRate < 4) return Interpolate_Linear(Term, Rate, nRate, TargetTerm);

	if (Term[0] >= TargetTerm)
	{
		// Extrapolation
		//hi = Term[1] - Term[0];
		//xp = TargetTerm - Term[0];
		//a = Rate[0];
		//b = (Rate[1] - Rate[0])/ hi - hi * (2.0 * C_Array[0] + C_Array[1]) / 3.0;
		//d = (C_Array[1] - C_Array[0]) / (3.0 * hi);
		//y = a + b * xp + C_Array[0] * xp * xp + d * xp * xp * xp;
		y = Rate[0];
	}
	else if (Term[nRate - 1] <= TargetTerm)
	{
		// Extrapolation
		//hi = Term[nRate - 1] - Term[nRate - 2];
		//xp = TargetTerm - Term[nRate - 2];
		//a = Rate[nRate - 2];
		//b = (Rate[nRate - 1] - Rate[nRate - 2]) / hi - hi * (2.0 * C_Array[nRate - 2] + C_Array[nRate - 1]) / 3.0;
		//d = (C_Array[nRate - 1] - C_Array[nRate - 2]) / (3.0 * hi);
		//y = a + b * xp + C_Array[nRate - 2] * xp * xp + d * xp * xp * xp;
		y = Rate[nRate - 1];
	}
	else
	{
		for (i = 1; i < nRate; i++)
		{
			if (Term[i] > TargetTerm) {
				hi = Term[i] - Term[i - 1];
				xp = TargetTerm - Term[i - 1];
				a = Rate[i - 1];
				b = (Rate[i] - Rate[i - 1]) / hi - hi * (2.0 * C_Array[i - 1] + C_Array[i]) / 3.0;
				d = (C_Array[i] - C_Array[i - 1]) / (3.0 * hi);
				y = a + b * xp + C_Array[i - 1] * xp * xp + d * xp * xp * xp;
				break;
			}
		}
	}
	return y;
}

