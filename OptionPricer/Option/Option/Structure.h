#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef STRUCTURE
#define STRUCTURE 1
#endif

double dw_over_dt(double** w, long n_parity, long n_term, long p, long q, double dt);
double dw_over_dy(double** w, long n_parity, long n_term, long p, long q, double* k, double dk);
double dwdw_over_dydy(double** w, long n_parity, long n_term, long p, long q, double* k, double dk);
double d1(double S, double K, double T, double r, double div, double Vol);
double dImvol_over_dK(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dK);
double dImvoldImvol_over_dKdK(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dK);
double dImvol_over_dT(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dT);
void fillna_Interpolate(double** Matrix, long N_Index, long N_Column);

//커브정보
class curveinfo {
private:
	long N_Term;
	long dynamicflag; // Implied vol dynamic flag 0: 할당되지 않음  1: 완전 Copy 할당
public:
	double* Term;
	double* Rate;

	curveinfo() { dynamicflag = 0; }

	void initialize(long length_term, double* termarray, double* ratearray)
	{
		if (dynamicflag == 0)
		{
			long i;
			dynamicflag = 1;
			N_Term = length_term;
			Term = (double*)malloc(sizeof(double) * length_term);
			Rate = (double*)malloc(sizeof(double) * length_term);
			for (i = 0; i < length_term; i++)
			{
				Term[i] = termarray[i];
				Rate[i] = ratearray[i];
			}
		}
	}

	curveinfo(long length_term, double* termarray, double* ratearray)
	{
		dynamicflag = 0;
		initialize(length_term, termarray, ratearray);
	}

	~curveinfo()
	{
		if (dynamicflag != 0)
		{
			free(Term);
			free(Rate);
		}
	}

	long nterm()
	{
		return N_Term;
	}

	double Interpolated_Rate(double targetx, long extrapolateflag)
	{
		long i;
		double result = 0.0;
		if (dynamicflag != 0)
		{
			if (N_Term == 1 || targetx == Term[0])
				return Rate[0];
			else if (targetx == Term[N_Term - 1])
				return Rate[N_Term - 1];


			if (targetx < Term[0])
			{
				if (extrapolateflag == 0) return Rate[0];
				else return (Rate[1] - Rate[0]) / (Term[1] - Term[0]) * (targetx - Term[0]) + Rate[0];
			}
			else if (targetx > Term[N_Term - 1])
			{
				if (extrapolateflag == 0) return Rate[N_Term - 1];
				else return (Rate[N_Term - 1] - Rate[N_Term - 2]) / (Term[N_Term - 1] - Term[N_Term - 2]) * (targetx - Term[N_Term - 1]) + Rate[N_Term - 1];
			}
			else
			{
				for (i = 1; i < N_Term; i++)
				{
					if (targetx < Term[i])
					{
						result = (Rate[i] - Rate[i - 1]) / (Term[i] - Term[i - 1]) * (targetx - Term[i - 1]) + Rate[i - 1];
						break;
					}
				}
				return result;
			}
		}
		else return 0.0;
	}

	double Interpolated_Rate(double targetx)
	{
		return Interpolated_Rate(targetx, 0);
	}

	void Rate_Up(double percentpoint)
	{
		long i, nt;
		if (dynamicflag != 0)
		{
			nt = N_Term;
			for (i = 0; i < nt; i++)
			{
				Rate[i] = Rate[i] + percentpoint;
			}
		}
	}
};

//변동성정보
class volinfo {
private:
	long N_Parity;
	long N_Term;
	long dynamicflag = 0;          // Implied vol dynamic flag 0: 할당되지 않음 1:포인터 위치로 할당 2: 완전 Copy
	long localvol_dynamicflag = 0; // local vol dynamic flag 0: 계산되지 않음 1:계산되어 할당됨

public:
	double** Vol_Matrix;       //ImpliedVolatily
	double** LocalVolMat;      //LocalVolatility
	double* Term;
	double* Parity;
	double* Term_Locvol;
	double* Parity_Locvol;     //Parity 보정 가능하게 하기 위해 Public

	void copyclass(long N_Parity2, long N_Term2, double** Vol_Matrix2, double** LocalVolMat2, double* Term2, double* Parity2, double* Term_Locvol2, double* Parity_Locvol2)
	{
		long i, j;
		N_Parity = N_Parity2;
		N_Term = N_Term2;
		dynamicflag = 2;
		localvol_dynamicflag = 1;
		Vol_Matrix = (double**)malloc(sizeof(double*) * N_Parity);
		LocalVolMat = (double**)malloc(sizeof(double*) * N_Parity);
		for (i = 0; i < N_Parity; i++)
		{
			Vol_Matrix[i] = (double*)malloc(sizeof(double) * N_Term);
			LocalVolMat[i] = (double*)malloc(sizeof(double) * N_Term);
			for (j = 0; j < N_Term; j++)
			{
				Vol_Matrix[i][j] = Vol_Matrix2[i][j];
				LocalVolMat[i][j] = LocalVolMat2[i][j];
			}
		}
		Term = (double*)malloc(sizeof(double) * N_Term);
		for (i = 0; i < N_Term; i++) Term[i] = Term2[i];

		Parity = (double*)malloc(sizeof(double) * N_Parity);
		for (i = 0; i < N_Parity; i++) Parity[i] = Parity2[i];

		Term_Locvol = (double*)malloc(sizeof(double) * N_Term);
		for (i = 0; i < N_Term; i++) Term_Locvol[i] = Term_Locvol2[i];

		Parity_Locvol = (double*)malloc(sizeof(double) * N_Parity);
		for (i = 0; i < N_Parity; i++) Parity_Locvol[i] = Parity_Locvol2[i];

	}

	volinfo()
	{
		dynamicflag = 0;
		localvol_dynamicflag = 0;
	}

	// 포인터로 할당
	void softcopy(long length_parity, double* parityarray, long length_term, double* termarray, double* reshapedvol)
	{
		if (dynamicflag == 0)
		{
			long i;
			dynamicflag = 1;
			N_Parity = length_parity;
			N_Term = length_term;
			Term = termarray;
			Parity = parityarray;
			Vol_Matrix = (double**)malloc(sizeof(double*) * N_Parity);
			for (i = 0; i < N_Parity; i++)
			{
				Vol_Matrix[i] = reshapedvol + i * N_Term;
			}
		}
	}

	// 완전 Copy
	void hardcopy(long length_parity, double* parityarray, long length_term, double* termarray, double* reshapedvol)
	{
		if (dynamicflag == 0)
		{
			long i, j, k;
			dynamicflag = 2;
			N_Parity = length_parity;
			N_Term = length_term;
			Term = (double*)malloc(sizeof(double) * N_Term);
			for (i = 0; i < N_Term; i++)
			{
				Term[i] = termarray[i];
			}
			Parity = (double*)malloc(sizeof(double) * N_Parity);
			for (i = 0; i < N_Parity; i++)
			{
				Parity[i] = parityarray[i];
			}
			Vol_Matrix = (double**)malloc(sizeof(double*) * N_Parity);
			k = 0;
			for (i = 0; i < N_Parity; i++)
			{
				Vol_Matrix[i] = (double*)malloc(sizeof(double) * N_Term);
				for (j = 0; j < N_Term; j++)
				{
					Vol_Matrix[i][j] = reshapedvol[k];
					k++;
				}
			}
		}
	}

	void hardcopyUp(long length_parity, double* parityarray, long length_term, double* termarray, double* reshapedvol, double percentpoint)
	{
		if (dynamicflag == 0)
		{
			long i, j, k;
			dynamicflag = 2;
			N_Parity = length_parity;
			N_Term = length_term;
			Term = (double*)malloc(sizeof(double) * N_Term);
			for (i = 0; i < N_Term; i++)
			{
				Term[i] = termarray[i];
			}
			Parity = (double*)malloc(sizeof(double) * N_Parity);
			for (i = 0; i < N_Parity; i++)
			{
				Parity[i] = parityarray[i];
			}
			Vol_Matrix = (double**)malloc(sizeof(double*) * N_Parity);
			k = 0;
			for (i = 0; i < N_Parity; i++)
			{
				Vol_Matrix[i] = (double*)malloc(sizeof(double) * N_Term);
				for (j = 0; j < N_Term; j++)
				{
					Vol_Matrix[i][j] = max(0.0001, reshapedvol[k] + percentpoint);
					k++;
				}
			}
		}
	}

	//변동성정보 생성자 (변동성 매핑: Parity길이, Parity, Term길이, Term, reshapedvol , copyflag)
	volinfo(long length_parity, double* parityarray, long length_term, double* termarray, double* reshapedvol, long copyflag)
	{
		if (copyflag == 1)
		{
			softcopy(length_parity, parityarray, length_term, termarray, reshapedvol);
			dynamicflag = 1;
		}
		else
		{
			hardcopy(length_parity, parityarray, length_term, termarray, reshapedvol);
			dynamicflag = 2;
		}
	}

	//변동성정보 생성자 (변동성 매핑: Parity길이, Parity, Term길이, Term, reshapedvol )
	volinfo(long length_parity, double* parityarray, long length_term, double* termarray, double* reshapedvol)
	{
		hardcopy(length_parity, parityarray, length_term, termarray, reshapedvol);
		dynamicflag = 2;
	}

	~volinfo()
	{
		long i;
		if (dynamicflag == 1)
		{
			if (Vol_Matrix) free(Vol_Matrix);
		}
		else if (dynamicflag == 0)
		{

		}
		else
		{
			for (i = 0; i < N_Parity; i++)
			{
				if (Vol_Matrix[i]) free(Vol_Matrix[i]);
			}
			if (Term) free(Term);
			if (Parity) free(Parity);
			if (Vol_Matrix) free(Vol_Matrix);
			dynamicflag = 0;
		}

		if (localvol_dynamicflag == 1)
		{
			for (i = 0; i < N_Parity; i++)
			{
				if (LocalVolMat[i]) free(LocalVolMat[i]);
			}
			if (LocalVolMat) free(LocalVolMat);
			if (Term_Locvol) free(Term_Locvol);
			if (Parity_Locvol) free(Parity_Locvol);
		}
	}
	long size(long dimension)
	{
		if (dynamicflag != 0)
		{
			if (dimension == 0)	return N_Parity;
			else return N_Term;
		}
		else
			return 0;
	}
	long size()
	{
		return size(0);
	}

	//Implied Data와 ,S, T를 넣으면 Interpolated Implied Vol을 계산해줌. (느린계산)
	double Calc_Implied_Volatility(double T, double S)
	{
		long node_T;
		long node_S;
		double slope_T;
		double slope_S;
		double result_vol;

		if (dynamicflag != 0)
		{
			if (T <= Term[0])
			{
				node_T = 0;
				if (S <= Parity[0])                      // T 가 minT 보다 작고, S가 minS 보다 작음
				{
					node_S = 0;
					result_vol = Vol_Matrix[node_S][node_T];
				}
				else if (S >= Parity[N_Parity - 1])          // T 가 minT 보다 작고, S가 maxS 보다 큼
				{
					node_S = N_Parity - 1;
					result_vol = Vol_Matrix[node_S][node_T];
				}
				else                                    //  T 가 minT 보다 작고, S가 minS ~ maxS 사이에 있음
				{
					for (node_S = 0; node_S < N_Parity - 1; node_S++)
						if (S >= Parity[node_S] && S < Parity[node_S + 1]) break;
					slope_S = (Vol_Matrix[node_S + 1][node_T] * Vol_Matrix[node_S + 1][node_T] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Parity[node_S + 1] - Parity[node_S]);
					// V_s * dS + V(S,t)
					result_vol = sqrt(slope_S * (S - Parity[node_S]) + Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]);
				}
			}
			else if (T >= Term[N_Term - 1])
			{
				node_T = N_Term - 1;
				if (S <= Parity[0])                     // T 가 maxT 보다 크고, S가 minS 보다 작음
				{
					node_S = 0;
					result_vol = Vol_Matrix[node_S][node_T];
				}
				else if (S >= Parity[N_Parity - 1])        // T 가 maxT 보다 크고, S가 maxS 보다 큼
				{
					node_S = N_Parity - 1;
					result_vol = Vol_Matrix[node_S][node_T];
				}
				else                                   // T 가 maxT 보다 크고 S가 minS ~ maxS 사이에 있음
				{
					for (node_S = 0; node_S < N_Parity - 1; node_S++)
						if (S >= Parity[node_S] && S < Parity[node_S + 1]) break;
					slope_S = (Vol_Matrix[node_S + 1][node_T] * Vol_Matrix[node_S + 1][node_T] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Parity[node_S + 1] - Parity[node_S]);
					// V_s * dS + V(S,t)
					result_vol = sqrt(slope_S * (S - Parity[node_S]) + Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]);
				}
			}
			else {
				for (node_T = 0; node_T < N_Term - 1; node_T++)
					if (T >= Term[node_T] && T < Term[node_T + 1]) break;
				if (S <= Parity[0])                   // T가 minT~ maxT 사이에 있고, S가 minS 보다 작음
				{
					node_S = 0;
					slope_T = (Vol_Matrix[node_S][node_T + 1] * Vol_Matrix[node_S][node_T + 1] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Term[node_T + 1] - Term[node_T]);
					// V_t * dt + V(S,t)
					result_vol = sqrt(slope_T * (T - Term[node_T]) + Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]);
				}
				else if (S >= Parity[N_Parity - 1])       // T가 minT~ maxT 사이에 있고, S가 maxS 보다 큼
				{
					node_S = N_Parity - 1;
					slope_T = (Vol_Matrix[node_S][node_T + 1] * Vol_Matrix[node_S][node_T + 1] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Term[node_T + 1] - Term[node_T]);
					// V_t * dt + V(S,t)
					result_vol = sqrt(slope_T * (T - Term[node_T]) + Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]);
				}
				else                                  // T가 minT~ maxT 사이에 있고, S가 minS~maxS 사이에 있음
				{
					for (node_S = 0; node_S < N_Parity - 1; node_S++)
						if (S >= Parity[node_S] && S < Parity[node_S + 1]) break;
					slope_S = (Vol_Matrix[node_S + 1][node_T] * Vol_Matrix[node_S + 1][node_T] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Parity[node_S + 1] - Parity[node_S]);
					slope_T = (Vol_Matrix[node_S][node_T + 1] * Vol_Matrix[node_S][node_T + 1] - Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]) / (Term[node_T + 1] - Term[node_T]);
					// V_t * dt + V_s * dS + V(S,t)
					result_vol = sqrt(slope_S * (S - Parity[node_S]) + slope_T * (T - Term[node_T]) + Vol_Matrix[node_S][node_T] * Vol_Matrix[node_S][node_T]);
				}
			}
			return result_vol;
		}
		else
		{
			return -1.0;
		}
	}

	// LocalVol Data와 ,S, T를 넣으면 Interpolated Local Vol을 계산해줌. 
	// (그때그때 가까운 node 위치를 기억하여 빠른계산)
	double Calc_Volatility_for_Simulation(double T, double S, long* Prev_NodeT, long* Prev_NodeS)
	{

		long Current_Node_T = 0;
		long Current_Node_S = 0;
		double Slope_T;
		double Slope_S;
		double Result_Vol = 0.0;

		if (T <= Term_Locvol[0])
		{
			Current_Node_T = 0;
			Prev_NodeT[0] = 0;
			if (S <= Parity_Locvol[0])
			{
				Current_Node_S = 0;
				Prev_NodeS[0] = 0;
				Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
				return Result_Vol;
			}
			else if (S >= Parity_Locvol[N_Parity - 1])
			{
				Current_Node_S = N_Parity - 1;
				Prev_NodeS[0] = N_Parity - 1;
				Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
				return Result_Vol;
			}
			else
			{
				if (Parity_Locvol[Prev_NodeS[0]] <= S)
				{
					for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
					{
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
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
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
						{
							Prev_NodeS[0] = Current_Node_S;
							break;
						}
					}
				}
				Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity_Locvol[Prev_NodeS[0] + 1] - Parity_Locvol[Prev_NodeS[0]]);
				// V_s * dS + V(S,t)
				Result_Vol = (Slope_S * (S - Parity_Locvol[Prev_NodeS[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
				return Result_Vol;
			}
		}
		else if (T >= Term_Locvol[N_Term - 1])
		{
			Current_Node_T = N_Term - 1;
			Prev_NodeT[0] = N_Term - 1;
			if (S <= Parity_Locvol[0])
			{
				Current_Node_S = 0;
				Prev_NodeS[0] = 0;
				Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
				return Result_Vol;
			}
			else if (S >= Parity_Locvol[N_Parity - 1])
			{
				Current_Node_S = N_Parity - 1;
				Prev_NodeS[0] = N_Parity - 1;
				Result_Vol = LocalVolMat[Current_Node_S][Current_Node_T];
				return Result_Vol;
			}
			else
			{
				if (Parity_Locvol[Prev_NodeS[0]] <= S)
				{
					for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
					{
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
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
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
						{
							Prev_NodeS[0] = Current_Node_S;
							break;
						}
					}
				}
				Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity_Locvol[Prev_NodeS[0] + 1] - Parity_Locvol[Prev_NodeS[0]]);
				// V_s * dS + V(S,t)
				Result_Vol = (Slope_S * (S - Parity_Locvol[Prev_NodeS[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
				return Result_Vol;
			}
		}
		else
		{
			for (Current_Node_T = Prev_NodeT[0]; Current_Node_T < N_Term - 1; Current_Node_T++)
				if (T >= Term_Locvol[Current_Node_T] && T < Term_Locvol[Current_Node_T + 1])
				{
					Prev_NodeT[0] = Current_Node_T;
					break;
				}

			if (S <= Parity_Locvol[0])                   // T가 minT~ maxT 사이에 있고, S가 minS 보다 작음
			{
				Current_Node_S = 0;
				Prev_NodeS[0] = 0;
				Result_Vol = LocalVolMat[Current_Node_S][Prev_NodeT[0]];
				return Result_Vol;
			}
			else if (S >= Parity_Locvol[N_Parity - 1])
			{
				Current_Node_S = N_Parity - 1;
				Prev_NodeS[0] = N_Parity - 1;
				Result_Vol = LocalVolMat[Current_Node_S][Prev_NodeT[0]];
				return Result_Vol;
			}
			else
			{
				if (Parity_Locvol[Prev_NodeS[0]] <= S)
				{
					for (Current_Node_S = Prev_NodeS[0]; Current_Node_S < N_Parity - 1; Current_Node_S++)
					{
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
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
						if (S >= Parity_Locvol[Current_Node_S] && S < Parity_Locvol[Current_Node_S + 1])
						{
							Prev_NodeS[0] = Current_Node_S;
							break;
						}
					}
				}
				Slope_S = (LocalVolMat[Prev_NodeS[0] + 1][Prev_NodeT[0]] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Parity_Locvol[Prev_NodeS[0] + 1] - Parity_Locvol[Prev_NodeS[0]]);
				Slope_T = (LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0] + 1] - LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]) / (Term_Locvol[Prev_NodeT[0] + 1] - Term_Locvol[Prev_NodeT[0]]);
				// V_t * dt + V_s * dS + V(S,t)
				Result_Vol = (Slope_S * (S - Parity_Locvol[Prev_NodeS[0]]) + Slope_T * (T - Term_Locvol[Prev_NodeT[0]]) + LocalVolMat[Prev_NodeS[0]][Prev_NodeT[0]]);
				return Result_Vol;
			}
		}
		return Result_Vol;
	}

	void set_localvol(curveinfo* rf_curve, curveinfo* div_curve, double MAXVOL, double MINVOL)
	{
		if (dynamicflag != 0)
		{
			localvol_dynamicflag = 1;

			long i;
			long j;
			double A, B;
			double S = 1.0;
			double K;
			double T;
			double Rf = 0.0;
			double Div = 0.0;
			double dT = (Term[N_Term - 1] - Term[0]) / ((double)N_Term - 1.0);
			double dK = (Parity[N_Parity - 1] - Parity[0]) / ((double)N_Parity - 1.0);
			double* RiskFree;
			double* Dividend;
			double** ImvolMat;
			double v, dv_dt, dv_dk, dvdv_dkdk;
			double d_1;
			double LocVar;

			Term_Locvol = (double*)calloc(N_Term, sizeof(double));                    // 리턴용 할당
			for (i = 0; i < N_Term; i++) Term_Locvol[i] = Term[0] + (double)i * dT;

			Parity_Locvol = (double*)calloc(N_Parity, sizeof(double));                // 리턴용
			for (i = 0; i < N_Parity; i++) Parity_Locvol[i] = Parity[0] + (double)i * dK;

			RiskFree = (double*)calloc(N_Term, sizeof(double));                       // 할당 1
			Dividend = (double*)calloc(N_Term, sizeof(double));                       // 할당 2
			for (i = 0; i < N_Term; i++)
			{
				RiskFree[i] = rf_curve->Interpolated_Rate(Term_Locvol[i]);
				Dividend[i] = div_curve->Interpolated_Rate(Term_Locvol[i]);
			}

			ImvolMat = (double**)malloc((N_Parity) * sizeof(double*));                // 할당 3       2차원
			LocalVolMat = (double**)malloc((N_Parity) * sizeof(double*));             // 리턴용 할당

			for (i = 0; i < N_Parity; i++)
			{
				ImvolMat[i] = (double*)calloc((N_Term), sizeof(double));
				LocalVolMat[i] = (double*)calloc((N_Term), sizeof(double));

				for (j = 0; j < N_Term; j++) ImvolMat[i][j] = Calc_Implied_Volatility(Term_Locvol[j], Parity_Locvol[i]);
			}

			for (i = 0; i < N_Parity; i++)
			{
				for (j = 0; j < N_Term; j++)
				{
					Rf = RiskFree[j];
					Div = Dividend[j];
					K = Parity_Locvol[i];
					T = Term_Locvol[j];
					v = ImvolMat[i][j];
					dv_dt = dImvol_over_dT(ImvolMat, N_Parity, N_Term, i, j, dT);
					dv_dk = dImvol_over_dK(ImvolMat, N_Parity, N_Term, i, j, dK);
					dvdv_dkdk = dImvoldImvol_over_dKdK(ImvolMat, N_Parity, N_Term, i, j, dK);
					d_1 = d1(S, K, T, Rf, Div, v);
					A = v * v + 2.0 * v * T * (dv_dt + (Rf - Div) * K * dv_dk);
					B = (1.0 + K * d_1 * dv_dk * sqrt(T)) * (1.0 + K * d_1 * dv_dk * sqrt(T)) + v * K * K * T * (dvdv_dkdk - d_1 * dv_dk * dv_dk * sqrt(T));
					LocVar = A / B;
					if (LocVar > 0)
						LocalVolMat[i][j] = max(min(sqrt(LocVar), MAXVOL), MINVOL);
					else
						LocalVolMat[i][j] = -9999999.99;
				}
			}

			fillna_Interpolate(LocalVolMat, N_Parity, N_Term);

			if (RiskFree) free(RiskFree);        // 할당 1
			if (Dividend) free(Dividend);        // 할당 2
			for (i = 0; i < N_Parity; i++)
			{
				if (ImvolMat[i]) free(ImvolMat[i]);     // 할당 3      2차원
			}
			if (ImvolMat) free(ImvolMat);
		}
		else
		{
			Term_Locvol = Term;
			Parity_Locvol = Parity;
			LocalVolMat = Vol_Matrix;
		}

	}

	void set_localvol(curveinfo* rf_curve, curveinfo* div_curve)
	{
		set_localvol(rf_curve, div_curve, 3.0, 0.00001);
	}

	void delete_localvol()
	{
		long i;
		if (localvol_dynamicflag == 1)
		{
			free(Term_Locvol);
			free(Parity_Locvol);
			for (i = 0; i < N_Parity; i++)
			{
				free(LocalVolMat[i]);
			}
			free(LocalVolMat);
			localvol_dynamicflag = 0;
		}
	}
};

double dw_over_dt(double** w, long n_parity, long n_term, long p, long q, double dt)
{
	double w_t;
	if (q == 0)
	{
		w_t = (w[p][q + 1] - w[p][q]) / dt;   // Forward Diff
	}
	else if (q == n_term - 1)
	{
		w_t = (w[p][q] - w[p][q - 1]) / dt;   // Backward Diff
	}
	else
	{
		w_t = (w[p][q + 1] - w[p][q - 1]) / (2.0 * dt);  // Central Diff
	}
	return w_t;
}

double dw_over_dy(double** w, long n_parity, long n_term, long p, long q, double* k, double dk)
{
	double w_y;
	if (p == 0)
	{
		w_y = (w[p + 1][q] - w[p][q]) / (dk)*k[p];   // Forward Diff
	}
	else if (p == n_parity - 1)
	{
		w_y = (w[p][q] - w[p - 1][q]) / (dk)*k[p];  // Backward Diff
	}
	else
	{
		w_y = (w[p + 1][q] - w[p - 1][q]) / (2.0 * dk) * k[p]; // Central Diff
	}
	return w_y;
}

double dwdw_over_dydy(double** w, long n_parity, long n_term, long p, long q, double* k, double dk)
{
	double w_yy;
	if (p == 0)
	{
		w_yy = (w[p + 2][q] + w[p][q] - 2.0 * w[p + 1][q]) / (dk * dk) * k[p + 1] * k[p + 1]; // 다음 노드의 감마를 대체로 사용
	}
	else if (p == n_parity - 1)
	{
		w_yy = (w[p][q] + w[p - 2][q] - 2.0 * w[p - 1][q]) / (dk * dk) * k[p - 1] * k[p - 1]; // 이전 노드의 감마를 대체로 사용
	}
	else
	{
		w_yy = (w[p + 1][q] + w[p - 1][q] - 2.0 * w[p][q]) / (dk * dk) * k[p] * k[p];
	}
	return w_yy;
}

// Matrix 안의 NaN값을 Interpolate한다.
// [0.5,   -999.99,     -999.99,  0.52]      [0.5,    0.51,       0.515,    0.52]
// [0.43,  -999.99,     0.43,     0.43]   >> [0.43,   0.39,       0.43,     0.43]
// [0.20,  -999.99,     0.22,     0.21]      [0.20,   0.27,       0.22,     0.21]
// [0.15,  -999.99,     0.15,     0.16]      [0.15,   0.15,       0.15,     0.16]
void fillna_Interpolate(double** Matrix, long N_Index, long N_Column)
{
	long i, j;
	long idx_i1;
	long idx_i2;
	long idx_j1;
	long idx_j2;
	double value;
	double value0;
	double value00;
	double value01;
	double value1;
	double value10;
	double value11;
	for (i = 0; i < N_Index; i++)
		for (j = 0; j < N_Column; j++)
			if (Matrix[i][j] < 0.0)
			{
				value0 = 0.0;
				value00 = 0.0;
				value01 = 0.0;
				value1 = 0.0;
				value10 = 0.0;
				value11 = 0.0;
				idx_i1 = i;
				idx_i2 = i;
				idx_j1 = j;
				idx_j2 = j;
				for (idx_i1 = i; idx_i1 >= 0; idx_i1--)
				{
					if (Matrix[idx_i1][j] > 0.0)
					{
						value0 = Matrix[idx_i1][j];
						break;
					}
				}

				if (idx_i1 < 0)
				{
					idx_i1 = 0;
					for (idx_j1 = j; idx_j1 >= 0; idx_j1--)
					{
						if (Matrix[idx_i1][idx_j1] > 0.0)
						{
							value00 = Matrix[idx_i1][idx_j1];
							break;
						}
					}

					for (idx_j2 = j; idx_j2 < N_Column; idx_j2++)
					{
						if (Matrix[idx_i1][idx_j2] > 0.0)
						{
							value01 = Matrix[idx_i1][idx_j2];
							break;
						}
					}
					value0 = 0.5 * value00 + 0.5 * value01;
				}

				for (idx_i2 = i; idx_i2 < N_Index; idx_i2++)
				{
					if (Matrix[idx_i2][j] > 0.0)
					{
						value1 = Matrix[idx_i2][j];
						break;
					}
				}

				if (idx_i2 == N_Index)
				{
					idx_i2 = N_Index - 1;
					for (idx_j1 = j; idx_j1 >= 0; idx_j1--)
					{
						if (Matrix[idx_i2][idx_j1] > 0.0)
						{
							value10 = Matrix[idx_i2][idx_j1];
							break;
						}
					}

					for (idx_j2 = j; idx_j2 < N_Column; idx_j2++)
					{
						if (Matrix[idx_i2][idx_j2] > 0.0)
						{
							value11 = Matrix[idx_i2][idx_j2];
							break;
						}
					}

					value1 = value10 * 0.5 + value11 * 0.5;
				}
				if (idx_i2 != idx_i1)
				{
					value = (value1 - value0) / ((double)idx_i2 - (double)idx_i1) * ((double)i - (double)idx_i1) + value0;
				}
				else
				{
					value = -99999.99;
				}
				Matrix[i][j] = value;

			}
}

double d1(double S, double K, double T, double r, double div, double Vol)
{
	double d1;
	d1 = (log(S / K) + (r - div + 0.5 * Vol * Vol) * T) / (Vol * sqrt(T));
	return d1;
}

double dImvol_over_dK(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dK)
{
	double dImvol;
	if (idx_parity == 0)
	{
		dImvol = ImVol[1][idx_term] - ImVol[0][idx_term];
		return dImvol / dK;
	}
	else if (idx_parity == NParity - 1)
	{
		dImvol = ImVol[NParity - 1][idx_term] - ImVol[NParity - 2][idx_term];
		return dImvol / dK;
	}
	else
	{
		dImvol = ImVol[idx_parity + 1][idx_term] - ImVol[idx_parity - 1][idx_term];
		return dImvol / (2.0 * dK);
	}
}

double dImvoldImvol_over_dKdK(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dK)
{
	double dImvoldImvol;
	if (idx_parity == 0)
	{
		dImvoldImvol = ImVol[2][idx_term] + ImVol[0][idx_term] - 2.0 * ImVol[1][idx_term];
		return dImvoldImvol / (dK * dK);
	}
	else if (idx_parity == NParity - 1)
	{
		dImvoldImvol = ImVol[NParity - 1][idx_term] + ImVol[NParity - 3][idx_term] - 2.0 * ImVol[NParity - 2][idx_term];
		return dImvoldImvol / (dK * dK);
	}
	else
	{
		dImvoldImvol = ImVol[idx_parity + 1][idx_term] + ImVol[idx_parity - 1][idx_term] - 2.0 * ImVol[idx_parity][idx_term];
		return dImvoldImvol / (dK * dK);
	}
}

double dImvol_over_dT(double** ImVol, long NParity, long NTerm, long idx_parity, long idx_term, double dT)
{
	double dImvol;
	if (idx_term == 0)
	{
		dImvol = ImVol[idx_parity][1] - ImVol[idx_parity][0];
		return dImvol / dT;
	}
	else if (idx_term == NTerm - 1)
	{
		dImvol = ImVol[idx_parity][NTerm - 1] - ImVol[idx_parity][NTerm - 2];
		return dImvol / dT;
	}
	else
	{
		dImvol = ImVol[idx_parity][idx_term + 1] - ImVol[idx_parity][idx_term - 1];
		return dImvol / (2.0 * dT);
	}
}