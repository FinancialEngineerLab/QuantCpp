#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

#include "CalcDate.h"

#ifndef UTILITY
#include "Util.h"
#endif

#ifndef STRUCTURE
#include "Structure.h"
#endif

#ifndef NULL
#define NULL 0
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

#include <crtdbg.h>

using namespace std;

double BSOption(double S, double K, double T, double sigma, double r, double div, long optType)
{
	double d1 = (log(S / K) + (r - div + 0.5 * sigma * sigma) * T)/(sigma * sqrt(T));
	double d2 = d1 - sigma * sqrt(T);
	if (optType == 0) return S * exp(-div * T) * CDF_N(d1) - K * exp(-r * T) * CDF_N(d2);
	else return K * exp(-r * T) * CDF_N(-d2) * S * exp(-div * T) * CDF_N(-d1);

}

double trapezoidalMethod(vector<double> x, vector<double> y)
{
	long i;
	long n = x.size();
	double answer = 0.0;
	double add = 0.0;
	for (i = 1; i < n; i++)
	{
		add = 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
		answer += add;
		if (i >= 500 && add < answer * 0.00001)
			break;
	}
	return answer;
}

double trapezoidalMethod(double* x, double* y, long n)
{
	long i;
	double answer = 0.0;
	double add = 0.0;
	for (i = 1; i < n; i++)
	{
		add = 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
		answer += add;
		if (i >= 500 && add < answer * 0.00001)
			break;
	}
	return answer;
}

//This function computes the value of a European option using the
//Heston model for stochastic volatility.
//
//S: Spot price
//K: Strike price
//r: Risk-free interest rate
//div: Dividend yield
//init_variance: Initial variance
//tau: Time to maturity (years)
//LR_variance: Long-run volatility
//kappa: Mean-reversion rate for volatility
//volofvol: Volatility of volatility
//rho: Price-volatility correlation
//gamma: Risk-aversion parameter
//optType: Option type (call or put)



double Heston_integrand(complex<double> K, double X, double init_variance,
	double tau, double LR_variance, double kappa, double volofvol,
	double rho, double gamma)
{
	complex<double> thetaadj;
	double omega = kappa * LR_variance;
	double ksi = volofvol;
	double theta = kappa;
	complex<double> t((ksi * ksi) * tau / 2.0, 0.0);
	complex<double> a((2.0 * omega) / (ksi * ksi), 0.0);
	if (gamma == 1.0) thetaadj = complex<double>(theta, 0.0);
	else thetaadj = complex<double>((1.0 - gamma) * rho * ksi + sqrt(theta * theta - gamma * (1.0 - gamma) * ksi * ksi), 0.0);
	complex<double> im(0.0, 1.0);
	complex<double> re(1.0, 0.0);
	complex<double> b = (2.0 / (ksi * ksi)) * (thetaadj);
	complex<double> c = (K * K - im * K) / (ksi * ksi);
	complex<double> d = sqrt(b * b + 4.0 * c);
	complex<double> g = (b + d) / 2.0;
	complex<double> h = (b + d) / (b - d);
	complex<double> f1 = a * (t * g - log((1.0 - h * exp(t * d)) / (1.0 - h)));
	complex<double> f2 = g * ((1.0 - exp(t * d)) / (1.0 - h * exp(t * d)));
	complex<double> H = exp(f1 + f2 * init_variance);
	//function to be integrated
	complex<double> integrand = exp(-im * K * X) * (H / (K * K - im * K));
	return real(integrand);
}

double HestonPrice(double S, double K, double r, double div,
	double init_variance, double tau,
	double LR_variance, double kappa, double volofvol, double rho,
	double gamma, long optType,
	long kmax, double* int_x, double* int_y)
{
	double ki = 0.5;
	double price;
	complex<double> pass_phi;
	double omega = kappa * LR_variance;
	double ksi = volofvol;
	double theta = kappa;
	long n = kmax * 5;
	//int kmax = 1000;
	//vector<double> int_x(kmax * 5);
	//vector<double> int_y(kmax * 5);
	double X = log(S / K) + (r - div) * tau;
	long count = 0;
	//setting up the numerical integration
	for (double phi = 0.000001; phi < kmax; phi += 0.2)
	{
		int_x[count] = phi;
		pass_phi = complex<double>(phi, ki);
		int_y[count] = Heston_integrand(pass_phi, X, init_variance, tau, LR_variance, kappa, volofvol, rho, gamma);
		count += 1;
	}
	//computing the price
	double callPrice = S * exp(-div * tau) - (1.0 / (PI)) * K * exp(-r * tau) * trapezoidalMethod(int_x, int_y, n);
	double putPrice = callPrice + K * exp(-r * tau) - S * exp(-div * tau);
	if (optType == 0) price = callPrice;
	else price = putPrice;
	return price;
}

int main()
{
	long i;
	long optType = 0;
	double S, K, T, r, div;
	S = 1.0;
	K = 1.0;
	T = 0.25;
	r = 0.02;
	div = 0.01;
	double init_variance, LR_variance, kappa, volofvol, rho, gamma = 0, init_vol;
	kappa = 1.0;
	init_vol = 0.2;
	init_variance = init_vol * init_vol;
	LR_variance = 0.2 * 0.2;
	volofvol = 0.01;
	rho = -0.1;

	const long k_max = 1000 ;
	double int_x[k_max * 5] = { 0.0, };
	double int_y[k_max * 5] = { 0.0, };

	double C = HestonPrice(S, K, r, div, init_variance, T, LR_variance, kappa, volofvol, rho, gamma, optType, k_max, int_x, int_y);
	double C2 = BSOption(S, K, T, init_vol, r, div, optType);
	
}