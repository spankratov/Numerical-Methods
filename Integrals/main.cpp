#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include "armadillo"

using namespace std;
using namespace arma;

const double SOLVED_INTERGAL = 18.6029; // by wolfram alpha
const double accuracy = 1e-12;

const double a_first = 1.1, b_first = 2.5, alpha = 0.4, beta = 0;
const double M0(double a, double b) { return 5.0 / 3.0 * (pow(b - a_first, 3.0 / 5.0) - pow(a - a_first, 3.0 / 5.0)); } // not true when beta != 0
const double M1(double a, double b) { return 5.0 / 8.0 * pow(b - a_first, 8.0 / 5.0) + 11.0 / 6.0 * pow(b - a_first, 3.0 / 5.0) - 5.0 / 8.0 * pow(a - a_first, 8.0 / 5.0) - 11.0 / 6.0 * pow(a - a_first, 3.0 / 5.0); } // not true when beta != 0
const double M2(double a, double b) { return 5.0 / 13.0 * pow(b - a_first, 13.0 / 5.0) + 11.0 / 8.0 * pow(b - a_first, 8.0 / 5.0) + 121.0 / 60.0*pow(b - a_first, 3.0 / 5.0) - 5.0 / 13.0 * pow(a - a_first, 13.0 / 5.0) - 11.0 / 8.0 * pow(a - a_first, 8.0 / 5.0) - 121.0 / 60.0*pow(a - a_first, 3.0 / 5.0); } // not true when beta != 0
const double M3(double a, double b) { return (pow(5 * b - 11.0 / 2.0, 3.0 / 5.0)*(416 * b*b*b + 528 * b*b + 726 * b + 1331) - pow(5 * a - 11.0 / 2.0, 3.0 / 5.0)*(416 * a*a*a + 528 * a*a + 726 * a + 1331))* pow(5.0, 2.0 / 5.0) / 7488.0; } // not true when beta != 0
const double M4(double a, double b) { return (pow(5 * b - 11.0 / 2.0, 3.0 / 5.0) * (3744 * b*b*b*b + 4576 * b*b*b + 5808 * b*b + 7986 * b + 14641) - pow(5 * a - 11.0 / 2.0, 3.0 / 5.0) * (3744 * a*a*a*a + 4576 * a*a*a + 5808 * a*a + 7986 * a + 14641)) * pow(5.0, 2.0 / 5.0) / 86112.0; }
const double M5(double a, double b) { return (pow(5 * b - 11.0 / 2.0, 3.0 / 5.0)*(172224 * pow(b, 5) + 205920 * b*b*b*b + 251680 * b*b*b + 319440 * b*b + 439230 * b + 805255) - pow(5 * a - 11.0 / 2.0, 3.0 / 5.0)*(172224 * pow(a, 5) + 205920 * a*a*a*a + 251680 * a*a*a + 319440 * a*a + 439230 * a + 805255)) * pow(5.0, 2.0 / 5.0) / 4822272.0; }

double p(double x)  { return pow(x - a_first, -alpha) * pow(b_first - x, -beta); }
double f(double x)  { return 0.5 * cos(2 * x) * exp(0.4 * x) + 2.4 * sin(1.5 * x) * exp(-6 * x) + 6 * x; }

vec Cardano(double a, double b, double c)
{
	double Q = (a*a - 3.0 * b) / 9.0, R = (2.0 * a*a*a - 9.0 * a*b + 27.0 * c) / 54.0;
	double S = Q*Q*Q - R*R;
	if (S > 0)
	{
		vec result;
		double F = 1.0 / 3.0 * acos(R / pow(Q, 1.5));
		result << -2.0 * sqrt(Q) * cos(F) - a / 3.0 << -2.0 * sqrt(Q) * cos(F + 2.0/3.0 * M_PI) - a / 3.0 << -2.0 * sqrt(Q) * cos(F - 2.0 / 3.0 * M_PI) - a / 3.0;
		return result;
	}
	else if (S == 0)
	{
		vec result;
		result << -2.0*pow(R, 1.0 / 3.0) - a / 3.0;
		result << pow(R, 1.0 / 3.0) - a / 3.0;
		result << pow(R, 1.0 / 3.0) - a / 3.0;
		return result;
	}
	else
		throw "Can not solve cubic function";
}

double NewtonCotes(double a, double b)
{
	double x1 = a, x2 = (a + b) / 2, x3 = b;
	mat X;
	X << 1 << 1 << 1 << endr
		<< x1 << x2 << x3 << endr
		<< x1*x1 << x2*x2 << x3*x3 << endr;
	vec M; M << M0(a,b) << M1(a,b) << M2(a,b);
	vec A = solve(X, M);
	double result = A[0] * f(x1) + A[1] * f(x2) + A[2] * f(x3);
	return result;
}

double Gauss(double a, double b)
{
	mat At;
	vec M;
	M << M0(a, b) << M1(a, b) << M2(a, b) << M3(a, b) << M4(a, b) << M5(a, b);
	At << M[0] << M[1] << M[2] << endr
		<< M[1] << M[2] << M[3] << endr
		<< M[2] << M[3] << M[4] << endr;
	vec m; m << -M[3] << -M[4] << -M[5];
	vec ai = solve(At, m);
	vec x = Cardano(ai[2], ai[1], ai[0]);
	mat X;
	X << 1 << 1 << 1 << endr
		<< x[0] << x[1] << x[2] << endr
		<< x[0]*x[0] << x[1]*x[1] << x[2]*x[2] << endr;
	vec Mt; Mt << M[0] << M[1] << M[2];
	vec A = solve(X, Mt);
	double result = A[0] * f(x[0]) + A[1] * f(x[1]) + A[2] * f(x[2]);
	return result;
}

double solveIntegral(bool NC, bool Aitken, double& Aitken_m)
{
	double interval = b_first - a_first;
	unsigned int intervalsNumber;
	double prevprevresult, prevresult = 0, result = 0;
	double error = 0;
	Aitken_m = 0;
	double m = (NC ? 3.0 : 6.0);

	if (Aitken)
	{
		for (intervalsNumber = 1; intervalsNumber != 8; intervalsNumber *= 2)
		{
			double step = interval / intervalsNumber;
			prevprevresult = prevresult; prevresult = result; result = 0;
			double a = a_first, b = a_first + step;
			if (NC)
				for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
					result += NewtonCotes(a, b);
			else
				for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
					result += Gauss(a, b);
			error = (result - prevresult) / (pow(2, m) - 1.0);
			cout << "result: " << result << " error: " << error << endl;
		}
		if (Aitken)
		{
			Aitken_m = -log((result - prevresult) / (prevresult - prevprevresult)) / log(2);
			cout << "Aitken process, m = " << Aitken_m << endl << endl;
		}
	}
	else
		intervalsNumber = 1;

	if ((error > accuracy) || (error == 0))
	do
	{
		double step = interval / intervalsNumber;
		prevprevresult = prevresult; prevresult = result; result = 0;
		double a = a_first, b = a_first + step;
		if (NC)
			for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
				result += NewtonCotes(a, b);
		else
			for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
				result += Gauss(a, b);
		error = (result - prevresult) / (pow(2, m) - 1.0);
		cout << "result: " << result << " error: " << error << endl;
		if (Aitken && error > 0)
		{
			Aitken_m = -log((result - prevresult) / (prevresult - prevprevresult)) / log(2);
			cout << "Aitken process, m = " << Aitken_m << endl << endl;
		}
		intervalsNumber *= 2;
	//} while (abs(result - SOLVED_INTERGAL) > accuracy);
	} while (error > accuracy);
	if (Aitken)
		cout << "Aitken process converges to m = " << Aitken_m << endl;
	return result;
}

double solveIntegralUntilADPHelper(unsigned int firstIntervalsNumber = 1, bool NC = true, bool ownADP = false, double ADP = 3.0)
{
	double interval = b_first - a_first;
	unsigned int intervalsNumber;
	double prevprevresult, prevresult = 0, result = 0;
	double error = 0, Cm = 0;
	if (!ownADP)
		ADP = (NC ? 3.0 : 6.0);

	double step = interval / firstIntervalsNumber;
	prevprevresult = prevresult; prevresult = result; result = 0;
	double a = a_first, b = a_first + step;
	if (NC)
		for (unsigned int i = 0; i < firstIntervalsNumber; ++i, a += step, b += step)
			result += NewtonCotes(a, b);
	else
		for (unsigned int i = 0; i < firstIntervalsNumber; ++i, a += step, b += step)
			result += Gauss(a, b);
	error = (result - prevresult) / (pow(2, ADP) - 1.0);
	cout << "result: " << result << " error: " << error << endl;

	step = interval / (firstIntervalsNumber * 2.0);
	prevprevresult = prevresult; prevresult = result; result = 0;
	a = a_first; b = a_first + step;
	if (NC)
		for (unsigned int i = 0; i < firstIntervalsNumber*2; ++i, a += step, b += step)
			result += NewtonCotes(a, b);
	else
		for (unsigned int i = 0; i < firstIntervalsNumber*2; ++i, a += step, b += step)
			result += Gauss(a, b);
	error = (result - prevresult) / (pow(2, ADP) - 1.0);
	cout << "result: " << result << " error: " << error << endl;

	double H0 = (interval / firstIntervalsNumber) * pow(accuracy * (1.0 - pow(2.0, -ADP)) / abs(result - prevresult), 1.0 / ADP);
	intervalsNumber = interval / H0;
	if (intervalsNumber % 2 == 1)
		++intervalsNumber;
	cout << "Optimal number of intervals: " << intervalsNumber << endl;

	step = (interval / intervalsNumber) * 2.0;
	prevprevresult = prevresult; prevresult = result; result = 0;
	a = a_first; b = a_first + step;
	if (NC)
		for (unsigned int i = 0; i < (intervalsNumber / 2); ++i, a += step, b += step)
			result += NewtonCotes(a, b);
	else
		for (unsigned int i = 0; i < (intervalsNumber / 2); ++i, a += step, b += step)
			result += Gauss(a, b);
	error = (result - prevresult) / (pow(2, ADP) - 1.0);
	cout << "result: " << result << " error: " << error << endl;

	step = interval / intervalsNumber;
	prevprevresult = prevresult; prevresult = result; result = 0;
	a = a_first; b = a_first + step;
	if (NC)
		for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
			result += NewtonCotes(a, b);
	else
		for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
			result += Gauss(a, b);
	error = (result - prevresult) / (pow(2, ADP) - 1.0);
	cout << "result: " << result << " error: " << error << endl;

	/*for (intervalsNumber = firstIntervalsNumber; intervalsNumber != 8 * firstIntervalsNumber; intervalsNumber *= 2)
	{
		double step = interval / intervalsNumber;
		prevprevresult = prevresult; prevresult = result; result = 0;
		double a = a_first, b = a_first + step;
		if (NC)
			for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
				result += NewtonCotes(a, b);
		else
			for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
				result += Gauss(a, b);
		error = (result - prevresult) / (pow(2, intervalsNumber) - 1.0);
		cout << "result: " << result << " error: " << error << endl;
	}
	Cm = (result - prevresult) / (prevresult - prevprevresult);

	if (Cm > powADP)
		do
		{
			double step = interval / intervalsNumber;
			prevprevresult = prevresult; prevresult = result; result = 0;
			double a = a_first, b = a_first + step;
			if (NC)
				for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
					result += NewtonCotes(a, b);
			else
				for (unsigned int i = 0; i < intervalsNumber; ++i, a += step, b += step)
					result += Gauss(a, b);
			error = (result - prevresult) / (pow(2, intervalsNumber) - 1.0);
			cout << "result: " << result << " error: " << error << endl;
			Cm = (result - prevresult) / (prevresult - prevprevresult);
			intervalsNumber *= 2;
		//} while (abs(result - SOLVED_INTERGAL) > accuracy);
		} while (Cm > powADP);*/

	return result;
}

double solveIntegralUntilADP(bool NC = true)
{
	cout << "With step == interval:\n";
	solveIntegralUntilADPHelper(1, NC);
	cout << "With step == interval / 3:\n";
	return solveIntegralUntilADPHelper(3, NC);
}

double solveIntergalUntilAitkenM(bool NC, double Aitken_m)
{
	cout << "With step == interval:\n";
	solveIntegralUntilADPHelper(1, NC, true, Aitken_m);
	cout << "With step == interval / 3:\n";
	solveIntegralUntilADPHelper(3, NC, true, Aitken_m);
	cout << "With step == interval / 5:\n";
	return solveIntegralUntilADPHelper(5, NC, true, Aitken_m);
}

int main()
{
	try
	{
		double Aitken_m_NC, Aitken_m_G;
		cout << "4) Solving integral with Newton-Cotes method:\n"; 
		solveIntegral(true, false, Aitken_m_NC);
		cout << "\n\n4) Solving integral with Gauss method:\n";
		solveIntegral(false, false, Aitken_m_G);
		cout << "\n\n5) Aitken process:\n";
		solveIntegral(true, true, Aitken_m_NC);
		cout << "\n\n6) Solving integral with ADP as rate of convergence:\n";
		solveIntegralUntilADP(true);
		cout << "\n\n7) Solving integral with m from Aitken process as rate of convergence:\n";
		solveIntergalUntilAitkenM(true, Aitken_m_NC);
	}
	catch (const char* error)
	{
		cout << "Program terminated with error:\n" << error << endl;
	}
	system("pause");
	return 0;
}