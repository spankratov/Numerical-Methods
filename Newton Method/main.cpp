#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>

#include "armadillo"
#include "functions.h"
#include "LU.h"

using namespace std;
using namespace arma;

mat calculateJacobyMat(const vector< vector<funn> > &J, const vec& x)
{
	unsigned int n = J.size();
	mat result(n, n);
	for (unsigned int i = 0; i < n; ++i)
		for (unsigned int j = 0; j < n; ++j)
			result(i, j) = J[i][j](x);
	return result;
}

vec calculateF(const vector<funn>& F, const vec& x)
{
	unsigned int n = F.size();
	vec result(n);
	for (unsigned int i = 0; i < n; ++i)
		result(i) = F[i](x);
	return result;
}

vec ModifiedNewtonMethod(const vector<funn>& F, const vector< vector<funn> > &J, const vec& x0, double e, unsigned int& iterations, unsigned int& operations, unsigned int& steps)
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	vec x = x0, delta;
	mat A = calculateJacobyMat(J, x);
	mat LU;
	vector<unsigned int> P;
	LUP(A, LU, P, iterations, operations);
	do
	{
		//delta = solve(calculateJacobyMat(J, x), -calculateF(F, x));
		
		vec b = -calculateF(F, x);
		unsigned int n = b.n_elem;
		vec bnext(n);
		for (unsigned int i = 0; i < n; ++i)
			bnext(i) = b(P[i]);
		delta = SLAU(LU, bnext, iterations, operations);
		x += delta;
		++steps;
		//cout << x << "\n\n";
	} while (norm(delta, "inf") > e);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Modified Newton Method. Iterations: " << steps << " Operations: " << operations << " Time: " << elapsed_seconds.count() << operations << endl;
	return x;
}

vec NewtonMethod(const vector<funn>& F, const vector< vector<funn> > &J, const vec& x0, double e, unsigned int& iterations, unsigned int& operations, unsigned int& steps, int k = -1)
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	
	unsigned int i = 0;
	vec xnext = x0, x;
	do
	{
		x = xnext;
		if (i >= k)
			return ModifiedNewtonMethod(F, J, x, e, iterations, operations, steps);
		mat inv;
		inverseMatrix(calculateJacobyMat(J, x), inv, iterations, operations);
		xnext = x - inv * calculateF(F, x);
		++i;
		++steps;
	} while (norm(xnext - x, "inf") > e);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Newton Method. Iterations: " << steps << " Operations: " << operations << " Time: " << elapsed_seconds.count() << endl;
	return xnext;
}

vec HybridNewtonMethod(const vector<funn>& F, const vector< vector<funn> > &J, const vec& x0, double e, unsigned int& iterations, unsigned int& operations, unsigned int& steps, int k = 1)
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	unsigned int i = 0;
	vec xnext = x0, x;
	mat inv;
	//cout << xnext << "\n\n\n";
	do
	{
		x = xnext;
		if (i % k == 0)
			inverseMatrix(calculateJacobyMat(J, x), inv, iterations, operations);
		xnext = x - inv * calculateF(F, x);
		++i;
		++steps;
		//xnext = x - calculateJacobyMat(J, x).i() * calculateF(F, x);
		//cout << xnext << "\n\n\n";
	} while (norm(xnext - x, "inf") > e);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Newton Method. Iterations: " << steps << " Operations: " << operations << " Time: " << elapsed_seconds.count() << endl;
	return xnext;
}

void firstTask()
{
	const double e = 1e-12;
	unsigned int iterations = 0, operations = 0, steps = 0;
	vec x0;
	x0 << 0 << 0;
	vec solution;
	cout << "Newton method:\n";
	solution = NewtonMethod(firstTaskVec, firstTaskJacoby, x0, e, iterations, operations, steps);
	cout << "First Task solution:\n" << solution;
	cout << "FirstTask1(solution) = " << firstTask1(solution) << endl;
	cout << "FirstTask2(solution) = " << firstTask2(solution) << endl;
	cout << "\n\nModified Newton method:\n";
	solution = ModifiedNewtonMethod(firstTaskVec, firstTaskJacoby, x0, e, iterations, operations, steps);
	cout << "First Task solution:\n" << solution;
	cout << "FirstTask1(solution) = " << firstTask1(solution) << endl;
	cout << "FirstTask2(solution) = " << firstTask2(solution) << endl;
}

void secondTask(double x5)
{
	const double e = 1e-12;
	unsigned int iterations = 0, operations = 0, steps = 0;
	vec x0;
	x0 << 0.5 << 0.5 << 1.5 << (-1.0) << x5 << 1.5 << 0.5 << (-0.5) << 1.5 << (-1.5);
	vec solution;
	cout << "Second Task with x5 = " << x5 << ":\n\n";
	cout << "a)\n";
	solution = NewtonMethod(secondTaskVec, secondTaskJacoby, x0, e, iterations, operations, steps);
	cout << "Second Task solution:\n" << solution;
	cout << "SecondTask1(solution) = " << secondTask1(solution) << endl;
	cout << "SecondTask2(solution) = " << secondTask2(solution) << endl;
	cout << "SecondTask3(solution) = " << secondTask3(solution) << endl;
	cout << "SecondTask4(solution) = " << secondTask4(solution) << endl;
	cout << "SecondTask5(solution) = " << secondTask5(solution) << endl;
	cout << "SecondTask6(solution) = " << secondTask6(solution) << endl;
	cout << "SecondTask7(solution) = " << secondTask7(solution) << endl;
	cout << "SecondTask8(solution) = " << secondTask8(solution) << endl;
	cout << "SecondTask9(solution) = " << secondTask9(solution) << endl;
	cout << "SecondTask10(solution) = " << secondTask10(solution) << endl;

	cout << "\n\nb)\n";
	iterations = 0; operations = 0; steps = 0;
	solution = ModifiedNewtonMethod(secondTaskVec, secondTaskJacoby, x0, e, iterations, operations, steps);
	cout << "Second Task solution:\n" << solution;
	cout << "SecondTask1(solution) = " << secondTask1(solution) << endl;
	cout << "SecondTask2(solution) = " << secondTask2(solution) << endl;
	cout << "SecondTask3(solution) = " << secondTask3(solution) << endl;
	cout << "SecondTask4(solution) = " << secondTask4(solution) << endl;
	cout << "SecondTask5(solution) = " << secondTask5(solution) << endl;
	cout << "SecondTask6(solution) = " << secondTask6(solution) << endl;
	cout << "SecondTask7(solution) = " << secondTask7(solution) << endl;
	cout << "SecondTask8(solution) = " << secondTask8(solution) << endl;
	cout << "SecondTask9(solution) = " << secondTask9(solution) << endl;
	cout << "SecondTask10(solution) = " << secondTask10(solution) << endl;

	cout << "\n\nc)From Newton method to modified Newton method with 3 iterations\n";
	iterations = 0; operations = 0; steps = 0;
	solution = NewtonMethod(secondTaskVec, secondTaskJacoby, x0, e, iterations, operations, steps, 7);
	cout << "Second Task solution:\n" << solution;
	cout << "SecondTask1(solution) = " << secondTask1(solution) << endl;
	cout << "SecondTask2(solution) = " << secondTask2(solution) << endl;
	cout << "SecondTask3(solution) = " << secondTask3(solution) << endl;
	cout << "SecondTask4(solution) = " << secondTask4(solution) << endl;
	cout << "SecondTask5(solution) = " << secondTask5(solution) << endl;
	cout << "SecondTask6(solution) = " << secondTask6(solution) << endl;
	cout << "SecondTask7(solution) = " << secondTask7(solution) << endl;
	cout << "SecondTask8(solution) = " << secondTask8(solution) << endl;
	cout << "SecondTask9(solution) = " << secondTask9(solution) << endl;
	cout << "SecondTask10(solution) = " << secondTask10(solution) << endl;
	
	cout << "\n\nd)Hybrid Newton method with calculating Jacoby matrix each 3 times\n";
	iterations = 0; operations = 0; steps = 0;
	solution = HybridNewtonMethod(secondTaskVec, secondTaskJacoby, x0, e, iterations, operations, steps, 4);
	cout << "Second Task solution:\n" << solution;
	cout << "SecondTask1(solution) = " << secondTask1(solution) << endl;
	cout << "SecondTask2(solution) = " << secondTask2(solution) << endl;
	cout << "SecondTask3(solution) = " << secondTask3(solution) << endl;
	cout << "SecondTask4(solution) = " << secondTask4(solution) << endl;
	cout << "SecondTask5(solution) = " << secondTask5(solution) << endl;
	cout << "SecondTask6(solution) = " << secondTask6(solution) << endl;
	cout << "SecondTask7(solution) = " << secondTask7(solution) << endl;
	cout << "SecondTask8(solution) = " << secondTask8(solution) << endl;
	cout << "SecondTask9(solution) = " << secondTask9(solution) << endl;
	cout << "SecondTask10(solution) = " << secondTask10(solution) << endl;
	cout << endl << endl;
}

void secondTaskWrapper()
{
	secondTask(-0.5);
	secondTask(-0.2);
}

int main()
{
	//firstTask();
	secondTaskWrapper();
	system("pause");
	return 0;
}
