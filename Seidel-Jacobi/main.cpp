#include "armadillo"
#include <time.h>
#include <stdlib.h>
#include <cmath>

using namespace std;
using namespace arma;

const double e = 0.00000001;

double matrixNorm(const mat& A)
{
	double norm = 0;
	for (unsigned int i = 0; i < A.n_rows; ++i)
	{
		double sum = 0;
		for (unsigned int j = 0; j < A.n_cols; ++j)
			sum += abs(A(i, j));
		if (sum > norm)
			norm = sum;
	}
	return norm;
}

mat generateGiagonallyDominantMatrix(unsigned int n)
{
	srand(time(NULL));
	mat result(n, n);
	for (unsigned int i = 0; i < n; ++i)
	{
		double diag = 0.0;
		for (unsigned int j = 0; j < i; ++j)
		{
			double temp = (((rand() % 2) == 0) ? 1 : -1) * (rand() % 11);
			result(i, j) = temp;
			diag += abs(temp);
		}
		for (unsigned int j = i + 1; j < n; ++j)
		{
			double temp = (((rand() % 2) == 0) ? 1 : -1) * (rand() % 11);
			result(i, j) = temp;
			diag += abs(temp);
		}
		result(i, i) = (((rand() % 2) == 0) ? 1 : -1) * (diag + 1 + rand() % 11);
	}
	return result;
}

vec seidel(const mat& A, const vec& b, unsigned int& iterations)
{
	unsigned int n = A.n_cols;
	iterations = 0;
	mat C(n, n);
	vec d(n);
	for (unsigned int i = 0; i < n; ++i)
	{
		for (unsigned int j = 0; j < i; ++j)
			C(i, j) = -A(i, j) / A(i, i);
		C(i, i) = 0;
		for (unsigned int j = i + 1; j < n; ++j)
			C(i, j) = -A(i, j) / A(i, i);
		d(i) = b(i) / A(i, i);
	}
	vec x = d;
	while (matrixNorm(A*x - b) > e)
	{
		for (unsigned int i = 0; i < n; ++i)
		{
			double temp = d(i);
			for (unsigned int j = 0; j < n; ++j)
				temp += C(i, j) * x(j);
			x(i) = temp;
			++iterations;
		}
		//cout << "\n\nx:\n" << x;
		//cout << "\nA*x - b:\n" << A*x - b;
	}
	return x;
}

vec jacobi(const mat& A, const vec& b, unsigned int& iterations)
{
	unsigned int n = A.n_cols;
	iterations = 0;
	mat C(n, n);
	vec d(n);
	for (unsigned int i = 0; i < n; ++i)
	{
		for (unsigned int j = 0; j < i; ++j)
			C(i, j) = -A(i, j) / A(i, i);
		C(i, i) = 0;
		for (unsigned int j = i + 1; j < n; ++j)
			C(i, j) = -A(i, j) / A(i, i);
		d(i) = b(i) / A(i, i);
	}
	vec x = d;
	vec xk(n);
	while (matrixNorm(A*x - b) > e)
	{
		for (unsigned int i = 0; i < n; ++i)
		{
			double temp = d(i);
			for (unsigned int j = 0; j < n; ++j)
				temp += C(i, j) * x(j);
			xk(i) = temp;
			++iterations;
		}
		x = xk;
		//cout << "\n\nx:\n" << x;
		//cout << "\nA*x - b:\n" << A*x - b;
	}
	return x;
}

int main()
{
	const unsigned int n = 5;
	unsigned int iterations1, iterations2;
	mat A = generateGiagonallyDominantMatrix(n);
	cout << "A:" << A << "\n\n";
	vec b = randi<vec>(n, distr_param(-10, 10));
	vec x = seidel(A, b, iterations1);
	cout << "Seidel:\n";
	cout << "\n\nx:\n" << x;
	cout << "\nA*x - b:\n" << A*x - b << "\n";
	cout << "number of iterations: " << iterations1 << "\n\n";
	vec y = jacobi(A, b, iterations2);
	cout << "Jacobi:\n";
	cout << "\n\nx:\n" << y;
	cout << "\nA*x - b:\n" << A*y - b;
	cout << "\nnumber of iterations: " << iterations2 << "\n\n";
	system("pause");
	return 0;
}