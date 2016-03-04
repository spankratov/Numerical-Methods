#include <iostream>
#include <algorithm>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

void QR(const mat& A, mat& Q, mat& R)
{
	unsigned int m = A.n_rows;
	unsigned int n = A.n_cols;
	Q = eye<mat>(n, m);
	R = A;
	for (unsigned int i = 0; i <= min(m - 1, n - 1); ++i)
	{
		for (unsigned int j = i + 1; j < n; ++j)
		{
			mat temp = eye<mat>(n, m);
			double c = R(i, i) / sqrt(R(i, i) * R(i, i) + R(j, i) * R(j, i));
			double s = -R(j, i) / sqrt(R(i, i) * R(i, i) + R(j, i) * R(j, i));
			temp(i, i) = c;
			temp(i, j) = -s;
			temp(j, i) = s;
			temp(j, j) = c;
			R = temp * R;
			Q = temp * Q;
			//cout << temp << "\n\n";
			//cout << R << "\n\n";
		}
	}
	Q = Q.t();
}

int main()
{
	mat A;
	bool status = A.load("slau.txt", raw_ascii);
	if (status == false)
	{
		cout << "Error! Can not open the file slau.txt or file is in incorrect format.\n";
		return 1;
	}
	unsigned int n = A.n_rows;
	cout << "SLAU:\n" << A << "\n\n";
	mat Q, R;
	vec b = A.col(n);
	A.shed_col(n);
	QR(A, Q, R);
	cout << "Q:\n" << Q << "\n\n";
	cout << "R:\n" << R << "\n\n";
	cout << "Q * Q.t" << Q * Q.t() << "n\n";
	cout << "Q*R - A:" << Q * R - A << "\n\n";

	vec bt = Q.t() * b;
	vec x(n);

	for (int i = n - 1; i >= 0; --i)
	{
		double minus = 0;
		for (int j = i + 1; j < n; ++j)
			minus += x[j] * R(i, j);
		x[i] = (bt(i) - minus) / R(i,i);
	}

	cout << "x:" << x << "\n\n";
	//cout << "A:" << A << "\n\n";
	//cout << "b:" << b << "\n\n";
	cout << "Ax - b" << A*x - b << "\n\n";

	system("pause");
	return 0;
}