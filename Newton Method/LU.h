/* TODO:
 * доделать СЛАУ
 */

#ifndef __LU_CPP__
#define __LU_CPP__

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include "armadillo"

using namespace std;
using namespace arma; // пространство имён библиотеки armadillo

// старая функция поиска строк с ненулевым элементом (когда ведущий элемент оказался нулевым)
/*unsigned int swappingRowsWhenZero(fmat& A, unsigned int first_row, unsigned int first_col)
{
    unsigned int i = first_row + 1;
    unsigned int n = A.n_rows;
    while((i < n) && (A(i,first_col) == 0))
        ++i;
    return i;
}*/

// найти максимальный элемент в столбце
unsigned int findMaxInCol(mat& A, unsigned int first_row, unsigned int first_col)
{
    unsigned int n = A.n_rows;
    unsigned int maxi = first_row;
    for(unsigned int i = first_row + 1; i < n; ++i)
        if(abs(A(i, first_col)) > abs(A(maxi, first_col)))
            maxi = i;
    return maxi;
}

// p-норма с p = бесконечность
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

unsigned int matrixRank(const mat& A)
{
	unsigned int n = A.n_rows;
	unsigned int m = A.n_cols;
	mat T = A;
	unsigned int i = 0, c = 0;
	/*for (; c < n; ++i, ++c)
	{
		if (i == m)
			return n;
		unsigned int maxIndex = findMaxInCol(T, i, i);
		if (T(maxIndex, i) == 0.0)
			break;
		if (maxIndex != i)
			T.swap_rows(i, maxIndex);

		double refelem = T(i, i);
		for (unsigned int j = i + 1; j < n; ++j)
		{
			double multiplier = T(j, i) / refelem;
			for (unsigned int k = i; k < m; ++k)
				T(j, k) -= T(i, k) * multiplier;
		}
		cout << T << endl;
	}*/
	for (; (i < n) && (c < m); ++i, ++c)
	{
		if (i == n)
			return n;

		unsigned int maxIndex = findMaxInCol(T, i, c);
		if (T(maxIndex, c) == 0.0)
		{
			while ((c < (m - 1)) && (T(maxIndex, c) == 0.0))
			{
				++c;
				maxIndex = findMaxInCol(T, i, c);
			}
			if ((c == m - 1) && (T(maxIndex, c) == 0.0))
				return i;
		}
		if (maxIndex != i)
			T.swap_rows(i, maxIndex);

		double refelem = T(i, c);

		for (unsigned int j = i + 1; j < n; ++j)
		{
			double multiplier = T(j, c) / refelem;
			for (unsigned int k = c; k < m; ++k)
				T(j, k) -= T(i, k) * multiplier;
		}
	}
	/*for (unsigned int j = 0; j < m; ++j)
		if (T(i,j) != 0)
			return n;
	return n - 1;*/
	return i;
}

bool LUP(const mat& A, mat& LU, vector<unsigned int>& P, unsigned int& iterations, unsigned int& operations)
{
    unsigned int n = A.n_rows;

    if(P.empty())
        for(unsigned int i = 0; i < n; ++i)
            P.push_back(i);

    LU = A;
	unsigned int i = 0, c = 0;
    for(; c < n; ++i, ++c)
    {
		++iterations;
        unsigned int maxIndex = findMaxInCol(LU, i, c);
		if (LU(maxIndex, c) == 0.0)
		{
			if (c == (n - 1))
				break;
			++c;
			maxIndex = findMaxInCol(LU, i, c);
			while ((c < n) && (LU(maxIndex, c) == 0.0))
			{
				++c;
				++iterations;
				maxIndex = findMaxInCol(LU, i, c);
			}
			if (c == n)
				break;
		}
        if(maxIndex != i)
        {
            LU.swap_rows(i, maxIndex);
			++operations;
            swap(P[i], P[maxIndex]);
        }

        double refelem = LU(i,c);

        for(unsigned int j = i + 1; j < n; ++j)
        {
			++iterations;
			++operations;
            double multiplier = LU(j,c) / refelem;
			for (unsigned int k = c; k < n; ++k)
			{
				LU(j, k) -= LU(i, k) * multiplier;
				++operations;
				++iterations;
			}
            LU(j,i) = multiplier;
			++operations;
        }
    }
	if (i < n - 1)
	{
		for (unsigned int j = i; j < n - 1; ++j)
		{
			++iterations;
			for (unsigned int k = j + 1; k < n; ++k)
			{
				LU(j, k) = 0.0;
				++iterations;
			}
		}
	}
    return true;
}


vec SLAU(const mat& LU, const vec& b, unsigned int& iterations, unsigned int& operations)
{
	unsigned int n = LU.n_rows;
	vec y(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		++iterations;
		double minus = 0;
		for (unsigned int j = 0; j < i; ++j)
		{
			minus += y[j] * LU(i, j);
			++iterations;
			++operations;
		}
		y[i] = b(i) - minus;
		++operations;
	}
	//cout << "y:\n" << y << "\n\n";
	vec x(n);
	vector<bool> x_filled(n, false);

	int k = n - 1;
	bool found = false;
	for (; k >= 0; --k)
	{
		for (int i = n - 1; i >= k; --i)
		//for (int i = k; i <= n - 1; ++i)
		{
			if (LU(i, k) != 0.0)
			{
				found = true; break;
			}
		}
		if (found)
			break;
	}
	if (k < 0)
	{
		x = zeros<vec>(n);
		return x;
	}
	//int c = k;
	for (int i = k; i >= 0; --i)
	{
		++iterations;
		int c = i;
		while (LU(i, c) == 0.0)
		{
			++c;
			++iterations;
		}
		for (int j = c + 1; j < n; ++j)
		{
			++iterations;
			if (x_filled[j] == false)
			{
				x[j] = 1;
				x_filled[j] = true;
			}
		}
		double minus = 0;
		for (int j = n - 1; j > c; --j)
		{
			minus += x[j] * LU(i, j);
			++iterations;
			++operations;
		}
		x[c] = (y(i) - minus) / LU(i, c);
		++operations;
		x_filled[c] = true;
	}

	return x;
}


int parityOfPermutation(const vector<unsigned int>& v)
{
    unsigned int n = v.size();
    vector<bool> checked(n, false);
    bool even = true;
    for(unsigned int i = 0; i < n; ++i)
    {
        if(!checked[i])
        {
            checked[i] = true;
            unsigned int j = v[i];
            unsigned int cycleLength = 1;
            while(i != j)
            {
                checked[j] = true;
                j = v[j];
                ++cycleLength;
            }
            if(cycleLength % 2 == 0)
                even = !even;
        }
    }
    if(even)
        return 1;
    else
        return -1;
}

bool inverseMatrix(const mat& A, mat& inv, unsigned int& iterations, unsigned int& operations)
{
	unsigned int n = A.n_rows;
	inv.zeros(0, 0);
	mat LU(n, n);
	vector<unsigned int> Pvec;
	if (LUP(A, LU, Pvec, iterations, operations) == false)
		return false;

	for (unsigned int i = 0; i < n; ++i)
	{
		vec b = zeros<vec>(n);
		b(distance(Pvec.begin(), find(Pvec.begin(), Pvec.end(), i))) = 1;
		inv.insert_cols(i, SLAU(LU, b, iterations, operations));
		++iterations;
	}
	return true;
}


void showRank()
{
	mat A;
	bool status = A.load("matrix.txt", raw_ascii);
	if (status == false)
	{
		cout << "Error! Can not open the file matrix.txt or file is in incorrect format.\n";
		return;
	}
	cout << "A =\n" << A << "\nRankA = " << matrixRank(A) << "\nand by armadillo: rank(A) = " << arma::rank(A) << endl;
}

#endif