/* TODO:
 * доделать СЛАУ
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include "armadillo"

using namespace std;
using namespace arma; // пространство имён библиотеки armadillo

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
	return i;
}

bool LUP(const mat& A, mat& LU, vector<unsigned int>& P)
{
    unsigned int n = A.n_rows;

    if(P.empty())
        for(unsigned int i = 0; i < n; ++i)
            P.push_back(i);

    LU = A;
	unsigned int i = 0, c = 0;
    for(; c < n; ++i, ++c)
    {

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
				maxIndex = findMaxInCol(LU, i, c);
			}
			if (c == n)
				break;
		}
        if(maxIndex != i)
        {
            LU.swap_rows(i, maxIndex);
            swap(P[i], P[maxIndex]);
        }

        double refelem = LU(i,c);

        for(unsigned int j = i + 1; j < n; ++j)
        {
            double multiplier = LU(j,c) / refelem;
            for(unsigned int k = c; k < n; ++k)
                LU(j,k) -= LU(i,k) * multiplier;
            LU(j,i) = multiplier;
        }
    }
	if (i < n - 1)
	{
		for (unsigned int j = i; j < n - 1; ++j)
			for (unsigned int k = j + 1; k < n; ++k)
				LU(j, k) = 0.0;
	}
    return true;
}

void showLUP()
{
    mat A;
    bool status = A.load("matrix.txt", raw_ascii);
    if(status == false)
    {
        cout << "Error! Can not open the file matrix.txt or file is in incorrect format.\n";
        return;
    }
    if(A.n_rows != A.n_cols)
    {
        cout << "Error! Matrix is not square.";
        return;
    }
    unsigned int n = A.n_rows;
    mat LU(n,n);
    vector<unsigned int> Pvec;
    if(LUP(A, LU, Pvec) == false)
    {
        cout << "Can not decompose the matrix.\n";
        return;
    }

    mat L = eye<mat>(n,n);
    mat U = LU;
    mat P(n,n,fill::zeros);
    P(0, Pvec[0]) = 1.0;
    for(unsigned int i = 1; i < n; ++i)
    {
        P(i, Pvec[i]) = 1;
        for(unsigned int j = 0; j < i; ++j)
        {
            L(i,j) = LU(i,j);
            U(i,j) = 0;
        }
    }

    cout << "A:\n" << A << "\nL:\n" << L << "\nU:\n" << U
         << "\nP:" << P << "\nL * U - P * A:\n" << L * U - P * A << '\n';
    /*A.raw_print(cout, "A:");
    L.raw_print(cout, "L:");
    U.raw_print(cout, "U:");
    P.raw_print(cout, "A:");
    mat test = L * U - P * A;
    test.raw_print(cout, "L * U - P * A:");*/
}

vec SLAU(const mat& LU, const vec& b)
{
	unsigned int n = LU.n_rows;
	vec y(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		double minus = 0;
		for (unsigned int j = 0; j < i; ++j)
			minus += y[j] * LU(i, j);
		y[i] = b(i) - minus;
	}
	cout << "y:\n" << y << "\n\n";
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
	for (int i = k; i >= 0; --i)
	{
		int c = i;
		while (LU(i, c) == 0.0)
			++c;
		for (int j = c + 1; j < n; ++j)
			if (x_filled[j] == false)
			{
				x[j] = 1;
				x_filled[j] = true;
			}
		double minus = 0;
		for (int j = n - 1; j > c; --j)
			minus += x[j] * LU(i, j);
		x[c] = (y(i) - minus) / LU(i, c);
		x_filled[c] = true;
	}

	return x;
}

void showSLAU()
{
    mat A;
    bool status = A.load("slau.txt", raw_ascii);
    if(status == false)
    {
        cout << "Error! Can not open the file slau.txt or file is in incorrect format.\n";
        return;
    }
    cout << "SLAU:\n" << A << '\n';
    unsigned int n = A.n_rows;
    vec b = A.col(n);
	unsigned int rankAb = matrixRank(A);
    A.shed_col(n);
	unsigned int rankA = matrixRank(A);
	if (rankA != rankAb)
	{
		cout << "System has not solution.\n";
		return;
	}
    mat LU(n,n);
    vector<unsigned int> Pvec;
    if(LUP(A, LU, Pvec) == false)
    {
        cout << "Can not decompose the matrix.\n";
        return;
    }
    vec bnext(n);
    for(unsigned int i = 0; i < n; ++i)
        bnext(i) = b(Pvec[i]);

	vec x = SLAU(LU, bnext);

    cout << "Solution:\n" << x << "\n";

    cout << "Ax - b:\n" << A * x - b << '\n';
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

void showDeterminant()
{
    mat A;
    bool status = A.load("matrix.txt", raw_ascii);
    if(status == false)
    {
        cout << "Error! Can not open the file matrix.txt or file is in incorrect format.\n";
        return;
    }
    if(A.n_rows != A.n_cols)
    {
		cout << "Error! Matrix is not square.";
		return;
    }
    unsigned int n = A.n_rows;
    mat LU(n,n);
    vector<unsigned int> Pvec;
    if(LUP(A, LU, Pvec) == false)
    {
		cout << "A:\n" << A << "\ndetA = 0" << "\nand by armadillo: detA = "
			<< det(A) << endl;
		return;
    }

    double result = 1;

    for(unsigned int i = 0; i < n; ++i)
        result *= LU(i,i);

    cout << "A:\n" << A << "\ndetA = " << parityOfPermutation(Pvec) * result << "\nand by armadillo: detA = "
         << det(A) << endl;
}

bool inverseMatrix(const mat& A, mat& inv)
{
	unsigned int n = A.n_rows;
	mat LU(n, n);
	vector<unsigned int> Pvec;
	if (LUP(A, LU, Pvec) == false)
		return false;

	for (unsigned int i = 0; i < n; ++i)
	{
		vec b = zeros<vec>(n);
		b(distance(Pvec.begin(), find(Pvec.begin(), Pvec.end(), i))) = 1;
		inv.insert_cols(i, SLAU(LU, b));
	}
	return true;
}

void showInverseMatrix()
{
	mat A;
	bool status = A.load("matrix.txt", raw_ascii);
	if (status == false)
	{
		cout << "Error! Can not open the file matrix.txt or file is in incorrect format.\n";
		return;
	}
	if (A.n_rows != A.n_cols)
	{
		cout << "Error! Matrix is not square.";
		return;
	}

	mat inv;
	if (inverseMatrix(A, inv) == false)
	{
		cout << "Can not decompose the matrix.\n";
		return;
	}

	cout << "A:\n" << A << "\nInverse Matrix:\n" << inv << "\nTest:\nA * Inverse:\n" << A * inv
		<< "\ninv * A:\n" << inv * A << endl;
}

double condNumber(const mat& A)
{
	mat Ainv;
	if (inverseMatrix(A, Ainv) == false)
		return 0;
	return matrixNorm(A) * matrixNorm(Ainv);
}

void showCondNumber()
{
	mat A;
	bool status = A.load("matrix.txt", raw_ascii);
	if (status == false)
	{
		cout << "Error! Can not open the file matrix.txt or file is in incorrect format.\n";
		return;
	}
	if (A.n_rows != A.n_cols)
	{
		cout << "Error! Matrix is not square.";
		return;
	}

	double condA = condNumber(A);
	if (condA == 0)
		cout << "Can not decompose the matrix.\n";
	else
		cout << "A:\n" << A << "\ncond(A) = " << condA << endl;
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

int main()
{
    string command;
    do
    {
        cout << "Enter the command (\"help\" for showing existing commands):\n";
        getline(cin, command);
        if(command == "help")
        {
            cout << "Existing commands:\n";
            cout << "help: for help.\n";
            cout << "LUP: for LUP-decomposition for matrix in matrix.txt.\n";
            cout << "SLAU: for solving SLAU.\n";
            cout << "det: for determinant of matrix in matrix.txt";
			cout << "inv: for inverse matrix for matrix in matrix.txt";
			cout << "cond: for condition number of matrix in matrix.txt";
			cout << "rank: for rank of matrix in matrix.txt";
            cout << "end: for exiting the program.\n";
        }
        else if(command == "LUP")
            showLUP();
        else if(command == "SLAU")
            showSLAU();
		else if (command == "det")
			showDeterminant();
		else if (command == "inv")
			showInverseMatrix();
		else if (command == "cond")
			showCondNumber();
		else if (command == "rank")
			showRank();
        else if(command != "end")
            cout << "Incorrect command.\n";
        cout << "\n\n";
    }while(command != "end");
    return 0;
}

