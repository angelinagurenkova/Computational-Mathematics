#include <iostream>
#include <vector>
#include<cmath>

using namespace std;

void show(vector <vector <double>> A, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			printf("%6f \t", A[i][j]);
		}
		cout << endl;
	}
}


void SWAP_STR(vector <vector <double>>& A, int ii, int jj, int n, int m) {
	if (ii != jj) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if ((i == ii) && (j == jj) && (i != j)) {
					for (int a = 0; a < m; a++) {
						double tmp1 = A[i][a];
						A[i][a] = A[j][a];
						A[j][a] = tmp1;
					}
				}
			}
		}
	}
}

void SWAP_STB(vector <vector <double>>& A, int ii, int jj, int n, int m) {
	if (ii != jj) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if ((i == ii) && (j == jj) && (i != j)) {
					for (int a = 0; a < n; a++) {
						double tmp1 = A[a][i];
						A[a][i] = A[a][j];
						A[a][j] = tmp1;
					}
				}
			}
		}
	}
}

void treugol(vector <vector <double>>& A, vector <vector <double>>& R, int n, int m, int& rang) {
	R = A;
	for (int k = 0; k < n; k++) {
		double max = abs(R[k][k]);
		int ii = k;
		int jj = k;
		for (int i = k; i < n; i++) { 
			for (int j = k; j < m; j++) {
				if (abs(R[i][j]) > max) {
					max = abs(R[i][j]);
					ii = i;
					jj = j;
				}
			}
		}


		if (max > 1e-15) {
			SWAP_STR(R, k, ii, n, m);
			SWAP_STB(R, k, jj, n, m);

			rang++;
			for (int i = k + 1; i < n; i++) {
				double w = R[i][k] / R[k][k];
				for (int j = k; j < m; j++) {
					R[i][j] -= w * R[k][j];
				}
			}
			cout << "obnulili stolbik k=" << k << endl;

			show(R, n, m);
			cout << endl;
		}
		else break;
	}

}



void LU(vector <vector <double>> A, vector <vector <double>>& L,
	vector <vector <double>>& U, int n, vector <vector <double>>& P)
{
	U = A;
	for (int k = 0; k < n; k++) {
		int pi = 0;
		int pj = 0;
		for (int i = k; i < n; i++) {
			if (abs(U[i][k]) > abs(U[k][k])) {
				pi = i;
				pj = k;
			}
			else {
				pi = k;
				pj = k;
			}
		}
		SWAP_STR(U, pi, pj, n, n);
		SWAP_STR(P, pi, pj, n, n);
	}

	for (int i = 0; i < n; i++)
		if (U[i][i] != 0) {
			for (int j = i; j < n; j++)
				L[j][i] = U[j][i] / U[i][i];
		}

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			if (U[i][i] != 0) {
				for (int j = i; j < n; j++)
					L[j][i] = U[j][i] / U[i][i];
			}

		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}
}

void proisv(vector <vector <double>> A, vector <vector <double>> B,
	vector <vector <double>>& R, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}


int main()
{
	int n = 4;
	vector <vector <double>> A(n), R(n), Rasshir_A(n), Treug_Rasshir_a(n);
	for (int i = 0; i < n; i++) {
		int k = (rand() % 6);
		for (int j = 0; j < n; j++)
		{
			Rasshir_A[i].push_back(0);
			Treug_Rasshir_a[i].push_back(0);
			R[i].push_back(0);
			A[i].push_back(rand() % 45);
		}
	}

	for (int i = 0; i < n; i++)
	{
		A[i][2] = A[i][0] + A[i][1];
		A[i][3] = A[i][0] - A[i][1];
	}

	cout << "-------------------- matrix A --------------------" << endl;
	show(A, n, n);
	int rang_A = 0;

	treugol(A, R, n, n, rang_A);

	cout << "-------------------- matrix treugol --------------------" << endl;
	show(R, n, n);

	cout << "-------------------- rang=" << rang_A << "--------------------" << endl;

	vector<double> b;
	Rasshir_A = A;
	for (int i = 0; i < n; i++) {
		b.push_back(A[i][0] + A[i][3]);
		cout << "b=" << b[i] << endl;
		Rasshir_A[i].push_back(b[i]);
		Treug_Rasshir_a[i].push_back(0);
	}
	cout << "-------------------- matrix Rasshir_A --------------------" << endl;
	show(Rasshir_A, n, n + 1);
	int rang_Rasshir_A = 0;
	treugol(Rasshir_A, Treug_Rasshir_a, n, n + 1, rang_Rasshir_A);

	cout << "-------------------- matrix Treug_Rasshir_a --------------------" << endl;
	show(Treug_Rasshir_a, n, n + 1);

	cout << "-------------------- rang=" << rang_Rasshir_A << "--------------------" << endl;

	if (rang_A == rang_Rasshir_A) {
		cout << "Matrix is sovmestnaya" << endl;


		vector <vector <double>> L(n), U(n), Rez(n), P(n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				L[i].push_back(0);
				U[i].push_back(0);
				Rez[i].push_back(0);
				if (i == j) {
					P[i].push_back(1);
				}
				else if (i != j) {
					P[i].push_back(0);
				}
			}
		}
		LU(A, L, U, n, P);
		cout << "-------------------- Fisrt matrix --------------------" << endl;
		show(A, n, n);
		cout << "-------------------- U matrix --------------------" << endl;
		show(U, n, n);
		cout << "-------------------- L matrix --------------------" << endl;
		show(L, n, n);
		proisv(L, U, Rez, n);
		cout << "-------------------- L*U matrix --------------------" << endl;
		show(R, n, n);
		cout << "-------------------- P matrix --------------------" << endl;
		show(P, n, n);

		vector<double> y, x, Ax;
		for (int i = 0; i < n; i++) {
			x.push_back(0);
			y.push_back(0);
			Ax.push_back(0);
			cout << "b = " << b[i] << endl;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				y[i] += P[i][j] * b[j];
			}
		}

		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= (i - 1); j++) {
				y[i - 1] -= L[i - 1][j - 1] * y[j - 1];
			}
			x[i - 1] = y[i - 1];
		}

		for (int i = n; i >= 1; i--) {
			if (U[i - 1][i - 1] != 0) {
				for (int j = (i + 1); j <= n; j++) {
					x[i - 1] -= U[i - 1][j - 1] * x[j - 1];
				}
				x[i - 1] /= U[i - 1][i - 1];
			}
		}
		for (int i = 0; i < n; i++) {
			cout << "x = " << x[i] << endl;
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Ax[i] += A[i][j] * x[j];
			}
			cout << "A*x = " << Ax[i] << "  b =" << b[i] << endl;
		}

	}
	else {
		cout << "Matrix is NE sovmestnaya" << endl;
	}
	system("pause");
	return 0;
}