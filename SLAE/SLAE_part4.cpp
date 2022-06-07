#include <math.h>
#include<iostream>
#include <iomanip>
#define eps 0.000001 
using namespace std;

double norm(double** A, int N) {
	double tmp_norm = 0;
	double max_norm = 0;
	for (int i = 0; i < N; ++i)
		max_norm += fabs(A[0][i]);

	for (int i = 1; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
			tmp_norm += fabs(A[i][j]);
		if (tmp_norm > max_norm)
			max_norm = tmp_norm;
		tmp_norm = 0;
	}

	return max_norm;
}

double Jacobi(int N, double** A, double* F, double* X) {

	double* TempX = new double[N];

	double norm, n;

	for (int k = 0; k < N; k++)
		TempX[k] = X[k];
	int cnt = 0;
	do {
		for (int i = 0; i < N; i++)
		{
			TempX[i] = F[i];
			for (int g = 0; g < N; g++)
				if (i != g)
					TempX[i] -= A[i][g] * X[g];
			TempX[i] /= A[i][i];
		}
		norm = abs(X[0] - TempX[0]);
		for (int h = 0; h < N; h++)
		{
			if (abs(X[h] - TempX[h]) > norm)
				norm = abs(X[h] - TempX[h]);
			X[h] = TempX[h];
		}
		cnt++;
		n = eps * fabs(norm - 1) / norm;
	} while (norm >= eps * 0.01);

	delete[] TempX;
	return cnt;

}



double razn(double* x, double* x_tmp, int n) {
	double raz = 0;
	for (int i = 0; i < n; i++) {
		if (abs(x[i] - x_tmp[i]) > raz) {
			raz = abs(x[i] - x_tmp[i]);
		}
	}
	return raz;
}

int Zeidel(int N, double** A, double* B, double* X, double& a) {

	double** C = new double* [N];
	for (int i = 0; i < N; i++) {
		C[i] = new double[N];
		for (int j = 0; j < N; j++) {
			C[i][j] = 0;
		}
	}

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (i != j)
				C[i][j] = -A[i][j] / A[i][i];

	double* d = new double[N];
	double norm_d = 0;

	for (int i = 0; i < N; ++i) {
		d[i] = B[i] / A[i][i];
		if (d[i] > norm_d) {
			norm_d = d[i];
		}
	}

	double norm_c = norm(C, N);



	a = (log((1 - norm_c) * eps / norm_d) / log(norm_c) - 1) / 2;

	double* x_tmp = new double[N];
	for (int i = 0; i < N; i++) {
		x_tmp[i] = 0;
	}

	for (int iter = 0; iter < 1000; ++iter)
	{
		for (int i = 0; i < N; ++i)
		{
			X[i] = 0;
			for (int j = 0; j < i; ++j)
				X[i] += C[i][j] * X[j];
			for (int j = i; j < N; ++j)
				X[i] += C[i][j] * x_tmp[j];
			X[i] += d[i];

		}

		if (razn(X, x_tmp, N) < eps * 0.1 * ((1 - norm_c) / norm_c))
			return iter + 1;

		for (int i = 0; i < N; i++) {
			x_tmp[i] = X[i];
		}
	}
	return -1;
}


int main() {
	const int n = 4;

	double apr;

	double** A = new double* [n];
	double* b = new double[n];
	double* x = new double[n];
	double* Ax = new double[n];

	for (int i = 0; i < n; i++) {
		b[i] = rand() % 40;
		Ax[i] = 0;
		x[i] = 1;
		A[i] = new double[n];
		for (int j = 0; j < n; j++) {

			A[i][j] = rand() % 100;
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][i] += A[i][j];
		}
	}


	for (int i = 0; i < n; i++) {
		cout << "b = " << b[i] << endl;
	}
	cout << "-------------------- matrix --------------------" << endl;;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	cout << "err = " << eps << endl;
	cout << endl;
	//ЯКОБИ
	cout << "-------------------- JAKOBI --------------------" << endl;
	int J = Jacobi(n, A, b, x);
	cout << "iter Jakobi: " << J << endl;

	for (int i = 0; i < n; i++) {
		cout << "x = " << x[i] << endl;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Ax[i] += (A[i][j] * x[j]);
		}
	}

	for (int i = 0; i < n; i++) {
		cout << fixed << setprecision(6) << "Ax = b = " << Ax[i] << endl;
		x[i] = 0;
	}

	//ЗЕЙДЕЛЬ
	cout << "-------------------- ZEIDEL --------------------`" << endl;
	int Z = Zeidel(n, A, b, x, apr);
	cout << "apriornay ozenka Ziedel`= " << apr << endl;
	cout << "iter Zeidel: " << Z << endl;

	for (int i = 0; i < n; i++) {
		cout << "x = " << x[i] << endl;
		Ax[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Ax[i] += (A[i][j] * x[j]);
		}
	}

	for (int i = 0; i < n; i++) {
		cout << fixed << setprecision(6) << "Ax = b = " << Ax[i] << endl;
	}


	system("pause");
	return 0;
}