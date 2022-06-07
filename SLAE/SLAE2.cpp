#include <iostream>
#include <vector>
#include<cmath>

using namespace std;

int chislo_perk = 0;

void show(vector <vector <double>> A, int n) // вывод матрицы
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%6f \t", A[i][j]);
		}
		cout << endl;
	}
}

void SWAP(vector <vector <double>>& U, int ii, int jj, int n) {
	if (ii != jj) {
		chislo_perk++;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if ((i == ii) && (j == jj) && (i != j)) {
					for (int a = 0; a < n; a++) {
						double tmp1 = U[i][a];
						U[i][a] = U[j][a];
						U[j][a] = tmp1;
					}
				}
			}
		}
	}

}

void LU(vector <vector <double>> A, vector <vector <double>>& L,
	vector <vector <double>>& U, int n, vector <vector <double>>& P)
{
	U = A;
	//в U записать строчки с максимальными элементами из столбца
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
		SWAP(U, pi, pj, n);
		SWAP(P, pi, pj, n);
	}


	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			L[j][i] = U[j][i] / U[i][i];

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++)
				L[j][i] = U[j][i] / U[i][i];

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
	int n = 5;
	vector <vector <double>> A(n), L(n), U(n), R(n), P(n), M(n), obratA(n), Otv_obr1(n), Otv_obr2(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i].push_back(rand() % 20 - 10);
			L[i].push_back(0);
			U[i].push_back(0);
			R[i].push_back(0);
			M[i].push_back(0);
			Otv_obr1[i].push_back(0);
			Otv_obr2[i].push_back(0);
			obratA[i].push_back(0);
			if (i == j) {
				P[i].push_back(1);
			}
			else if (i != j) {
				P[i].push_back(0);
			}
		}
	}
	LU(A, L, U, n, P);
	cout << "------------------- Fisrt matrix -------------------" << endl;
	show(A, n);
	cout << "------------------- U matrix -------------------" << endl;
	show(U, n);
	cout << "------------------- L matrix -------------------" << endl;
	show(L, n);
	proisv(L, U, R, n);
	cout << "------------------- L*U matrix -------------------" << endl;
	show(R, n);
	cout << "------------------- P matrix -------------------" << endl;
	show(P, n);
	proisv(P, A, M, n);
	cout << "------------------- P*A matrix -------------------" << endl;
	show(M, n);


	double opredelitel = 1;// определитель находим как произведение диагональных элементов
	for (int i = 0; i < n; i++) {
		opredelitel *= U[i][i];
	}
	cout << "------------------- opredelitel A = " << opredelitel * pow((-1), chislo_perk / 2) << "-------------------" << endl;

	vector<double> b, y, x, Ax;
	for (int i = 0; i < n; i++) {
		b.push_back(rand() % 20 - 10);
		x.push_back(0);
		y.push_back(0);
		Ax.push_back(0);
		cout << "b = " << b[i] << endl;
	}
	for (int i = 0; i < n; i++) {//y = U*x
		for (int j = 0; j < n; j++) {
			y[i] += P[i][j] * b[j];
		}
	}

	for (int i = 1; i <= n; i++) { //L*y=0
		for (int j = 1; j <= (i - 1); j++) {
			y[i - 1] -= L[i - 1][j - 1] * y[j - 1];
		}
		x[i - 1] = y[i - 1];
	}

	for (int i = n; i >= 1; i--) { //U*x=y
		for (int j = (i + 1); j <= n; j++) {
			x[i - 1] -= U[i - 1][j - 1] * x[j - 1];
		}
		x[i - 1] /= U[i - 1][i - 1];
		cout << "x = " << x[i - 1] << endl;
	}

	for (int i = 0; i < n; i++) { // Проверяем и выводим
		for (int j = 0; j < n; j++) {
			Ax[i] += A[i][j] * x[j];
		}
		cout << "A*x = " << Ax[i] << "  b =" << b[i] << endl;
	}

	for (int it = 0; it < n; it++) {
		vector<double> b, y, x, Ax;
		for (int i = 0; i < n; i++) {
			if (i == it) {
				b.push_back(1);
			}
			if (i != it) {
				b.push_back(0);
			}
			x.push_back(0);
			y.push_back(0);
			Ax.push_back(0);
			//cout << "bi = " << b[i] << "it = " << it << endl;
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
			for (int j = (i + 1); j <= n; j++) {
				x[i - 1] -= U[i - 1][j - 1] * x[j - 1];
			}
			x[i - 1] /= U[i - 1][i - 1];
			//cout << "x = " << x[i - 1] << endl;
		}
		for (int a = 0; a < n; a++) {
			obratA[a][it] = x[a];
		}
	}
	cout << "------------------- obratnaya matrix -------------------" << endl;
	show(obratA, n);

	proisv(A, obratA, Otv_obr1, n);
	cout << "------------------- A*obrat matrix -------------------" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << abs(round(Otv_obr1[i][j])) << "\t";
		}
		cout << endl;
	}

	proisv(obratA, A, Otv_obr2, n);
	cout << "------------------- obrat*A matrix -------------------" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << abs(round(Otv_obr2[i][j])) << "\t";
		}
		cout << endl;
	}

	double normaA = 0, normaObratA = 0;
	for (int i = 0; i < n; i++) {
		double sumAij = 0, sumObratAij = 0;
		for (int j = 0; j < n; j++) {
			sumAij += abs(A[i][j]);
			sumObratAij += abs(obratA[i][j]);
		}
		if (sumAij > normaA) {
			normaA = sumAij;
		}
		if (sumObratAij > normaObratA) {
			normaObratA = sumObratAij;
		}
	}
	cout << "------------------- CHislo OBUSLOVLENOSTI = " << normaA * normaObratA << "-------------------" << endl;


	system("pause");

	return 0;
}