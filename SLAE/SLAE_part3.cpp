#include <iostream>
#include <vector>

using namespace std;

double scal(vector <double> a, int N) {
    double m = 0;
    for (int i = 0; i < N; i++) {
        m += (a[i] * a[i]);
    }
    return m;
}

bool equal(vector <vector <double>> A, vector <vector <double>> B, double err, int n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (std::abs(A[i][j] - B[i][j]) > err)
                return 0;
    return 1;
}

void proisv(vector <vector <double>> A, vector <vector <double>> B,
    vector <vector <double>>& R, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                R[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void show(vector <vector <double>> A, int n)
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

void trans(vector <vector <double>> Q, vector <vector <double>>& Q_trans, int n) {
    Q_trans = Q;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double tmp = Q_trans[i][j];
            Q_trans[i][j] = Q_trans[j][i];
            Q_trans[j][i] = tmp;

        }
    }
}

void QR_dec(vector <vector <double>> A, vector <vector <double>>& Q,
    vector <vector <double>>& R, int N) {

    vector <vector <double>> mulp(N), H(N), E(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mulp[i].push_back(0);
            H[i].push_back(0);
            if (i == j) {
                E[i].push_back(1);
                Q[i][j] = 1;
            }
            else {
                E[i].push_back(0);
            }
        }
    }

    R = A;
    for (int iter = 0; iter < N; ++iter)
    {
        vector <double> y;
        for (int i = 0; i < N; i++) {
            y.push_back(0);
        }
        for (int i = iter; i < N; ++i)
            y[i] = R[i][iter];

        double norm_y = std::sqrt(scal(y, N));

        vector <double> e;
        for (int i = 0; i < N; i++) {
            e.push_back(0);
        }
        int sign;
        if (R[iter][iter] > 0)
            sign = 1;
        else
            sign = -1;

        norm_y *= sign;
        e[iter] = norm_y;
        vector <double> w;
        for (int i = 0; i < N; i++) {
            w.push_back(0);
            w[i] = y[i] + e[i];
        }

        double k = scal(w, N);


        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                mulp[i][j] = w[i] * w[j];
            }
        }


        for (int i = iter; i < N; ++i)
            for (int j = iter; j < N; ++j)
                mulp[i][j] /= k * 0.5;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                H[i][j] = E[i][j] - mulp[i][j];
            }
        }
        proisv(H, R, R, N);
        proisv(Q, H, Q, N);
    }
}

void sys_lin(vector <vector <double>> Q, vector <vector <double>> R,
    vector <double>& x, vector <double> b, int N)
{
    vector <vector <double>> Qtb(N), A(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Qtb[i].push_back(0);
            A[i].push_back(0);
        }
    }

    vector<double> y, a, c, ax;
    for (int i = 0; i < N; i++) {
        y.push_back(0);
        ax.push_back(0);
    }

    trans(Q, Qtb, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            y[i] += Qtb[i][j] * b[j];
        }
    }

    for (int i = N; i >= 1; i--) {
        x[i - 1] = y[i - 1];
        for (int j = N; j > i; --j) {
            x[i - 1] -= x[j - 1] * R[i - 1][j - 1];
        }
        x[i - 1] /= R[i - 1][i - 1];
    }


}


int main() {
    int n = 5;
    vector <vector <double>> A(n), Q(n), R(n), Ed(n), Qtr(n), QR(n);
    vector <double> b, x, Ax;
    for (int i = 0; i < n; i++) {
        b.push_back(rand() % 100);
        x.push_back(0);
        Ax.push_back(0);
        for (int j = 0; j < n; j++) {
            A[i].push_back(rand() % 100);
            Q[i].push_back(0);
            R[i].push_back(0);
            Ed[i].push_back(0);
            Qtr[i].push_back(0);
            QR[i].push_back(0);
        }
    }

    QR_dec(A, Q, R, n);
    cout << "-------------------- matrix A --------------------" << endl;
    show(A, n);

    cout << "-------------------- matrix Q --------------------" << endl;
    show(Q, n);

    cout << "-------------------- matrix R --------------------" << endl;
    show(R, n);

    cout << "-------------------- matrix Q*Q_tr --------------------" << endl;
    trans(Q, Qtr, n);
    proisv(Q, Qtr, Ed, n);
    show(Ed, n);

    cout << "-------------------- matrix Q*R --------------------" << endl;
    proisv(Q, R, Ed, n);
    show(Ed, n);

    for (int i = 0; i < n; i++) {
        cout << "b = " << b[i] << endl;
    }
    sys_lin(Q, R, x, b, n);

    proisv(Q, R, QR, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Ax[i] += QR[i][j] * x[j];
        }
    }

    for (int i = 0; i < n; i++) {
        cout << "Ax = " << Ax[i] << endl;
    }
}