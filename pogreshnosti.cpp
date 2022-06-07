#include<iostream>
#include <math.h>
#include <iomanip> 

const double err = 0.0000000001;
using namespace std;

int fact(int x) {
    if (x < 0)
        return 0;
    if (x == 0)
        return 1;
    else
        return x * fact(x - 1);
}

double my_sqrt(double x) {
    double a = (1.0 / 2.0) * (15.0 + x / 15.0);
    for (int i = 0; i < 50; i++) {
        double a_n = (1.0 / 2.0) * (a + x / a);
        a = a_n;
    }
    return a;
}

double my_sin(double x) {
    double n = x;
    double sum = 0.0;
    int i = 3;
    do {
        sum += n;
        n *= -1.0 * x * x / fact(i);
        i += 2;
    } while (fabs(n) > err);
    return sum;
}

double my_cos(double x) {
    double n = 1.0;
    double sum = 0.0;
    int i = 2;
    do {
        sum += n;
        n *= -1.0 * x * x / fact(i);
        i += 2;
    } while (fabs(n) > err);
    return sum;
}
double f_approx(double x) {
    return my_sqrt(1 + x * x) * (my_sin(3 * x + 0.1) + my_cos(2 * x + 0.3));
}
double f_exact(double x) {
    return sqrt(1 + x * x) * (sin(3 * x + 0.1) + cos(2 * x + 0.3));
}
int main() {
    for (double x = 0.2; x <= 0.3; x += 0.01) {
        double error = f_exact(x) - f_approx(x);

        std::cout << setprecision(5) << x << '\t' << f_exact(x) << '\t' << f_approx(x) << '\t';
        std::cout << error << '\n';
    }
    return 0;
}