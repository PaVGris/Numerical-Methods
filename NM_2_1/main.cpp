#include <iostream>
#include <math.h>
using namespace std;

class newtonMethod{
public:
    double phi_x(const double &x) {
        return pow(M_E, x) - 2 * x - 2;
    }

    double dx_phi_x(const double &x) {
        return  pow(M_E, x) - 2;
    }

    double d2x_phi_x(const double &x) {
        return  pow(M_E, x);
    }

    bool check_precision(const double &x_0, const double &x_1, const double &eps) {
        return abs(x_1 - x_0) > eps;
    }
};

class simpleIterations {

    double dx_phi_x(const double &x) {
        return 1/(x + 1);
    }

public:
    double get_q(const double &a, const double &b) {
        return max(abs(dx_phi_x(a)), abs(dx_phi_x(b)));
    }

    double phi_x(const double &x) {
        return log(2 * x + 2);
    }

    bool check_precision(const double &q, const double &x_0, const double &x_1, const double &eps) {
        return q/(1-q) * abs(x_1 - x_0) > eps;
    }
};
int main() {
    double a, b, eps;
    simpleIterations simple;
    a = 1.2; b = 2, eps = 0.0001;
    double q = simple.get_q(a, b);
    double x_0 = (a + b) / 2;
    double x_1 = simple.phi_x(x_0);
    int k = 0;
    while (simple.check_precision(q, x_0, x_1, eps)) {
        x_0 = x_1;
        x_1 = simple.phi_x(x_0);
        k++;
    }
    cout << "Simple iterations:\n" << x_1 << "\nIterations: "<< k <<"\n\n";


    //newton
    newtonMethod newton;
    x_0 = b;
    double f_x = newton.phi_x(x_0);
    double df_x = newton.dx_phi_x(x_0);
    double d2f_x = newton.d2x_phi_x(x_0);

    if (f_x * d2f_x < 0) {
        cout << "Can't apply method or choose different interval\n";
        return -1;
    }

    x_1 = x_0 - f_x/df_x;
    k = 0;
    while(newton.check_precision(x_0, x_1, eps)) {
        x_0 = x_1;
        f_x = newton.phi_x(x_0);
        df_x = newton.dx_phi_x(x_0);
        x_1 = x_0 - f_x/df_x;
        k++;
    }

    cout << "Newton's method:\n" << x_1 << "\nIterations: "<< k <<"\n\n";;
    return 0;
}
