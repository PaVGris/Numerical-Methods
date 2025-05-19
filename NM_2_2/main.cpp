#include <iostream>
#include <math.h>
#include <vector>

using namespace std;


class newtonMethod{
public:

    double get_det(const vector<vector<double>> &A) {
        return A[0][0] * A[1][1] - A[1][0]*A[0][1];
    }

    vector<double> phi_x(const vector<double> &current) {
        double x_1 = current[0];
        double x_2 = current[1];
        vector<double> res(2);
        res[0] = x_1 - cos(x_2) - 1;
        res[1] = x_2 - log10(x_1 + 1) - 3;
        return res;
    }

    vector<vector<double>> get_J(const vector<double> &current) {
        double x_1 = current[0];
        double x_2 = current[1];
        vector<vector<double>> res = {{1, sin(x_2)},
                                      { -1 / (log(10) * (x_1 + 1)) , 1}};
        return res;
    }

    vector<double> get_next(const vector<double> &current,
                            double d_J, double d_A1, double d_A2) {
        vector<double> res (2);
        res[0] = current[0] - d_A1/d_J;
        res[1] = current[1] - d_A2/d_J;
        return res;
    }

    bool check_precision(const vector<double> &current, const vector<double> &next, const double &eps) {
        return max(abs(next[0] - current[0]), abs(next[1] - current[1])) > eps;
    }
};

class simpleIterations {

public:
    double get_q(const vector<double> &G) {
        vector<double> dx_phi_x(2);
        dx_phi_x[0] = -sin(max(G[2],G[3]));
        dx_phi_x[1] = 1/(min(G[0],G[1]) + 1)/ log(10);
        return max(abs(dx_phi_x[0]), abs(dx_phi_x[1]));
    }

    vector<double> phi_x(const vector<double> &current) {
        double x_1 = current[0];
        double x_2 = current[1];
        vector<double> res (2);
        res[0] = 1 + cos(x_2);
        res[1] = 3 + log10(x_1 + 1);
        return res;
    }

    bool check_precision(const double &q, const vector<double> &current, const vector<double> &next, const double &eps) {
        return q/(1-q) * max(abs(next[0] - current[0]), abs(next[1] - current[1])) > eps;
    }
};

int main() {
    double eps = 0.000001;
    vector<double> G = {0, 0.5, 3, 3.5};
    simpleIterations simple;
    double q = simple.get_q(G);
    vector<double> current (2);
    vector<double> next (2);
    current[0] = (G[0] + G[1]) / 2; current[1] = (G[2] + G[3]) / 2;

    int k = 0;
    while (simple.check_precision(q, current, next, eps)) {
        current = std::move(next);
        next = simple.phi_x(current);
        k++;
    }
    cout << "Simple iterations:\n" << next[0] << "; " << next[1] << "\nIterations: " << k <<"\n\n";


    //newton
    newtonMethod newton;
    current[0] = (G[0] + G[1]) / 2; current[1] = (G[2] + G[3]) / 2;
    vector<vector<double>> J = newton.get_J(current);
    vector<double> f = newton.phi_x(current);
    vector<vector<double>> A_1 = {{f[0] , J[0][1]},
                                  {f[1], J[1][1]}};
    vector<vector<double>> A_2 = {{J[0][0] , f[0]},
                                  {J[1][0], f[1]}};

    double d_J = newton.get_det(J);
    double d_A1 = newton.get_det(A_1);
    double d_A2 = newton.get_det(A_2);
    next = newton.get_next(current, d_J, d_A1, d_A2);

    k = 0;
    while(newton.check_precision(current, next, eps)) {
        current = std::move(next);
        f = newton.phi_x(current);
        J = newton.get_J(current);
        A_1 = {{f[0] , J[0][1]},
               {f[1], J[1][1]}};
        A_2 = {{J[0][0] , f[0]},
               {J[1][0], f[1]}};
        d_J = newton.get_det(J);
        d_A1 = newton.get_det(A_1);
        d_A2 = newton.get_det(A_2);
        next = newton.get_next(current, d_J, d_A1, d_A2);
        k++;
    }
    cout << "Newton's method:\n" << next[0] << "; " << next[1] << "\nIterations: " << k;
    return 0;
}
