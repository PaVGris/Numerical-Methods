#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

void LU(vector <vector <double>> A, vector <vector <double>> &L,
        vector <vector <double>> &U, vector <double> &b, int n) {
    U=std::move(A);

    for(int k = 1; k < n; k++) {
        for(int i = k-1; i < n; i++)
            for(int j = i; j < n; j++){
                L[j][i]=U[j][i]/U[i][i];
            }

        for(int i = k; i < n; i++)
            for(int j = k-1; j < n; j++) {
                U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j];
            }
    }

}

void showVector(const vector<double> &x, ofstream& outfile) {
    cout << "\n";
    cout << "P(x) = " << x[0] << " + ";
    for(int i = 0; i < x.size(); i++) {
        outfile << x[i] << " ";
        if(i != 0 && i != x.size()-1){
            cout << x[i] << " * x^" << i << " + ";
        }
        else if(i == x.size()-1) {
            cout << x[i] << " * x^" << i;
        }
    }
    cout << "\n";
    outfile << "\n";
}

vector<double> solve(const vector <vector <double>> &L, const vector <vector <double>> &U,
                     const vector <double> &b, int n) {
    vector<double> z(n);
    z[0] = b[0];
    for (int i = 1; i < n; i++) {
        z[i]=b[i];
        for (int j = 0; j < i; j++) {
            z[i] -= z[j]*L[i][j];
        }
    }
    vector<double> x(n);
    x[n-1] = z[n-1]/U[n-1][n-1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = z[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= x[j]*U[i][j];
        }
        x[i] = x[i]/U[i][i];
    }
    return x;
}

vector<double> solveSystem (const vector<vector<double>> a, vector<double> b) {
    const int dimension = (int)a.size();

    vector <vector <double>> L(dimension);
    vector <vector <double>> U(dimension);
    vector <vector <double>> R(dimension);
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            L[i].push_back(0);
            U[i].push_back(0);
            R[i].push_back(0);
        }
    }
    LU(a, L, U, b,dimension);
    vector<double> res = solve(L, U, b, dimension);
    cout << "Solution: ";
    return res;
}

void matrixInit (vector<vector<double>> &a) {
    int d = a.size();
    for (int i = 0; i < d; i++) {
        for(int j = 0; j < d; j ++) {
            a[i].push_back(0);
        }
    }
}

double findPolynomValues (vector<double> &a, double x, int dimension) {
    double res = 0;
    for (int i = 0; i < dimension; i++) {
        res += a[i] * pow(x, i);
    }
    return res;
}

int main() {
    std::ofstream outfile (R"(C:\Users\pavel\Desktop\DA\NM_3_3\test.txt)");
    vector<double> val = {-3, -2, -1, 0, 1, 2};
    vector<double> f = {0.04979, 0.13534, 0.36788, 1.0, 2.7183, 7.3891};
    int max_x_pow = 5;
    int max_y_pow = 3;
    vector<double> pow_x_sum(max_x_pow);
    vector<double> pow_y_sum(max_y_pow);
    for (int i = 0; i < max_x_pow; i++) {
        pow_x_sum[i] = 0;
        for (int j = 0; j < val.size(); j++) {
            pow_x_sum[i] += pow(val[j], i);
        }
    }

    for (int i = 0; i < max_y_pow; i++) {
        pow_y_sum[i] = 0;
        for (int j = 0; j < val.size(); j++) {
            pow_y_sum[i] += (f[j] * pow(val[j], i));
        }
    }

    vector<vector<double>> first (2);
    matrixInit(first);
    for (int i = 0; i < first.size(); i++) {
        for (int j = 0; j < first.size(); j++) {
            first[i][j] = pow_x_sum[i + j];
        }
    }
    vector<double> coef_first = solveSystem(first, pow_y_sum);
    showVector(coef_first, outfile);
    double phi_1 = 0;
    for (int i = 0; i < val.size(); i++) {
        phi_1 += pow(findPolynomValues(coef_first, val[i], first.size()) - f[i], 2);
    }
    cout << "Error Sum of Squares: " << phi_1 << "\n\n";



    vector<vector<double>> second (3);
    matrixInit(second);
    for (int i = 0; i < second.size(); i++) {
        for (int j = 0; j < second.size(); j++) {
            second[i][j] = pow_x_sum[i + j];
        }
    }
    vector<double> coef_second = solveSystem(second, pow_y_sum);
    showVector(coef_second, outfile);
    double phi_2 = 0;
    for (int i = 0; i < val.size(); i++) {
        phi_2 += pow(findPolynomValues(coef_second, val[i], second.size()) - f[i], 2);
    }
    cout << "Error Sum of Squares: " << phi_2 << "\n\n";
    outfile.close();
    system( R"(python C:\Users\pavel\Desktop\DA\NM_3_3\plot.py)");
    return 0;
}
