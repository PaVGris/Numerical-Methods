#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;

void showLagrangeInterpolation(const vector<double> &a, const vector<double> &values) {
    cout << "L(x) = ";
    int k = 0;

    while (k != 4) {
        string coef_sign = "+";
        if (a[k] < 0) coef_sign = "-";
        cout << coef_sign << abs(a[k]);
        for (int i = 0; i < values.size(); i++) {
            string sign = "-";
            if (k != i) {
                if (values[i] == 0) cout << "x";
                else {
                    if (values[i] < 0) sign = "+";
                    cout << "(x" << sign << abs(values[i]) << ")";
                }
            }
        }
        k++;
    }

    cout << "\n";
}


void showNewtonInterpolation(const vector<double> &a, const vector<double> &values) {
    cout << "P(x) = ";
    int k = 0;
    string tmp;
    for (double el: a) {
        string coef_sign = "+";
        if (el < 0) coef_sign = "-";
        cout << coef_sign << abs(el) << tmp;
        string sign = "-";
        if (values[k] == 0) tmp += "x";
        else {
            if (values[k] < 0) sign = "+";
            tmp += "(x" + sign;
            stringstream stream;
            stream << fixed << setprecision(1) << abs(values[k]);
            tmp += stream.str();
            tmp += ")";
        }
        k++;
    }
    cout << "\n";
}


vector<double> find_w(const vector<double>& values) {
    vector<double> res(values.size());
    for(int i = 0; i < values.size(); i++) {
        double w = 1;
        for (int j = 0; j < values.size(); j++) {
            if (i == j) continue;
            else w *= (values[i] - values[j]);
        }
        res[i] = w;
    }
    return res;
}

vector<double> func_values(const vector<double>& values) {
    vector<double> res (values.size());
    for (int i = 0; i < values.size(); i++) {
        res[i] = exp(values[i]);
    }
    return res;
}

vector<double> find_coef(const vector<double>& values) {
    vector<double> res (values.size());
    vector<double> w = find_w(values);
    vector<double> f = func_values(values);
    for (int i = 0; i < values.size(); i++) {
        res[i] = f[i]/w[i];
    }
    return res;
}

double interpolation(const vector<double>& coef, const vector<double>& values, double x) {
    double res = 0;
    for (int i = 0; i < coef.size(); i++) {
        double sum = coef[i];
        for (int j = 0; j < coef.size(); j++) {
            if (i == j) continue;
            else {
                sum *= (x - values[j]);
            }
        }
        res += sum;
    }
    return res;
}


void coef_interpolation(vector<double> &values, vector<double>& coef) {
    int x_size = (int)values.size();
    int k = (1 + x_size) * x_size / 2 - 1;
    int l = 0;
    coef.push_back(exp(values[0]));
    for (int i = x_size; i <= k; i += x_size - l) {
        int k = 0;
        for (int j = x_size - l; j > 1; j--) {
            double first = i == x_size ? exp(values[i - j]) : values[i - j];
            double second = i == x_size ? exp(values[i - j + 1]) : values[i - j + 1];
            double diff = values[k] - values[k + l + 1];
            values.push_back((first - second) / diff);
            if (j == x_size - l) coef.push_back((first - second) / diff);
            k++;
        }
        l++;
    }
}

double newton_interpolation(const vector<double> &x_1, double x) {
    int k = 4;
    int l = 3;
    double res = exp(x_1[0]);
    for (int i = k; i < x_1.size(); i += k) {
        double sum = x_1[i];
        for (int j = 0; j <= 3 - l; j++) {
            sum *= (x - x_1[j]);
        }
        res += sum;
        k--;
        l--;
    }
    return res;
}


int main() {
    vector<double> x_1 = {-2, -1, 0, 1};
    vector<double> x_2 = {-2, -1, 0.2, 1};
    double x = -0.5;

    vector<double> coef_1 = find_coef(x_1);
    double f_1 = interpolation(coef_1, x_1, x);

    vector<double> coef_2 = find_coef(x_2);
    double f_2 = interpolation(coef_2, x_2, x);
    cout << "Lagrange interpolation\nVariant a:\n";
    showLagrangeInterpolation(coef_1, x_1);
    cout << "Interpolation for x_1: "<< f_1 << "\n" << "Difference: " << abs(exp(x) - f_1) << "\n\n";

    cout << "Variant b:\n";
    showLagrangeInterpolation(coef_2, x_2);
    cout << "Interpolation for x_2: "<< f_2 << "\n" << "Difference: " << abs(exp(x) - f_2) << "\n\n";

    cout << "Newton interpolation\nVariant a:\n";
    vector<double> first_coef;
    coef_interpolation(x_1, first_coef);
    showNewtonInterpolation(first_coef, x_1);
    double res_1 = newton_interpolation(x_1, x);
    cout << "Interpolation for x_1: "<< res_1 << "\n" << "Difference: " << abs(exp(x) - res_1) << "\n\n";

    cout << "Variant b:\n";
    vector<double> second_coef;
    coef_interpolation(x_2, second_coef);
    showNewtonInterpolation(second_coef, x_2);
    double res_2 = newton_interpolation(x_2, x);
    cout << "Interpolation for x_2: "<< res_2 << "\n" << "Difference: " << abs(exp(x) - res_2) << "\n\n";
    return 0;
}
