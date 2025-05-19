#include <iostream>
#include <valarray>

using namespace std;

double f(double x) {
    return x / ((2*x+7) * (3*x+4));
}

double richardsonAbs(double res1, double res2, int k) {
    double coef = (res2 - res1)/(pow(2,k) - 1);
    double absolute = abs(coef);
    return absolute;
}

double richardson(double res1, double res2, int k) {
    double coef = (res2 - res1)/(pow(2,k) - 1);
    return res2 + coef;
}

double rectangleMethod(double h, double x_0, double x_k) {
    int k = (int)((x_k - x_0)/h);
    double res = 0;
    for (int i = 0; i < k; i++) {
        res += h * (f(x_0 + h * i + 0.5 * h));
    }
    return res;
}

double trapezoidMethod(double h, double x_0, double x_k) {
    int k = (int)((x_k - x_0)/h);
    double res = 0;
    for (int i = 0; i < k; i++) {
        res += 0.5 * h * (f(x_0 + (i+1) * h) + f(x_0 + h * i));
    }
    return res;
}

double simpsonRule(double h, double x_0, double x_k) {
    int k = (int)((x_k - x_0)/h);
    double res = f(x_0);
    for (int i = 1; i < k; i++) {
        if (i % 2 != 0) res += 4 * f(x_0 + h * i);
        else res += 2 * f(x_0 + h * i);
    }
    res += f(x_0 + h * k);
    return h/3 * res;
}

int main() {
    double exact = -0.0413027217305138;
    double h1 = 0.5;
    double h2 = 0.25;
    double x_0 = -1;
    double x_k = 1;


    cout << "Rectangle Method: \n";
    double rec1 = rectangleMethod(h1, x_0, x_k);
    double rec2 = rectangleMethod(h2, x_0, x_k);
    cout << "For h1: " << rec1 << "\n" << "For h2: " << rec2 << "\n";
    cout << "Clarification: " << richardson(rec1, rec2, 2) << "\n";
    cout << "Absolute error: " << richardsonAbs(rec1, rec2, 2) << "\n\n";

    cout << "Trapezoid Method: \n";
    double trap1 = trapezoidMethod(h1, x_0, x_k);
    double trap2 = trapezoidMethod(h2, x_0, x_k);
    cout << "For h1: " << trap1 << "\n" << "For h2: " << trap2 << "\n";
    cout << "Clarification: " << richardson(trap1, trap2, 2) << "\n";
    cout << "Absolute error: " << richardsonAbs(rec1, rec2, 2) << "\n\n";

    cout << "Simpson Rule: \n";
    double simp1 = simpsonRule(h1, x_0, x_k);
    double simp2 = simpsonRule(h2, x_0, x_k);
    cout << "For h1: " << simp1 << "\n" << "For h2: " << simp2 << "\n";
    cout << "Clarification: " << richardson(simp1, simp2, 4) << "\n";
    cout << "Absolute error: " << richardsonAbs(simp1, simp2, 4) << "\n\n";

    cout.precision(5);
    cout << "Exact value: " << exact << "\n";
}
