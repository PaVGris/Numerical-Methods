#include <iostream>
#include <vector>

using namespace std;

void showInterpolation(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, int pos) {
    cout << "f(x) = " << a[pos - 1];
    int k = 0;
    while (k < 3) {
        k++;
        string sign = "+";
        if (b[pos - 1] < 0) sign = "";
        if (k == 1) cout <<  sign << b[pos-1];

        sign = "+";
        if (c[pos - 1] < 0) sign = "";
        if (k == 2) cout <<  sign << c[pos-1];

        sign = "+";
        if (d[pos - 1] < 0) sign = "";
        if (k == 3) cout <<  sign << d[pos-1];

        cout << "*x^" << k;
    }
    cout << "\n";
}

void showVector(const vector<double> &a) {
    cout << "{";
    for (auto el : a) {
        cout << el << "; ";
    }
    cout << "}\n";
}

int findPosition(double x, vector<double> &val) {
    for (int i = 1; i < val.size(); i++) {
        if (x > val[i-1] && val[i] >= x) {
            return i;
        }
    }
    return -1;
}

vector<double> shuttle(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d) {
    const int dimension = 3;
    vector<double> P(dimension);
    vector<double> Q(dimension);

    P[0] = -c[0]/b[0];
    Q[0] = d[0]/b[0];

    for(int i = 1; i < dimension; i++) {
        double sum = b[i]+a[i]*P[i-1];
        P[i] = -c[i]/sum;
        Q[i] = (d[i] - a[i]*Q[i-1])/sum;
    }

    vector<double> res(dimension);
    res[dimension-1] = Q[dimension-1];
    for(int i = dimension - 2; i >= 0; i--) {
        res[i] = P[i] * res[i+1] + Q[i];
    }
    return res;
}


int main() {
    double x = -0.5;
    vector<double> val = {-2, -1, 0, 1, 2};
    int pos = findPosition(x, val);
    cout << "Interval: x = " << x << " in [" << val[pos-1] << ";" << val[pos] << "]\n";
    vector<double> f = {0.13534, 0.36788, 1.0, 2.7183, 7.3891};
    vector<double> h (val.size() - 1);
    for (int i = 1; i < h.size() + 1; i++) {
        h[i - 1] = val[i] - val[i-1];
    }

    vector<double> a = {0, h[1], h[2]};;
    vector<double> b (a.size());
    vector<double> c = {h[1], h[2], 0};

    for (int i = 0; i < b.size(); i++) {
        b[i] = 2 * (h[i] + h[i+1]);
    }
    vector<double> d (a.size());
    for (int i = 0; i < d.size(); i++) {
        d[i] = 3 * (((f[i + 2] - f[i + 1])/h[i + 1]) - ((f[i + 1]-f[i])/h[i]));
    }
    vector<double> c_coef = shuttle(a, b, c, d);;
    c_coef.insert(c_coef.begin(), 0);


    vector<double> a_coef;
    for (int i = 0; i < f.size() - 1; i++) {
        a_coef.push_back(f[i]);
    }

    vector<double> b_coef;
    for (int i = 0; i < f.size() - 1; i++) {
        b_coef.push_back((f[i + 1] - f[i])/h[i] - h[i] * (2 * c_coef[i] + c_coef[i+1]) / 3);
    }
    b.push_back((f[4] - f[3])/h[4] - 2*h[4]*c_coef[4]/3);

    vector<double> d_coef;
    for (int i = 0; i < f.size() - 2; i++) {
        d_coef.push_back((c_coef[i+1] - c_coef[i])/h[i]/3);
    }
    d_coef.push_back(-c_coef[c_coef.size()-1]/h[h.size()-1]/3);

    cout << "Coef a: ";
    showVector(a_coef);
    cout << "Coef b: ";
    showVector(b_coef);
    cout << "Coef c: ";
    showVector(c_coef);
    cout << "Coef d: ";
    showVector(d_coef);
    double res = a_coef[pos-1] + b_coef[pos-1] * (x + 1) + c_coef[pos-1] * (x+1) * (x+1) + d_coef[pos-1] * (x+1) * (x+1) * (x+1);
    showInterpolation(a_coef,b_coef,c_coef,d_coef,pos);
    cout << "f(" << x << ") = "<< res;
    return 0;
}
