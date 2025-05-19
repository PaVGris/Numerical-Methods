#include <iostream>
#include <vector>

using namespace std;

void showVector(vector<double> x) {
    cout << "\n";
    for(int i = 0; i < x.size(); i++) {
        cout << "x[" << i+1 << "] = "<< x[i] << "\n";
    }
}

bool check(vector<double> &a, vector<double> &b, vector<double> &c, int n) {
    for (int i = 0; i < n; i++) {
        if (abs(b[i]) < abs(a[i]) + abs(c[i])) return false;
    }
    return true;
}

int main() {
    const int dimension = 5;
    vector<double> a = {0, -6, 9, 8, 6};
    vector<double> b = {6, 16, -17, 22, -13};
    vector<double> c = {-5, 9, -3, -8, 0};

    if (!check(a,b,c,dimension)) {
        return -1;
    }

    cout << "Success! Conditions for stability and correctness are met\n";
    vector<double> d = {-58, 161, -114, -90, -55};
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

    showVector(res);
}
