#include <iostream>
#include <vector>
#include <valarray>

using namespace std;

void showMatrix(const vector <vector <double>> &A, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n - 1; j++) {
            cout << A[i][j] << "  ";
        }
        cout << A[i][n-1] << "\n";
    }
}

void showEigenValues(const vector<vector<double>> &a, vector<double> &complex) {
    cout << "Eigen values:\n";
    for(int k = 0; k < a.size(); k++) {
        if(complex[k] == 0) cout << "A[" << k << "] = " << a[k][k] << " " << "\n";
        else {
            cout << "A[" << k << "] = " << complex[k] << " + " << complex[k+1] << " * i" << "\n";
            k++;
            cout << "A[" << k << "] = " << complex[k] << " - " << complex[k+1] << " * i" << "\n";
        }
    }
    cout << "\n";
}

void zeroing(vector<vector<double>> &a, int pos) {
    for (int i = pos + 1; i < a.size(); i++) {
        a[i][pos] = 0;
    }
}

int sign(double x) {
    if (x < 0) return -1;
    else if (x == 0) return 0;
    else return 1;
}

double normalization (const vector<vector<double>> &a, int pos) {
    double sum = 0;
    for (int i = pos; i < a.size(); i++) {
        sum += (a[i][pos] * a[i][pos]);
    }
    return sqrt(sum);
}

void multiNumber(vector<vector<double>> &a, double number) {
    for (int i = 0; i < a.size(); i++){
        for (int j = 0; j < a.size(); j++) {
            a[i][j] = a[i][j] * number;
        }
    }
}

void matrixDiff (vector<vector<double>> &a, vector<vector<double>> &b) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a.size(); j++) {
            a[i][j] = a[i][j] - b[i][j];
        }
    }
}


vector<double> solveRoots(double a, double b, double c) {
    vector<double> res = {0, 0};
    res[0] = -b/(2*a);
    res[1] = sqrt(-(b*b - 4*a*c))/(2*a);
    return res;
}


vector<vector<double>> expression(const vector<vector<double>> &tmp, int iter) {
    vector<vector<double>> V(tmp.size() + 1);
    for (int i = 0; i < tmp.size() + 1; i++){
        for (int j = 0; j < tmp.size() + 1; j++) {
            V[i].push_back(0);
            V[i][j] = tmp[iter][i] * tmp[iter][j];
        }
    }
    double v = 0;
    for (int i = 0; i < tmp.size() + 1; i++) {
        v += (tmp[iter][i] * tmp[iter][i]);
    }
    v = 2 / v;
    multiNumber(V, v);
    return V;
}
vector<vector<double>> householderMatrix(vector<vector<double>> &tmp, int iter) {
    vector<vector<double>> E = {{ 1, 0, 0},
                                { 0, 1, 0},
                                { 0, 0, 1}};
    vector<vector<double>> ex = expression(tmp, iter);
    matrixDiff(E, ex);
    return E;
}

vector <vector <double>> multiply(vector <vector <double>> &A, vector <vector <double>> &B, int n) {
    vector <vector <double>> R(n);
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            R[i].push_back(0);
            for(int k = 0; k < n; k++)
                R[i][j] += (A[i][k] * B[k][j]);
        }
    return R;
}


vector<vector<double>> householderTransform(vector<vector<double>> &a) {
    vector<vector<double>> Q = {{ 1, 0, 0},
                                { 0, 1, 0},
                                { 0, 0, 1}};
    vector<vector<double>> tmp (a.size() - 1);
    for (int i = 0; i < a.size() - 1; i++) {
        for (int j = 0; j < a.size(); j++) {
            tmp[i].push_back(0);
            if (i == j) {
                double norm = normalization(a, i);
                tmp[i][i] = a[i][i] + sign(a[i][i]) * norm;
            }
            else if (i < j) tmp[i][j] = a[j][i];
        }
        vector<vector<double>> H = householderMatrix(tmp, i);
        a = multiply(H, a, H.size());
        zeroing(a, i);
        Q = multiply(Q, H, Q.size());
    }
    a = multiply(a, Q, a.size());
    return a;
}

bool precision (const vector<vector<double>> &a, int n, vector<double> &complex,
                vector<vector<double>> &prev_roots, double eps) {
    vector<double> lambda(n);
    bool stop = true;
    for(int i = 0; i < n; i++) {
        bool zerosUnder = true;
        for (int j = i + 1; j < n; j++) {
            if (abs(a[j][i]) > eps) {
                zerosUnder = false;
                break;
            }
        }
        if(zerosUnder) lambda[i] = a[i][i];

        bool complexCheck = true;
        if(i < n - 1) {
            vector<double> roots = solveRoots(1, -a[i][i] - a[i + 1][i + 1], a[i][i] * a[i + 1][i + 1] - a[i][i + 1] * a[i + 1][i]);
            if(roots[0] - prev_roots[i][0] < eps && roots[1] - prev_roots[i][1] < eps) {
                complex[i] = roots[0];
                complex[i + 1] = roots[1];
            } else complexCheck = false;
            prev_roots[i][0] = roots[0];
            prev_roots[i][1] = roots[1];
        }
        if(!zerosUnder && !complexCheck) {
            stop = false;
            break;
        }
    }

    return stop;
}


int main() {
    const int dimension = 3;
    double eps = 0.00001;
    cout.precision(4);
    vector<double> complex(dimension);

    vector<vector<double>> prev_roots(dimension);
    for (int i = 0; i < dimension; i++) {
        prev_roots[i].push_back(0);
        prev_roots[i].push_back(0);
    }

    vector<vector<double>> a = {{ 8, -1, -3},
                                {-5,  9, -8},
                                { 4, -5,  7}};


    while (!precision(a, dimension, complex, prev_roots, eps)) {
        a = householderTransform(a);
    }

    showMatrix(a, a.size());
    showEigenValues(a, complex);
}
