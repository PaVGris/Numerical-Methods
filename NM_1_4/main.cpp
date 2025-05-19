#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

pair<int, int> findAbsMax(const vector<vector<double>> &a, const int n) {
    pair<int, int> tmp;
    double maximum = 0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (maximum < abs(a[i][j])) {
                maximum = abs(a[i][j]);
                tmp.first = i;
                tmp.second = j;
            }
        }
    }
    return tmp;
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

vector<vector<double>> createRotationMatrix(const vector <vector <double>> &a, int n, pair<int, int> pos) {
    vector<vector<double>> U(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) U[i].push_back(1);
            else U[i].push_back(0);
        }
    }
    int row = pos.first;
    int column = pos.second;
    double angle = atan((2 * a[row][column]) / (a[row][row] - a[column][column])) / 2;
    U[row][row] = U[column][column] = cos(angle);
    U[column][row] = sin(angle);
    U[row][column] = -U[column][row];

    return U;
}

void matrixTranspose (vector<vector<double>> &U) {
    for (int i = 0; i < U.size() - 1; i++) {
        for (int j = i + 1; j < U.size(); j++) {
            double tmp = U[i][j];
            U[i][j] = U[j][i];
            U[j][i] = tmp;
        }
    }
}

double precision (const vector<vector<double>> &a, int n) {
    double tmp = 0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            tmp += (a[i][j] * a[i][j]);
        }
    }
    return sqrt(tmp);
}
void showEigenValues(const vector<vector<double>> &a) {
    cout << "Eigen values:\n";
    for (int i = 0; i < a.size(); i++) {
        cout << "A[" << i+1 << "] = " << a[i][i] << "\n";
    }
    cout << "\n";
}
void showEigenVectors(const vector<vector<double>> &a) {
    cout << "Eigen vectors:\n";
    for(int i = 0; i < a.size(); i++) {
        cout << "x[" << i + 1 << "] = (";
        for(int j = 0; j < a.size() - 1; j++) {
            cout << a[j][i] << ";" << " ";
        }
        cout << a[a.size()-1][i] << ")" << "\n";
    }
}

int main() {
    const int dimension = 3;
    double eps = 0.0001;
    cout.precision(3);
    vector<vector<double>> a = {{ 9,  2, -7},
                                { 2, -4, -1},
                                {-7, -1,  1}};
    vector<vector<double>> eigenMatrix = {{1, 0, 0},
                                          {0, 1, 0},
                                          {0, 0, 1}};
    int k = 0;
    while (precision(a, dimension) - eps > 0) {
        k++;
        pair<int, int> position = findAbsMax(a, dimension);
        vector<vector<double>> U = createRotationMatrix(a, dimension, position);
        eigenMatrix = multiply(eigenMatrix, U, dimension);
        a = multiply(a, U, dimension);
        matrixTranspose(U);
        a = multiply(U, a, dimension);
    }
    showEigenValues(a);
    showEigenVectors(eigenMatrix);
}
