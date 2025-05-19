#include <iostream>
#include <vector>

using namespace std;

void showMatrix(const vector <vector <double>> &A, int n) {
    for(int i = 0; i < n; i++) {
        cout << A[i][0] << "\t";
        for(int j = 1; j < n; j++) {
            cout << A[i][j] << "\t";
        }
        cout << "\n";
    }
}

void showVector(const vector<double> &x) {
    for(int i = 0; i < x.size(); i++) {
        cout << "x[" << i+1 << "] = "<< x[i] << "\n";
    }
}

vector<double> multiplyByVector(vector <vector <double>> A, vector <double> B) {
    vector<double> R(B.size());
    for(int i = 0; i < B.size(); i++)
        for(int j = 0; j < B.size(); j++)
            R[i] += A[i][j] * B[j];
    return R;
}

vector<double> add(const vector<double> &R, const vector<double> &B) {
    vector<double> res(R.size());
    for (int i = 0; i < R.size(); i++) {
        res[i] = R[i] + B[i];
    }
    return res;
}

bool check(vector<vector<double>> &a) {
    for (int i = 0; i < a.size(); i++) {
        double diag_elem = abs(a[i][i]);
        double sum = 0;
        for (int j = 0; j < a.size(); j++) {
            if (i == j) continue;
            sum += abs(a[i][j]);
        }
        if (diag_elem <= sum) return true;
    }
    return false;
}

double maxDiff (const vector<double> &prev, const vector<double> &curr) {
    double maximum = abs(curr[0] - prev[0]);
    for (int i = 1; i < prev.size(); i++) {
        maximum = std::max(maximum, abs(curr[i] - prev[i]));
    }
    return maximum;
}

void seidelMethod (const vector<double>& B, const vector<vector<double>>& A, vector<double> &tmp) {
    for (int i = 0; i < tmp.size(); i++) {
        tmp[i] = B[i];
        for (int j = 0; j < tmp.size(); j++) {
            if (i != j) tmp[i] += A[i][j] * tmp[j];
        }
    }
}

int main() {
    const int dimension = 4;
    cout.precision(4);
    vector<vector<double>> a = {{23, -6, -5, 9},
                               {8, 22, -2, 5},
                               {7, -6, 18, -1},
                               {3, 5, 5, -19}};
    vector<double> b = {232, -82, 202, -57};

    vector<double> B(dimension);

    if (check(a)) {
        return -1;
    }

    cout << "Success! Diagonal elements are less than the sum of non-diagonal elements\n\n";

    for (int i = 0; i < dimension; i++) {
        B[i] = b[i]/a[i][i];
    }

    vector<double> sum(dimension);
    vector<vector<double>> A(dimension);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            A[i].push_back(-a[i][j]/a[i][i]);
            sum[i] += abs(A[i][j]);
        }
        sum[i] -= abs(A[i][i]);
        A[i][i] = 0;
    }

    cout << "Modified matrix a: \n";
    showMatrix(A, A.size());

    cout << "Modified vector b: \n";
    showVector(B);

    double eps = 0.001;
    vector<double> res = B;
    double currentEps = 10;
    int k = 0;
    cout << "\nSimple method:\n";
    while (currentEps - eps > 0) {
        vector<double> tmp  = multiplyByVector(A, res);
        tmp = add(tmp, B);
        currentEps = maxDiff(res, tmp);
        res = std::move(tmp);
        k++;
    }
    showVector(res);
    cout << "Number of iteration: "<< k << "\n\n";
    k = 0;
    res = B;
    currentEps = 10;
    cout << "Seidel's method:\n";
    while (currentEps - eps > 0) {
        vector<double> tmp = res;
        seidelMethod(B, A, res);
        currentEps = maxDiff(tmp, res);
        k++;
    }
    showVector(res);
    cout << "Number of iterations: " << k << "\n";
}
