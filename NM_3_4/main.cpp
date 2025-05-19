#include <iostream>
#include <vector>

using namespace std;

int findPosition(double x, vector<double> &val) {
    for (int i = 1; i < val.size(); i++) {
        if (x > val[i-1] && val[i] >= x) {
            return i;
        }
    }
    return -1;
}

int main() {
    vector<double> val = {-0.2, 0, 0.2, 0.4, 0.6};
    vector<double> y = {-0.20136, 0, 0.20136, 0.41152, 0.64350};
    double x = 0.2;
    int pos = findPosition(x, val);
    double  l = (y[pos] - y[pos-1])/(val[pos] - val[pos-1]);
    double r = (y[pos+1] - y[pos])/(val[pos+1] - val[pos]);
    double tmp = (r-l)/(val[pos+1] - val[pos-1]);
    double c = l + tmp * (2 * x - val[pos-1] - val[pos]);
    double d2x = 2*tmp;
    cout << "Left = " << l << "\n";
    cout << "Right = " << r << "\n";
    cout << "Central = "<< c << "\n";
    cout << "d2x = " << d2x;
    return 0;
}
