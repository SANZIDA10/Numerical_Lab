#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream in("D:\\Numerical project\\Least Square Regression Linear Equation\\input.txt");
    ofstream out("D:\\Numerical project\\Least Square Regression Linear Equation\\output.txt");

    if (!in) {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int n;
    in >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++)
        in >> x[i] >> y[i];

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }

    double b = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
    double a = (sumY - b * sumX) / n;

    out << fixed << setprecision(4);
    out << "Linear Equation: y = " << a << " + " << b << "x\n";

    in.close();
    out.close();
}
