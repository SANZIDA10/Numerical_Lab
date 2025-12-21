#include <bits/stdc++.h>
using namespace std;

int main() {

    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in) {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int n;
    in >> n;

    vector<double> x(n), y(n);


    for(int i = 0; i < n; i++)
        in >> x[i] >> y[i];


    vector<double> fx(n);
    for(int i = 0; i < n; i++)
        fx[i] = exp(x[i] / 4.0);

    double sumFx = 0, sumY = 0, sumFxY = 0, sumFx2 = 0;


    for(int i = 0; i < n; i++) {
        sumFx  += fx[i];
        sumY   += y[i];
        sumFxY += fx[i] * y[i];
        sumFx2 += fx[i] * fx[i];
    }

    double b = (n * sumFxY - sumFx * sumY) / (n * sumFx2 - sumFx * sumFx);
    double a = (sumY - b * sumFx) / n;

    out << fixed << setprecision(6);
    out << "Regression Equation: T = " 
        << a << " + " << b << " * e^(t/4)\n";

    
    double t = 6;
    double T = a + b * exp(t / 4.0);

    out << "T when t = 6: " << T << endl;

    in.close();
    out.close();
    return 0;
}
