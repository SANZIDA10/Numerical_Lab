#include <bits/stdc++.h>
using namespace std;

const double E = 1e-6;
const int max_it = 1000;

double f(double x) {
    return 3*x - cos(x) - 1;
}

void secant(double x, double x1, double &r, ofstream &out) {
    double fx = f(x);
    double fx1 = f(x1);

    r = x1 - ((fx1 * (x1 - x)) / (fx1 - fx));
    int counts = 0;

    while (fabs(x1 - x) >= E && counts < max_it) {
        counts++;
        fx = f(x);
        fx1 = f(x1);

        r = x1 - ((fx1 * (x1 - x)) / (fx1 - fx));
        x = x1;
        x1 = r;
    }

    if (counts >= max_it) {
        out << "No real root is found in this interval.\n";
        return;
    }

    out << "The root " << r << " is found after " << counts << " iterations.\n";
}

int main() {
    ifstream input("D:\\Numerical project\\secant method\\input.txt");
    ofstream out("D:\\Numerical project\\secant method\\output.txt");

    if (!input) {
        cout << "ERROR: input.txt not found!\n";
        return 0;
    }

    double start, ends;
    input >> start >> ends;
    input.close();

    double step = 0.5;
    double a = start, b, r;

    out << fixed << setprecision(6);
    out << "Secant Method\n";
    out << "Checking intervals for possible roots...\n";

    while (a < ends) {
        b = a + step;

        if (f(a) * f(b) < 0) {
            out << "\nInterval [" << a << ", " << b << "] seems to contain a root.\n";
            secant(a, b, r, out);
        }

        a = b;
    }

    out.close();
    return 0;
}
