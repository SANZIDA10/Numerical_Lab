#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

double f(double x, double y) {
    return x + y * y;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double x0, y0, xn, h;

    fin >> x0;
    fin >> y0;
    fin >> xn;
    fin >> h;

    double x = x0;
    double y = y0;
    int iteration = 0;

    while (x < xn) {
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2, y + k1/2);
        double k3 = h * f(x + h/2, y + k2/2);
        double k4 = h * f(x + h, y + k3);

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        x = x + h;
        iteration++;
    }

    fout << fixed << setprecision(5);
    fout << "Value of y at x = " << xn << " is: " << y << endl;
    fout << "Total iterations = " << iteration << endl;

    fin.close();
    fout.close();

    return 0;
}
