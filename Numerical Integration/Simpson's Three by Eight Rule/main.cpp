#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f(double x) {
    return sin(x);
}

int main() {
    ifstream fin("input2.txt");
    ofstream fout("output2.txt");

    double a, b;
    int n;

    fin >> a >> b >> n;

    if (n % 3 != 0) {
        fout << "Error: Number of intervals must be a multiple of 3." << endl;
        return 0;
    }

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 3 == 0)
            sum += 2 * f(x);
        else
            sum += 3 * f(x);
    }

    double I = (3 * h / 8.0) * sum;

    fout << "Value of integral : " << I << endl;

    fin.close();
    fout.close();

    return 0;
}
