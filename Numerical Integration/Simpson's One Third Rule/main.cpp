#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f(double x) {
    return sin(x);
}

int main() {
    ifstream fin("input1.txt");
    ofstream fout("output1.txt");

    double a, b;
    int n;

    fin >> a >> b >> n;

    if (n % 2 != 0) {
        fout << "Error: Number of intervals must be even for Simpson's 1/3 Rule." << endl;
        return 0;
    }

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0)
            sum += 2 * f(x);  
        else
            sum += 4 * f(x);   
    }

    double I = (h / 3.0) * sum;

    fout << "Value of integral : " << I << endl;

    fin.close();
    fout.close();

    return 0;
}
