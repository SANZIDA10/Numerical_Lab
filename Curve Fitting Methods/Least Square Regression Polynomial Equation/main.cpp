#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main() {
    ifstream fin("input5.txt");
    ofstream fout("output5.txt");

    int n, m;
    fin >> n;
    fin >> m;

    double x[20], y[20];

    for (int i = 0; i < n; i++)
        fin >> x[i];

    for (int i = 0; i < n; i++)
        fin >> y[i];

    double A[20][21] = {0};

    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= m; j++) {
            for (int k = 0; k < n; k++)
                A[i][j] += pow(x[k], i + j);
        }
        for (int k = 0; k < n; k++)
            A[i][m + 1] += y[k] * pow(x[k], i);
    }

    for (int i = 0; i <= m; i++) {
        for (int j = i + 1; j <= m; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = 0; k <= m + 1; k++)
                A[j][k] -= ratio * A[i][k];
        }
    }

    double a[20];
    for (int i = m; i >= 0; i--) {
        a[i] = A[i][m + 1];
        for (int j = i + 1; j <= m; j++)
            a[i] -= A[i][j] * a[j];
        a[i] /= A[i][i];
    }

    fout << fixed << setprecision(4);
    fout << "y = ";
    for (int i = 0; i <= m; i++) {
        fout << a[i];
        if (i > 0) fout << "x^" << i;
        if (i != m) fout << " + ";
    }
    fout << endl;

    fin.close();
    fout.close();

    return 0;
}
