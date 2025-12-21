#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    int a[10][11], x[10];

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            int f = a[k][i] / a[i][i];
            for (int j = i; j <= n; j++)
                a[k][j] -= f * a[i][j];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }

    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << x[i] << endl;

    fout << "\nFinal Upper Triangular Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++)
            fout << setw(5) << a[i][j] << " ";
        fout << endl;
    }

    fin.close();
    fout.close();

    return 0;
}
