#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    int a[10][11];

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n; i++) {
        int p = a[i][i];
        for (int j = 0; j <= n; j++)
            a[i][j] /= p; 

        for (int k = 0; k < n; k++) {
            if (k != i) {
                int f = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= f * a[i][j]; 
            }
        }
    }

    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << a[i][n] << endl;

    fout << "\nFinal Reduced Matrix (RREF):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++)
            fout << setw(5) << a[i][j] << " ";
        fout << endl;
    }

    fin.close();
    fout.close();

    return 0;
}
