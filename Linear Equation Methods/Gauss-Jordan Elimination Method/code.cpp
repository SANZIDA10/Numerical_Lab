#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    double a[100][101]; 


    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            fin >> a[i][j];
        }
    }

    
    for (int i = 0; i < n; i++) {

        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[i][i])) {
                for (int j = 0; j <= n; j++)
                    swap(a[i][j], a[k][j]);
            }
        }


        double pivot = a[i][i];
        for (int j = 0; j <= n; j++)
            a[i][j] /= pivot;

        
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }
        }
    }

    
    fout << "Solution:\n";
    fout << fixed << setprecision(6);
    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << a[i][n] << "\n";

    fin.close();
    fout.close();

    return 0;
}
