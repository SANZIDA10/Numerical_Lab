#include <iostream>
#include <fstream>

using namespace std;

int main() {
    ifstream fin("input4.txt");
    ofstream fout("output4.txt");

    int n;
    fin >> n;

    float a[10][11];

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n; i++) {
        float p = a[i][i];
        for (int j = 0; j <= n; j++)
            a[i][j] /= p;

        for (int k = 0; k < n; k++) {
            if (k != i) {
                float f = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= f * a[i][j];
            }
        }
    }

    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << a[i][n] << endl;

    fin.close();
    fout.close();

    return 0;
}
