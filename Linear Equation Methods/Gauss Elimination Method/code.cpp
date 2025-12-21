#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("input.txt");  
    ofstream fout("output.txt"); 

    if(!fin) {
        cout << "Error: input.txt not found!" << endl;
        return 0;
    }

    int n;
    fin >> n;  

    vector<vector<float>> a(n, vector<float>(n + 1));
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> a[i][j];

    
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            float factor = a[k][i] / a[i][i];
            for (int j = i; j <= n; j++) {
                a[k][j] -= factor * a[i][j];
            }
        }
    }

    
    vector<float> x(n, 1);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }

    
    fout << fixed << setprecision(2);
    fout << "Solution:\n";
    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << x[i] << "\n";

    fin.close();
    fout.close();

    return 0;
}
