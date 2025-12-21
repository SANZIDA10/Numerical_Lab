#include <bits/stdc++.h>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> B, ofstream &out) {
    int n = A.size();
    
    out << "\n--- Gaussian Elimination Steps ---\n";

    
    out << "\nInitial Augmented Matrix (A | B):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            out << setw(12) << A[i][j] << " ";
        out << " | " << setw(12) << B[i] << "\n";
    }
    out << "\n";


    for (int i = 0; i < n; i++) {

        double pivot = A[i][i];
        for (int j = 0; j < n; j++)
            A[i][j] /= pivot;
        B[i] /= pivot;

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = 0; j < n; j++)
                A[k][j] -= factor * A[i][j];
            B[k] -= factor * B[i];
        }

        
        out << "After step " << i+1 << ":\n";
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++)
                out << setw(12) << A[r][c] << " ";
            out << " | " << setw(12) << B[r] << "\n";
        }
        out << "\n";
    }

    
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
    }
    return x;
}

int main() {

    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in) {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int n, m;
    in >> n >> m;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++)
        in >> x[i] >> y[i];

    vector<double> S(2*m + 1, 0.0);
    vector<double> T(m + 1, 0.0);

    for (int i = 0; i < n; i++) {
        double xp = 1;
        for (int k = 0; k <= 2*m; k++) {
            S[k] += xp;
            xp *= x[i];
        }

        xp = 1;
        for (int k = 0; k <= m; k++) {
            T[k] += xp * y[i];
            xp *= x[i];
        }
    }

    vector<vector<double>> A(m+1, vector<double>(m+1));
    vector<double> B(m+1);

    for (int row = 0; row <= m; row++) {
        for (int col = 0; col <= m; col++)
            A[row][col] = S[row + col];

        B[row] = T[row];
    }

    out << fixed << setprecision(2);

    
    vector<double> a = gaussianElimination(A, B, out);

    out << "\nPolynomial Regression Equation:\n";
    out << "y = ";

    for (int i = 0; i <= m; i++) {
        out << a[i] << "x^" << i;
        if (i < m) out << " + ";
    }
    out << endl;

    in.close();
    out.close();
    return 0;
}
