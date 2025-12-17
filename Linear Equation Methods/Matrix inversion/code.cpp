#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

double det(double A[10][10], int n) {
  if (n == 1)
    return A[0][0];

  if (n == 2)
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];

  double temp[10][10];
  double d = 0;

  for (int p = 0; p < n; p++) {
    int h = 0, k = 0;

    for (int i = 1; i < n; i++) {
      k = 0;
      for (int j = 0; j < n; j++) {
        if (j == p)
          continue;
        temp[h][k++] = A[i][j];
      }
      h++;
    }

    d += pow(-1, p) * A[0][p] * det(temp, n - 1);
  }
  return d;
}

void adjoint(double A[10][10], double adj[10][10], int n) {
  if (n == 1) {
    adj[0][0] = 1;
    return;
  }

  double temp[10][10];

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {

      int h = 0, k = 0;
      for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
          if (r == i || c == j)
            continue;
          temp[h][k++] = A[r][c];
          if (k == n - 1) {
            h++;
            k = 0;
          }
        }
      }

      adj[j][i] = pow(-1, i + j) * det(temp, n - 1);
    }
  }
}

int main() {
  ifstream fin("D:\\Numerical project\\Linear Equation Methods\\Matrix "
               "inversion\\input.txt");
  ofstream fout("D:\\Numerical project\\Linear Equation Methods\\Matrix "
                "inversion\\output.txt");

  int n;
  fin >> n;

  double aug[10][11], A[10][10], B[10];

  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++)
      fin >> aug[i][j];

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      A[i][j] = aug[i][j];
    B[i] = aug[i][n];
  }

  double D = det(A, n);
  fout << "Determinant = " << D << endl;

  if (D == 0) {
    fout << "System has NO unique solution\n";
    return 0;
  }

  double adj[10][10], inv[10][10];
  adjoint(A, adj, n);

  fout << "\nInverse Matrix:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      inv[i][j] = adj[i][j] / D;
      fout << setw(10) << inv[i][j] << " ";
    }
    fout << endl;
  }

  double X[10] = {0};
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      X[i] += inv[i][j] * B[j];

  fout << "\nSolution Vector X:\n";
  for (int i = 0; i < n; i++)
    fout << "x" << i + 1 << " = " << X[i] << endl;

  fin.close();
  fout.close();
  return 0;
}
