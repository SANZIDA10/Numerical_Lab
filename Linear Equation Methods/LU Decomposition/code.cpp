#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int main() {
  ifstream fin("input.txt");
  ofstream fout("output.txt");

  int n;
  fin >> n;

  double A[20][20], L[20][20] = {0}, U[20][20] = {0};
  double B[20], Y[20], X[20];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      fin >> A[i][j];

  for (int i = 0; i < n; i++)
    fin >> B[i];

  for (int i = 0; i < n; i++) {

    for (int k = i; k < n; k++) {
      double sum = 0;
      for (int j = 0; j < i; j++)
        sum += L[i][j] * U[j][k];
      U[i][k] = A[i][k] - sum;
    }

    for (int k = i; k < n; k++) {
      if (i == k)
        L[i][i] = 1;
      else {
        double sum = 0;
        for (int j = 0; j < i; j++)
          sum += L[k][j] * U[j][i];
        L[k][i] = (A[k][i] - sum) / U[i][i];
      }
    }
  }

  fout << "L Matrix:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fout << setw(10) << L[i][j] << " ";
    fout << endl;
  }

  fout << "\nU Matrix:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      fout << setw(10) << U[i][j] << " ";
    fout << endl;
  }

  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < i; j++)
      sum += L[i][j] * Y[j];
    Y[i] = B[i] - sum;
  }

  for (int i = n - 1; i >= 0; i--) {
    double sum = 0;
    for (int j = i + 1; j < n; j++)
      sum += U[i][j] * X[j];
    X[i] = (Y[i] - sum) / U[i][i];
  }

  fout << "\nSolution Vector X:\n";
  for (int i = 0; i < n; i++)
    fout << X[i] << " ";

  fin.close();
  fout.close();

  return 0;
}
