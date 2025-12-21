#include <bits/stdc++.h>
using namespace std;

const double E = 1e-6;
const int max_it = 1000;

int degree;
vector<double> coeff;

double f(double x) {
  double res = 0;
  for (int i = 0; i <= degree; i++) {
    res += coeff[i] * pow(x, degree - i);
  }
  return res;
}

void secant(double x, double x1, double &r, ofstream &out) {
  double fx = f(x);
  double fx1 = f(x1);
  int counts = 0;

  while (fabs(x1 - x) >= E && counts < max_it) {
    fx = f(x);
    fx1 = f(x1);

    if (fx1 - fx == 0) {
      out << "Division by zero encountered.\n";
      return;
    }

    r = x1 - (fx1 * (x1 - x)) / (fx1 - fx);
    x = x1;
    x1 = r;
    counts++;
  }

  if (counts >= max_it) {
    out << "No real root found in this interval.\n";
    return;
  }

  out << "The root = " << r << " found after " << counts << " iterations.\n";
}

int main() {
  ifstream input("D:\\Numerical project\\Non-Linear Equation Methods\\secant "
                 "method\\input.txt");
  ofstream out("D:\\Numerical project\\Non-Linear Equation Methods\\secant "
               "method\\output.txt");

  if (!input) {
    cout << "ERROR: input.txt not found!\n";
    return 0;
  }

  input >> degree;

  coeff.resize(degree + 1);
  for (int i = 0; i <= degree; i++) {
    input >> coeff[i];
  }

  double start, ends;
  input >> start >> ends;
  input.close();

  out << fixed << setprecision(6);

  out << "Equation:\n";
  for (int i = 0; i <= degree; i++) {
    if (coeff[i] == 0)
      continue;

    if (i != 0 && coeff[i] > 0)
      out << " + ";
    if (coeff[i] < 0)
      out << " - ";

    out << fabs(coeff[i]);
    int power = degree - i;
    if (power > 0)
      out << "x";
    if (power > 1)
      out << "^" << power;
  }
  out << "\n\n";

  out << "Secant Method\n";
  out << "Searching interval: [" << start << ", " << ends << "]\n";

  double step = 0.5;
  double a = start, b, r;

  while (a < ends) {
    b = a + step;
    if (b > ends)
      b = ends;

    if (f(a) * f(b) < 0) {
      out << "\nInterval [" << a << ", " << b << "] contains a root.\n";
      secant(a, b, r, out);
    }

    a = b;
  }

  out.close();
  return 0;
}
