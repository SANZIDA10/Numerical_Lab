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

double diff(double x) {
  double res = 0;
  for (int i = 0; i < degree; i++) {
    res += coeff[i] * (degree - i) * pow(x, degree - i - 1);
  }
  return res;
}

void newton(double x, double &r, ofstream &out) {
  int counts = 0;

  while (counts < max_it) {
    double fx = f(x);
    double fpx = diff(x);

    if (fabs(fpx) < 1e-12) {
      x += 0.5;
      continue;
    }

    r = x - fx / fpx;

    if (fabs(r - x) < E) {
      out << "The root = " << r << " found after " << counts
          << " iterations.\n\n";
      return;
    }

    x = r;
    counts++;
  }

  out << "No real root found starting from this initial guess.\n\n";
}

int main() {
  ifstream input("D:\\Numerical project\\Non-Linear Equation Methods\\Newton "
                 "Raphson method\\input.txt");
  ofstream out("D:\\Numerical project\\Non-Linear Equation Methods\\Newton "
               "Raphson method\\output.txt");

  if (!input) {
    cout << "Error: input.txt not found!\n";
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

  out << "Newton-Raphson Method\n";
  out << "Search interval: " << start << " to " << ends << "\n\n";

  double step = 0.5;
  double a = start, b, r;
  bool found = false;

  while (a < ends) {
    b = a + step;
    if (b > ends)
      b = ends;

    if (f(a) * f(b) < 0) {
      out << "Interval [" << a << ", " << b << "] contains a root.\n";
      newton(a, r, out);
      found = true;
    }

    a = b;
  }

  out.close();
  return 0;
}
