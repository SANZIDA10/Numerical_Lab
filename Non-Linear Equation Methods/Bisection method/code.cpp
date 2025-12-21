#include <bits/stdc++.h>
using namespace std;

int degree;
vector<double> coeff;

double f(double x) {
  double res = 0;
  for (int i = 0; i <= degree; i++) {
    res += coeff[i] * pow(x, degree - i);
  }
  return res;
}

double bisection(double a, double b, double E, ofstream &out) {
  double mid;
  int count = 0;

  double fa = f(a);

  while (true) {
    mid = (a + b) / 2.0;
    count++;

    double fmid = f(mid);

    if (fabs(fmid) < E) {
      out << "Iterations : " << count << "\n";
      return mid;
    }

    if (fa * fmid < 0) {
      b = mid;
    } else {
      a = mid;
      fa = fmid;
    }
  }
}

int main() {
  ifstream input("D:\\Numerical project\\Non-Linear Equation "
                 "Methods\\Bisection method\\input.txt");
  ofstream out("D:\\Numerical project\\Non-Linear Equation Methods\\Bisection "
               "method\\output.txt");

  input >> degree;

  coeff.resize(degree + 1);
  for (int i = 0; i <= degree; i++) {
    input >> coeff[i];
  }

  double start, ends;
  input >> start >> ends;

  out << fixed << setprecision(6);

  out << "The equation is:\n";
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

  out << "The search interval is: " << start << " to " << ends << "\n\n";

  double E = 1e-4;
  double step = 0.1;

  double a = start;
  while (a < ends) {
    double b = a + step;
    if (b > ends)
      b = ends;

    if (fabs(f(a)) < E) {
      out << "The root is at x = " << a << "\n";
    }

    if (f(a) * f(b) < 0) {
      double root = bisection(a, b, E, out);
      out << "Root lies between " << a << " and " << b << "\nRoot = " << root
          << "\n\n";
    }

    a = b;
  }

  input.close();
  out.close();
  return 0;
}
