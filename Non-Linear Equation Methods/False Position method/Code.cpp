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

double falsePosition(double a, double b, double E, ofstream &out) {
  double c;
  int count = 0;

  double fa = f(a);
  double fb = f(b);

  while (true) {
    c = a - fa * (b - a) / (fb - fa);
    count++;

    double fc = f(c);

    if (fabs(fc) < E) {
      out << "Iterations : " << count << "\n";
      return c;
    }

    if (fa * fc < 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }
}

int main() {
  ifstream input("D:\\Numerical project\\Non-Linear Equation Methods\\False "
                 "Position method\\input.txt");
  ofstream out("D:\\Numerical project\\Non-Linear Equation Methods\\False "
               "Position method\\output.txt");

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

  out << "False Position Method\n";
  out << "The search interval is: " << start << " to " << ends << "\n\n";

  double E = 1e-4;
  double step = 0.1;

  double a = start, b;
  while (a < ends) {
    b = a + step;
    if (b > ends)
      b = ends;

    if (fabs(f(a)) < E) {
      out << "Root found at x = " << a << "\n";
    }

    if (f(a) * f(b) < 0) {
      out << "Root lies between [" << a << ", " << b << "]\n";
      double root = falsePosition(a, b, E, out);
      out << "Root = " << root << "\n\n";
    }

    a = b;
  }

  out.close();
  return 0;
}
