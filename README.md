# Numerical Methods Laboratory

A comprehensive collection of numerical methods implementations in C++ for solving linear and non-linear equations, interpolation, differentiation, curve fitting, and numerical integration.

---

# Table of Contents
- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton-Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion Method](#matrix-inversion-method)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)
- [Differential Equations](#differential-equations)
  - [Runge-Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)
- [Interpolation](#interpolation)
  - [Newton's Forward Interpolation](#newtons-forward-interpolation)
    - [Theory](#newtons-forward-theory)
    - [Code](#newtons-forward-code)
    - [Input](#newtons-forward-input)
    - [Output](#newtons-forward-output)
  - [Newton's Backward Interpolation](#newtons-backward-interpolation)
    - [Theory](#newtons-backward-theory)
    - [Code](#newtons-backward-code)
    - [Input](#newtons-backward-input)
    - [Output](#newtons-backward-output)
  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)
- [Numerical Differentiation](#numerical-differentiation)
  - [Newton's Forward Differentiation](#newtons-forward-differentiation)
    - [Theory](#newtons-forward-differentiation-theory)
    - [Code](#newtons-forward-differentiation-code)
    - [Input](#newtons-forward-differentiation-input)
    - [Output](#newtons-forward-differentiation-output)
  - [Newton's Backward Differentiation](#newtons-backward-differentiation)
    - [Theory](#newtons-backward-differentiation-theory)
    - [Code](#newtons-backward-differentiation-code)
    - [Input](#newtons-backward-differentiation-input)
    - [Output](#newtons-backward-differentiation-output)
- [Curve Fitting](#curve-fitting)
  - [Linear Equation](#curve-fitting-linear)
    - [Theory](#linear-theory)
    - [Code](#linear-code)
    - [Input](#linear-input)
    - [Output](#linear-output)
  - [Transcendental Equation](#curve-fitting-transcendental)
    - [Theory](#transcendental-theory)
    - [Code](#transcendental-code)
    - [Input](#transcendental-input)
    - [Output](#transcendental-output)
  - [Polynomial Equation](#curve-fitting-polynomial)
    - [Theory](#polynomial-theory)
    - [Code](#polynomial-code)
    - [Input](#polynomial-input)
    - [Output](#polynomial-output)
- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpsons-13-theory)
    - [Code](#simpsons-13-code)
    - [Input](#simpsons-13-input)
    - [Output](#simpsons-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpsons-38-theory)
    - [Code](#simpsons-38-code)
    - [Input](#simpsons-38-input)
    - [Output](#simpsons-38-output)

---

## Solution of Non-Linear Equations

A non-linear equation is one in which the power of the variable is greater than 1 or the variables appear in non-linear forms (squared, cubed, exponential, trigonometric, logarithmic, etc.).

**Examples:**
- e^x - 3x = 0
- x^2 - 4x - 10 = 0
- sin(x) = x/2

---

### Bisection Method

#### Bisection Theory
Bisection Method is a numerical method used to find a root of a continuous function. It works by repeatedly dividing an interval [a,b][a, b][a,b] into two halves where the function values at the ends have opposite signs (f(a)⋅f(b)<0)(f(a)\cdot f(b) < 0)(f(a)⋅f(b)<0), and then selecting the half in which the root lies. This process continues until the root is found with the desired accuracy.

- Binary chopping or half-interval method

- If f(x) is real and continuous in the interval a < x < b, and f(a) and f(b) are of opposite sign, that is, f(a) f(b) < 0, then there is at least one real root in the interval between a and b.

- Let x1 = a and x2 = b. Define x0 to be the midpoint between a and b, that is,

X0 =(x1 + x2)/2

- Now there exist the following 3 conditions:
- (i) If f(x0) = 0, then the root is x0.
- (ii) If f(x0) * f(x1) < 0, then root is between x0 and x1.
- (iii) If f(x0) * f(x2) < 0, then root is between x0 and x2.

#### Bisection Code
```cpp
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

```

#### Bisection Input
```
4
1 -3 2 6 0
-5 5


```

#### Bisection Output
```
The equation is:
1.00x^4 - 3.00x^3 + 2.00x^2 + 6.00x

The search interval is: -5.00 to 5.00

The root is at x = -1.00
Iterations : 14
Root lies between -1.00 and -0.90
Root = -1.00

The root is at x = -0.00
Iterations : 13
Root lies between -0.00 and 0.10
Root = 0.00


```

---

### False Position Method

#### False Position Theory

False Position Method (Regula Falsi) is a numerical method used to find a root of a continuous function.

- Linear interpolation method

- If f(x) is real and continuous in the interval a < x < b, and f(a) and f(b) are of opposite sign, that is, f(a)f(b) < 0, then there is at least one real root in the interval between a and b.

- Let x1 = a and x2 = b. Define x0 to be the midpoint between a and b, that is,

X0 = x1 - f(x1) * (x2 – x1)/ (f(x2) - f(x1))

- Now there exist the following 3 conditions:
- (i) If f(x0) = 0, then the root is x0.
- (ii) If f(x0) * f(x1) < 0, then root is between x0 and x1.
- (iii) If f(x0) * f(x2) < 0, then root is between x0 and x2.

#### False Position Code
```cpp
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

```

#### False Position Input
```
4
1 -3 2 6 0
-2.24 2.24

```

#### False Position Output
```
Equation:
1.000000x^4 - 3.000000x^3 + 2.000000x^2 + 6.000000x

False Position Method
The search interval is: -2.240000 to 2.240000

Root lies between [-1.040000, -0.940000]
Iterations : 4
Root = -0.999999

Root lies between [-0.040000, 0.060000]
Iterations : 2
Root = -0.000014


```

---

### Newton-Raphson Method

#### Newton-Raphson Theory
Newton–Raphson Method is a numerical method used to find a root of a real-valued function. It starts with an initial guess and repeatedly improves it using the function and its derivative.

- If f(x) is a real and continuously differentiable function, the Newton–Raphson method is used to find a root of the equation f(x) = 0.

- Let the initial approximation to the root be x0.

- A better approximation x1 is obtained using the formula: x1 = x0 – (f(x0)/f'(x0)).

- This process is repeated iteratively as:

X(n+1) = xn – (f(xn)/f'(xn))

- The iteration continues until the difference between two successive values of x becomes sufficiently small.

#### Newton-Raphson Code
```cpp
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

```

#### Newton-Raphson Input
```
3
1 0 -1 0
-2.23 2.23

```

#### Newton-Raphson Output
```
Equation:
1.000000x^3 - 1.000000x

Newton-Raphson Method
Search interval: -2.230000 to 2.230000

Interval [-1.230000, -0.730000] contains a root.
The root = -1.000000 found after 4 iterations.

Interval [-0.230000, 0.270000] contains a root.
The root = 0.000000 found after 3 iterations.

Interval [0.770000, 1.270000] contains a root.
The root = 1.000000 found after 5 iterations.


```

---

### Secant Method

#### Secant Theory

Secant method is a recursive method for finding the root of a polynomial by successive approximation. In this method, the neighborhood roots are approximated by secant line (A line that intersects the curve at two distinct points) to the function f(x).

Formula for Secant Method: xk + 1 = xk – (((xk – 1) – xk)/ (f(xk-1) – f(xk))*f(xk)

The secant method will find a root of a function f if the starting points x0 and x1 are close enough to the actual root and if the function behaves well. If the function is smooth and has a simple root (i.e., the root only occurs once), the method converges to the root at a rate related to the golden ratio: Φ = 1.618

This means the method improves quickly but not as fast as some other methods.

#### Secant Code
```cpp
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

```

#### Secant Input
```
3
3 0 -1 -1
0 2

```

#### Secant Output
```
Equation:
3.000000x^3 - 1.000000x - 1.000000

Secant Method
Searching interval: [0.000000, 2.000000]

Interval [0.500000, 1.000000] contains a root.
The root = 0.851383 found after 6 iterations.

```

---

## Solution of Linear Equations

A linear equation is an equation in which the highest power (degree) of the variable is 1.

**Example:**
- 3x - 2 = 0

---

### Gauss Elimination Method

#### Gauss Elimination Theory

The Gauss Elimination Method is used to solve a system of linear equations by converting it into an upper triangular matrix and then solving it by back substitution.

Steps:
1. Write the augmented matrix of the system of equations.
2. Forward Elimination:
   - Eliminate elements below the main diagonal using row operations.
   - Convert the matrix into an upper triangular form.
3. Back Substitution:
   - Start from the last row and calculate the value of the last variable.
   - Substitute this value into the rows above to find the remaining variables.
4. Result:
   - All variables are calculated step by step.

#### Gauss Elimination Code
```cpp
#include <iostream>
#include <fstream>

using namespace std;

int main() {
    ifstream fin("input.txt");     
    ofstream fout("output.txt");   

    int n;
    fin >> n;

    float a[10][11], x[10];

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            float f = a[k][i] / a[i][i];
            for (int j = i; j <= n; j++)
                a[k][j] -= f * a[i][j];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }

    for (int i = 0; i < n; i++)
        fout << "x" << i + 1 << " = " << x[i] << endl;

    fin.close();
    fout.close();

    return 0;
}
```

#### Gauss Elimination Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

#### Gauss Elimination Output
```
x1 = 2
x2 = 3
x3 = -1
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

The Gauss–Jordan Elimination Method is used to solve a system of linear equations by reducing the augmented matrix to diagonal form (or reduced row echelon form), so that the solution can be read directly without back substitution.

Steps:
1. Write the augmented matrix of the system of equations.
2. Forward Elimination:
   - Eliminate elements below the main diagonal to form an upper triangular matrix.
3. Backward Elimination:
   - Eliminate elements above the main diagonal to form a diagonal matrix.
4. Normalization:
   - Divide each row by its pivot element so that all diagonal elements become 1.
5. Result:
   - The solution of the system is directly obtained from the last column of the matrix.

#### Gauss Jordan Code
```cpp
#include <iostream>
#include <fstream>

using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

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
```

#### Gauss Jordan Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

#### Gauss Jordan Output
```
x1 = 2
x2 = 3
x3 = -1
```

---

### LU Decomposition Method

#### LU Decomposition Theory
LU Decomposition Method is a matrix factorization technique in which a square matrix A is expressed as the product of a lower triangular matrix L and an upper triangular matrix U, such that
                A=LU.
It is mainly used to solve systems of linear equations efficiently by first decomposing the matrix and then applying forward and backward substitution.

To factor any square matrix into two triangular matrices, i.e., one is a lower triangular matrix and the other is an upper triangular matrix, we can use the following steps.

**Steps for LU Decomposition:**

Start with a square matrix A: Given a square matrix A of size n ×n, the goal is to factor it into the product of two matrices: A = L×U, where:

- L is a lower triangular matrix with 1s on the diagonal.

- U is an upper triangular matrix.

Apply Gaussian elimination to convert matrix A into upper triangular form U. This step involves row operations to eliminate elements below the diagonal, resulting in an upper triangular matrix.

As we perform row operations, keep track of the multipliers used to eliminate the elements below the diagonal. These multipliers form the entries of the lower triangular matrix L.

- The entries of L will be the factors used during the elimination steps.

- The diagonal of L will consist of 1s.

#### LU Decomposition Code
```cpp
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
```

#### LU Decomposition Input
```
3
2 -1 -2
-4 6 3
-4 -2 8
-2
9
-5
```

#### LU Decomposition Output
```
L Matrix:
         1          0          0 
        -2          1          0 
        -2         -1          1 

U Matrix:
         2         -1         -2 
         0          4         -1 
         0          0          3 

Solution Vector X:
-1.875 0.916667 -1.33333 
```

---

### Matrix Inversion Method

#### Matrix Inversion Theory
Matrix Inversion is the process of finding a matrix A^-1 for a given square matrix A such that
                  AA^-1 = A^-1A= 1
where I is the identity matrix. The inverse exists only if the matrix is non-singular (has a non-zero determinant).

AX = B

A^(-1)AX = A^(-1)B

IX = A^(-1)B

X = A^(-1)B

**Steps to find the inverse of a matrix:**

1. Find the determinant.
2. Find the cofactor matrix.
3. Find the adjoint (Transpose of cofactor matrix).
4. Compute the inverse: A^(-1) = (1 / (|A|)) * adj(A).

#### Matrix Inversion Code
```cpp
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
  ifstream fin("input.txt");
  ofstream fout("output.txt");

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
```

#### Matrix Inversion Input
```
3
2 -1 -2 -2
-4 6 3 9
-4 -2 8 -5
```

#### Matrix Inversion Output
```
Determinant = 24

Inverse Matrix:
      2.25        0.5      0.375 
  0.833333   0.333333  0.0833333 
   1.33333   0.333333   0.333333 

Solution Vector X:
x1 = -1.875
x2 = 0.916667
x3 = -1.33333
```

---

## Differential Equations

---

### Runge-Kutta Method

#### Runge-Kutta Theory

The Runge-Kutta 4th Order Method (RK4) is a numerical technique to solve ordinary differential equations (ODEs) of the form:
    dy / dx = f(x , y) ,  y(x0) = y0
It is more accurate than simple methods because it estimates the slope at multiple points within a step.

Steps:
1. Start with initial values:
   x = x0 ,   y = y0
   and choose a step size h.
2. Calculate intermediate slopes:
   k1 = h ⋅ f( x, y )
   k2 = h ⋅ f( x + h / 2 , y + k1 / 2)
   k3 = h ⋅ f( x + h / 2 , y + k2 / 2)
   k4 = h ⋅ f( x + h , y + k3 ) 
3. Update the solution:
   y(next) = y + ( k1 + 2*k2 + 2*k3 + k4 ) / 6
   x(next) = x + h 
4. Repeat the process until reaching the desired value of x = xn .

#### Runge-Kutta Code
```cpp
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

double f(double x, double y) {
    return x + y * y;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double x0, y0, xn, h;

    fin >> x0;
    fin >> y0;
    fin >> xn;
    fin >> h;

    double x = x0;
    double y = y0;
    int iteration = 0;

    while (x < xn) {
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2, y + k1/2);
        double k3 = h * f(x + h/2, y + k2/2);
        double k4 = h * f(x + h, y + k3);

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        x = x + h;
        iteration++;
    }

    fout << fixed << setprecision(5);
    fout << "Value of y at x = " << xn << " is: " << y << endl;
    fout << "Total iterations = " << iteration << endl;

    fin.close();
    fout.close();

    return 0;
}
```

#### Runge-Kutta Input
```
0
1
0.2
0.1
```

#### Runge-Kutta Output
```
Value of y at x = 0.20000 is: 1.27356
Total iterations = 2
```

---

## Interpolation

Interpolation is the process of estimating the value of a function at a point that lies between two known data points. It uses known values of a function to construct a new function (often a polynomial) that approximates the original function and allows us to compute missing values inside the data interval.

---

### Newton's Forward Interpolation

#### Newtons Forward Theory

- Newton's forward interpolation method is used to estimate values near the beginning of a table when data points (x is data point) are equally spaced.

- Let x₀, x₁, x₂, … , xₙ₋₁, xₙ be a set of equally spaced values of the independent variable x.

- So x₁ − x₀ = x₂ − x₁ = x₃ − x₂ = … = xₙ − xₙ₋₁ = h

- Let, u = (x − x₀) / h

- The Newton's forward difference interpolation formula for equal intervals is:
    y = y₀ + u Δy₀ + [u(u − 1) / 2!] Δ²y₀+ [u(u − 1)(u − 2) / 3!] Δ³y₀+ … + [u(u − 1)(u − 2) … (u − n - 1) / n!] Δⁿy₀

#### Newtons Forward Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double facto(int n){
    double f=1;
    for(int i=2;i<=n;i++)f*=i;
    return f;
}

int main(){
    int n;cin>>n;
    vector<double>x(n),y(n);

    for(int i=0;i<n;i++)
        cin>>x[i]>>y[i];

    double x_val;cin>>x_val;

    vector<vector<double>>diff(n,vector<double>(n));

    for(int i=0;i<n;i++)
        diff[i][0]=y[i];

    for(int j=1;j<n;j++){
        for(int i=0;i<n-j;i++){
            diff[i][j]=diff[i+1][j-1]-diff[i][j-1];
            diff[i][j]=round(diff[i][j]*1e6)/1e6;
        }
    }
    cout<<"Forward Difference Table:\n";
    for(int i=0;i<n;i++){
        cout<<setw(6)<<x[i];
        for(int j=0;j<n-i;j++)
            cout<<setw(10)<<diff[i][j];
        cout<<endl;
    }
    double h=x[1]-x[0],u=(x_val-x[0])/h,ans=diff[0][0],u_prod=1;
    for(int i=1;i<n;i++){
        u_prod*=(u-(i-1));
        ans+=(u_prod*diff[0][i])/facto(i);
    }
    cout<<"\nInterpolated Value at x="<<x_val<<" is: "<<ans<<endl;
    return 0;
}
```

#### Newtons Forward Input
```
5
1 1
2 8
3 27
4 64
5 125
2.5
```

#### Newtons Forward Output
```
Forward Difference Table:
     1         1         7        12         6         0
     2         8        19        18         6
     3        27        37        24
     4        64        61
     5       125

Interpolated Value at x=2.5 is: 15.625
```

---

### Newton's Backward Interpolation

#### Newtons Backward Theory

- Backward interpolation is used when the value of x lies near the end of the table.

- Let x₀, x₁, x₂, … , xₙ₋₁, xₙ be a set of equally spaced values of the independent variable x.

- The interval between successive values is constant, i.e., x₁ − x₀ = x₂ − x₁ = x₃ − x₂ = … = xₙ − xₙ₋₁ = h

- Let, u = (x − xₙ) / h

- The Newton's backward difference interpolation formula for equal intervals is:
  y = yₙ + u∇yₙ+ [u(u + 1) / 2!] ∇²yₙ+ [u(u + 1)(u + 2) / 3!] ∇³yₙ+ … + [u(u + 1)(u + 2)…(u + n − 1) / n!] ∇ⁿyₙ

#### Newtons Backward Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double facto(int n){
    double f=1;
    for(int i=2;i<=n;i++)f*=i;
    return f;
}

int main(){
    int n;cin>>n;
    vector<double>x(n),y(n);
    for(int i=0;i<n;i++)cin>>x[i]>>y[i];
    double x_val;cin>>x_val;

    vector<vector<double>>diff(n,vector<double>(n));
    for(int i=0;i<n;i++)diff[i][0]=y[i];

    for(int j=1;j<n;j++)
        for(int i=n-1;i>=j;i--){
            diff[i][j]=diff[i][j-1]-diff[i-1][j-1];
            diff[i][j]=round(diff[i][j]*1e6)/1e6;
        }

    cout<<"Backward Difference Table:\n";
    for(int i=0;i<n;i++){
        cout<<setw(6)<<x[i];
        for(int j=0;j<=i;j++)cout<<setw(10)<<diff[i][j];
        cout<<endl;
    }

    double h=x[1]-x[0];
    double u=(x_val-x[n-1])/h;
    double ans=diff[n-1][0],u_prod=1;
    for(int i=1;i<n;i++){
        u_prod*=(u+(i-1));
        ans+=(u_prod*diff[n-1][i])/facto(i);
    }

    cout<<"\nInterpolated Value at x="<<x_val<<" is: "<<ans<<endl;
    return 0;
}
```

#### Newtons Backward Input
```
5
1 1
2 8
3 27
4 64
5 125
4.5
```

#### Newtons Backward Output
```
Backward Difference Table:
     1         1
     2         8         7
     3        27        19        12
     4        64        37        18         6
     5       125        61        24         6         0

Interpolated Value at x=4.5 is: 91.125
```

---

### Divided Difference Method

#### Divided Difference Theory

- The Newton's Divided Difference Interpolation Method is used to estimate the value of a function when the given data points are not equally spaced.

- The interpolation polynomial is given by:

f(xₙ)= f(x₀) + (x - x₀) f [x₁, x₀] + (x - x₀)(x - x₁) f[x₂, x₁, x₀] +...+ (x - x₀)(x - x₁)...(x − xₙ₋₁) f[xₙ, xₙ₋₁,..., x₁, x₀]...(1)

Here divided differences are defined as:

- First order divided difference:
f [ xᵢ, xⱼ ] = [ f(xᵢ) - f(xⱼ) ] / xᵢ - xⱼ …..(2)

- Second order divided difference:
f [ xᵢ, xⱼ, xₖ ] = [ f(xᵢ, xⱼ) - f(xⱼ, xₖ) ] / xᵢ - xₖ …..(3)

#### Divided Difference Code
```cpp
#include <iostream>
#include <iomanip>
using namespace std;

int main(){
    int n;
    cin>>n;

    double x[50],y[50][50],X;

    for(int i=0;i<n;i++)
        cin>>x[i]>>y[i][0];

    cin>>X;

    for(int j=1;j<n;j++)
        for(int i=0;i<n-j;i++)
            y[i][j]=(y[i+1][j-1]-y[i][j-1])/(x[i+j]-x[i]);

    cout<<"Divided Difference Table:\n";
    for(int i=0;i<n;i++){
        cout<<setw(6)<<x[i];
        for(int j=0;j<n-i;j++)
            cout<<setw(12)<<y[i][j];
        cout<<endl;
    }

    double ans=y[0][0],term=1;
    for(int i=1;i<n;i++){
        term*=(X-x[i-1]);
        ans+=term*y[0][i];
    }

    cout<<"\nInterpolated value: "<<ans<<endl;
    return 0;
}
```

#### Divided Difference Input
```
4
1 1
2 4
4 16
7 49
3
```

#### Divided Difference Output
```
Divided Difference Table:
     1           1           3           1           0
     2           4           6           1
     4          16          11
     7          49

Interpolated value: 9
```

---

## Numerical Differentiation

---

### Newton's Forward Differentiation

#### Newtons Forward Differentiation Theory

[Add your Newton's Forward Differentiation theory here]

#### Newtons Forward Differentiation Code
```cpp
// Add your Newton's Forward Differentiation code here
```

#### Newtons Forward Differentiation Input
```
Add your Newton's Forward Differentiation input here
```

#### Newtons Forward Differentiation Output
```
Add your Newton's Forward Differentiation output here
```

---

### Newton's Backward Differentiation

#### Newtons Backward Differentiation Theory

[Add your Newton's Backward Differentiation theory here]

#### Newtons Backward Differentiation Code
```cpp
// Add your Newton's Backward Differentiation code here
```

#### Newtons Backward Differentiation Input
```
Add your Newton's Backward Differentiation input here
```

#### Newtons Backward Differentiation Output
```
Add your Newton's Backward Differentiation output here
```

---

## Curve Fitting

---

### Curve Fitting: Linear Equation

#### Linear Theory

[Add your Linear Curve Fitting theory here]

#### Linear Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin>>n;

    vector<double>x(n),y(n);
    for(int i=0;i<n;i++)
        cin>>x[i]>>y[i];

    double sumx =0,sumy=0,sumxy=0, sumx2=0;

    for (int i=0;i<n;i++) {
        sumx  += x[i];
        sumy  += y[i];
        sumxy += x[i]*y[i];
        sumx2 += x[i]*x[i];
    }

    double b=(n*sumxy-sumx* sumy)/(n*sumx2-sumx*sumx);
    double a=(sumy-b*sumx)/n;

    cout<<fixed<<setprecision(4);
    cout<<"a(intercept) = "<<a<<"\n";
    cout<<"b(slope) = "<<b<<"\n\n";

    cout<<"Predicted y values:\n";
    for(int i =0; i<n; i++) {
        double y_calc=a+b*x[i];
        cout<<"x = "<<x[i]<<", y = "<<y_calc<<"\n";
    }
    return 0;
}
```

#### Linear Input
```
5
1 2
2 3
3 5
4 4
5 6
```

#### Linear Output
```
a (intercept) = 1.3000
b (slope) = 0.9000

Predicted y values:
x = 1 , y = 2.2000
x = 2 , y = 3.1000
x = 3 , y = 4.0000
x = 4 , y = 4.9000
x = 5 , y = 5.8000
```

---

### Curve Fitting: Transcendental Equation

#### Transcendental Theory

[Add your Transcendental Curve Fitting theory here]

#### Transcendental Code
```cpp
#include <bits/stdc++.h>//y=ae^bx
using namespace std;

int main() {

    int n;
    cin>>n;

    vector<double>x(n),y(n);
    for (int i=0;i<n;i++)
        cin>>x[i]>>y[i];

    double sumx =0, sumY =0, sumxY =0, sumx2 =0;

    for (int i=0; i <n; i++) {
        double Y = log(y[i]);   // Y = ln(y)
        sumx  +=x[i];
        sumY  +=Y;
        sumxY +=x[i]*Y;
        sumx2 +=x[i]*x[i];
    }

    double B =(n * sumxY -sumx *sumY)/(n * sumx2 -sumx *sumx);
    double A =(sumY - B * sumx)/n;

    double a =exp(A);
    double b =B;

    cout<<fixed<<setprecision(6);
    cout<<"a = "<<a<<"\n";
    cout<<"b = "<<b<<"\n";
    cout<<"Model: y = "<<a<<" * e^("<<b<<"x)\n\n";

    cout<<"Calculated y values:\n";
    for (int i= 0;i<n;i++) {
        double y_calc=a*exp(b* x[i]);
        cout<<"x = "<< x[i]<<", y = "<<y_calc<<"\n";
    }
    return 0;
}
```

#### Transcendental Input
```
5
1 2.7
2 4.1
3 6.2
4 9.1
5 13.5
```

#### Transcendental Output
```
a = 1.827779
b = 0.401616
Model: y = 1.827779 * e^(0.401616x)

Calculated y values:
x = 1, y = 2.731136
x = 2, y = 4.080968
x = 3, y = 6.097936
x = 4, y = 9.111768
x = 5, y = 13.615149
```

---

### Curve Fitting: Polynomial Equation

#### Polynomial Theory

Curve fitting is the process of finding a mathematical function that best approximates a set of data points. In polynomial curve fitting, we fit a polynomial of degree m to n data points (xi, yi).

Steps:
1. Collect data points:
   - n pairs of (xi, yi).
2. Choose the degree m of the polynomial:
   - The polynomial has the form:
     y = a0 + a1x + a2x^2 + ⋯ + amx^m
3. Form the normal equations:
   - Use the least squares method to minimize the sum of squared differences between actual and predicted values:
     ∑(yi – polynomial)^2
   - This gives a system of linear equations in the unknown coefficients a0, a1 ,..., am.
4. Solve the system using Gaussian Elimination:
   - Find the coefficients a0, a1, ..., am.
5. Result:
   - The polynomial function y = a0 + a1x +...+ amx^m approximates the given data.

#### Polynomial Code
```cpp
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n, m;
    fin >> n;
    fin >> m;

    double x[20], y[20];

    for (int i = 0; i < n; i++)
        fin >> x[i];

    for (int i = 0; i < n; i++)
        fin >> y[i];

    double A[20][21] = {0};

    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= m; j++) {
            for (int k = 0; k < n; k++)
                A[i][j] += pow(x[k], i + j);
        }
        for (int k = 0; k < n; k++)
            A[i][m + 1] += y[k] * pow(x[k], i);
    }

    for (int i = 0; i <= m; i++) {
        for (int j = i + 1; j <= m; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = 0; k <= m + 1; k++)
                A[j][k] -= ratio * A[i][k];
        }
    }

    double a[20];
    for (int i = m; i >= 0; i--) {
        a[i] = A[i][m + 1];
        for (int j = i + 1; j <= m; j++)
            a[i] -= A[i][j] * a[j];
        a[i] /= A[i][i];
    }

    fout << fixed << setprecision(4);
    fout << "y = ";
    for (int i = 0; i <= m; i++) {
        fout << a[i];
        if (i > 0) fout << "x^" << i;
        if (i != m) fout << " + ";
    }
    fout << endl;

    fin.close();
    fout.close();

    return 0;
}
```

#### Polynomial Input
```
5
2
1 2 3 4 5
2 4 5 4 5
```

#### Polynomial Output
```
y = 0.3714 + 1.3143x^1 + 0.0857x^2
```

---

## Numerical Integration

---

### Simpson's 1/3 Rule

#### Simpsons 13 Theory

Simpson's 1/3 Rule is a numerical method used to approximate the definite integral of a function when the exact integration is difficult. It approximates the area under a curve using parabolic (quadratic) segments.

Steps:
1. Divide the interval [a,b] into n sub-intervals, where n must be even.
2. Calculate the step size:
   h = (b−a) / n
3. Apply the formula:
   - Multiply odd-indexed points by 4 and even-indexed points by 2.
4. Sum up the results to get the approximate value of the integral.

#### Simpsons 13 Code
```cpp
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f(double x) {
    return sin(x);
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double a, b;
    int n;

    fin >> a >> b >> n;

    if (n % 2 != 0) {
        fout << "Error: Number of intervals must be even for Simpson's 1/3 Rule." << endl;
        return 0;
    }

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0)
            sum += 2 * f(x);  
        else
            sum += 4 * f(x);   
    }

    double I = (h / 3.0) * sum;

    fout << "Value of integral : " << I << endl;

    fin.close();
    fout.close();

    return 0;
}
```

#### Simpsons 13 Input
```
0
3.1416
6
```

#### Simpsons 13 Output
```
Value of integral : 2.00086
```

---

### Simpson's 3/8 Rule

#### Simpsons 38 Theory

Simpson's 3/8 Rule is a numerical method used to approximate the definite integral of a function when finding the exact integral is difficult. It is a type of formula that approximates the area under a curve using cubic polynomials.

Steps:
1. Divide the interval [a,b] into n sub-intervals, where n must be a multiple of 3.
2. Calculate the step size:
   h=(b−a)/n
3. Apply the formula:
   - Multiply function values at every third point by 2.
   - Multiply function values at other points by 3.
4. Sum up the results to get the approximate value of the integral.

#### Simpsons 38 Code
```cpp
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double f(double x) {
    return sin(x);
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double a, b;
    int n;

    fin >> a >> b >> n;

    if (n % 3 != 0) {
        fout << "Error: Number of intervals must be a multiple of 3." << endl;
        return 0;
    }

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 3 == 0)
            sum += 2 * f(x);
        else
            sum += 3 * f(x);
    }

    double I = (3 * h / 8.0) * sum;

    fout << "Value of integral : " << I << endl;

    fin.close();
    fout.close();

    return 0;
}
```

#### Simpsons 38 Input
```
0
3.1416
6
```

#### Simpsons 38 Output
```
Value of integral : 2.00201
```

---

