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
  - [Theory](#numerical-differentiation-theory)
  - [Code](#numerical-differentiation-code)
  - [Input](#numerical-differentiation-input)
  - [Output](#numerical-differentiation-output)
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

double f(double x){
    return pow(x, 4) - 3*(pow(x, 3)) + 2*x*x + 6*x;
}

double bisection(double a, double b, double E, ofstream &out){
    double mid;
    int count = 0;

    double fa = f(a);
    double fb = f(b);

    while(true){
        mid = (a + b)/2;
        count++;
        
        double fmid = f(mid);

        if(fabs(fmid) < E){
            out << "iterations : " << count << "\n";
            return mid;
        }

        if(fa * fmid < 0){
            b = mid;
            fb = fmid;
        }
        else{
            a = mid;
            fa = fmid;  
        }
    }

    out << "iterations : " << count << "\n";
    return mid;
}

int main(){
    ifstream input("input.txt");
    ofstream out("output.txt");

    double start, ends;
    input >> start >> ends;

    out << fixed << setprecision(6);
    out << "The search interval is: " << start << " to " << ends << "\n\n";

    double E = 1e-4;
    double a = start, b;
    double step = 0.1;  

    while(a < ends){
        b = a + step;
        if(b > ends) b = ends;

        if(fabs(f(a)) < E){
            out << "The root is at = " << a << "\n";
        }

        if(f(a) * f(b) < 0){
            double root = bisection(a, b, E, out);
            out << "The roots are in between " << a << " and " << b
                << " the root is = " << root << "\n\n";
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
-2.24  2.24
```

#### Bisection Output
```
The search interval is: -2.240000 to 2.240000

iterations : 13
The roots are in between -1.040000 and -0.940000 the root is = -0.999998

iterations : 11
The roots are in between -0.040000 and 0.060000 the root is = -0.000010
```

---

### False Position Method

#### False Position Theory

- Regular Falsi method in Latin

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

double f(double x){
    return pow(x, 4) - 3*(pow(x, 3)) + 2*x*x + 6*x;
}

double falsePosition(double a, double b, double E, ofstream &out){
    double c;
    int count = 0;

    double fa = f(a);
    double fb = f(b);

    while(true){
        c = a - fa * (b - a)/(fb - fa);
        count++;

        if(fabs(f(c)) < E){
            out << "iterations : " << count << "\n";
            return c;
        }

        double fc = f(c);
        
        if(fa * fc < 0){
            b = c;
            fb = fc;
        }
        else{
            a = c;
            fa = fc;
        }
    }

    out << "iterations : " << count << "\n";
    return c;
}

int main(){
    ifstream input("input.txt");
    ofstream out("output.txt");

    double start, ends;
    input >> start >> ends;

    out << fixed << setprecision(6);
    out << "The search interval is: " << start << " to " << ends << "\n\n";

    double E = 1e-4;
    double a = start, b;
    double step = 0.1;

    while(a < ends){
        b = a + step;
        if(b > ends) b = ends;

        if(fabs(f(a)) < E){
            out << "The root is at = " << a << "\n";
        }

        if(f(a) * f(b) < 0){
            double root = falsePosition(a, b, E, out);
            out << "The roots are in between " << a << " and " << b
                << " the root is = " << root << "\n\n";
        }

        a = b;
    }

    input.close();
    out.close();
    return 0;
}
```

#### False Position Input
```
-2.24 2.24
```

#### False Position Output
```
The search interval is: -2.240000 to 2.240000

iterations : 4
The roots are in between -1.040000 and -0.940000 the root is = -0.999999

iterations : 2
The roots are in between -0.040000 and 0.060000 the root is = -0.000014
```

---

### Newton-Raphson Method

#### Newton-Raphson Theory

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


double f(double x){
    
    return x*x*x - x;
}

double diff(double x){
    return 3*x*x - 1;
}


void newton(double x, double &r, ofstream &out) {
    double fx = f(x);
    double f_px = diff(x);

    if (f_px == 0) { 
        x = x + 0.5;
        f_px = diff(x);
    }

    r = x - (fx / f_px);
    int counts = 0;

    while (fabs(r - x) >= E && counts < max_it) {
        x = r;
        counts++;
        fx = f(x);
        f_px = diff(x);

        if (f_px == 0) {
            x = x + 0.5;
            f_px = diff(x);
        }

        r = x - (fx / f_px);
    }

    if (counts >= max_it) {
        out << "No real root is found starting from " << x << ".\n";
        return;
    }

    out << "The root " << r << " is found after " << counts << " iterations.\n\n";
}

int main() {
    ifstream input("input.txt");
    ofstream out("output.txt");

    if(!input){
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    double start, ends;
    input >> start >> ends;
    input.close();

    out << fixed << setprecision(6);
    out << "Newton-Raphson Method\n";
    out << "Search interval: " << start << " to " << ends << "\n\n";

    double step = 0.5;
    double a = start, b, r;
    bool found = false;

    while(a < ends){
        b = a + step;
        if(b > ends) b = ends;

        if(f(a) * f(b) < 0){
            out << "Interval [" << a << ", " << b << "] seems to contain a root.\n";
            newton(a, r, out);
            found = true;
        }

        a = b;
    }

    if(!found){
        out << "No sign change detected automatically.\n";
        out << "Please enter a guessed value of x in the input file for manual root search.\n";
    }

    out.close();
    return 0;
}
```

#### Newton-Raphson Input
```
-2.23 2.23
```

#### Newton-Raphson Output
```
Newton-Raphson Method
Search interval: -2.230000 to 2.230000

Interval [-1.230000, -0.730000] seems to contain a root.
The root -1.000000 is found after 4 iterations.

Interval [-0.230000, 0.270000] seems to contain a root.
The root 0.000000 is found after 3 iterations.

Interval [0.770000, 1.270000] seems to contain a root.
The root 1.000000 is found after 5 iterations.
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
// Add your Secant Method code here
```

#### Secant Input
```
Add your Secant Method input here
```

#### Secant Output
```
Add your Secant Method output here
```

---

## Solution of Linear Equations

A linear equation is an equation in which the highest power (degree) of the variable is 1.

**Example:**
- 3x - 2 = 0

---

### Gauss Elimination Method

#### Gauss Elimination Theory

[Add your Gauss Elimination theory here]

#### Gauss Elimination Code
```cpp
// Add your Gauss Elimination code here
```

#### Gauss Elimination Input
```
Add your Gauss Elimination input here
```

#### Gauss Elimination Output
```
Add your Gauss Elimination output here
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

[Add your Gauss Jordan theory here]

#### Gauss Jordan Code
```cpp
// Add your Gauss Jordan code here
```

#### Gauss Jordan Input
```
Add your Gauss Jordan input here
```

#### Gauss Jordan Output
```
Add your Gauss Jordan output here
```

---

### LU Decomposition Method

#### LU Decomposition Theory

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

[Add your Runge-Kutta theory here]

#### Runge-Kutta Code
```cpp
// Add your Runge-Kutta code here
```

#### Runge-Kutta Input
```
Add your Runge-Kutta input here
```

#### Runge-Kutta Output
```
Add your Runge-Kutta output here
```

---

## Interpolation
Interpolation is the process of estimating the value of a function at a point that lies between two known data points.It uses known values of a function to construct a new function (often a polynomial) that approximates the original function and allows us to compute missing values inside the data interval.

---

### Newton's Forward Interpolation

#### Newtons Forward Theory


-Newton’s forward interpolation method is used to estimate values near the beginning of a table when data points(x is data point) are equally spaced.

-Let x₀, x₁, x₂, … , xₙ₋₁, xₙ  be a set of equally spaced values of the independent variable x.

-So x₁ − x₀ = x₂ − x₁ = x₃ − x₂ = … = xₙ − xₙ₋₁ = h

-Let, u = (x − x₀) / h

-The Newton’s forward difference interpolation formula for equal intervals is:
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

 -Backward interpolation is used when the value of x lies near the end of the table.
 
 -Let x₀, x₁, x₂, … , xₙ₋₁, xₙ be a set of equally spaced values of the independent variable x.
 -The interval between successive values is constant, i.e., x₁ − x₀ = x₂ − x₁ = x₃ − x₂ = … = xₙ − xₙ₋₁ = h
 -Let, u = (x − xₙ) / h
 -The Newton’s backward difference interpolation formula for equal intervals is:
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

-The Newton’s Divided Difference Interpolation Method is used to estimate the value of a function when the given data points are not equally spaced.

-The interpolation polynomial is given by:

f(xₙ)= f(x₀) + (x - x₀) f [x₁, x₀] + (x - x₀)(x - x₁) f[x₂, x₁, x₀] +...+ (x - x₀)(x - x₁)...(x − xₙ₋₁) f[xₙ, xₙ₋₁,..., x₁, x₀]...(1)

Here divided differences are defined as:

-First order divided difference:
f [ xᵢ, xⱼ ] = [ f(xᵢ) - f(xⱼ) ] / xᵢ - xⱼ …..(2)

-Second order divided difference:
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

#### Numerical Differentiation Theory

-This formula is the Newton’s Forward Differentiation Formula, used to find the first derivative of a function when the tabulated values of x are equally spaced.

-Let, x₀, x₁, x₂, … , xₙ be equally spaced values with h = x₁ − x₀

-Let, u = (x − x₀) / h measures how far the point x is from the starting value x0​ in terms of step size.

-So, the first derivative formula for forward differentiation is:
y’ = (1/h) [ Δy₀ + (2u − 1)/2! · Δ²y₀ + (3u² − 6u + 2)/3! · Δ³y₀ + (4u³ − 18u² + 22u - 6)/4! · Δ⁴y₀ + … ]

-Δy0,Δ^2y0,Δ^3y0​,… are forward differences of the function values.

#### Numerical Differentiation Code
```cpp
#include <iostream>
#include <iomanip>
using namespace std;

double facto(int n){
    double f=1;
    for(int i=2;i<=n;i++)f*=i;
    return f;
}

int main(){
    int n;cin>>n;
    double x[50],y[50][50],X;

    for(int i=0;i<n;i++)cin>>x[i]>>y[i][0];
    cin>>X;

    for(int j=1;j<n;j++)
        for(int i=0;i<n-j;i++)
            y[i][j]=y[i+1][j-1]-y[i][j-1];

    cout<<"Forward Difference Table:\n";
    for(int i=0;i<n;i++){
        cout<<setw(6)<<x[i];
        for(int j=0;j<n-i;j++)
            cout<<setw(10)<<y[i][j];
        cout<<endl;
    }

    double h=x[1]-x[0];
    double u=(X-x[0])/h;

    double dydx=y[0][1]+(2*u-1)*y[0][2]/2 +(3*u*u-6*u+2)*y[0][3]/6;

    dydx/=h;

    cout<<"\nFirst derivative at x="<<X<<" is: "<<dydx<<endl;
    return 0;
}

```

#### Numerical Differentiation Input
```
5
1 1
2 8
3 27
4 64
5 125
2
```

#### Numerical Differentiation Output
```
Forward Difference Table:
     1         1         7        12         6         0
     2         8        19        18         6
     3        27        37        24
     4        64        61
     5       125

First derivative at x=2 is: 12

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

[Add your Polynomial Curve Fitting theory here]

#### Polynomial Code
```cpp
// Add your Polynomial Curve Fitting code here
```

#### Polynomial Input
```
Add your Polynomial Curve Fitting input here
```

#### Polynomial Output
```
Add your Polynomial Curve Fitting output here
```

---

## Numerical Integration

---

### Simpson's 1/3 Rule

#### Simpsons 13 Theory

[Add your Simpson's 1/3 Rule theory here]

#### Simpsons 13 Code
```cpp
// Add your Simpson's 1/3 Rule code here
```

#### Simpsons 13 Input
```
Add your Simpson's 1/3 Rule input here
```

#### Simpsons 13 Output
```
Add your Simpson's 1/3 Rule output here
```

---

### Simpson's 3/8 Rule

#### Simpsons 38 Theory

[Add your Simpson's 3/8 Rule theory here]

#### Simpsons 38 Code
```cpp
// Add your Simpson's 3/8 Rule code here
```

#### Simpsons 38 Input
```
Add your Simpson's 3/8 Rule input here
```

#### Simpsons 38 Output
```
Add your Simpson's 3/8 Rule output here
```

---
