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
    ifstream input("D:\\Numerical project\\Newton Raphson method\\input.txt");
    ofstream out("D:\\Numerical project\\Newton Raphson method\\output.txt");

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
