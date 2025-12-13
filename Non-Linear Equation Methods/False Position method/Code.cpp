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
    ifstream input("D:\\Numerical project\\False Position method\\input.txt");
    ofstream out("D:\\Numerical project\\False Position method\\output.txt");

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