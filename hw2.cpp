#include <iostream>
#include <cmath>
#include <iomanip>
#define EPSILON 0.0000001

using namespace std;

double f(double x)
{
    return x * x * x - 2 * x - 2;
}

double g(double x)
{
    return powf64(2 * x + 2, 1 / 3.);
}

int main()
{
    double x1 = 0.5;
    double x2 = g(x1);
    while ((x2 - x1) / x2 > EPSILON)
    {
        x1 = x2;
        x2 = g(x1);
    }
    cout << setprecision(15) << x2 << '\n';
}
