#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#define EPSILON 0.0000000001

using namespace std;

// double f(double x)
// {
//     return x * x * x - 2 * x - 2;
// }

// double g(double x)
// {
//     return powf64(2 * x + 2, 1 / 3.);
// }

// double f(double x)
// {
//     return x * x * x + 4 * x * x - 10;
// }

// double g(double x)
// {
//     return powf64((-x * x * x + 10) / 4, 1 / 2.);
// }

double f(double x)
{
    assert((x + x * x) != 0);
    return (3 - 2 * x - x * x) / (x + x * x) - 3.06;
}

double g(double x)
{
    return powf64(-3.06 * x * x - 5.06 * x + 3, 1 / 3.);
}

int main()
{
    cout << fixed << setprecision(6);
    for (int i = 0; i < 100; i++)
    {
        double x1 = 0.4524 + 0.000001 * i;
        double x2 = g(x1);
        double _x1 = x1;
        // cout << "x1, x2 = " << x1 << ", " << x2 << '\n';
        while ((x2 - x1) / x2 > EPSILON)
        {
            x1 = x2;
            x2 = g(x1);
        }
        if (!isnan(x2))
        {
            cout << x1 << ' ';
            cout << x2 << ' ';
            cout << f(x2) << '\n';
        }
    }

    cout << f(0.4385) << '\n';

    return 0;
}
