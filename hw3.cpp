#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#define EPSILON 0.0000001

using namespace std;

struct Question
{
    const char *expr;
    function<double(double)> f;
    function<double(double)> g;
};

int main()
{
    Question qs[] = {
        Question{
            .expr = "x^5 + x = 1",
            .f = [](double x)
            { return x * x * x * x * x + x - 1; },
            .g = [](double x)
            { return 5 * x * x * x * x + 1; },
        },
        Question{
            .expr = "ln x + x^2 = 3",
            .f = [](double x)
            { return log(x) + x * x - 3; },
            .g = [](double x)
            { return 2 * x + 1 / x; },
        },
        Question{
            .expr = "sin x = 6x + 5",
            .f = [](double x)
            { return sin(x) - 6 * x - 5; },
            .g = [](double x)
            { return cos(x) - 6; },
        },
    };
    cout << setprecision(4);
    for (size_t i = 0; i < 3; i++)
    {
        const auto &q = qs[i];
        double x = 0.5;
        do
        {
            x = x - q.f(x) / q.g(x);
        } while (fabs(q.f(x)) > EPSILON);
        cout << q.expr << "; root = " << x << '\n';
    }

    return 0;
}
