#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>

#define EPSILON 0.000000001
#define DEBUG 1

template <typename T>
class Matrix2D
{
public:
    using vec = typename std::vector<T>;
    using iter = typename vec::iterator;

    Matrix2D(
        size_t row,
        size_t col) : row_(row),
                      col_(col),
                      mat_(row * col)
    {
    }

    Matrix2D(std::vector<std::vector<T>> mat)
    {
        size_t col = mat[0].size();
        for (const auto &row : mat)
        {
            assert(row.size() == col);
        }
        size_t row = mat.size();

        this->row_ = row;
        this->col_ = col;
        this->mat_.reserve(row * col);
        for (auto &&row : mat)
        {
            this->mat_.insert(this->mat_.end(), row.begin(), row.end());
        }
    }

    Matrix2D(std::initializer_list<std::initializer_list<T>> mat) : Matrix2D(std::vector<std::vector<T>>(mat.begin(), mat.end())) {}

    static Matrix2D<T> identity(size_t n)
    {
        return Matrix2D::identity(n, n);
    }

    static Matrix2D<T> identity(size_t row, size_t col)
    {
        Matrix2D<T> ret(row, col);
        for (size_t i = 0, n = std::min(row, col); i < n; i++)
        {
            ret.get(i, i) = 1;
        }
        return ret;
    }

    T &get(size_t r, size_t c)
    {
        assert(r < this->row());
        assert(c < this->col());
        return this->mat_.at(r * this->col() + c);
    }

    const T &get(size_t r, size_t c) const
    {
        assert(r < this->row());
        assert(c < this->col());
        return this->mat_.at(r * this->col() + c);
    }

    iter get_it(size_t r, size_t c)
    {
        assert(r < this->row());
        assert(c < this->col());
        auto b = this->mat_.begin();
        return b + (r * this->col()) + c;
    }

    size_t row() const
    {
        return this->row_;
    }

    size_t col() const
    {
        return this->col_;
    }

    Matrix2D<T> mul(const Matrix2D<T> &other) const
    {
        assert(this->col() == other.row());

        Matrix2D<T> ret(this->row(), other.col());
        for (size_t r = 0; r < ret.row(); r++)
        {
            for (size_t c = 0; c < ret.col(); c++)
            {
                T acc = 0;
                for (size_t i = 0; i < this->row(); i++)
                    acc += (this->get(r, i) * other.get(i, c));
                ret.get(r, c) = acc;
            }
        }

        return ret;
    }

    bool operator==(const Matrix2D<T> &other) const
    {
        if (this->row() != other.row() || this->col() != other.col())
            return false;
        return this->mat_ == other.mat_;
    }

    Matrix2D<T> operator+(const Matrix2D<T> &other) const
    {
        assert(this->row() == other.row());
        assert(this->col() == other.col());
        Matrix2D<T> ret(this->row(), this->col());
        for (size_t i = 0; i < this->row(); i++)
        {
            for (size_t j = 0; j < this->col(); j++)
            {
                ret.get(i, j) = this->get(i, j) + other.get(i, j);
            }
        }

        return ret;
    }

    Matrix2D<T> operator-(const Matrix2D<T> &other) const
    {
        assert(this->row() == other.row());
        assert(this->col() == other.col());
        Matrix2D<T> ret(this->row(), this->col());
        for (size_t i = 0; i < this->row(); i++)
        {
            for (size_t j = 0; j < this->col(); j++)
            {
                ret.get(i, j) = this->get(i, j) - other.get(i, j);
            }
        }

        return ret;
    }

    // Elementary row operations

    // Swap two rows
    // R_i <-> R_j
    void swap(size_t i, size_t j)
    {
        assert(i < this->row());
        assert(j < this->row());
        for (size_t k = 0; k < this->col(); k++)
            std::swap(this->get(i, k), this->get(j, k));
    }

    // Multiply any row by a constant
    // R_i *= c
    void mul_scalar(size_t row, double c)
    {
        assert(row < this->row());
        for (size_t i = 0; i < this->col(); i++)
        {
            this->get(row, i) *= c;
        }
    }

    // Add one row to another
    // R_i += c * R_j
    void add_row(size_t dst, size_t src, double c)
    {
        assert(dst < this->row());
        assert(src < this->row());

        for (size_t i = 0; i < this->col(); i++)
        {
            this->get(dst, i) += this->get(src, i) * c;
        }
    }

private:
    size_t row_, col_;
    std::vector<T> mat_;
};

template <typename T>
void print(const Matrix2D<T> mat)
{
    for (size_t r = 0; r < mat.row(); r++)
    {
        for (size_t c = 0; c < mat.col(); c++)
        {
            std::cout << std::setw(10) << mat.get(r, c) << ' ';
        }
        std::cout << '\n';
    }
}

// Find P, L, U so that PA = LU
std::tuple<Matrix2D<double>, Matrix2D<double>, Matrix2D<double>> plu_factorization(const Matrix2D<double> &A)
{
    Matrix2D<double> U = A;
    auto P = Matrix2D<double>::identity(A.row(), A.col());
    auto L = Matrix2D<double>(A.row(), A.col());

    for (size_t c = 0; c < U.col(); c++)
    {
        if (U.get(c, c) == 0)
        {
            size_t k = c + 1;
            for (; k < U.row(); k++)
            {
                if (U.get(k, c) != 0)
                {
                    U.swap(k, c);
                    P.swap(k, c);
                    L.swap(k, c);
                    break;
                }
            }
            assert(k != U.row());
        }
        auto d = U.get(c, c);
        assert(d != 0);
        for (size_t r = c + 1; r < U.row(); r++)
        {
            auto s = -U.get(r, c) / d;
            U.add_row(r, c, s);
            L.get(r, c) = -s;
        }
    }

    for (size_t i = 0, n = std::min(L.row(), L.col()); i < n; i++)
    {
        L.get(i, i) = 1;
    }

    return std::make_tuple(P, L, U);
}

Matrix2D<double> jacobi(Matrix2D<double> a, Matrix2D<double> b)
{
    assert(a.col() == b.row());
    Matrix2D<double> x = b;
    for (size_t i = 0; i < x.row(); i++)
    {
        x.get(i, 0) = rand() % 16;
    }
    Matrix2D<double> D(a.row(), a.col()), D_inv(a.row(), a.col());
    for (size_t i = 0; i < a.col(); i++)
    {
        D.get(i, i) = a.get(i, i);
        D_inv.get(i, i) = 1 / a.get(i, i);
    }
    Matrix2D<double> L = a, U = a;
    for (size_t i = 0; i < a.row(); i++)
    {
        for (size_t j = 0; j < a.col(); j++)
        {
            if (i <= j)
                L.get(i, j) = 0;
            if (i >= j)
                U.get(i, j) = 0;
        }
    }
    auto LU = L + U;

    bool done = false;
    while (!done)
    {
        auto _x = x;
        x = D_inv.mul(b - LU.mul(x));
        auto diff = x - _x;
        bool _done = true;
        for (size_t i = 0; i < x.row(); i++)
        {
            if (diff.get(i, 0) > EPSILON)
            {
                _done = false;
                break;
            }
        }
        done = _done;
    }
    return x;
}

class func
{
public:
    func(double x, double y) : x_(x), y_(y)
    {
        this->f_ = [](double x, double y, double h, double k, double r)
        {
            return (x - h) * (x - h) + (y - k) * (y - k) - r * r;
        };
    }

    func(double x, double y, std::function<double(double, double, double, double, double)> f) : x_(x), y_(y), f_(f) {}

    double eval(double h, double k, double r)
    {
        return this->f_(this->x_, this->y_, h, k, r);
    }

    // HACK: hard-coded
    func dh()
    {
        return func(
            this->x_,
            this->y_,
            [](double x, double y, double h, double k, double r)
            {
                return -2 * (x - h);
            });
    }
    func dk()
    {
        return func(
            this->x_,
            this->y_,
            [](double x, double y, double h, double k, double r)
            {
                return -2 * (y - k);
            });
    }
    func dr()
    {
        return func(
            this->x_,
            this->y_,
            [](double x, double y, double h, double k, double r)
            {
                return -2 * r;
            });
    }

private:
    double x_, y_;
    std::function<double(double, double, double, double, double)> f_;
};

class DF
{
public:
    DF(func f1, func f2, func f3)
    {
        std::vector<func> df1 = {f1.dh(), f1.dk(), f1.dr()};
        std::vector<func> df2 = {f2.dh(), f2.dk(), f2.dr()};
        std::vector<func> df3 = {f3.dh(), f3.dk(), f3.dr()};

        this->funcs.push_back(df1);
        this->funcs.push_back(df2);
        this->funcs.push_back(df3);
    }

    DF(std::vector<std::vector<func>> _funcs) : funcs(_funcs) {}

    Matrix2D<double> eval(double h, double k, double r)
    {
        Matrix2D<double> ret(3, 3);
        for (size_t i = 0; i < 3; i++)
        {
            for (size_t j = 0; j < 3; j++)
            {
                assert(i < ret.row() && j < ret.col());
                ret.get(i, j) = this->funcs[i][j].eval(h, k, r);
            }
        }

        return ret;
    }

private:
    // Should be 3x3
    std::vector<std::vector<func>> funcs;
};

class F
{
public:
    F(std::initializer_list<func> funcs) : funcs_(funcs) {}
    F(std::vector<func> funcs) : funcs_(funcs) {}

    Matrix2D<double> eval(double h, double k, double r)
    {
        Matrix2D<double> ret(this->funcs_.size(), 1);
        for (size_t i = 0; i < ret.row(); i++)
        {
            ret.get(i, 0) = this->funcs_[i].eval(h, k, r);
        }
        return ret;
    }

    DF gen_df()
    {
        return DF(
            this->funcs_[0],
            this->funcs_[1],
            this->funcs_[2]);
    }

private:
    std::vector<func> funcs_;
};

template <typename T>
std::vector<T> gaussian_elimination(const Matrix2D<T> &A)
{
#if DEBUG > 1
    std::cout << ">>gaussian_elimination:\n";
    print(A);
#endif

    Matrix2D<T> a = A;
    for (size_t c = 0; c < a.col(); c++)
    {
        for (size_t r = c + 1; r < a.row(); r++)
        {
            assert(a.get(c, c) != 0);
            a.add_row(r, c, -a.get(r, c) / a.get(c, c));
        }
    }
    for (size_t r = 0; r < a.row(); r++)
    {
        assert(a.get(r, r) != 0);
        a.mul_scalar(r, 1 / a.get(r, r));
    }

    std::vector<T> ans(a.row());

    for (size_t r = a.row() - 1;; r--)
    {
        auto new_ans = a.get(r, a.col() - 1);
        for (size_t c = a.col() - 2; c >= r; c--)
        {
            new_ans -= a.get(r, c) * ans[c];
            if (r == 0 && c == 0)
                break;
        }
        ans[r] = new_ans;
        if (r == 0)
            break;
    }

#if DEBUG > 1
    for (auto v : ans)
        std::cout << v << ' ';
    std::cout << '\n';
    std::cout << "<<gaussian_elimination:\n";
#endif

    return ans;
}

void solve_p7(std::vector<std::pair<double, double>> inputs)
{
    std::vector<func> funcs;
    for (auto [x, y] : inputs)
    {
        funcs.push_back(func(x, y));
    }

    F _F(funcs);
    DF _DF = _F.gen_df();

    Matrix2D<double> x = {
        {(double)(rand() % 4 + 1)},
        {(double)(rand() % 4 + 1)},
        {(double)(rand() % 4 + 1)},
    };
    assert(x.row() == 3 && x.col() == 1);

#if DEBUG
    std::cout << "Init x =";
    print(x);
#endif

    int _try = 100000;
    while (_try--)
    {
        auto A = _DF.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0));
        auto b = _F.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0));
        for (size_t i = 0; i < b.row(); i++)
        {
            b.get(i, 0) = -b.get(i, 0);
        }
        Matrix2D<double> s(A.row(), A.col() + 1);
        for (size_t i = 0; i < s.row(); i++)
        {
            for (size_t j = 0; j < s.col() - 1; j++)
            {
                s.get(i, j) = A.get(i, j);
            }
            s.get(i, s.col() - 1) = b.get(i, 0);
        }

        auto s_ = gaussian_elimination(s);
        Matrix2D<double> s__(s_.size(), 1);
        for (size_t i = 0; i < s_.size(); i++)
        {
            s__.get(i, 0) = s_[i];
        };

        // double sum = 0;
        // for (size_t i = 0; i < s.row(); i++)
        // {
        //     sum += abs(s__.get(i, 0));
        // }
        // if (sum < EPSILON)
        //     break;
        x = x + s__;
    }

#if DEBUG
    std::cout << "try = " << _try << '\n';
#endif

    std::cout << "x = \n";
    print(x);

#if DEBUG
    std::cout << "F(x) = \n";
    for (auto fn : funcs)
    {
        std::cout << fn.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0)) << '\n';
    }
#endif
}

void solve_mid_e()
{

    auto f1 = [](double p, double q, double x, double y, double z)
    {
        return sin(x) + y * y + log(z) - 7;
    };
    // Maybe -sin
    auto f1dx = [](double p, double q, double x, double y, double z)
    {
        return cos(x);
    };
    auto f1dy = [](double p, double q, double x, double y, double z)
    {
        return 2 * y;
    };
    auto f1dz = [](double p, double q, double x, double y, double z)
    {
        return 1 / z;
    };
    auto f2 = [](double p, double q, double x, double y, double z)
    {
        return 3 * x + 2 * y - z * z * z + 1;
    };
    auto f2dx = [](double p, double q, double x, double y, double z)
    {
        return 3;
    };
    auto f2dy = [](double p, double q, double x, double y, double z)
    {
        return 2;
    };
    auto f2dz = [](double p, double q, double x, double y, double z)
    {
        return -3 * z * z;
    };
    auto f3 = [](double p, double q, double x, double y, double z)
    {
        return x + y + z - 5;
    };
    auto f3dx = [](double p, double q, double x, double y, double z)
    {
        return 1;
    };
    auto f3dy = [](double p, double q, double x, double y, double z)
    {
        return 1;
    };
    auto f3dz = [](double p, double q, double x, double y, double z)
    {
        return 1;
    };

    std::vector<func> funcs{
        func(0, 0, f1),
        func(0, 0, f2),
        func(0, 0, f3)};

    std::vector<std::vector<func>> funcs2{
        {func(0, 0, f1dx), func(0, 0, f1dy), func(0, 0, f1dz)},
        {func(0, 0, f2dx), func(0, 0, f2dy), func(0, 0, f2dz)},
        {func(0, 0, f3dx), func(0, 0, f3dy), func(0, 0, f3dz)},
    };

    F _F(funcs);
    DF _DF(funcs2);

    Matrix2D<double> x = {
        {0},
        {2},
        {2},
    };
    assert(x.row() == 3 && x.col() == 1);

#if DEBUG
    std::cout << "Init x =";
    print(x);
#endif

    int _try = 100000;
    while (_try--)
    {
        auto A = _DF.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0));
        auto b = _F.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0));
        for (size_t i = 0; i < b.row(); i++)
        {
            b.get(i, 0) = -b.get(i, 0);
        }
        Matrix2D<double> s(A.row(), A.col() + 1);
        for (size_t i = 0; i < s.row(); i++)
        {
            for (size_t j = 0; j < s.col() - 1; j++)
            {
                s.get(i, j) = A.get(i, j);
            }
            s.get(i, s.col() - 1) = b.get(i, 0);
        }

        auto s_ = gaussian_elimination(s);
        Matrix2D<double> s__(s_.size(), 1);
        for (size_t i = 0; i < s_.size(); i++)
        {
            s__.get(i, 0) = s_[i];
        };

        // double sum = 0;
        // for (size_t i = 0; i < s.row(); i++)
        // {
        //     sum += abs(s__.get(i, 0));
        // }
        // if (sum < EPSILON)
        //     break;
        x = x + s__;
    }

#if DEBUG
    std::cout << "try = " << _try << '\n';
#endif

    std::cout << "x = \n";
    print(x);

#if DEBUG
    std::cout << "F(x) = \n";
    for (auto fn : funcs)
    {
        std::cout << fn.eval(x.get(0, 0), x.get(1, 0), x.get(2, 0)) << '\n';
    }
#endif
}

int main()
{
    srand(time(NULL));

    // std::vector<std::pair<double, double>> inputs = {
    //     {-8, -4},
    //     {6, 9},
    //     {4, -9},
    // };
    // solve_p7(inputs);

    // std::vector<std::pair<double, double>> inputs2 = {
    //     {-1, 6},
    //     {-2, -6},
    //     {5, 0},
    // };
    // solve_p7(inputs2);
    solve_mid_e();

    return 0;
}
