#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <ctime>

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

Matrix2D<double> input_a()
{
    Matrix2D<double> ret(10, 10);
    for (size_t i = 0; i < 10; i++)
    {
        for (size_t j = 0; j < 10; j++)
        {
            if (i == j)
                ret.get(i, j) = 3;
            else if (abs(i - j) == 1)
                ret.get(i, j) = -1;
        }
    }

    return ret;
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

    print(D_inv);

    bool done = false;
    while (!done)
    {
        auto _x = x;
        x = D_inv.mul(b - LU.mul(x));
        auto diff = x - _x;
        bool _done = true;
        for (size_t i = 0; i < x.row(); i++)
        {
            if (diff.get(i, 0) != 0)
            {
                _done = false;
                break;
            }
        }
        done = _done;
    }
    return x;
}

int main()
{
    srand(time(NULL));

    auto A = input_a();
    // print(A);
    // std::cout << "===\n";
    std::vector<std::vector<double>> raw_b{
        {2},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {2},
    };
    Matrix2D<double> b(raw_b);
    auto x = jacobi(A, b);

    auto b_ = A.mul(x);
    // print(b_);
    // std::cout << "===\n";
    // print(b);
    // std::cout << "===\n";
    assert(b == b_);
    print(x);

    return 0;
}
