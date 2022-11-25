#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

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

    bool operator==(const Matrix2D<T> other) const
    {
        if (this->row() != other.row() || this->col() != other.col())
            return false;
        return this->mat_ == other.mat_;
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
            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << mat.get(r, c) << ' ';
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

std::tuple<std::vector<double>, Matrix2D<double>, Matrix2D<double>> lu_factorization(
    const Matrix2D<double> &A,
    const std::vector<double> &B)
{
    assert(A.row() == B.size());

    Matrix2D<double> U = A;
    Matrix2D<double> L(A.row(), A.col());

    // Make L an identity matrix
    for (size_t r = 0; r < U.row(); r++)
        L.get(r, r) = 1;
    for (size_t c = 0; c < U.col(); c++)
    {
        for (size_t r = c + 1; r < U.row(); r++)
        {
            assert(U.get(c, c) != 0);
            auto s = -U.get(r, c) / U.get(c, c);
            U.add_row(r, c, s);
            L.get(r, c) = -s;
            // std::cout << "XXX\n";
            // print(U);
        }
    }

    // Solve Lc = B
    auto x = B;
    for (size_t r = 0; r < L.row(); r++)
    {
        for (size_t c = 0; c < r; c++)
        {
            x[r] -= x[c] * L.get(r, c);
        }
        x[r] /= L.get(r, r);
    }
    // Solve Ux = c
    for (size_t r = U.row() - 1;; r--)
    {
        for (size_t c = r + 1; c < U.row(); c++)
        {
            x[r] -= x[c] * U.get(r, c);
        }
        x[r] /= U.get(r, r);

        if (r == 0)
            break;
    }

    return std::make_tuple(x, L, U);
}

int main()
{
    // std::vector<std::vector<double>> raw_A = {
    //     {0, 0, -1, 1},
    //     {1, 1, -1, 2},
    //     {-1, -1, 2, 0},
    //     {1, 2, 0, 2},
    // };
    std::vector<std::vector<double>> raw_A = {
        {1.19, 2.11, -100, 1},
        {14.2, +0.122, 12.2, -1},
        {0, 100, -99.9, 1},
        {15.3, 0.110, -13.1, -1},
    };
    Matrix2D<double> A(raw_A);
    std::vector<double> b{1.12, 3.44, 2.15, 4.16};

    auto [x, L, U] = lu_factorization(A, b);

    std::cout << "L = \n";
    print(L);
    std::cout << "U = \n";
    print(U);
    std::cout << "LU = \n";
    print(L.mul(U));
    std::cout << "A = \n";
    print(A);

    for (auto v : x)
        std::cout << v << ' ';
    std::cout << '\n';

    // const auto [P, L, U] = plu_factorization(A);
    // std::cout << "P = \n";
    // print(P);
    // std::cout << "L = \n";
    // print(L);
    // std::cout << "U = \n";
    // print(U);

    // auto PA = P.mul(A);
    // auto LU = L.mul(U);
    // std::cout << "PA = \n";
    // print(PA);
    // std::cout << "LU = \n";
    // print(LU);
    // std::cout << (PA == LU ? "PA == LU" : "PA != LU") << '\n';

    return 0;
}
