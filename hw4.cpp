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

    // Elementary row operations

    // Swap two rows
    // R_i <-> R_j
    void swap(size_t i, size_t j)
    {
        assert(i < this->row());
        assert(j < this->row());
        // Copy R_i
        std::vector<T> i_copy(this->get_it(i, 0), this->get_it(i, this->col()));
        assert(i_copy.size() == this->col());
        // R_i = R_j
        this->mat_.assign(
            this->get_it(i, 0),
            this->get_it(j, 0),
            this->get_it(j, this->col()));

        this->mat_.assign(
            this->get_it(j, 0),
            i_copy.begin(),
            i_copy.end());
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

// Solve Ax = B
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

template <typename T>
std::vector<T> gaussian_elimination(const Matrix2D<T> &A)
{
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

    return ans;
}

int main()
{
    std::vector<std::vector<double>> raw_mat = {
        {2, 1, 0, 0, 1},
        {0, 1, 2, 0, 1},
        {2, 4, 5, 1, 2},
        {8, 5, 0, 3, 0},
    };
    Matrix2D<double> mat(raw_mat);

    const auto ans = gaussian_elimination(mat);

    // Check
    assert(ans.size() == 4);
    for (size_t r = 0; r < mat.row(); r++)
    {
        double y = 0;
        for (size_t c = 0; c < mat.col() - 1; c++)
        {
            y += mat.get(r, c) * ans[c];
        }
        const double epsilon = 0.000001;
        assert(abs(y - mat.get(r, mat.col() - 1)) < epsilon);
    }

    std::cout << "Gaussian elimination\n";
    std::cout << "[x y z w] = [";
    for (const auto v : ans)
        std::cout << std::setw(4) << v << ' ';
    std::cout << "]\n";

    std::vector<std::vector<double>> raw_A = {
        {2, 1, 0, 0},
        {0, 1, 2, 0},
        {2, 4, 5, 1},
        {8, 5, 0, 3},
    };
    Matrix2D<double> A(raw_A);
    std::vector<double> B{1, 1, 2, 0};

    const auto [x, L, U] = lu_factorization(A, B);
    std::cout << "LU Factorization\n";
    std::cout << "[x y z w] = [";
    for (const auto v : x)
        std::cout << std::setw(4) << v << ' ';
    std::cout << "]\n";
    std::cout << "L = \n";
    print(L);
    std::cout << "U = \n";
    print(U);

    return 0;
}
