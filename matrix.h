template<typename T>
class Matrix {
private:
    std::vector<std::vector<T>> mat;

public:
    Matrix(const std::vector<std::vector<T>>& ini) : mat(ini) {}

    template<typename T1>
    Matrix(const Matrix<T1>& another) {
        auto[n, m] = another.size();
        mat = std::vector<std::vector<T>>(n, std::vector<T>(m));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                mat[i][j] = another[i][j];
            }
        }
    }

    class iterator {
    private:
        size_t i;
        size_t j;
        const Matrix* p;

    public:
        iterator(std::pair<size_t, size_t> a, const Matrix* b) : i(a.first), j(a.second), p(b) {}

        iterator(size_t a, size_t b, const Matrix* c) : i(a), j(b), p(c) {}

        iterator(const iterator& a) : i(a.geti()), j(a.getj()), p(a.getptr()) {}

        size_t geti() const {
            return i;
        }

        size_t getj() const {
            return j;
        }

        const Matrix* getptr() const {
            return p;
        }

        bool operator == (const iterator& a) const {
            return i == a.geti() && j == a.getj() && p == a.getptr();
        }

        bool operator != (const iterator& a) const {
            return !(*this == a);
        }

        iterator& operator++ () {
            auto[n, m] = p->size();
            j += 1;
            if (j == m) {
                i += 1;
                if (i == n) {
                    j = m;
                } else {
                    j = 0;
                }
            }

            return *this;
        }

        iterator operator++ (int) {
            iterator ret(*this);
            ++(*this);

            return ret;
        }

        const T& operator * () const {
            return p->mat[i][j];
        }
    };

    const iterator begin() const {
        return iterator(0, 0, this);
    }

    const iterator end() const {
        return iterator(size(), this);
    }

    std::pair<size_t, size_t> size() const {
        return {mat.size(), (mat.empty() ? 0 : mat[0].size())};
    }

    std::vector<T>& operator[] (const int i) {
        return mat[i];
    }

    const std::vector<T>& operator[] (const int i) const {
        return mat[i];
    }

    Matrix& operator += (const Matrix& another) {
        auto[n, m] = size();

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                mat[i][j] += another[i][j];
            }
        }

        return *this;
    }

    Matrix operator + (const Matrix& another) const {
        Matrix ret(*this);
        ret += another;
        return ret;
    }

    Matrix& operator -= (const Matrix& another) {
        auto[n, m] = size();

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                mat[i][j] -= another[i][j];
            }
        }

        return *this;
    }

    Matrix operator - (const Matrix& another) const {
        Matrix ret(*this);
        ret -= another;
        return ret;
    }

    Matrix& operator *= (const Matrix& another) {
        auto[n, m] = size();
        auto l = another.size().second;
        std::vector<std::vector<T>> vec(n, std::vector<T>(l));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t k = 0; k < l; ++k) {
                    vec[i][k] += mat[i][j] * another[j][k];
                }
            }
        }

        mat = vec;
        return *this;
    }

    template<typename TN>
    Matrix& operator *= (const TN& k) {
        auto[n, m] = size();

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                mat[i][j] *= k;
            }
        }

        return *this;
    }

    template<typename TN>
    Matrix operator * (const TN& k) const {
        Matrix ret(*this);
        ret *= k;
        return ret;
    }

    Matrix& transpose() {
        auto[n, m] = size();
        std::vector<std::vector<T>> vec(m, std::vector<T>(n));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                vec[j][i] = mat[i][j];
            }
        }

        this->mat = vec;
        return *this;
    }

    Matrix transposed() const {
        auto[n, m] = size();
        std::vector<std::vector<T>> vec(m, std::vector<T>(n));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                vec[j][i] = mat[i][j];
            }
        }

        Matrix ret(vec);
        return ret;
    }

    template<typename U>
    std::vector<U> solve(const std::vector<U>& b) {
        size_t n = size().first;
        Matrix<U> a_matrix(*this);
        std::vector<int> where_is_x_number(n, -1);
        std::vector<bool> used(n);
        std::vector<U> b_vector(b);

        for (size_t col = 0, row = 0; row < n && col < n; ++col) {
            size_t cur_row = 0;

            while (used[cur_row]) {
                ++cur_row;
            }

            for (size_t i = 0; i < n; ++i) {
                if (used[i]) {
                    continue;
                }
                if (abs(a_matrix[i][col]) > abs(a_matrix[cur_row][col])) {
                    cur_row = i;
                }
            }

            used[cur_row] = true;
            where_is_x_number[col] = cur_row;

            for (size_t i = 0; i < n; ++i) {
                if (i != cur_row) {
                    long double c = a_matrix[i][col] / a_matrix[cur_row][col];
                    for (size_t j = col; j < n; ++j) {
                        a_matrix[i][j] -= c * a_matrix[cur_row][j];
                    }
                    b_vector[i] -= c * b_vector[cur_row];
                }
            }

            ++row;
        }

        std::vector<U> ans(n);
        for (size_t col = 0; col < n; ++col) {
            if (a_matrix[where_is_x_number[col]][col] == 0) {
                ans[col] = 0;
            } else {
                ans[col] = b_vector[where_is_x_number[col]] /
                           a_matrix[where_is_x_number[col]][col];
            }
        }

        return ans;
    }
};

template<typename T>
Matrix<T> operator * (const Matrix<T>& a, const Matrix<T>& b) {
    Matrix<T> ret(a);
    return ret *= b;
}

template<typename T, typename TN>
Matrix<T> operator * (const TN& k, const Matrix<T>& m) {
    Matrix<T> ret(m);
    return ret *= k;
}

template<typename T>
int tr(const Matrix<T>& ma) {
    auto[n, m] = ma.size();
    int ans = 0;

    for (size_t i = 0; i < std::min(n, m); ++i) {
        ans += ma[i][i];
    }

    return ans;
}

template <typename T>
std::ostream& operator << (std::ostream& out, const Matrix<T>& Mat) {
    auto[n, m] = Mat.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            out << Mat[i][j];
            if (j != m - 1) {
                out << '\t';
            }
        }
        if (i != n - 1) {
            out << '\n';
        }
    }

    return out;
}