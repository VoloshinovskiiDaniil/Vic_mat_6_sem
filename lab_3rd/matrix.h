#pragma once

#include <array>
#include <vector>

#ifndef MATRIX
#define MATRIX

template <typename T, unsigned int N> class Matrix 
{
private:
    std::array<std::array<T, 3>, N> element;

public:
    Matrix(const std::array<std::array<T, 3>, N>& data) : element{ data } {}
    T operator()(const int i, const int j) const 
    {
        if (j == i - 1) 
        {
            return element[i][0];
        }
        else if (j == i) 
        {
            return element[i][1];
        }
        else if (j == i + 1) 
        {
            return element[i][2];
        }
        else 
        {
            return 0;
        }
    }

    std::array<T, N> solve(const std::array<T, N>& column) const 
    {
        int n = N - 1;
        std::vector<T> p{ 0 }, q{ 0 };
        p.reserve(N);
        q.reserve(N);
        for (int i = 0; i < n + 1; ++i) 
        {
            p.push_back(-1 * (*this)(i, i + 1) / ((*this)(i, i - 1) * p[i] + (*this)(i, i)));
            q.push_back((column[i] - (*this)(i, i - 1) * q[i]) / ((*this)(i, i - 1) * p[i] + (*this)(i, i)));
        }
        std::array<T, N> solution;
        solution[n] = (column[n] - (*this)(n, n - 1) * q[n]) / ((*this)(n, n - 1) * p[n] + (*this)(n, n));
        for (int i = n - 1; i >= 0; i--) 
        {
            solution[i] = p[i + 1] * solution[i + 1] + q[i + 1];
        }
        return solution;
    }
};

#endif
