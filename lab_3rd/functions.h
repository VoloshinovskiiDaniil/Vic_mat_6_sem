#pragma once

#include "matrix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>

#ifndef FUNCTIONS
#define FUNCTIONS

template <typename T, unsigned int NX> struct Mesh 
{
    const double LEFT_B;
    const double RIGHT_B;
    const double h;
    std::array<T, NX> points;

    Mesh(const double x_left, const double x_right)
        : LEFT_B(x_left), RIGHT_B(x_right), h((RIGHT_B - LEFT_B) / (NX - 1)) {}

    T& operator[](std::size_t ind) 
    {
        return points[ind]; 
    }

    const T& operator[](std::size_t ind) const 
    { 
        return points[ind]; 
    }

    friend std::ostream& operator<<(std::ostream& os,
        const Mesh<T, NX>& grid) 
    {
        for (std::size_t i = 0; i < NX; i++) 
        {
            double x = grid.LEFT_B + i * grid.h;
            os << x << " " << grid[i] << std::endl;
        }
        return os;
    }
};

struct Equation 
{
    double a;
    std::function<double(double, double)> func;
    std::array<double, 2> left_b_cond;
    std::array<double, 2> right_b_cond;
    std::function<double(double)> left_u;
    std::function<double(double)> right_u;
};

template <typename T, unsigned int NX>
void solve(const Equation& eq, const Mesh<T, NX>& grid,
    const double T0, const double tau, std::ostream& file) 
{
    Mesh<T, NX> curr_grid = grid;
    Mesh<T, NX> next_grid = grid;
    double curr_time = 0;
    double Co = eq.a * tau / grid.h / grid.h / 2.;
    std::array<std::array<T, 3>, NX> data;
    for (int i = 1; i < NX - 1; i++) 
    {
        data[i][0] = -Co;
        data[i][1] = 1 + 2 * Co;
        data[i][2] = -Co;
    }
    data[0][0] = 0;
    data[0][1] = eq.left_b_cond[0] - eq.left_b_cond[1] / grid.h;
    data[0][2] = eq.left_b_cond[1] / grid.h - eq.left_b_cond[1] / 2. / grid.h / Co;
    data[NX - 1][0] = -eq.right_b_cond[1] / grid.h + eq.right_b_cond[1] / 2. / grid.h / Co;
    data[NX - 1][1] = eq.right_b_cond[0] + eq.right_b_cond[1] / grid.h;
    data[NX - 1][2] = 0;

    Matrix<T, NX> matrix(data);
    std::array<double, NX> column;
    file << curr_time << std::endl;
    file << curr_grid << std::endl << std::endl;

    while (curr_time < T0) 
    {
        column[0] = eq.left_u(curr_time + tau) - (curr_grid[1] +
                Co * (curr_grid[0] - 2 * curr_grid[1] + curr_grid[2]) +
                (eq.func(curr_time, grid.LEFT_B + grid.h) +
                    eq.func(curr_time + tau, grid.LEFT_B + grid.h)) * tau / 2) *
            (eq.left_b_cond[1] / Co / 2. / grid.h);
        for (int i = 1; i < NX - 1; i++) 
        {
            column[i] = curr_grid[i] + Co * (curr_grid[i - 1] - 2 * curr_grid[i] + curr_grid[i + 1]) +
                (eq.func(curr_time, grid.LEFT_B + i * grid.h) +
                    eq.func(curr_time + tau, grid.LEFT_B + i * grid.h)) * tau / 2;
        }
        column[NX - 1] = eq.right_u(curr_time + tau) + (curr_grid[NX - 2] +
                Co * (curr_grid[NX - 1] - 2 * curr_grid[NX - 2] + curr_grid[NX - 3]) +
                (eq.func(curr_time, grid.RIGHT_B - grid.h) +
                    eq.func(curr_time + tau, grid.RIGHT_B - grid.h)) * tau / 2) *
            (eq.right_b_cond[1] / Co / 2. / grid.h);

        next_grid.points = matrix.solve(column);
        file << curr_time + tau << std::endl;
        file << next_grid << std::endl << std::endl;
        curr_grid.points = next_grid.points;
        curr_time += tau;
    }
}

#endif
