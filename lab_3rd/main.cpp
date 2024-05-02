#include "functions.h"
#include "matrix.h"

double func(double t, double x) 
{ 
    return t * (x + 1); 
}

double start_u(double x)
{
    return 0;
}
double left_u(double t) 
{
    return t * t; 
}
double right_u(double t) 
{
    return t * t; 
}

const int NX = 400;

int main() 
{
    std::ofstream data("data_400.txt");
    const double T0 = 10;
    const double tau = 0.01;
    Equation heat;
    heat.a = 1;
    heat.left_b_cond = { 0, 1 };
    heat.right_b_cond = { 1, 0 };
    heat.left_u = left_u;
    heat.right_u = right_u;
    heat.func = func;

    Mesh<double, NX> grid(0, 1);
    for (int i = 0; i < NX; i++) 
    {
        grid[i] = 0;
    }

    solve<double, NX>(heat, grid, T0, tau, data);
    data.close();

    return 0;
}

