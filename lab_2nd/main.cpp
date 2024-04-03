#include "functions.h"

int main() 
{
    const double k = 5.0 / 3;
    const double p_l = 13;
    const double v_l = 0;
    const double e_l = 116913;
    const double p_r = 1.3;
    const double v_r = 0;
    const double e_r = 116913;

    const double tau = 0.0001;
    const double t_s = tau;
    const double t_e = 0.02;
    const double L = 10;

    const double Co_max = 1;
    const unsigned int M = 100;

    DVE_1D<DVE, M> mesh(-L, L);
    for (unsigned int x = 0; x < M / 2; x++) 
    {
        mesh[x].density = p_l;
        mesh[x].velocity = v_l;
        mesh[x].energy = e_l;

        mesh[M - x - 1].density = p_r;
        mesh[M - x - 1].velocity = v_r;
        mesh[M - x - 1].energy = e_r;
    }
    std::ofstream data("data.txt");
    solve<100>(mesh, k, t_s, t_e, tau, Co_max, data);
    data.close();

    return 0;
}


