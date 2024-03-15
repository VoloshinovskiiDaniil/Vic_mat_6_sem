#pragma once
#include<iostream>
#include<vector>
#include<cmath>


struct mesh 
{
    double x;
    double u;
};

struct mesh_with_time
{
    std::vector<mesh> layer;
    double t;
};

template<typename Lable>
mesh_with_time Prototype(const double L, const double h, const Lable& u_0) 
{
    std::vector<mesh> data;
    const unsigned int N = 20;

    mesh m;
    m.x = 0;
    m.u = u_0(0);
    data.push_back(m);
    for (unsigned int i = 1; i < N; i++) 
    {
        mesh m;
        m.x = data[i - 1].x + h;
        m.u = u_0(m.x);
        data.push_back(m);
    }

    return { data, 0.0 };
}


mesh_with_time L_angle_method(mesh_with_time input_data, const double a, const double tau, const double h)
{
    mesh_with_time data = input_data;
    unsigned int N = (input_data.layer).size();
    data.t += tau;
    for (unsigned int i = 1; i < (data.layer).size() - 1; i++) 
    {
        data.layer[i].u -= a * tau / h * (input_data.layer[i].u - input_data.layer[i - 1].u);
    }

    data.layer[0].u -= a * tau / h * (input_data.layer[0].u - input_data.layer[N - 1].u);
    data.layer[N - 1].u -= a * tau / h * (input_data.layer[N - 1].u - input_data.layer[N - 2].u);

    return data;
}

mesh_with_time LW_method(mesh_with_time input_data, const double a, const double tau, const double h)
{
    mesh_with_time data = input_data;
    unsigned int N = (input_data.layer).size();
    data.t += tau;
    const double Co = a * tau / h;
    for (unsigned int i = 1; i < (data.layer).size() - 1; i++) {
        data.layer[i].u = (Co * Co + Co) * input_data.layer[i - 1].u / 2.0 +
            (1 - Co * Co) * input_data.layer[i].u + (Co * Co - Co) * input_data.layer[i + 1].u / 2.0;
    }

    data.layer[0].u = (Co * Co + Co) * input_data.layer[N - 1].u / 2.0 +
        (1 - Co * Co) * input_data.layer[0].u + (Co * Co - Co) * input_data.layer[1].u / 2.0;

    data.layer[N - 1].u = (Co * Co + Co) * input_data.layer[N - 2].u / 2.0 +
        (1 - Co * Co) * input_data.layer[N - 1].u + (Co * Co - Co) * input_data.layer[0].u / 2.0;

    return data;
}