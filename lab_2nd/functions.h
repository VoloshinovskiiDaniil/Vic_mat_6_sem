#pragma once

#include <iostream>
#include <cmath>
#include <eigen-3.4.0\Eigen\Dense>
#include <fstream>
#include <array>

struct DVE 
{
    double density;
    double velocity;
    double energy;

    DVE() : density(0), energy(0), velocity(0) {}

    DVE(const Eigen::Matrix<double, 3, 1>& vec) : density(vec(0)), velocity(vec(1) / vec(0)), energy(vec(2) / vec(0)) {}

    operator Eigen::Vector3d() const 
    {
        return Eigen::Vector3d(density, velocity * density, energy * density);
    }
};

template <typename T, unsigned int NX> struct DVE_1D 
{
    const double X_L;
    const double X_R;
    const double h;
    std::array<T, NX> points;

    DVE_1D(const double x_left, const double x_right)
        : X_L(x_left), X_R(x_right), h((X_R - X_L) / (NX - 1)) {}

   T& operator[](std::size_t ind) 
    { 
        return points[ind];
    }

    const T& operator[](std::size_t ind) const
    { 
        return points[ind]; 
    }

    friend std::ostream& operator<<(std::ostream& os,
        const DVE_1D<DVE, NX>& mesh) 
    {
        for (unsigned int i = 0; i < NX; i++) 
        {
            double x = mesh.X_L + i * mesh.h;
            os << x << " " << mesh[i].density << " " << mesh[i].velocity << " "
                << mesh[i].energy << std::endl;
        }
        return os;
    }
};

template <unsigned int NX>
double get_max_abs_eigenvalue(const DVE_1D<DVE, NX>& mesh, const double c) 
{
    double maxAbsEigenvalue = 0;
    for (const auto& element : mesh.points) 
    {
        double abs_V_p_C = std::abs(element.velocity + c);
        double abs_V_m_C = std::abs(element.velocity - c);

        maxAbsEigenvalue = std::max(maxAbsEigenvalue, abs_V_p_C);
        maxAbsEigenvalue = std::max(maxAbsEigenvalue, abs_V_m_C);
    }
    return maxAbsEigenvalue;
}

template <unsigned int NX>
void solve(DVE_1D<DVE, NX> mesh, const double gamma,
    const double t_s, const double t_e,
    double t_increment, const double Co_max,
    std::ostream& data) 
{
    double t = t_s;
    double speed_of_sound = 0;
    double tau = t_increment;
    DVE_1D<DVE, NX> mesh_next = mesh;
    Eigen::Matrix<double, 3, 3> Omega_T{ {0, 0, gamma - 1},
                                        {0, 0, gamma - 1},
                                        {0, 0, gamma - 1} };
    Eigen::Matrix<double, 3, 3> Lambda = Eigen::Matrix<double, 3, 3>::Zero();
    Eigen::Matrix<double, 3, 3> A;

    while (t < t_e) 
    {
        double max_abs_eigenvalue =
            get_max_abs_eigenvalue<NX>(mesh, speed_of_sound);
        if (tau > mesh.h / max_abs_eigenvalue) 
        {
            tau = mesh.h * Co_max / max_abs_eigenvalue;
        }

        for (unsigned int x = 1; x < NX - 1; x++) 
        {
            speed_of_sound = sqrt(gamma * (gamma - 1) * mesh[x].energy);

            Omega_T(0, 0) = -mesh[x].velocity * speed_of_sound;
            Omega_T(0, 1) = speed_of_sound;
            Omega_T(1, 0) = -speed_of_sound * speed_of_sound;
            Omega_T(2, 0) = mesh[x].velocity * speed_of_sound;
            Omega_T(2, 1) = -speed_of_sound;

            Lambda(0, 0) = mesh[x].velocity + speed_of_sound;
            Lambda(1, 1) = mesh[x].velocity;
            Lambda(2, 2) = mesh[x].velocity - speed_of_sound;

            A = Omega_T.inverse() * Lambda * Omega_T;

            mesh_next[x] = DVE(Eigen::Vector3d(mesh[x]) - tau * A *
                (Eigen::Vector3d(mesh[x + 1]) -
                    Eigen::Vector3d(mesh[x - 1])) / (2 * mesh.h) +
                tau * (Omega_T.inverse() * Lambda.array().abs().matrix() * Omega_T) *
                (Eigen::Vector3d(mesh[x + 1]) - 2 * Eigen::Vector3d(mesh[x]) +
                    Eigen::Vector3d(mesh[x - 1])) / (2 * mesh.h));
        }
        mesh_next[0] = mesh_next[1];
        mesh_next[NX - 1] = mesh_next[NX - 2];
        mesh.points = mesh_next.points;
        data << "T: " << t << std::endl;
        data << mesh;
        data << std::endl;
        t += tau;
        tau = std::max(t_increment, tau);
    }
}
