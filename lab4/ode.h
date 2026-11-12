#pragma once

#include <functional>
#include <vector>

// Система двух ОДУ первого порядка, разрешённая относительно производных.
struct OdeSystem {
    std::function<double(double, double, double)> fy;  // правая часть для y'
    std::function<double(double, double, double)> fz;  // правая часть для z'
};

// Результат численного решения на равномерной сетке
struct OdeSolution {
    double h = 0.0;
    std::vector<double> x, y, z;        // узлы сетки и компоненты решения
    std::vector<double> dy, dz;         // приращения на каждом шаге
    std::vector<double> theta1, theta2; // контроль шага
    std::vector<double> g;              // значения fz(x_k, y_k, z_k) в узлах
};

OdeSolution solveEuler(const OdeSystem& sys, double x0, double y0, double z0,
                       double h, int steps);

OdeSolution solveRK4(const OdeSystem& sys, double x0, double y0, double z0,
                     double h, int steps, bool withTheta = false);

OdeSolution solveAdams4(const OdeSystem& sys, double x0, double y0, double z0,
                        double h, int steps);
