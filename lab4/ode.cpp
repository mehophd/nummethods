#include "ode.h"

#include <algorithm>
#include <cmath>

// Явный метод Эйлера (4.2): y_{k+1} = y_k + h*fy(x_k, y_k, z_k) и аналогично для z. Первый порядок точности
OdeSolution solveEuler(const OdeSystem& sys, double x0, double y0, double z0,
                       double h, int steps) {
    OdeSolution s;
    s.h = h;
    s.x.reserve(steps + 1);
    s.y.reserve(steps + 1);
    s.z.reserve(steps + 1);
    s.x.push_back(x0);
    s.y.push_back(y0);
    s.z.push_back(z0);
    for (int k = 0; k < steps; ++k) {
        const double x = s.x.back(), y = s.y.back(), z = s.z.back();
        const double dy = h * sys.fy(x, y, z);
        const double dz = h * sys.fz(x, y, z);
        s.dy.push_back(dy);
        s.dz.push_back(dz);
        s.x.push_back(x + h);
        s.y.push_back(y + dy);
        s.z.push_back(z + dz);
    }
    return s;
}

// Метод Рунге-Кутты 4-го порядка для системы
//   K_i = h*fy(x + a_i h, y + ..., z + ...),  L_i = h*fz(...),
//   dy = (K1 + 2K2 + 2K3 + K4)/6,  dz = (L1 + 2L2 + 2L3 + L4)/6.
OdeSolution solveRK4(const OdeSystem& sys, double x0, double y0, double z0,
                     double h, int steps, bool withTheta) {
    OdeSolution s;
    s.h = h;
    s.x.reserve(steps + 1);
    s.y.reserve(steps + 1);
    s.z.reserve(steps + 1);
    s.x.push_back(x0);
    s.y.push_back(y0);
    s.z.push_back(z0);
    for (int k = 0; k < steps; ++k) {
        const double x = s.x.back(), y = s.y.back(), z = s.z.back();
        const double K1 = h * sys.fy(x, y, z);
        const double L1 = h * sys.fz(x, y, z);
        const double K2 = h * sys.fy(x + h / 2, y + K1 / 2, z + L1 / 2);
        const double L2 = h * sys.fz(x + h / 2, y + K1 / 2, z + L1 / 2);
        const double K3 = h * sys.fy(x + h / 2, y + K2 / 2, z + L2 / 2);
        const double L3 = h * sys.fz(x + h / 2, y + K2 / 2, z + L2 / 2);
        const double K4 = h * sys.fy(x + h, y + K3, z + L3);
        const double L4 = h * sys.fz(x + h, y + K3, z + L3);
        const double dy = (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
        const double dz = (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;
        s.dy.push_back(dy);
        s.dz.push_back(dz);
        if (withTheta) {
            // theta = (K3-K2)/(K2-K1); если theta
            auto theta = [](double a, double b, double c) {
                return std::fabs(b - a) > 1e-14 ? (c - b) / (b - a) : 0.0;
            };
            s.theta1.push_back(theta(K1, K2, K3));
            s.theta2.push_back(theta(L1, L2, L3));
        }
        s.x.push_back(x + h);
        s.y.push_back(y + dy);
        s.z.push_back(z + dz);
    }
    return s;
}

// Метод Адамса 4-го порядка
//   y_{k+1} = y_k + h/24 * (55 f_k - 59 f_{k-1} + 37 f_{k-2} - 9 f_{k-3})
OdeSolution solveAdams4(const OdeSystem& sys, double x0, double y0, double z0,
                        double h, int steps) {
    const int start = std::min(3, steps);
    OdeSolution s = solveRK4(sys, x0, y0, z0, h, start, false);
    s.dy.reserve(steps);
    s.dz.reserve(steps);
    for (int k = start; k < steps; ++k) {
        const double gy = s.z[k];  // fy(x_k, y_k, z_k) = z_k
        const double gk = sys.fz(s.x[k], s.y[k], s.z[k]);
        const double gk1 = sys.fz(s.x[k - 1], s.y[k - 1], s.z[k - 1]);
        const double gk2 = sys.fz(s.x[k - 2], s.y[k - 2], s.z[k - 2]);
        const double gk3 = sys.fz(s.x[k - 3], s.y[k - 3], s.z[k - 3]);
        const double dy =
            h / 24.0 * (55.0 * gy - 59.0 * s.z[k - 1] + 37.0 * s.z[k - 2] -
                        9.0 * s.z[k - 3]);
        const double dz = h / 24.0 * (55.0 * gk - 59.0 * gk1 + 37.0 * gk2 -
                                      9.0 * gk3);
        s.dy.push_back(dy);
        s.dz.push_back(dz);
        s.x.push_back(s.x[k] + h);
        s.y.push_back(s.y[k] + dy);
        s.z.push_back(s.z[k] + dz);
    }
    // Значения правой части z' = fz 
    s.g.resize(s.x.size());
    for (std::size_t i = 0; i < s.x.size(); ++i)
        s.g[i] = sys.fz(s.x[i], s.y[i], s.z[i]);
    return s;
}
