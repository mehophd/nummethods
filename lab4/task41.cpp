#include "task41.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"
#include "ode.h"

namespace {

//   y'' - (1/sqrt(x)) y' + (x + sqrt(x) - 8)/(4 x^2) y = 0,
//   y(1) = 2e, y'(1) = 2e, x из [1; 2].
// Точное решение: y = (x^2 + 1/x) e^{sqrt(x)}
const double kA = 1.0;  // левый конец отрезка
const double kB = 2.0;  // правый конец отрезка
const double kY0 = 2.0 * std::exp(1.0);  // y(1) = 2e
const double kZ0 = 2.0 * std::exp(1.0);  // y'(1) = 2e

// Точное решение y = (x^2 + 1/x) e^{sqrt(x)}.
double yExact(double x) {
    return (x * x + 1.0 / x) * std::exp(std::sqrt(x));
}

// Точная производная z = y':
// y' = e^{sqrt x} [ 2x - 1/x^2 + (x^2 + 1/x)/(2 sqrt x) ]
double zExact(double x) {
    const double s = std::sqrt(x);
    return std::exp(s) *
           (2.0 * x - 1.0 / (x * x) + (x * x + 1.0 / x) / (2.0 * s));
}

// Система после замены z = y':
//   y' = z,   z' = (1/sqrt x) z - (x + sqrt x - 8)/(4 x^2) y
OdeSystem makeSystem() {
    OdeSystem sys;
    sys.fy = [](double, double, double z) { return z; };
    sys.fz = [](double x, double y, double z) {
        const double s = std::sqrt(x);
        return z / s - (x + s - 8.0) / (4.0 * x * x) * y;
    };
    return sys;
}

// Максимальная фактическая погрешность компоненты y
double maxErrorY(const OdeSolution& s) {
    double m = 0.0;
    for (std::size_t i = 0; i < s.x.size(); ++i)
        m = std::max(m, std::fabs(s.y[i] - yExact(s.x[i])));
    return m;
}

// Таблица метода Эйлера
void printEulerTable(const OdeSolution& s) {
    printSubHeader("Метод Эйлера");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5},  {"y_k", 14}, {"z_k", 14}, {"Δy_k", 14},
        {"Δz_k", 14}, {"y_ист", 14}, {"ε_y", 11}, {"ε_z", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (std::size_t k = 0; k < s.x.size(); ++k) {
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(s.x[k]) << "| "
            << fCell(s.y[k]) << "| " << fCell(s.z[k]) << "| ";
        if (k + 1 < s.x.size()) {
            row << fCell(s.dy[k]) << "| " << fCell(s.dz[k]) << "| ";
        } else {
            row << std::string(14, ' ') << "| " << std::string(14, ' ') << "| ";
        }
        row << fCell(yExact(s.x[k])) << "| "
            << eCell(s.y[k] - yExact(s.x[k])) << "| "
            << eCell(s.z[k] - zExact(s.x[k]));
        std::cout << row.str() << "\n";
    }
}

// Таблица метода Рунге-Кутты 4-го порядка 
void printRK4Table(const OdeSolution& s) {
    printSubHeader("Метод Рунге-Кутты 4-го порядка");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5},  {"y_k", 14}, {"z_k", 14}, {"Δy_k", 14},
        {"Δz_k", 14}, {"θ1", 9}, {"θ2", 9}, {"y_ист", 14}, {"ε_y", 11},
        {"ε_z", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (std::size_t k = 0; k < s.x.size(); ++k) {
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(s.x[k]) << "| "
            << fCell(s.y[k]) << "| " << fCell(s.z[k]) << "| ";
        if (k + 1 < s.x.size()) {
            row << fCell(s.dy[k]) << "| " << fCell(s.dz[k]) << "| ";
        } else {
            row << std::string(14, ' ') << "| " << std::string(14, ' ') << "| ";
        }
        if (k + 1 < s.x.size() && !s.theta1.empty() && k < s.theta1.size()) {
            std::ostringstream t1, t2;
            t1 << std::fixed << std::setprecision(5) << std::setw(9)
               << s.theta1[k];
            t2 << std::fixed << std::setprecision(5) << std::setw(9)
               << s.theta2[k];
            row << t1.str() << "| " << t2.str() << "| ";
        } else {
            row << std::string(9, ' ') << "| " << std::string(9, ' ') << "| ";
        }
        row << fCell(yExact(s.x[k])) << "| "
            << eCell(s.y[k] - yExact(s.x[k])) << "| "
            << eCell(s.z[k] - zExact(s.x[k]));
        std::cout << row.str() << "\n";
    }
}

// Таблица метода Адамса
void printAdamsTable(const OdeSolution& s) {
    printSubHeader("Метод Адамса 4-го порядка");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5}, {"y_k", 14}, {"z_k", 14}, {"g_k", 14},
        {"y_ист", 14}, {"ε_y", 11}, {"ε_z", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (std::size_t k = 0; k < s.x.size(); ++k) {
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(s.x[k]) << "| "
            << fCell(s.y[k]) << "| " << fCell(s.z[k]) << "| "
            << fCell(s.g[k]) << "| " << fCell(yExact(s.x[k])) << "| "
            << eCell(s.y[k] - yExact(s.x[k])) << "| "
            << eCell(s.z[k] - zExact(s.x[k]));
    }
}

// Сводная таблица решений трех методов.
void printSummary(const OdeSolution& eu, const OdeSolution& rk,
                  const OdeSolution& ad) {
    printSubHeader("Сводная таблица и сравнение с точным решением");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5}, {"y_ист", 14}, {"y_Эйлер", 14}, {"ε_Эйлер", 11},
        {"y_РК4", 14}, {"ε_РК4", 11}, {"y_Адамс", 14}, {"ε_Адамс", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (std::size_t k = 0; k < eu.x.size(); ++k) {
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(eu.x[k]) << "| "
            << fCell(yExact(eu.x[k])) << "| " << fCell(eu.y[k]) << "| "
            << eCell(eu.y[k] - yExact(eu.x[k])) << "| " << fCell(rk.y[k])
            << "| " << eCell(rk.y[k] - yExact(eu.x[k])) << "| "
            << fCell(ad.y[k]) << "| " << eCell(ad.y[k] - yExact(eu.x[k]));
        std::cout << row.str() << "\n";
    }
    std::cout << "Максимум фактической погрешности |y_k - y(x_k)|:\n"
              << "  метод Эйлера:            " << std::scientific
              << std::setprecision(3) << maxErrorY(eu) << "\n"
              << "  метод Рунге-Кутты 4:     " << maxErrorY(rk) << "\n"
              << "  метод Адамса 4:          " << maxErrorY(ad) << "\n"
              << std::fixed;
}

// Таблица метода Рунге-Ромберга для одного метода:
//   R_h = (y_2h - y_h)/(2^p - 1)  - оценка погрешности решения с шагом h,
//   y_уточн = y_h + R_h = (2^p y_h - y_2h)/(2^p - 1)
void printRR(const std::string& name, const OdeSolution& fine,
             const OdeSolution& coarse, int p) {
    const int n2 = static_cast<int>(coarse.x.size()) - 1;
    const double denom = std::pow(2.0, p) - 1.0;
    printSubHeader("Рунге-Ромберг: " + name + ", p = " + std::to_string(p) +
                   " (сравнение шагов h и 2h)");
    const std::vector<std::pair<std::string, int>> cols = {
        {"x_k", 5}, {"y_h", 14}, {"y_2h", 14}, {"R_h", 11}, {"y_уточн", 14},
        {"|y_h-y_ист|", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    double maxR = 0.0, maxT = 0.0;
    for (int i = 0; i <= n2; ++i) {
        const int j = 2 * i;  // общий узел сеток h и 2h
        const double R = (fine.y[j] - coarse.y[i]) / denom;
        const double ref = fine.y[j] + R;
        const double terr = std::fabs(fine.y[j] - yExact(fine.x[j]));
        maxR = std::max(maxR, std::fabs(R));
        maxT = std::max(maxT, terr);
        std::ostringstream row;
        row << " " << xCell(fine.x[j]) << "| " << fCell(fine.y[j]) << "| "
            << fCell(coarse.y[i]) << "| " << eCell(R) << "| " << fCell(ref)
            << "| " << eCell(terr);
        std::cout << row.str() << "\n";
    }
    std::cout << "max|R_h| = " << std::scientific << std::setprecision(3)
              << maxR << ",  max|y_h - y_ист| = " << maxT << std::fixed << "\n";
}

void writePlotData(const OdeSolution& eu, const OdeSolution& rk,
                   const OdeSolution& ad) {
    std::ofstream out("cauchy_plot_data.txt");
    out << "# x y_exact y_euler y_rk4 y_adams\n";
    for (std::size_t k = 0; k < eu.x.size(); ++k) {
        out << std::fixed << std::setprecision(6) << eu.x[k] << " "
            << yExact(eu.x[k]) << " " << eu.y[k] << " " << rk.y[k] << " "
            << ad.y[k] << "\n";
    }
}

}  // namespace

void task41(double hIn) {
    printHeader("Задание 4.1. Задача Коши для ОДУ 2-го порядка");
    const OdeSystem sys = makeSystem();

    int n = static_cast<int>(std::lround((kB - kA) / hIn));
    if (n < 1) n = 1;
    const double h = (kB - kA) / n;
    if (std::fabs(h - hIn) > 1e-12)
        std::cout << "(шаг скорректирован так, чтобы (b-a)/h было целым: h = "
                  << h << ")\n";
    std::cout << "Сетка: x_k = 1 + k*h, k = 0.." << n << ", h = " << h << "\n";

    // Решения с шагом h
    const OdeSolution eu = solveEuler(sys, kA, kY0, kZ0, h, n);
    const OdeSolution rk = solveRK4(sys, kA, kY0, kZ0, h, n, true);
    const OdeSolution ad = solveAdams4(sys, kA, kY0, kZ0, h, n);

    printEulerTable(eu);
    printRK4Table(rk);
    printAdamsTable(ad);
    printSummary(eu, rk, ad);

    // Решения с удвоенным шагом 2h для оценки погрешности по Рунге-Ромбергу.
    if (n % 2 == 0) {
        const int n2 = n / 2;
        const double h2 = 2.0 * h;
        const OdeSolution eu2 = solveEuler(sys, kA, kY0, kZ0, h2, n2);
        const OdeSolution rk2 = solveRK4(sys, kA, kY0, kZ0, h2, n2);
        const OdeSolution ad2 = solveAdams4(sys, kA, kY0, kZ0, h2, n2);
        printRR("метод Эйлера", eu, eu2, 1);
        printRR("метод Рунге-Кутты 4", rk, rk2, 4);
        printRR("метод Адамса 4", ad, ad2, 4);
    } else {
        std::cout << "\nДля оценки по Рунге-Ромбергу нужно четное число "
                     "шагов; таблицы Рунге-Ромберга пропущены.\n";
    }

    writePlotData(eu, rk, ad);
}
