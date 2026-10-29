#include "lsm.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"

namespace {

// Таблица функции y = arccos(x)
const std::vector<double> kX = {-0.7, -0.4, -0.1, 0.2, 0.5, 0.8};
const std::vector<double> kY = {2.3462, 1.9823, 1.6710, 1.3694, 1.0472, 0.6435};
const char* kPlotFile = "lsm_plot_data.txt";

std::vector<double> gaussSolve(std::vector<std::vector<double>> a,
                               std::vector<double> b) {
    const int n = static_cast<int>(b.size());
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row)
            if (std::fabs(a[row][col]) > std::fabs(a[pivot][col])) pivot = row;
        std::swap(a[col], a[pivot]);
        std::swap(b[col], b[pivot]);
        for (int row = col + 1; row < n; ++row) {
            const double factor = a[row][col] / a[col][col];
            for (int j = col; j < n; ++j) a[row][j] -= factor * a[col][j];
            b[row] -= factor * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) sum -= a[i][j] * x[j];
        x[i] = sum / a[i][i];
    }
    return x;
}

// coeffs[0] + coeffs[1]*x + ...
double polyEval(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i)
        result = result * x + coeffs[i];
    return result;
}

}  // namespace

void task33() {
    printHeader("Задание 3.3. Метод наименьших квадратов");
    std::cout << "Для таблично заданной функции найти путем решения нормальной\n"
                 "системы МНК приближающие многочлены 1-й и 2-й степени, вычислить\n"
                 "суммы квадратов ошибок и построить графики.\n";

    const int N = static_cast<int>(kX.size()) - 1;  // 5

    std::cout << "\nТаблица значений:\n i |   xi   |   yi\n---+--------+---------\n";
    for (int j = 0; j <= N; ++j)
        std::cout << " " << j << " | " << std::setw(6) << std::fixed
                  << std::setprecision(2) << kX[j] << " | " << std::setw(7)
                  << std::setprecision(4) << kY[j] << "\n";

    // Степенные суммы для нормальной системы
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
    double sy = 0.0, sxy = 0.0, sx2y = 0.0;
    for (int j = 0; j <= N; ++j) {
        const double x = kX[j], y = kY[j];
        s0 += 1.0;
        s1 += x;
        s2 += x * x;
        s3 += x * x * x;
        s4 += x * x * x * x;
        sy += y;
        sxy += x * y;
        sx2y += x * x * y;
    }

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\nСтепенные суммы:\n"
              << "  N = " << N + 1 << ",  Σx = " << s1 << ",  Σx^2 = " << s2
              << ",  Σx^3 = " << s3 << ",  Σx^4 = " << s4 << "\n"
              << "  Σy = " << sy << ",  Σxy = " << sxy << ",  Σx^2*y = " << sx2y << "\n";

    // F1(x) = a0 + a1*x
    std::cout << "\nа) Нормальная система МНК для многочлена 1-й степени:\n";
    std::cout << std::setprecision(2) << "  " << std::setw(5) << s0 << "*a0 + "
              << std::setw(6) << s1 << "*a1 = " << std::setprecision(5) << sy << "\n"
              << std::setprecision(2) << "  " << std::setw(5) << s1 << "*a0 + "
              << std::setw(6) << s2 << "*a1 = " << std::setprecision(5) << sxy << "\n";

    const std::vector<double> c1 = gaussSolve({{s0, s1}, {s1, s2}}, {sy, sxy});

    std::cout << "Решение: a0 = " << std::setprecision(5) << c1[0]
              << ",  a1 = " << c1[1] << "\n";
    std::cout << "F1(x) = " << c1[0] << " " << std::showpos << c1[1]
              << std::noshowpos << "*x\n";

    // F2(x) = a0 + a1*x + a2*x^2
    std::cout << "\nб) Нормальная система МНК для многочлена 2-й степени:\n";
    std::cout << std::setprecision(4) << "  " << std::setw(6) << s0 << "*a0 + "
              << std::setw(7) << s1 << "*a1 + " << std::setw(7) << s2 << "*a2 = "
              << std::setprecision(5) << sy << "\n"
              << std::setprecision(4) << "  " << std::setw(6) << s1 << "*a0 + "
              << std::setw(7) << s2 << "*a1 + " << std::setw(7) << s3 << "*a2 = "
              << std::setprecision(5) << sxy << "\n"
              << std::setprecision(4) << "  " << std::setw(6) << s2 << "*a0 + "
              << std::setw(7) << s3 << "*a1 + " << std::setw(7) << s4 << "*a2 = "
              << std::setprecision(5) << sx2y << "\n";

    const std::vector<double> c2 =
        gaussSolve({{s0, s1, s2}, {s1, s2, s3}, {s2, s3, s4}}, {sy, sxy, sx2y});

    std::cout << "Решение: a0 = " << std::setprecision(5) << c2[0]
              << ",  a1 = " << c2[1] << ",  a2 = " << c2[2] << "\n";
    std::cout << "F2(x) = " << c2[0] << " " << std::showpos << c2[1] << "*x "
              << c2[2] << std::noshowpos << "*x^2\n";

    // Значения в узлах и суммы квадратов ошибок
    std::cout << "\nЗначения приближающих многочленов в узлах:\n";
    std::cout << " j |   xj   |   yj    |  F1(xj)  |  F2(xj)\n";
    std::cout << "---+--------+---------+----------+---------\n";
    double phi1 = 0.0, phi2 = 0.0;
    for (int j = 0; j <= N; ++j) {
        const double f1 = polyEval(c1, kX[j]);
        const double f2 = polyEval(c2, kX[j]);
        phi1 += (f1 - kY[j]) * (f1 - kY[j]);
        phi2 += (f2 - kY[j]) * (f2 - kY[j]);
        std::cout << " " << j << " | " << std::setprecision(2) << std::setw(6)
                  << kX[j] << " | " << std::setprecision(4) << std::setw(7) << kY[j]
                  << " | " << std::setprecision(5) << std::setw(8) << f1 << " | "
                  << std::setw(8) << f2 << "\n";
    }

    std::cout << "\nСуммы квадратов ошибок:\n"
              << "  Φ1 = Σ(F1(xj) - yj)^2 = " << phi1 << "\n"
              << "  Φ2 = Σ(F2(xj) - yj)^2 = " << phi2 << "\n";

    std::ofstream out(kPlotFile);
    out << "# Задание 3.3: данные для построения графика\n";
    out << "# Блок 1: табличные точки\n# x y\n";
    for (int j = 0; j <= N; ++j) out << kX[j] << " " << kY[j] << "\n";
    out << "\n# Блок 2: приближающие многочлены на плотной сетке\n# x F1 F2\n";
    const double xmin = kX.front() - 0.15, xmax = kX.back() + 0.15;
    const int steps = 120;
    for (int s = 0; s <= steps; ++s) {
        const double x = xmin + (xmax - xmin) * s / steps;
        out << x << " " << polyEval(c1, x) << " " << polyEval(c2, x) << "\n";
    }
    out.close();
}
