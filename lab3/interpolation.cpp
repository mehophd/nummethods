#include "interpolation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common.h"

namespace {

const double kXStar = 0.1;  // точка вычисления погрешности интерполяции

// (x - узел) для печати многочленов.
std::string factorStr(double node) {
    if (std::fabs(node) < 1e-12) return "(x)";
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << std::fabs(node);
    return node < 0.0 ? "(x + " + out.str() + ")" : "(x - " + out.str() + ")";
}

// w'(x_i) = произведение (x_i - x_j) по всем j != i
double omegaDerivative(const std::vector<double>& xs, int i) {
    double p = 1.0;
    for (int j = 0; j < static_cast<int>(xs.size()); ++j)
        if (j != i) p *= (xs[i] - xs[j]);
    return p;
}

// Лагранж в точке x
double lagrangeValue(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
    double sum = 0.0;
    for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
        double term = ys[i];
        for (int j = 0; j < static_cast<int>(xs.size()); ++j)
            if (j != i) term *= (x - xs[j]) / (xs[i] - xs[j]);
        sum += term;
    }
    return sum;
}

// Таблица dd[k][i] = f[x_i, ..., x_{i+k}]
std::vector<std::vector<double>> dividedDiffTable(const std::vector<double>& xs,
                                                  const std::vector<double>& ys) {
    const int n = static_cast<int>(xs.size());
    std::vector<std::vector<double>> dd(n);
    dd[0] = ys;
    for (int k = 1; k < n; ++k) {
        dd[k].assign(n - k, 0.0);
        for (int i = 0; i + k < n; ++i)
            dd[k][i] = (dd[k - 1][i + 1] - dd[k - 1][i]) / (xs[i + k] - xs[i]);
    }
    return dd;
}

// Ньютон в точке x
double newtonValue(const std::vector<double>& xs,
                   const std::vector<std::vector<double>>& dd, double x) {
    double result = dd[0][0];
    double product = 1.0;
    for (int k = 1; k < static_cast<int>(xs.size()); ++k) {
        product *= (x - xs[k - 1]);
        result += dd[k][0] * product;
    }
    return result;
}

void printWorkTable(const std::vector<double>& xs, const std::vector<double>& ys) {
    std::cout << " i |    Xi    |  Yi=f(Xi) |  w4'(Xi)  | Yi/w4'(Xi) |  X*-Xi\n";
    std::cout << "---+----------+-----------+-----------+------------+---------\n";
    for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
        const double w = omegaDerivative(xs, i);
        std::cout << " " << i << " | " << std::setprecision(2) << std::setw(6)
                  << xs[i] << " | " << std::setprecision(6) << std::setw(9) << ys[i]
                  << " | " << std::setw(9) << w << " | " << std::setw(10)
                  << ys[i] / w << " | " << std::setprecision(2) << std::setw(6)
                  << (kXStar - xs[i]) << "\n";
    }
}

void printLagrangeForm(const std::vector<double>& xs, const std::vector<double>& ys) {
    const int n = static_cast<int>(xs.size());
    std::cout << "\nМногочлен Лагранжа:\n";
    std::cout << "L3(x) = ";
    for (int i = 0; i < n; ++i) {
        const double coef = ys[i] / omegaDerivative(xs, i);
        if (i > 0) std::cout << "\n        ";
        std::cout << std::showpos << std::fixed << std::setprecision(5) << coef
                  << std::noshowpos << "*";
        bool firstFactor = true;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (!firstFactor) std::cout << "*";
            std::cout << factorStr(xs[j]);
            firstFactor = false;
        }
    }
    std::cout << "\n";
}

void printDividedDiffTable(const std::vector<double>& xs,
                           const std::vector<std::vector<double>>& dd) {
    const int n = static_cast<int>(xs.size());
    std::cout << "\nТаблица разделённых разностей:\n";
    std::cout << " i |    Xi    |  f(Xi)   | I порядок | II порядок | III порядок\n";
    std::cout << "---+----------+----------+-----------+------------+------------\n";
    for (int i = 0; i < n; ++i) {
        std::cout << " " << i << " | " << std::setw(8) << xs[i] << " | "
                  << std::setw(8) << dd[0][i] << " |";
        for (int k = 1; k < n; ++k) {
            if (i + k < n)
                std::cout << std::setw(11) << dd[k][i] << "  |";
            else
                std::cout << "             |";
        }
        std::cout << "\n";
    }
}

void printNewtonForm(const std::vector<double>& xs,
                     const std::vector<std::vector<double>>& dd) {
    const int n = static_cast<int>(xs.size());
    std::cout << "\nМногочлен Ньютона:\n";
    std::cout << "P3(x) = " << std::fixed << std::setprecision(5) << dd[0][0] << "\n";
    for (int k = 1; k < n; ++k) {
        std::cout << "        " << std::showpos << dd[k][0] << std::noshowpos << "*";
        for (int j = 0; j < k; ++j)
            std::cout << (j > 0 ? "*" : "") << factorStr(xs[j]);
        std::cout << "\n";
    }
}

}  // namespace

void task31() {
    printHeader("Задание 3.1. Интерполяция: многочлены Лагранжа и Ньютона");
    std::cout << "Функция y = arccos(x) задана в четырёх узлах. Требуется построить\n"
                 "интерполяционные многочлены Лагранжа и Ньютона, проходящие через\n"
                 "точки {Xi, Yi}, и вычислить погрешность интерполяции в точке X* = 0.1.\n";

    const double fStar = std::acos(kXStar);

    struct Variant {
        const char* title;
        std::vector<double> xs;
    };
    const std::vector<Variant> variants = {
        {"а) узлы Xi = -0.4; -0.1; 0.2; 0.5", {-0.4, -0.1, 0.2, 0.5}},
        {"б) узлы Xi = -0.4;  0.0; 0.2; 0.5", {-0.4, 0.0, 0.2, 0.5}},
    };

    std::cout << std::fixed;

    for (const auto& var : variants) {
        printSubHeader(var.title);

        // Значения функции y = arccos(x) в узлах
        std::vector<double> ys;
        ys.reserve(var.xs.size());
        for (double x : var.xs) ys.push_back(std::acos(x));

        printWorkTable(var.xs, ys);
        printLagrangeForm(var.xs, ys);

        const auto dd = dividedDiffTable(var.xs, ys);
        printDividedDiffTable(var.xs, dd);
        printNewtonForm(var.xs, dd);

        // Значения многочленов и погрешность в точке X*
        const double lag = lagrangeValue(var.xs, ys, kXStar);
        const double newt = newtonValue(var.xs, dd, kXStar);
        const double error = std::fabs(lag - fStar);

        // R <= M1/4! * |w4(X*)|,  M1 = max |f'(t)| = max 1/sqrt(1-t^2) на [x0; x3]
        const double m1 =
            std::max(1.0 / std::sqrt(1.0 - var.xs.front() * var.xs.front()),
                     1.0 / std::sqrt(1.0 - var.xs.back() * var.xs.back()));
        double omegaStar = 1.0;
        for (double x : var.xs) omegaStar *= (kXStar - x);
        const double bound = m1 * std::fabs(omegaStar) / 24.0;

        std::cout << "\nЗначения в точке X* = " << kXStar << ":\n"
                  << std::setprecision(6) << "  L3(X*) = " << lag << "\n"
                  << "  P3(X*) = " << newt << "   \n"
                  << "  f(X*)  = arccos(0.1) = " << fStar << "\n"
                  << "  Абсолютная погрешность Δ = |L3(X*) - f(X*)| = " << error << "\n"
                  << "  M1 = max|f'(t)| = " << m1 << ",  |w4(X*)| = "
                  << std::fabs(omegaStar) << "\n"
                  << "  Априорная оценка (3.9): R <= M1/4!*|w4(X*)| = " << bound << "\n";
    }
}
