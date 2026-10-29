#include "differentiation.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "common.h"

namespace {

// Таблица задания 3.4.
const std::vector<double> kX = {-1.0, 0.0, 1.0, 2.0, 3.0};
const std::vector<double> kY = {2.3562, 1.5708, 0.7854, 0.46365, 0.32175};
const double kXStar = 1.0;

}  // namespace

void task34() {
    printHeader("Задание 3.4. Численное дифференцирование");
    std::cout << "Вычислить первую и вторую производные таблично заданной\n"
                 "функции y_i = f(x_i), i = 0..4, в точке X* = 1.0.\n";

    const int n = static_cast<int>(kX.size());

    std::cout << "\nТаблица значений:\n i |  xi  |    yi\n---+------+---------\n";
    for (int i = 0; i < n; ++i)
        std::cout << " " << i << " | " << std::setw(4) << std::fixed
                  << std::setprecision(1) << kX[i] << " | " << std::setw(8)
                  << std::setprecision(5) << kY[i] << "\n";

    // X* совпадает с внутренним узлом, используем x_{i-1}, x_i, x_{i+1}
    int i = -1;
    for (int k = 1; k + 1 < n; ++k)
        if (std::fabs(kX[k] - kXStar) < 1e-12) i = k;

    std::cout << "\nТочка X* = " << std::setprecision(1) << kXStar
              << " совпадает с внутренним узлом x" << i
              << ", поэтому используются узлы x" << i - 1 << ", x" << i << ", x"
              << i + 1 << ".\n";

    const double hl = kX[i] - kX[i - 1];  
    const double hr = kX[i + 1] - kX[i];   

    // Левосторонняя и правосторонняя производные 1-й порядок точности
    const double dLeft = (kY[i] - kY[i - 1]) / hl;
    const double dRight = (kY[i + 1] - kY[i]) / hr;

    // Производная со 2-м порядком точности в узле x_i
    const double dCentral = dLeft * (kX[i] - kX[i + 1]) / (kX[i - 1] - kX[i + 1]) +
                            dRight * (kX[i] - kX[i - 1]) / (kX[i + 1] - kX[i - 1]);

    // Вторая производная
    const double d2 = 2.0 * (dRight - dLeft) / (kX[i + 1] - kX[i - 1]);

    std::cout << std::fixed;

    std::cout << "\n1) Левосторонняя производная на отрезке [x" << i - 1 << "; x"
              << i << "], 1-й порядок точности:\n   y'(X*) = (y" << i << " - y" << i - 1
              << ")/(x" << i << " - x" << i - 1 << ") = (" << std::setprecision(5)
              << kY[i] << " - " << kY[i - 1] << ")/(" << std::setprecision(1) << kX[i]
              << " - " << kX[i - 1] << ") = " << std::setprecision(5) << dLeft << "\n";

    std::cout << "\n2) Правосторонняя производная на отрезке [x" << i << "; x"
              << i + 1 << "], 1-й порядок точности:\n   y'(X*) = (y" << i + 1 << " - y"
              << i << ")/(x" << i + 1 << " - x" << i << ") = (" << std::setprecision(5)
              << kY[i + 1] << " - " << kY[i] << ")/(" << std::setprecision(1) << kX[i + 1]
              << " - " << kX[i] << ") = " << std::setprecision(5) << dRight << "\n";

    std::cout << "\n3) Производная со 2-м порядком точности:\n"
              << "   y'(X*) = (y" << i << " - y" << i - 1 << ")/(x" << i << " - x"
              << i - 1 << ")·(x" << i << " - x" << i + 1 << ")/(x" << i - 1 << " - x"
              << i + 1 << ") + (y" << i + 1 << " - y" << i << ")/(x" << i + 1 << " - x"
              << i << ")·(x" << i << " - x" << i - 1 << ")/(x" << i + 1 << " - x"
              << i - 1 << ")\n"
              << "   = (" << std::setprecision(5) << dLeft << ")·(-1)/(-2) + ("
              << dRight << ")·(1)/(2) = " << dCentral
              << "\n   (на равномерной сетке — полусумма значений)\n";

    std::cout << "\n4) Вторая производная:\n"
              << "   y''(X*) = 2*[(y" << i + 1 << " - y" << i << ")/(x" << i + 1 << " - x"
              << i << ") - (y" << i << " - y" << i - 1 << ")/(x" << i << " - x" << i - 1
              << ")] / (x" << i + 1 << " - x" << i - 1 << ")\n"
              << "   = 2*(" << dRight << " - (" << dLeft << ")) / 2.0 = " << d2 << "\n";

    std::cout << "\nОтвет: y'(X*) = " << dCentral << ",  y''(X*) = " << d2 << "\n";
}
