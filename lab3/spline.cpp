#include "spline.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "common.h"

namespace {

// Таблица функции y = arccos(x)
const std::vector<double> kX = {-0.4, -0.1, 0.2, 0.5, 0.8};
const std::vector<double> kF = {1.9823, 1.6710, 1.3694, 1.0472, 0.6435};
const double kXStar = 0.1;

// Метод прогонки
// a поддиагональ (a[0] не используется), b главная диагональ
// c наддиагональ (c[n-1] не используется), d правая часть
std::vector<double> thomasSolve(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d) {
    const int n = static_cast<int>(b.size());
    std::vector<double> p(n, 0.0), q(n, 0.0), x(n, 0.0);
    p[0] = (n > 1) ? -c[0] / b[0] : 0.0;
    q[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        const double denom = b[i] + a[i] * p[i - 1];
        p[i] = (i < n - 1) ? -c[i] / denom : 0.0;
        q[i] = (d[i] - a[i] * q[i - 1]) / denom;
    }
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = p[i] * x[i + 1] + q[i];
    return x;
}

}  // namespace

void task32() {
    printHeader("Задание 3.2. Кубический сплайн (нулевая кривизна на концах)");
    std::cout << "Для таблично заданной функции построить кубический сплайн, считая\n"
                 "кривизну сплайна равной нулю при x = x0 и x = x4 (естественные\n"
                 "граничные условия), и вычислить значение функции в точке X* = 0.1.\n";

    const int n = static_cast<int>(kX.size()) - 1;  // число отрезков 4

    std::cout << "\nТаблица значений:\n i |   xi   |   fi\n---+--------+--------\n";
    for (int i = 0; i <= n; ++i)
        std::cout << " " << i << " | " << std::setw(6) << std::fixed
                  << std::setprecision(2) << kX[i] << " | " << std::setw(7)
                  << std::setprecision(4) << kF[i] << "\n";

    // h[i] = x[i] - x[i-1]
    std::vector<double> h(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) h[i] = kX[i] - kX[i - 1];

    //   h[i]*c[i] + 2*(h[i]+h[i+1])*c[i+1] + h[i+1]*c[i+2] = 3*((f[i+1]-f[i])/h[i+1] - (f[i]-f[i-1])/h[i]),  i = 1..n-1
    const int m = n - 1;  // число неизвестных
    std::vector<double> a(m, 0.0), b(m, 0.0), c(m, 0.0), d(m, 0.0);
    for (int i = 1; i <= n - 1; ++i) {
        const int row = i - 1;
        if (row > 0) a[row] = h[i];          // коэффициент при c[i] (для i >= 2)
        b[row] = 2.0 * (h[i] + h[i + 1]);    // коэффициент при c[i+1]
        if (row < m - 1) c[row] = h[i + 1];  // коэффициент при c[i+2]
        d[row] = 3.0 * ((kF[i + 1] - kF[i]) / h[i + 1] - (kF[i] - kF[i - 1]) / h[i]);
    }

    std::cout << "\nСистема относительно c2..c4 (c1 = c5 = 0):\n";
    for (int row = 0; row < m; ++row) {
        bool first = true;
        auto printTerm = [&first](double coef, int cIndex) {
            if (std::fabs(coef) < 1e-15) return;
            if (!first) std::cout << " + ";
            std::cout << std::fixed << std::setprecision(2) << coef << "*c" << cIndex;
            first = false;
        };
        if (row > 0) printTerm(a[row], row + 1);
        printTerm(b[row], row + 2);
        if (row < m - 1) printTerm(c[row], row + 3);
        std::cout << std::setprecision(5) << " = " << std::setw(9) << d[row] << "\n";
    }

    const std::vector<double> cSol = thomasSolve(a, b, c, d);  // c2, c3, c4

    std::cout << "\nРешение системы (метод прогонки): ";
    for (int i = 0; i < m; ++i)
        std::cout << "c" << i + 2 << " = " << std::setprecision(5) << cSol[i]
                  << (i + 1 < m ? ";  " : "\n");

    // Полный массив коэффициентов c_i, i = 1..n+1 (c1 = c_{n+1} = 0).
    std::vector<double> cFull(n + 2, 0.0);
    for (int i = 0; i < m; ++i) cFull[i + 2] = cSol[i];

    // Коэффициенты сплайна
    //   a_i = f_{i-1};
    //   b_i = (f_i - f_{i-1})/h_i - h_i*(2c_i + c_{i+1})/3;
    //   d_i = (c_{i+1} - c_i)/(3*h_i).
    std::vector<double> A(n + 1, 0.0), B(n + 1, 0.0), C(n + 1, 0.0), D(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) {
        A[i] = kF[i - 1];
        B[i] = (kF[i] - kF[i - 1]) / h[i] - h[i] * (2.0 * cFull[i] + cFull[i + 1]) / 3.0;
        C[i] = cFull[i];
        D[i] = (cFull[i + 1] - cFull[i]) / (3.0 * h[i]);
    }

    std::cout << "\nКоэффициенты сплайна, на [x_{i-1}; x_i]:\n"
              << "S_i(x) = a_i + b_i*(x-x_{i-1}) + c_i*(x-x_{i-1})^2 + d_i*(x-x_{i-1})^3\n";
    std::cout << " i |    отрезок    |    a_i   |    b_i   |    c_i   |    d_i\n";
    std::cout << "---+---------------+----------+----------+----------+----------\n";
    for (int i = 1; i <= n; ++i) {
        std::cout << " " << i << " | [" << std::setw(5) << std::setprecision(2)
                  << kX[i - 1] << "; " << std::setw(5) << kX[i] << "] | "
                  << std::setw(8) << std::setprecision(5) << A[i] << " | "
                  << std::setw(8) << B[i] << " | " << std::setw(8) << C[i] << " | "
                  << std::setw(8) << D[i] << "\n";
    }

    // Поиск отрезка, содержащего X*, и вычисление S(X*)
    int seg = 1;
    while (seg < n && kXStar > kX[seg]) ++seg;  // X* в (x_{seg-1}; x_seg]
    const double t = kXStar - kX[seg - 1];
    const double sVal = A[seg] + B[seg] * t + C[seg] * t * t + D[seg] * t * t * t;

    std::cout << "\nТочка X* = 0.1 принадлежит отрезку [x" << seg - 1 << "; x" << seg
              << "] = [" << std::setprecision(2) << kX[seg - 1] << "; " << kX[seg]
              << "] (i = " << seg << "),\n"
              << "t = X* - x" << seg - 1 << " = " << t << "\n";
    std::cout << "S_" << seg << "(X*) = " << std::setprecision(5) << A[seg] << " + ("
              << B[seg] << ")*" << std::setprecision(2) << t << " + ("
              << std::setprecision(5) << C[seg] << ")*" << std::setprecision(3) << t * t
              << " + (" << std::setprecision(5) << D[seg] << ")*" << std::setprecision(4)
              << t * t * t << " = " << std::setprecision(5) << sVal << "\n";

    const double fExact = std::acos(kXStar);
    std::cout << "Справочно: табличные данные соответствуют y = arccos(x):\n"
              << "  arccos(0.1) = " << std::setprecision(6) << fExact
              << ",  |S(X*) - f(X*)| = " << std::fabs(sVal - fExact) << "\n";
}
