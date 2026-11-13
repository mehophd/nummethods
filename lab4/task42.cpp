#include "task42.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"
#include "ode.h"

namespace {

//   x y'' - (2x + 1) y' + (x + 1) y = 0,
//   y'(0) = 1,   y'(1) - 2 y(1) = 0,   x из [0; 1].
// Точное решение: y = e^x (x^2 + 1)
const double kA = 0.0;  // левый конец отрезка
const double kB = 1.0;  // правый конец отрезок

double yExact(double x) { return std::exp(x) * (x * x + 1.0); }

// y' = e^x (x + 1)^2
double zExact(double x) {
    const double e = std::exp(x);
    return e * (x + 1.0) * (x + 1.0);
}

// В точке x = 0 уравнение вырождается в -y'(0) + y(0) = 0, т.е. y(0) = y'(0). 
// Поэтому условие y'(0) = 1 равносильно y(0) = 1, а
// 1/x снимается заменой v = (y' - y)/x:
//   y' = y + x v,    v' = v.
OdeSystem makeRegularizedSystem() {
    OdeSystem sys;
    // y' = y + x v ; v' = v
    sys.fy = [](double x, double y, double v) { return y + x * v; };
    sys.fz = [](double, double, double v) { return v; };
    return sys;
}

// Структура одной итерации стрельбы
struct ShotIter {
    double nu;   // параметр (начальное значение v(0))
    double y1;   // y(1) полученного решения задачи Коши
    double z1;   // y'(1) = y(1) + v(1)
    double phi;  // невязка Phi(nu) = y'(1) - 2 y(1)
};

// Решение задачи Коши с параметром nu и невязка
//   Phi(nu) = y'(1) - 2 y(1) = v(1) - y(1)
double shootResidual(const OdeSystem& sys, double nu, double h, int steps,
                     double* y1, double* z1) {
    const OdeSolution s = solveRK4(sys, kA, 1.0, nu, h, steps);
    *y1 = s.y.back();
    *z1 = s.y.back() + kB * s.z.back();  // y'(1) = y + x v при x = 1
    return s.z.back() - s.y.back();      // v(1) - y(1)
}

// Метод стрельбы поиск параметра nu методом секущих
std::vector<ShotIter> shooting(const OdeSystem& sys, double h, int steps,
                               double* nuStar) {
    const double kEps = 1e-10;
    const int kMaxIt = 30;
    std::vector<ShotIter> iters;

    // Начальные значения параметра
    double nu0 = 0.0, nu1 = 2.5;
    double y0 = 0.0, z0 = 0.0, y1 = 0.0, z1 = 0.0;
    double P0 = shootResidual(sys, nu0, h, steps, &y0, &z0);
    double P1 = shootResidual(sys, nu1, h, steps, &y1, &z1);
    iters.push_back({nu0, y0, z0, P0});
    iters.push_back({nu1, y1, z1, P1});
    *nuStar = nu1;

    // Секущие nu_{j+2} = nu_j - Phi_j (nu_{j+1} - nu_j)/(Phi_{j+1} - Phi_j)
    for (int j = 2; j <= kMaxIt && std::fabs(P1) >= kEps; ++j) {
        const double nu2 = nu0 - P0 * (nu1 - nu0) / (P1 - P0);
        double y2 = 0.0, z2 = 0.0;
        const double P2 = shootResidual(sys, nu2, h, steps, &y2, &z2);
        iters.push_back({nu2, y2, z2, P2});
        nu0 = nu1; P0 = P1;
        nu1 = nu2; P1 = P2;
        *nuStar = nu2;
    }
    return iters;
}

void printShooting(const OdeSystem& sys, double h, int steps, double* nuStar) {
    printSubHeader("Метод стрельбы: итерации уточнения параметра");
    const std::vector<ShotIter> iters = shooting(sys, h, steps, nuStar);
    const std::vector<std::pair<std::string, int>> cols = {
        {"j", 2}, {"nu_j", 15}, {"y(1)", 15}, {"y'(1)", 15}, {"Φ(ν_j)", 12}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (std::size_t j = 0; j < iters.size(); ++j) {
        std::ostringstream row;
        row << " " << std::setw(2) << j << "| " << fCell(iters[j].nu, 15, 10)
            << "| " << fCell(iters[j].y1, 15, 10) << "| "
            << fCell(iters[j].z1, 15, 10) << "| " << eCell(iters[j].phi, 12);
        std::cout << row.str() << "\n";
    }
    std::cout << "Итерации прекращены при |Φ(ν)| < 1e-10; найденный параметр\n"
              << "ν* = " << std::fixed << std::setprecision(10) << *nuStar
              << " (точное значение ν = 2: ему отвечает C2 = 0,\n"
                 "т.е. решение совпадает с e^x (x^2+1)).\n";
}

void printShootingSolution(const OdeSystem& sys, double nuStar, double h,
                           int steps) {
    printSubHeader("Решение краевой задачи методом стрельбы (ν = ν*)");
    const OdeSolution s = solveRK4(sys, kA, 1.0, nuStar, h, steps);
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5}, {"y_k", 14}, {"z_k=y'_k", 14}, {"y_ист", 14},
        {"ε_y", 11}, {"ε_z", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    double maxE = 0.0;
    for (std::size_t k = 0; k < s.x.size(); ++k) {
        const double z = s.y[k] + s.x[k] * s.z[k];  // y' = y + x v
        const double ey = s.y[k] - yExact(s.x[k]);
        const double ez = z - zExact(s.x[k]);
        maxE = std::max(maxE, std::fabs(ey));
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(s.x[k]) << "| "
            << fCell(s.y[k]) << "| " << fCell(z) << "| "
            << fCell(yExact(s.x[k])) << "| " << eCell(ey) << "| " << eCell(ez);
        std::cout << row.str() << "\n";
    }
    std::cout << "max |ε_y| = " << std::scientific << std::setprecision(3)
              << maxE << std::fixed << "\n";
}

// Коэффициенты уравнения, приведенного к виду y'' + p(x) y' + q(x) y = 0:
//   p(x) = -(2x+1)/x,   q(x) = (x+1)/x   (при x > 0)
double pCoef(double x) { return -(2.0 * x + 1.0) / x; }
double qCoef(double x) { return (x + 1.0) / x; }

// Система линейных уравнений конечно-разностной схемы:
//   внутренние узлы: (1 - h p_k/2) y_{k-1} + (h^2 q_k - 2) y_k +
//                    (1 + h p_k/2) y_{k+1} = 0,   k = 1..N-1;
//   левое ГУ y'(0) = 1 и правое ГУ y'(1) = 2 y(1) аппроксимируются
//   односторонними разностями 2-го порядка:
//     (-3 y_0 + 4 y_1 - y_2)/(2h) = 1,
//     (y_{N-2} - 4 y_{N-1} + 3 y_N)/(2h) = 2 y_N.
struct LinSystem {
    std::vector<std::vector<double>> a;
    std::vector<double> b;
};

LinSystem buildFdSystem(double h, int n) {
    LinSystem s;
    s.a.assign(n + 1, std::vector<double>(n + 1, 0.0));
    s.b.assign(n + 1, 0.0);
    // Левое краевое условие y'(0) = 1, разность 2-го порядка
    s.a[0][0] = -3.0;
    s.a[0][1] = 4.0;
    s.a[0][2] = -1.0;
    s.b[0] = 2.0 * h;
    // Внутренние узлы
    for (int k = 1; k <= n - 1; ++k) {
        const double x = k * h;
        const double p = pCoef(x);
        const double q = qCoef(x);
        s.a[k][k - 1] = 1.0 - h * p / 2.0;
        s.a[k][k] = h * h * q - 2.0;
        s.a[k][k + 1] = 1.0 + h * p / 2.0;
        s.b[k] = 0.0;
    }
    // Правое краевое условие y'(1) - 2 y(1) = 0, разность 2-го порядка
    s.a[n][n - 2] = 1.0;
    s.a[n][n - 1] = -4.0;
    s.a[n][n] = 3.0 - 4.0 * h;
    s.b[n] = 0.0;
    return s;
}

std::vector<double> gaussSolve(std::vector<std::vector<double>> a,
                               std::vector<double> b) {
    const int m = static_cast<int>(a.size());
    for (int col = 0; col < m; ++col) {
        int piv = col;
        for (int r = col + 1; r < m; ++r)
            if (std::fabs(a[r][col]) > std::fabs(a[piv][col])) piv = r;
        std::swap(a[col], a[piv]);
        std::swap(b[col], b[piv]);
        for (int r = col + 1; r < m; ++r) {
            const double f = a[r][col] / a[col][col];
            if (f != 0.0) {
                for (int c = col; c < m; ++c) a[r][c] -= f * a[col][c];
                b[r] -= f * b[col];
            }
        }
    }
    std::vector<double> y(m, 0.0);
    for (int r = m - 1; r >= 0; --r) {
        double t = b[r];
        for (int c = r + 1; c < m; ++c) t -= a[r][c] * y[c];
        y[r] = t / a[r][r];
    }
    return y;
}

void printCoeffTable(double h, int n) {
    printSubHeader("Коэффициенты уравнения во внутренних узлах");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 3}, {"x_k", 6}, {"p(x_k)", 12}, {"q(x_k)", 12}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (int k = 1; k <= n - 1; ++k) {
        const double x = k * h;
        std::ostringstream row;
        row << " " << std::setw(3) << k << "| " << xCell(x, 6) << "| "
            << fCell(pCoef(x), 12, 6) << "| " << fCell(qCoef(x), 12, 6);
        std::cout << row.str() << "\n";
    }
}

void printFdSystem(const LinSystem& s, int n) {
    printSubHeader("Система линейных уравнений (N = " + std::to_string(n) +
                   ")");
    std::vector<std::pair<std::string, int>> cols = {{"k", 3}};
    for (int c = 0; c <= n; ++c)
        cols.push_back({"y" + std::to_string(c), 9});
    cols.push_back({"RHS", 10});
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    for (int r = 0; r <= n; ++r) {
        std::ostringstream row;
        row << " " << std::setw(3) << r << "| ";
        for (int c = 0; c <= n; ++c)
            row << fCell(s.a[r][c], 9, 4) << "| ";
        row << fCell(s.b[r], 10, 4);
        std::cout << row.str() << "\n";
    }
}

void printFdSolution(const std::vector<double>& xs,
                     const std::vector<double>& ys) {
    printSubHeader("Решение конечно-разностным методом");
    const std::vector<std::pair<std::string, int>> cols = {
        {"k", 2}, {"x_k", 5}, {"y_k", 14}, {"y_ист", 14}, {"ε_y", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    double maxE = 0.0;
    for (std::size_t k = 0; k < xs.size(); ++k) {
        const double e = ys[k] - yExact(xs[k]);
        maxE = std::max(maxE, std::fabs(e));
        std::ostringstream row;
        row << " " << std::setw(2) << k << "| " << xCell(xs[k]) << "| "
            << fCell(ys[k]) << "| " << fCell(yExact(xs[k])) << "| "
            << eCell(e);
        std::cout << row.str() << "\n";
    }
    std::cout << "max |ε_y| = " << std::scientific << std::setprecision(3)
              << maxE << std::fixed << "\n";
}

// Таблица Рунге-Ромберга для наборов значений на сетках h и 2h
void printRRValues(const std::string& name,
                   const std::vector<double>& xf, const std::vector<double>& yf,
                   const std::vector<double>& yc, int p) {
    const int n2 = static_cast<int>(yc.size()) - 1;
    const double denom = std::pow(2.0, p) - 1.0;
    printSubHeader("Рунге-Ромберг: " + name + ", p = " + std::to_string(p) +
                   " (сравнение шагов h и 2h)");
    const std::vector<std::pair<std::string, int>> cols = {
        {"x_k", 5}, {"y_h", 14}, {"y_2h", 14}, {"R_h", 11}, {"y_уточн", 14},
        {"|y_h-y_ист|", 11}};
    std::cout << tableHeader(cols) << "\n" << tableSep(cols) << "\n";
    double maxR = 0.0, maxT = 0.0;
    for (int i = 0; i <= n2; ++i) {
        const int j = 2 * i;
        const double R = (yf[j] - yc[i]) / denom; 
        const double ref = yf[j] + R;
        const double terr = std::fabs(yf[j] - yExact(xf[j]));
        maxR = std::max(maxR, std::fabs(R));
        maxT = std::max(maxT, terr);
        std::ostringstream row;
        row << " " << xCell(xf[j]) << "| " << fCell(yf[j]) << "| "
            << fCell(yc[i]) << "| " << eCell(R) << "| " << fCell(ref) << "| "
            << eCell(terr);
        std::cout << row.str() << "\n";
    }
    std::cout << "max|R_h| = " << std::scientific << std::setprecision(3)
              << maxR << ",  max|y_h - y_ист| = " << maxT << std::fixed << "\n";
}

}  // namespace

void task42(double hIn) {
    printHeader("Задание 4.2. Краевая задача для ОДУ 2-го порядка");
    int n = static_cast<int>(std::lround((kB - kA) / hIn));
    if (n < 2) n = 2;
    const double h = (kB - kA) / n;
    if (std::fabs(h - hIn) > 1e-12)
        std::cout << "(шаг скорректирован так, чтобы (b-a)/h было целым: h = "
                  << h << ")\n";
    std::cout << "Сетка: x_k = k*h, k = 0.." << n << ", h = " << h << "\n";

    const OdeSystem reg = makeRegularizedSystem();

    double nuStar = 0.0;
    printShooting(reg, h, n, &nuStar);
    printShootingSolution(reg, nuStar, h, n);

    printCoeffTable(h, n);
    const LinSystem sysH = buildFdSystem(h, n);
    printFdSystem(sysH, n);
    const std::vector<double> yH = gaussSolve(sysH.a, sysH.b);
    std::vector<double> xH(n + 1);
    for (int k = 0; k <= n; ++k) xH[k] = k * h;
    printFdSolution(xH, yH);

    if (n % 2 == 0) {
        const int n2 = n / 2;
        const double h2 = 2.0 * h;
        // Стрельба на сетке h (то же решение, что в таблице выше)
        const OdeSolution sFine = solveRK4(reg, kA, 1.0, nuStar, h, n);
        std::vector<double> yShootH(n + 1);
        for (int k = 0; k <= n; ++k) yShootH[k] = sFine.y[k];

        // Стрельба на грубой сетке (со своим параметром ν)
        double nu2 = 0.0;
        shooting(reg, h2, n2, &nu2);
        const OdeSolution s2 = solveRK4(reg, kA, 1.0, nu2, h2, n2);
        std::vector<double> yShoot2(n2 + 1), x2(n2 + 1);
        for (int k = 0; k <= n2; ++k) {
            x2[k] = k * h2;
            yShoot2[k] = s2.y[k];
        }
        printRRValues("метод стрельбы (РК4)", xH, yShootH, yShoot2, 4);

        // КР метод на грубой сетке
        const LinSystem sys2 = buildFdSystem(h2, n2);
        const std::vector<double> yFd2 = gaussSolve(sys2.a, sys2.b);
        printRRValues("конечно-разностный метод", xH, yH, yFd2, 2);
    } else {
        std::cout << "\nДля оценки по Рунге-Ромбергу нужно четное число "
                     "шагов; таблицы Рунге-Ромберга пропущены.\n";
    }

    std::ofstream out("bvp_plot_data.txt");
    out << "# x y_exact y_shooting y_fd\n";
    const OdeSolution sFin = solveRK4(reg, kA, 1.0, nuStar, h, n);
    for (int k = 0; k <= n; ++k) {
        out << std::fixed << std::setprecision(6) << xH[k] << " "
            << yExact(xH[k]) << " " << sFin.y[k] << " " << yH[k] << "\n";
    }
}
