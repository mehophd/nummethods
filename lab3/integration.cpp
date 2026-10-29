#include "integration.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"

namespace {

const double kA = 0.0;   // нижний предел интегрирования
const double kB = 2.0;   // верхний предел интегрирования
const double kH1 = 0.5;  // шаг первой сетки
const double kH2 = 0.25; // шаг второй сетки

double func(double x) { return x / (x + 9.0); }

std::string padRight(const std::string& s, int width) {
    int visible = 0;
    for (unsigned char ch : s)
        if ((ch & 0xC0) != 0x80) ++visible;
    return s + std::string(width > visible ? width - visible : 0, ' ');
}

struct StepTables {
    std::vector<double> x, y;             // узлы и значения функции
    std::vector<double> rect, trap, simp; // накопленные интегралы
    double fRect = 0.0, fTrap = 0.0, fSimp = 0.0;  // итоговые значения F(h)
};

StepTables integrate(double h) {
    StepTables t;
    const int n = static_cast<int>(std::lround((kB - kA) / h));
    t.x.resize(n + 1);
    t.y.resize(n + 1);
    t.rect.assign(n + 1, 0.0);
    t.trap.assign(n + 1, 0.0);
    t.simp.assign(n + 1, 0.0);
    for (int i = 0; i <= n; ++i) {
        t.x[i] = kA + i * h;
        t.y[i] = func(t.x[i]);
        if (i == 0) continue;
        // Метод прямоугольников
        t.rect[i] = t.rect[i - 1] + h * func((t.x[i - 1] + t.x[i]) / 2.0);
        // Метод трапеций
        t.trap[i] = t.trap[i - 1] + h * (t.y[i - 1] + t.y[i]) / 2.0;
        // Метод Симпсона
        if (i >= 2)
            t.simp[i] =
                t.simp[i - 2] + h * (t.y[i - 2] + 4.0 * t.y[i - 1] + t.y[i]) / 3.0;
    }
    t.fRect = t.rect[n];
    t.fTrap = t.trap[n];
    t.fSimp = t.simp[n];
    return t;
}

void printStepTables(const StepTables& t, double h) {
    const int n = static_cast<int>(t.x.size()) - 1;
    std::cout << "\nТаблица значений и накопленных интегралов (шаг h = "
              << std::fixed << std::setprecision(2) << h << "):\n";
    std::cout << " i |  Xi   | Yi=f(Xi) | прямоугольников | трапеций | Симпсона\n";
    std::cout << "---+-------+----------+-----------------+----------+----------\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << " " << i << " | " << std::setw(5) << std::setprecision(2) << t.x[i]
                  << " | " << std::setw(8) << std::setprecision(5) << t.y[i] << " | "
                  << std::setw(15) << t.rect[i] << " | " << std::setw(8) << t.trap[i]
                  << " | ";
        if (i % 2 != 0)
            std::cout << "         -";
        else
            std::cout << std::setw(8) << t.simp[i];
        std::cout << "\n";
    }
}

}  // namespace

void task35() {
    printHeader("Задание 3.5. Численное интегрирование, метод Рунге-Ромберга");
    std::cout << "Вычислить F = INT[0;2] ( x/(x+9) ) dx методами прямоугольников,\n"
                 "трапеций и Симпсона с шагами h1 = 0.5 и h2 = 0.25;\n"
                 "оценить погрешность и уточнить результат методом Рунге-Ромберга.\n";

    for (double h : {kH1, kH2}) {
        const double intervals = (kB - kA) / h;
        if (std::fabs(intervals - std::lround(intervals)) > 1e-9 ||
            static_cast<int>(std::lround(intervals)) % 2 != 0) {
            std::cout << "Ошибка: шаг h = " << h
                      << " дает нечетное число отрезков для метода Симпсона.\n";
            return;
        }
    }

    const StepTables t1 = integrate(kH1);
    const StepTables t2 = integrate(kH2);
    printStepTables(t1, kH1);
    printStepTables(t2, kH2);

    // F = F(h2) + (F(h2) - F(h1)) / (k^p - 1),  k = h1/h2 = 2,
    // p = 2 для прямоугольников и трапеций, p = 4 для Симпсона
    const double k = kH1 / kH2;

    struct RResult {
        double est = 0.0;      // оценка погрешности F(h2)
        double refined = 0.0;  // уточненное значение
    };
    auto refine = [k](double finer, double coarser, double p) {
        RResult r;
        r.est = (finer - coarser) / (std::pow(k, p) - 1.0);
        r.refined = finer + r.est;
        return r;
    };
    const RResult rRect = refine(t2.fRect, t1.fRect, 2.0);
    const RResult rTrap = refine(t2.fTrap, t1.fTrap, 2.0);
    const RResult rSimp = refine(t2.fSimp, t1.fSimp, 4.0);

    // F = [x - 9*ln(x+9)] на [0; 2] = 2 - 9*ln(11/9).
    const double exact = (kB - kA) - 9.0 * std::log((kB + 9.0) / (kA + 9.0));

    std::cout << "\nМетод Рунге-Ромберга: F = F(h2) + (F(h2)-F(h1))/(k^p - 1),\n"
              << "k = h1/h2 = " << k
              << "; p = 2 (прямоугольники, трапеции), p = 4 (Симпсон).\n\n";
    std::cout << " метод       |  F(h1)   |  F(h2)   | оценка погр. F(h2) | уточненное значение\n";
    std::cout << "-------------+----------+----------+--------------------+--------------------\n";
    auto printRow = [](const char* name, double f1, double f2, double est,
                       double ref) {
        std::cout << " " << padRight(name, 12) << "| " << std::right
                  << std::fixed << std::setprecision(5) << std::setw(8) << f1 << " | "
                  << std::setw(8) << f2 << " | " << std::scientific << std::setprecision(2)
                  << std::setw(18) << est << " | " << std::fixed << std::setprecision(6)
                  << std::setw(18) << ref << "\n";
    };
    printRow("прямоугольников", t1.fRect, t2.fRect, rRect.est, rRect.refined);
    printRow("трапеций", t1.fTrap, t2.fTrap, rTrap.est, rTrap.refined);
    printRow("Симпсона", t1.fSimp, t2.fSimp, rSimp.est, rSimp.refined);

    std::cout << "\nТочное значение: F = [x - 9*ln(x+9)] на [0; 2] = 2 - 9*ln(11/9) = "
              << std::setprecision(6) << exact << "\n";
    std::cout << "Фактическая погрешность уточненных значений:\n"
              << "  прямоугольники: " << std::scientific << std::setprecision(2)
              << std::fabs(rRect.refined - exact) << "\n"
              << "  трапеции:       " << std::fabs(rTrap.refined - exact) << "\n"
              << "  Симпсон:        " << std::fabs(rSimp.refined - exact) << "\n"
              << std::fixed;
}
