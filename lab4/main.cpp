#include <cmath>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

void task41(double h);
void task42(double h);

namespace {

// Запрос шага сетки h - входного данных задания. Пустой ввод -> h = 0.1
double readStep() {
    std::cout << "Введите шаг сетки h (например, 0.1; Enter - 0.1): ";
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return 0.1;
    const double h = std::atof(line.c_str());
    if (h <= 0.0 || !std::isfinite(h)) {
        std::cout << "Некорректный шаг, используется h = 0.1.\n";
        return 0.1;
    }
    return h;
}

}  // namespace

int main() {
    std::cout << "############################################################\n"
              << "#   Лабораторная работа: численные методы решения          #\n"
              << "#   начальных и краевых задач для ОДУ (задания 4.1, 4.2)   #\n"
              << "############################################################\n";

    const double h = readStep();

    task41(h);
    task42(h);

    return 0;
}
