#include <iostream>
#include "solvemethods.h"
#include "matrixtools.h"
#include <cmath>
#include "Slae.h"

using namespace std;

void printSolution(const char* methodName, float* solution, int size, int iters = -1) {
    cout << "\n" << methodName << endl;
    if (!solution) {
        cout << "Метод не сошелся или неприменим." << endl;
        return;
    }
    for (int i = 0; i < size; ++i) {
        cout << "x" << i << ": " << solution[i] << endl;
    }
    if (iters >= 0) {
        cout << "Итераций: " << iters << endl;
    }
}

bool readMatrix(SLAE& problem, int size) {
    float value;
    cout << "Введите коэффициенты: " << endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << "c" << i << j << ":";
            if (!(cin >> value)) {
                cout << "Ошибка: коэффициенты матрицы должны быть заданы числом." << endl;
                return false;
            }
            problem.setCoef(i, j, value);
        }
    }
    return true;
}

bool readConstants(SLAE& problem, int size) {
    float value;
    cout << "Введите свободные члены: " << endl;
    for (int i = 0; i < size; ++i) {
        cout << "b" << i << ":";
        if (!(cin >> value)) {
            cout << "Ошибка: свободные члены должны быть заданы числом." << endl;
            return false;
        }
        problem.setConst(i, value);
    }
    return true;
}

void showMenu() {
    cout << "\n=== Выберите метод решения ===" << endl;
    cout << "1. Метод Гаусса (LU-разложение)" << endl;
    cout << "2. Метод прогонки" << endl;
    cout << "3. Метод простых итераций" << endl;
    cout << "4. Метод Зейделя" << endl;
    cout << "5. Все методы сразу" << endl;
    cout << "0. Выход" << endl;
    cout << "Ваш выбор: ";
}

void solveLU(SLAE& problem, int size) {
    float* x = methodLU(problem);
    printSolution("Метод Гаусса (LU-разложение)", x, size);
    delete[] x;
}

void solveTDMA(SLAE& problem, int size) {
    float* x = methodTDMA(problem);
    if (!x) {
        cout << "\nМетод прогонки" << endl;
        cout << "Не выполнено диагональное преобладание или матрица не трехдиагональная." << endl;
    } else {
        printSolution("Метод прогонки", x, size);
        delete[] x;
    }
}

void solveSI(SLAE& problem, int size, float eps, int limit) {
    int iters = 0;
    float* x = methodSI(problem, eps, limit, iters);
    printSolution("Метод простых итераций", x, size, iters);
    delete[] x;
}

void solveZ(SLAE& problem, int size, float eps, int limit) {
    int iters = 0;
    float* x = methodZ(problem, eps, limit, iters);
    printSolution("Метод Зейделя", x, size, iters);
    delete[] x;
}

int main() {
    int size, limit;
    float eps;
    
    cout << "Введите порядок системы: ";
    if (!(cin >> size) || size <= 0) {
        cout << "Ошибка: порядок системы должен быть положительным числом." << endl;
        return 1;
    }

    SLAE problem(size);
    
    if (!readMatrix(problem, size) || !readConstants(problem, size)) {
        return 1;
    }

    cout << "Введите эпсилон: ";
    if (!(cin >> eps)) {
        cout << "Ошибка: эпсилон должен быть задан числом." << endl;
        return 1;
    }

    cout << "Введите предел итераций: ";
    if (!(cin >> limit)) {
        cout << "Ошибка: предел итераций должен быть задан числом." << endl;
        return 1;
    }

    problem.printSLAE();
    
    float det = determinant(problem.getMatrixCopy(), size);
    cout << "\nДетерминант: " << det << endl;
    
    if (fabs(det) < 1e-10) {
        cout << "Матрица вырожденная. Решение невозможно." << endl;
        return 1;
    }
    
    float** inverse = inverseMatrix(problem.getMatrixCopy(), size);
    cout << "\nОбратная матрица: " << endl;
    printMatrix(inverse, size);

    int choice;
    while (true) {
        showMenu();
        if (!(cin >> choice)) {
            cout << "Ошибка ввода. Попробуйте снова." << endl;
            cin.clear();
            cin.ignore(10000, '\n');
            continue;
        }

        if (choice == 0) {
            cout << "Выход из программы." << endl;
            break;
        }

        switch (choice) {
            case 1:
                solveLU(problem, size);
                break;
            case 2:
                solveTDMA(problem, size);
                break;
            case 3:
                solveSI(problem, size, eps, limit);
                break;
            case 4:
                solveZ(problem, size, eps, limit);
                break;
            case 5:
                solveLU(problem, size);
                solveTDMA(problem, size);
                solveSI(problem, size, eps, limit);
                solveZ(problem, size, eps, limit);
                break;
            default:
                cout << "Неверный выбор. Попробуйте снова." << endl;
                break;
        }
    }

    return 0;
}