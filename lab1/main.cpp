#include <iostream>
#include "solvemethods.h"
#include "matrixtools.h"
#include <cmath>
#include "Slae.h"

using namespace std;

void printSolution(const char* methodName, float* solution, int size) {
    cout << "\n" << methodName << endl;
    if (!solution) {
        cout << "Метод не сошелся." << endl;
        return;
    }
    for (int i = 0; i < size; ++i) {
        cout << "x" << i << ": " << solution[i] << endl;
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
        cout << "Ошибка: эпсилон должен быть задан числом.\n" << endl;
        return 1;
    }

    cout << "Введите предел итераций: ";
    if (!(cin >> limit)) {
        cout << "Ошибка: предел итераций должен быть задан числом.\n" << endl;
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

    float* x = methodLU(problem);
    printSolution("Метод Гаусса LA: ", x, size);

    x = methodTDMA(problem);
    if (!x) {
        cout << "\nМетод прогонки" << endl;
        cout << "Не выполнено диагональное преобладание." << endl;
    } else {
        printSolution("Метод прогонки", x, size);
    }

    x = methodSI(problem, eps, limit);
    printSolution("Метод простых итераций", x, size);

    return 0;
}