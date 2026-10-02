#include <iostream>
#include "solvemethods.h"
#include "matrixtools.h"
#include <cmath>
#include "Slae.h"
#include "eigen.h"

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

void clearInput() {
    cin.clear();
    cin.ignore(10000, '\n');
}

int readSize() {
    int size;
    while (true) {
        cout << "Введите порядок системы: ";
        if (cin >> size && size > 0) {
            return size;
        }
        cout << "Ошибка: порядок системы должен быть положительным числом." << endl;
        clearInput();
    }
}

float readEps() {
    float eps;
    while (true) {
        cout << "Введите эпсилон: ";
        if (cin >> eps && eps > 0) {
            return eps;
        }
        cout << "Ошибка: эпсилон должен быть положительным числом." << endl;
        clearInput();
    }
}

int readLimit() {
    int limit;
    while (true) {
        cout << "Введите предел итераций: ";
        if (cin >> limit && limit > 0) {
            return limit;
        }
        cout << "Ошибка: предел итераций должен быть положительным числом." << endl;
        clearInput();
    }
}

bool readMatrix(SLAE& problem, int size) {
    float value;
    cout << "Введите коэффициенты: " << endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            while (true) {
                cout << "c" << i << j << ":";
                if (cin >> value) {
                    problem.setCoef(i, j, value);
                    break;
                }
                cout << "Ошибка: коэффициент должен быть числом." << endl;
                clearInput();
            }
        }
    }
    return true;
}

bool readConstants(SLAE& problem, int size) {
    float value;
    cout << "Введите свободные члены: " << endl;
    for (int i = 0; i < size; ++i) {
        while (true) {
            cout << "b" << i << ":";
            if (cin >> value) {
                problem.setConst(i, value);
                break;
            }
            cout << "Ошибка: свободный член должен быть числом." << endl;
            clearInput();
        }
    }
    return true;
}

void showMenu() {
    cout << "\n.----Выберите метод решения----" << endl;
    cout << "1) Метод Гаусса (LU-разложение)" << endl;
    cout << "2) Метод прогонки" << endl;
    cout << "3) Метод простых итераций" << endl;
    cout << "4) Метод Зейделя" << endl;
    cout << "5) Метод вращений Якоби" << endl;
    cout << "6) Метод QR (собственные значения)" << endl;
    cout << "7) Все методы сразу" << endl;
    cout << "0) Выход" << endl;
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
    float* x = methodSI(problem, eps, limit, &iters);
    printSolution("Метод простых итераций", x, size, iters);
    delete[] x;
}

void solveZ(SLAE& problem, int size, float eps, int limit) {
    int iters = 0;
    float* x = methodZ(problem, eps, limit, &iters);
    printSolution("Метод Зейделя", x, size, iters);
    delete[] x;
}

void solveJacobi(SLAE& problem, int size, float eps, int limit) {
    float** matrixCopy = problem.getMatrixCopy();
    float** eigenvectors = nullptr;
    int iters = 0;
    
    float* eigenvalues = methodJacobi(matrixCopy, eps, size, limit, &eigenvectors, &iters);
    
    cout << "\nМетод вращений Якоби" << endl;
    if (!eigenvalues || !eigenvectors) {
        cout << "Метод не сошелся или матрица не симметрична." << endl;
    } else {
        cout << "Собственные значения:" << endl;
        for (int i = 0; i < size; ++i) {
            cout << "λ" << i << ": " << eigenvalues[i] << endl;
        }

        cout << "\nСобственные векторы (по столбцам):" << endl;
        for (int j = 0; j < size; ++j) {
            cout << "x" << j << ": [";
            for (int i = 0; i < size; ++i) {
                cout << eigenvectors[i][j];
                if (i < size - 1) cout << ", ";
            }
            cout << "]" << endl;
        }
        cout << "Итераций: " << iters << endl;
    }

    delete[] eigenvalues;
    if (eigenvectors) freeMatrix(eigenvectors, size);
    freeMatrix(matrixCopy, size);
}

void solveQR(SLAE& problem, int size, float eps, int limit) {
    float** matrixCopy = problem.getMatrixCopy();
    int iters = 0;
    float* imaginaryParts = new float[size];
    
    cout << "\nQR-алгоритм (собственные значения)" << endl;
    
    float* eigenvalues = methodQR(matrixCopy, eps, size, limit, &iters, imaginaryParts);
    
    if (!eigenvalues) {
        cout << "Ошибка при вычислении собственных значений." << endl;
    } else {
        cout << "Собственные значения:" << endl;
        for (int i = 0; i < size; ++i) {
            if (imaginaryParts[i] == 0.0f) {
                cout << "λ" << i << ": " << eigenvalues[i] << endl;
            } else {
                cout << "λ" << i << ": " << eigenvalues[i] 
                     << (imaginaryParts[i] > 0 ? "+" : "") 
                     << imaginaryParts[i] << "i" << endl;
            }
        }
        cout << "Итераций: " << iters << endl;
        
        delete[] eigenvalues;
    }
    
    delete[] imaginaryParts;
    freeMatrix(matrixCopy, size);
}

int main() {
    int size = readSize();
    SLAE problem(size);
    
    readMatrix(problem, size);
    readConstants(problem, size);

    float eps = readEps();
    int limit = readLimit();

    problem.printSLAE();
    
    float det = determinant(problem.getMatrixCopy(), size);
    cout << "\nДетерминант: " << det << endl;
    
    if (fabs(det) < 1e-7) {
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
            clearInput();
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
                solveJacobi(problem, size, eps, limit);
                break;
            case 6:
                solveQR(problem, size, eps, limit);
                break;
            case 7:
                solveLU(problem, size);
                solveTDMA(problem, size);
                solveSI(problem, size, eps, limit);
                solveZ(problem, size, eps, limit);
                solveJacobi(problem, size, eps, limit);
                solveQR(problem, size, eps, limit);
                break;
            default:
                cout << "Неверный выбор. Попробуйте снова." << endl;
                break;
        }
    }

    return 0;
}