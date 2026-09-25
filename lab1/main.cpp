#include <iostream>
#include "solvemethods.h"
#include "matrixtools.h"
#include "Slae.h"

using namespace std;

int main() {
    int size;

    cout << "Введите порядок системы: ";
    if (!(cin >> size)) {
        cout << "Ошибка: порядок системы должен быть задан числом.";
        return 1;
    }

    if (size < 0) {
        cout << "Ошибка: порядок матрицы не может быть меньше 1." << endl;
        return 1;
    }

    SLAE problem(size);
    
    size = problem.getSize();
    cout << "Введите коэффициенты: " << endl;
    float value;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << "c" << i << j << ":";
            if (cin>>value) {
                problem.setCoef(i, j, value);
            } else {
                cout << "Ошибка: коэффициенты матрицы должны быть заданы числом." << endl;
                return 1;
            }
        }
    }

    cout << "Введите свободные члены: " << endl;
    for (int i = 0; i < size; ++i) {
        cout << "b" << i << ":";
        if (cin>>value) {
            problem.setConst(i, value);
        } else {
            cout << "Ошибка: свободные члены должны быть заданы числом." << endl;
            return 1;
        }
    }

    problem.printSLAE();
    
    cout << "\nОбратная матрица: " << endl;
    float** inverse = inverseMatrix(problem.getMatrixCopy(), size);
    printMatrix(inverse, size);

    cout << "\nДетерминант: " << determinant(problem.getMatrixCopy(), size) << endl;

    cout << "\nМетод Гаусса LA: " << endl;
    float* x = methodLU(problem);
    for (int i = 0; i < size; ++i) {
        cout << "x" << i << ": " << x[i] << endl;
    }

    cout << "Метод прогонки" << endl;
    x = methodTDMA(problem);
    if (!x) {
        cout << "Не выполнено диагональное преобладание." << endl;
    }
    for (int i = 0; i < size; ++i) {
        cout << "x" << i << ": " << x[i] << endl;
    }

    return 0;
}