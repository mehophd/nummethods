#include <iostream>
#include "solvemethods.h"
#include "Slae.h"

using namespace std;

int main() {
    int size;

    cout << "Введите порядок системы: ";
    if (!(cin >> size)) {
        cout << "Ошибка: размер матрицы должен быть задан числом.";
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
                cout << "Ошибка: коэффициенты матрицы должны быть заданы числом.";
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
            cout << "Ошибка: свободные члены должны быть заданы числом.";
            return 1;
        }
    }
    methodLA(problem);

    problem.printMatrix();
    return 0;
}