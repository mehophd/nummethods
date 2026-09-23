#include <cmath>
#include <algorithm>
#include <iostream>
#include "solvemethods.h"

using namespace std;

float* methodLA(SLAE &problem) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    int idxMax;
    for (int start = 0; start < size; ++start) {
        idxMax = start;
        for (int i = start + 1; i < size; ++i) {
            if (abs(matrix[i][start]) > abs(matrix[idxMax][start])) {
                idxMax = i;
            }
        }
        swap(constants[start], constants[idxMax]);
        swap(matrix[start], matrix[idxMax]);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cout << matrix[i][j] << " ";
            }
            cout << "| " << constants[i] << endl;
        }
        cout << "\n";

        for (int i = start + 1; i < size; ++i) {
            matrix[i][start] = matrix[i][start] / matrix[start][start];
        }

        for (int i = start + 1; i < size; ++i) {
            for (int j = start + 1; j < size; ++j) {
                matrix[i][j] = matrix[i][j] - matrix[start][j] * matrix[i][start];
            }
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << "| " << constants[i] << endl;
    }

    cout << "\nПоиск Y" << endl;
    float ytab[size];
    ytab[0] = constants[0];
    for (int i = 1; i < size; ++i) {
        float scalarProduct = 0;
        for (int j = 0; j < i; ++j)
            scalarProduct += ytab[j]*matrix[i][j];
        ytab[i] = constants[i] - scalarProduct;
    }

    cout << "\nВывод Y: [";
    for (int i = 0; i < size; ++i) {
        cout << ytab[i] << (i == size - 1 ? "" : ", ");
    }
    cout << "]" << endl;

    cout << "\nПоиск X" << endl;
    float* solve = new float[size];
    solve[size-1] = ytab[size-1] / matrix[size-1][size-1];
    for (int i = size-2; i >= 0; --i) {
        float scalarProduct = 0;
        for (int j = size-1; j > i; --j)
            scalarProduct += solve[j]*matrix[i][j];
        solve[i] = (ytab[i] - scalarProduct) / matrix[i][i];
    }
    
    cout << "\nВывод X: [";
    for (int i = 0; i < size; ++i) {
        cout << solve[i] << (i == size - 1 ? "" : ", ");
    }
    cout << "]" << endl;

    for (int i = 0; i < size; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] constants;

    return solve;
}