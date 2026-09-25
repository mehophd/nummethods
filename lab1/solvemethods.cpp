#include <cmath>
#include <algorithm>
#include <iostream>
#include "solvemethods.h"
#include "matrixtools.h"

using namespace std;

float* methodLU(SLAE &problem) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    int maxRow;
    for (int start = 0; start < size; ++start) {
        maxRow = start;
        for (int i = start + 1; i < size; ++i) {
            if (abs(matrix[i][start]) > abs(matrix[maxRow][start])) {
                maxRow = i;
            }
        }
        
        if (maxRow != start) {
            swap(constants[start], constants[maxRow]);
            swap(matrix[start], matrix[maxRow]);
        }

        if (abs(matrix[start][start]) < 1e-7f) {
            freeMatrix(matrix, size);
            delete[] constants;
            return nullptr;
        }

        for (int i = start + 1; i < size; ++i) {
            matrix[i][start] = matrix[i][start] / matrix[start][start];
        }

        for (int i = start + 1; i < size; ++i) {
            for (int j = start + 1; j < size; ++j) {
                matrix[i][j] = matrix[i][j] - matrix[start][j] * matrix[i][start];
            }
        }
    }

    // Ly=Pb Прямой ход
    // y[0] = b'[0]
    // L10*y[0] + y[1] = b'[1] 
    // L20*y[0] + L21*y[1] + y[2] = b'[2]
    // L30*y[0] + L31*y[1] + L32*y[2] + y[3] = b'[2]
    float* ytab = new float[size];
    ytab[0] = constants[0];
    for (int i = 1; i < size; ++i) {
        float scalarProduct = 0.0f;
        for (int j = 0; j < i; ++j)
            scalarProduct += ytab[j]*matrix[i][j];
        ytab[i] = constants[i] - scalarProduct;
    }

    // Ux=y Обратный ход
    // U00*x[0] + U01*x[1] + U02*x[2] + U03*x[3] = y[0]
    // U11*x[1] + U12*x[2] + U13*x[3] = y[1]
    //    U22*x[2] + U23*x[3] = y[2]
    //             U33*x[3] = y[3]
    float* solve = new float[size];
    solve[size-1] = ytab[size-1] / matrix[size-1][size-1];
    for (int i = size - 2; i >= 0; --i) {
        float scalarProduct = 0.0f;
        for (int j = size - 1; j > i; --j)
            scalarProduct += solve[j]*matrix[i][j];
        solve[i] = (ytab[i] - scalarProduct) / matrix[i][i];
    }

    delete[] ytab;
    freeMatrix(matrix, size);
    delete[] constants;

    return solve;
}

float* methodTDMA(SLAE &problem) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < i - 1; ++j)
            if (matrix[i][j] > 1e-7f) {
                freeMatrix(matrix, size);
                delete[] constants;

                return nullptr;
            }
        
        for (int j = i + 2; j < size; ++j)
            if (matrix[i][j] > 1e-7f) {
                freeMatrix(matrix, size);
                delete[] constants;

                return nullptr;
            }
    }

    float* P = new float[size];
    float* Q = new float[size];
    float* solve = new float[size];

    for (int i = 0; i < size; ++i) {
        float a = (i > 0) ? matrix[i][i - 1] : 0.0f;
        float b = matrix[i][i];
        float c = (i < size - 1) ? matrix[i][i + 1] : 0.0f;
        float d = constants[i];

        float denominator = b + a * ((i > 0) ? P[i - 1] : 0.0f);

        if (abs(denominator) < 1e-7f) {
            delete[] P;
            delete[] Q;
            delete[] solve;
            freeMatrix(matrix, size);
            delete[] constants;
            return nullptr;
        }

        P[i] = -c / denominator;
        Q[i] = (d - a * ((i > 0) ? Q[i - 1] : 0.0f)) / denominator;
    }

    solve[size - 1] = Q[size - 1];

    for (int i = size - 2; i >= 0; --i) {
        solve[i] = P[i] * solve[i + 1] + Q[i];
    }

    delete[] P;
    delete[] Q;
    freeMatrix(matrix, size);
    delete[] constants;

    return solve;
}

float* methodSI(SLAE &problem) {
}