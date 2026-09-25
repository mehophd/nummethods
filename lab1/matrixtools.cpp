#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

namespace {
    float** cloneMatrix(float** original, int size) {
        float** copy = new float*[size];
        for (int i = 0; i < size; ++i) {
            copy[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                copy[i][j] = original[i][j];
            }
        }
        return copy;
    }
}

void freeMatrix(float** matrix, int size) {
    if (matrix == nullptr) return;
    for (int i = 0; i < size; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

float determinant(float** original, int size) {
    float** matrix = cloneMatrix(original, size);
    float result = 1.0f;
    int sign = 1;

    for (int start = 0; start < size; ++start) {
        int maxRow = start;
        float maxVal = abs(matrix[start][start]);
        for (int i = start + 1; i < size; ++i) {
            if (abs(matrix[i][start]) > maxVal) {
                maxVal = abs(matrix[i][start]);
                maxRow = i;
            }
        }

        if (maxVal < 1e-7f) {
            freeMatrix(matrix, size);
            return 0.0f; 
        }

        if (maxRow != start) {
            swap(matrix[start], matrix[maxRow]);
            sign = -sign; 
        }

        float pivot = matrix[start][start];
        for (int i = start + 1; i < size; ++i) {
            float factor = matrix[i][start] / pivot;
            for (int j = start; j < size; ++j) {
                matrix[i][j] -= matrix[start][j] * factor;
            }
        }

        result *= matrix[start][start];
    }

    freeMatrix(matrix, size);
    return result * sign;
}

float** inverseMatrix(float** original, int size) {
    float** matrix = cloneMatrix(original, size);

    float** inverse = new float*[size];
    for (int i = 0; i < size; ++i) {
        inverse[i] = new float[size];
        for (int j = 0; j < size; ++j) {
            inverse[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }

    for (int start = 0; start < size; ++start) {
        int maxRow = start;
        float maxVal = abs(matrix[start][start]);
        for (int i = start + 1; i < size; ++i) {
            if (abs(matrix[i][start]) > maxVal) {
                maxVal = abs(matrix[i][start]);
                maxRow = i;
            }
        }

        if (maxVal < 1e-7f) {
            freeMatrix(matrix, size);
            freeMatrix(inverse, size);
            return nullptr;
        }

        if (maxRow != start) {
            swap(matrix[start], matrix[maxRow]);
            swap(inverse[start], inverse[maxRow]);
        }
   
        float pivot = matrix[start][start];
        for (int j = 0; j < size; ++j) {
            matrix[start][j]  /= pivot;
            inverse[start][j] /= pivot;
        }

        for (int i = 0; i < size; ++i) {
            if (i != start) {
                float factor = matrix[i][start];
                for (int j = 0; j < size; ++j) {
                    matrix[i][j]  -= factor * matrix[start][j];
                    inverse[i][j] -= factor * inverse[start][j];
                }
            }
        }
    }

    freeMatrix(matrix, size);

    return inverse;
}

void printMatrix(float** original, int size) {
    float** matrix = cloneMatrix(original, size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << "\n";
    }

    freeMatrix(matrix, size);
}