#include <cmath>
#include <algorithm>
#include "solvemethods.h"
#include "matrixtools.h"

using namespace std;

#include <cmath>
#include <algorithm>
#include "solvemethods.h"
#include "matrixtools.h"

using namespace std;

float* methodLU(SLAE &problem, float**& out_inverse, float& out_det) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    out_det = 1.0f;
    int swapCount = 0;

    // Инициализируем обратную матрицу как единичную
    float** inverse = new float*[size];
    for (int i = 0; i < size; ++i) {
        inverse[i] = new float[size]();
        inverse[i][i] = 1.0f;
    }

    for (int start = 0; start < size; ++start) {
        int maxRow = start;
        for (int i = start + 1; i < size; ++i) {
            if (fabs(matrix[i][start]) > fabs(matrix[maxRow][start])) {
                maxRow = i;
            }
        }

        if (maxRow != start) {
            swap(constants[start], constants[maxRow]);
            swap(matrix[start], matrix[maxRow]);
            swap(inverse[start], inverse[maxRow]); // Применяем перестановку и к обратной матрице
            swapCount++;
        }

        if (fabs(matrix[start][start]) < 1e-7f) {
            out_det = 0.0f;
            freeMatrix(matrix, size);
            delete[] constants;
            for (int i = 0; i < size; ++i) delete[] inverse[i];
            delete[] inverse;
            return nullptr;
        }

        out_det *= matrix[start][start];

        for (int i = start + 1; i < size; ++i) {
            float factor = matrix[i][start] / matrix[start][start];
            matrix[i][start] = factor;
            for (int j = start + 1; j < size; ++j) {
                matrix[i][j] -= factor * matrix[start][j];
            }
            for (int j = 0; j < size; ++j) {
                inverse[i][j] -= factor * inverse[start][j];
            }
        }
    }

    if (swapCount % 2 != 0) {
        out_det = -out_det;
    }

    // U * X = inverse, где inverse сейчас содержит Y = L^{-1} * P
    // В результате X станет равной A^{-1}
    for (int j = 0; j < size; ++j) {
        for (int i = size - 1; i >= 0; --i) {
            float sum = 0.0f;
            for (int k = i + 1; k < size; ++k) {
                sum += matrix[i][k] * inverse[k][j];
            }
            inverse[i][j] = (inverse[i][j] - sum) / matrix[i][i];
        }
    }

    out_inverse = inverse;

    // Ly=Pb Прямой ход
    // y[0] = b'[0]
    // L10*y[0] + y[1] = b'[1] 
    // L20*y[0] + L21*y[1] + y[2] = b'[2]
    // L30*y[0] + L31*y[1] + L32*y[2] + y[3] = b'[2]
    float* y = new float[size];
    y[0] = constants[0];
    for (int i = 1; i < size; ++i) {
        float sum = 0.0f;
        for (int j = 0; j < i; ++j) {
            sum += matrix[i][j] * y[j];
        }
        y[i] = constants[i] - sum;
    }

    // Ux=y Обратный ход
    // U00*x[0] + U01*x[1] + U02*x[2] + U03*x[3] = y[0]
    //            U11*x[1] + U12*x[2] + U13*x[3] = y[1]
    //                       U22*x[2] + U23*x[3] = y[2]
    //                                  U33*x[3] = y[3]
    float* x = new float[size];
    x[size - 1] = y[size - 1] / matrix[size - 1][size - 1];
    for (int i = size - 2; i >= 0; --i) {
        float sum = 0.0f;
        for (int j = size - 1; j > i; --j) {
            sum += matrix[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / matrix[i][i];
    }

    delete[] y;
    freeMatrix(matrix, size);
    delete[] constants;

    return x;
}

float* methodTDMA(SLAE &problem) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (abs(i - j) > 1 && fabs(matrix[i][j]) > 1e-7f) {
                freeMatrix(matrix, size);
                delete[] constants;
                return nullptr;
            }
        }
    }

    float* P = new float[size];
    float* Q = new float[size];
    
    for (int i = 0; i < size; ++i) {
        float a = (i > 0) ? matrix[i][i - 1] : 0.0f;
        float b = matrix[i][i];
        float c = (i < size - 1) ? matrix[i][i + 1] : 0.0f;
        float d = constants[i];

        float denom = b + a * ((i > 0) ? P[i - 1] : 0.0f);
        if (fabs(denom) < 1e-7f) {
            delete[] P;
            delete[] Q;
            freeMatrix(matrix, size);
            delete[] constants;
            return nullptr;
        }

        P[i] = -c / denom;
        Q[i] = (d - a * ((i > 0) ? Q[i - 1] : 0.0f)) / denom;
    }

    float* x = new float[size];
    x[size - 1] = Q[size - 1];
    for (int i = size - 2; i >= 0; --i) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }

    delete[] P;
    delete[] Q;
    freeMatrix(matrix, size);
    delete[] constants;

    return x;
}

float* methodSI(SLAE &problem, float eps, int limit, int* iterations) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    float** alpha = new float*[size];
    float* beta = new float[size];
    
    for (int i = 0; i < size; ++i) {
        alpha[i] = new float[size];
        if (fabs(matrix[i][i]) < 1e-7f) {
            for (int k = 0; k <= i; ++k) delete[] alpha[k];
            delete[] alpha;
            delete[] beta;
            freeMatrix(matrix, size);
            delete[] constants;
            return nullptr;
        }
        
        beta[i] = constants[i] / matrix[i][i];
        for (int j = 0; j < size; ++j) {
            alpha[i][j] = (i != j) ? -matrix[i][j] / matrix[i][i] : 0.0f;
        }
    }

    float* prev = new float[size];
    float* curr = new float[size];
    
    for (int j = 0; j < size; ++j) {
        prev[j] = beta[j];
    }

    float inaccuracy = eps + 1;
    int k = 0;
    
    while (k < limit && inaccuracy > eps) {
        for (int i = 0; i < size; ++i) {
            float sum = 0.0f;
            for (int j = 0; j < size; ++j) {
                sum += alpha[i][j] * prev[j];
            }
            curr[i] = beta[i] + sum;
        }

        float* diff = new float[size];
        for (int i = 0; i < size; ++i) {
            diff[i] = curr[i] - prev[i];
        }
        
        inaccuracy = normVector(diff, size);
        delete[] diff;

        for(int i = 0; i < size; ++i) {
            if (isnan(curr[i]) || isinf(curr[i])) {
                for (int k = 0; k < size; ++k) delete[] alpha[k];
                    delete[] alpha;
                    delete[] beta;
                    delete[] prev;
                    delete[] curr;
                    freeMatrix(matrix, size);
                    delete[] constants;
                    return nullptr;
                }
        }

        swap(prev, curr);
        ++k;
    }

    float* result = new float[size];
    for (int i = 0; i < size; ++i) {
        result[i] = prev[i];
    }

    for (int i = 0; i < size; ++i) {
        delete[] alpha[i];
    }
    delete[] alpha;
    delete[] beta;
    delete[] prev;
    delete[] curr;
    
    freeMatrix(matrix, size);
    delete[] constants;

    *iterations = k;
    
    return result;
}

float* methodZ(SLAE &problem, float eps, int limit, int* iterations) {
    int size = problem.getSize();
    float** matrix = problem.getMatrixCopy();
    float* constants = problem.getConstCopy();

    float** alpha = new float*[size];
    float* beta = new float[size];
    
    for (int i = 0; i < size; ++i) {
        alpha[i] = new float[size];
        if (fabs(matrix[i][i]) < 1e-7f) {
            for (int k = 0; k <= i; ++k) delete[] alpha[k];
            delete[] alpha;
            delete[] beta;
            freeMatrix(matrix, size);
            delete[] constants;
            return nullptr;
        }
        
        beta[i] = constants[i] / matrix[i][i];
        for (int j = 0; j < size; ++j) {
            alpha[i][j] = (i != j) ? -matrix[i][j] / matrix[i][i] : 0.0f;
        }
    }

    float* prev = new float[size];
    float* curr = new float[size];
    
    for (int j = 0; j < size; ++j) {
        prev[j] = beta[j];
    }

    float inaccuracy = eps + 1;
    int k = 0;
    
    while (k < limit && inaccuracy > eps) {
        for (int i = 0; i < size; ++i) {
            float sum = 0.0f;
            for (int j = 0; j < size; ++j) {
                if (j < i) {
                    sum += alpha[i][j] * curr[j];   
                } else {
                    sum += alpha[i][j] * prev[j]; 
                }
            }
            curr[i] = beta[i] + sum;
        }

        float* diff = new float[size];
        for (int i = 0; i < size; ++i) {
            diff[i] = curr[i] - prev[i];
        }
        
        inaccuracy = normVector(diff, size);
        delete[] diff;

        for(int i = 0; i < size; ++i) {
            if (isnan(curr[i]) || isinf(curr[i])) {
                for (int k = 0; k < size; ++k) delete[] alpha[k];
                    delete[] alpha;
                    delete[] beta;
                    delete[] prev;
                    delete[] curr;
                    freeMatrix(matrix, size);
                    delete[] constants;
                    return nullptr;
                }
        }

        swap(prev, curr);
        ++k;
    }

    float* result = new float[size];
    for (int i = 0; i < size; ++i) {
        result[i] = prev[i];
    }

    for (int i = 0; i < size; ++i) {
        delete[] alpha[i];
    }
    delete[] alpha;
    delete[] beta;
    delete[] prev;
    delete[] curr;
    
    freeMatrix(matrix, size);
    delete[] constants;

    *iterations = k;
    
    return result;
}