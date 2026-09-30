#include <cmath>
#include <algorithm>
#include "eigen.h"
#include "matrixtools.h"

using namespace std;

float* methodJacobi(float** original, float eps, int size, int limit, float*** eigenvectors, int* iterations) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (fabs(original[i][j] - original[j][i]) > 1e-7f) {
                return nullptr;
            }
        }
    }

    float** curr = cloneMatrix(original, size);
    float** Uresult = new float*[size];

    for (int i = 0; i < size; ++i) {
        Uresult[i] = new float[size];
        for (int j = 0; j < size; ++j) {
            Uresult[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }

    float tsum = eps + 1.0f;
    int k = 0;

    while (k < limit && tsum > eps) {
        float maxA = 0.0f;
        int maxRow = 0, maxCol = 1;

        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (fabs(curr[i][j]) > maxA) {
                    maxA = fabs(curr[i][j]);
                    maxRow = i;
                    maxCol = j;
                }
            }
        }

        if (maxA < 1e-10f) break;

        float angleF;
        if (fabs(curr[maxRow][maxRow] - curr[maxCol][maxCol]) < 1e-10f) {
            angleF = M_PI / 4.0f;
        } else {
            angleF = atan(2 * curr[maxRow][maxCol] / 
                         (curr[maxRow][maxRow] - curr[maxCol][maxCol])) / 2.0f;
        }

        float c = cos(angleF);
        float s = sin(angleF);

        float** U = new float*[size];
        for (int i = 0; i < size; ++i) {
            U[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                U[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
        U[maxRow][maxCol] = -s;
        U[maxCol][maxRow] = s;
        U[maxRow][maxRow] = c;
        U[maxCol][maxCol] = c;

        float** temp = new float*[size];
        for (int i = 0; i < size; ++i) {
            temp[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                temp[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    temp[i][j] += U[m][i] * curr[m][j];
                }
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                curr[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    curr[i][j] += temp[i][m] * U[m][j];
                }
            }
        }

        float** tempU = new float*[size];
        for (int i = 0; i < size; ++i) {
            tempU[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                tempU[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    tempU[i][j] += Uresult[i][m] * U[m][j];
                }
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Uresult[i][j] = tempU[i][j];
            }
        }

        tsum = 0.0f;
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                tsum += curr[i][j] * curr[i][j];
            }
        }
        tsum = sqrt(tsum);

        freeMatrix(temp, size);
        freeMatrix(tempU, size);
        freeMatrix(U, size);
        ++k;
    }

    float* eigenvalues = new float[size];
    for (int i = 0; i < size; ++i) {
        eigenvalues[i] = curr[i][i];
    }

    *eigenvectors = Uresult;

    freeMatrix(curr, size);

    *iterations = k;

    return eigenvalues;
}