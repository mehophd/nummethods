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

        if (maxA < 1e-7f) break;

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

bool qrDecomposition(float** original, int size, float*** Q, float ***R) {
    float** A = cloneMatrix(original, size);
    float** Hresult = new float*[size];

    for (int i = 0; i < size; ++i) {
        Hresult[i] = new float[size];
        for (int j = 0; j < size; ++j) {
           Hresult[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }


    for (int k = 0; k < size; ++k) {
        float sign = (A[k][k] >= 0) ? 1.0f : -1.0f;
        float norm = 0.0f;

        for (int i = k; i < size; ++i) {
            norm += A[i][k] * A[i][k];
        }
        norm = sqrtf(norm);

        if (norm < 1e-7f) continue;

        float* vec = new float[size];
        for (int i = 0; i < k; ++i) vec[i] = 0.0f;

        vec[k] = A[k][k] + sign * norm;
        for (int i = k + 1; i < size; ++i) {
            vec[i] = A[i][k];
        }

        float scalardown = 0.0f;
        for (int i = 0; i < size; ++i) {
            scalardown += vec[i] * vec[i];
        }
        scalardown = 2.0f / scalardown;

        float** H = new float*[size];
        for (int i = 0; i < size; ++i) {
            H[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                H[i][j] = (i == j) ? 1.0f : 0.0f;
                H[i][j] -= vec[i] * vec[j] * scalardown;
            }
        }

        // A[i] = H[i] * A[i - 1]
        float** temp = new float*[size];
        for (int i = 0; i < size; ++i) {
            temp[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                temp[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    temp[i][j] += H[i][m] * A[m][j];
                }
            }
        }

        // Q[i] = Q[i-1] * H[i]
        float** tempH = new float*[size];
        for (int i = 0; i < size; ++i) {
            tempH[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                tempH[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    tempH[i][j] += Hresult[i][m] * H[m][j];
                }
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Hresult[i][j] = tempH[i][j];
            }
        }

        freeMatrix(A, size);
        A = temp;
        freeMatrix(tempH, size);
        freeMatrix(H, size);
        delete[] vec;
    }

    *R = A;
    *Q = Hresult;
    
    return true;
}

float* methodQR(float** original, float eps, int size, int limit, int* iterations, float* imaginaryParts) {
    float** A = cloneMatrix(original, size);
    int k = 0;
    bool converged = false;

    while (k < limit && !converged) {
        float** Q = nullptr;
        float** R = nullptr;
        
        qrDecomposition(A, size, &Q, &R);

        // Новая матрица A_new = R * Q
        float** A_new = new float*[size];
        for (int i = 0; i < size; ++i) {
            A_new[i] = new float[size];
            for (int j = 0; j < size; ++j) {
                A_new[i][j] = 0.0f;
                for (int m = 0; m < size; ++m) {
                    A_new[i][j] += R[i][m] * Q[m][j];
                }
            }
        }

        converged = true;
        for (int i = 1; i < size; ++i) {
            if (fabs(A_new[i][i - 1]) > eps) {
                converged = false;
                break;
            }
        }

        freeMatrix(A, size);
        A = A_new;

        freeMatrix(Q, size);
        freeMatrix(R, size);

        ++k;
    }

    *iterations = k;

    float* eigenvalues = new float[size];
    int i = 0;
    
    while (i < size) {
        // Если это последний элемент или поддиагональный элемент близок к 0
        if (i == size - 1 || fabs(A[i + 1][i]) <= eps) {
            eigenvalues[i] = A[i][i];
            imaginaryParts[i] = 0.0f;
            ++i;
        } else {

            float trace = A[i][i] + A[i + 1][i + 1];
            float det = A[i][i] * A[i + 1][i + 1] - A[i][i + 1] * A[i + 1][i];
            float discriminant = trace * trace - 4.0f * det;

            if (discriminant >= 0.0f) {
                eigenvalues[i] = (trace + sqrtf(discriminant)) / 2.0f;
                eigenvalues[i + 1] = (trace - sqrtf(discriminant)) / 2.0f;
                imaginaryParts[i] = 0.0f;
                imaginaryParts[i + 1] = 0.0f;
            } else {
                float realPart = trace / 2.0f;
                float imagPart = sqrtf(-discriminant) / 2.0f;
                eigenvalues[i] = realPart;
                eigenvalues[i + 1] = realPart;
                imaginaryParts[i] = imagPart;
                imaginaryParts[i + 1] = -imagPart;
            }
            i += 2;
        }
    }

    freeMatrix(A, size);
    return eigenvalues;
}