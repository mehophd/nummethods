#ifndef SLAE_H
#define SLAE_H

#include <iostream>

using namespace std;

class SLAE {
    private:
        int size;               // порядок системы
        float** matrix;         // матрица коэффициентов
        float* constants;       // столбец с константами

    public:
        SLAE(int n) {
            size = n;
            matrix = new float*[n];
            for (int i = 0; i < n; ++i)
                matrix[i] = new float[n];
            constants = new float[n];
        }

        ~SLAE() {
            for (int i = 0; i < size; ++i)
                delete[] matrix[i];
            delete[] matrix;
        }

        int getSize() {
            return size;
        }

        void setCoef(int i, int j, float value) {
            matrix[i][j] = value;
        }

        void setConst(int i, float value) {
            constants[i] = value;
        }

        float getCoef(int i, int j) {
            return matrix[i][j];
        }

        float getConst(int i) {
            return constants[i];
        }

        float** getMatrixCopy() {
            float** copy = new float*[size];
            for (int i = 0; i < size; ++i) {
                copy[i] = new float[size];
                for (int j = 0; j < size; ++j) {
                    copy[i][j] = matrix[i][j];
                }
            }
            return copy; 
        }

        float* getConstCopy() {
            float* copy = new float[size];
            for (int i = 0; i < size; ++i) {
                copy[i] = constants[i];
            }
            return copy; 
        }

        void printSLAE() {
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    cout << matrix[i][j] << " ";
                }
                cout << "| " << constants[i] << endl;
            }
        }
};

#endif