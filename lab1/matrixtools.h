#ifndef MATRIXTOOLS_H
#define MATRIXTOOLS_H

void freeMatrix(float** matrix, int size);

float determinant(float** matrix, int size);

float** inverseMatrix(float** matrix, int size);

void printMatrix(float** original, int size);

float normVector(float* vec, int size);

float normMatrix(float** matrix, int size);


#endif