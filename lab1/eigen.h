#ifndef EIGEN_H
#define EIGEN_H

float* methodJacobi(float** original, float eps, int size, int limit, float*** eigenvectors, int* iterations);

float* methodQR(float** original, float eps, int size, int limit, int* iterations, float* imaginaryPart);

#endif