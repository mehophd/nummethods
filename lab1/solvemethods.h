#ifndef SOLVEMETHODS_H
#define SOLVEMETHODS_H

#include "Slae.h"

float* methodLU(SLAE &problem);

float* methodTDMA(SLAE &problem);

float* methodSI(SLAE &problem, float eps, int limit, int* iterations);

float* methodZ(SLAE &problem, float eps, int limit, int* iterations);

#endif