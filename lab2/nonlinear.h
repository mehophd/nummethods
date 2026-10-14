#ifndef NONLINEAR_H
#define NONLINEAR_H

#include <vector>
#include <cmath>

using namespace std;

struct IterRecord1D {
    int k;
    float x;
    float fx;
    float error;
};

struct Vec2 {
    float x1;
    float x2;
};

struct IterRecord2D {
    int k;
    Vec2 x;
    float f1;
    float f2;
    float error;
};

float f_task1(float x);
float df_task1(float x);
float phi_task1(float x);

vector<IterRecord1D> solveNewton1D(float x0, float eps, int limit);
vector<IterRecord1D> solveSimpleIter1D(float x0, float eps, int limit);

void task1();

Vec2 F_task2(Vec2 x);
void getJacobian_task2(Vec2 x, float J[2][2]);
Vec2 phi_task2(Vec2 x);

vector<IterRecord2D> solveNewton2D(Vec2 x0, float eps, int limit);
vector<IterRecord2D> solveSimpleIter2D(Vec2 x0, float eps, int limit);

void task2();

#endif