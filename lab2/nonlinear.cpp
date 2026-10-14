#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "nonlinear.h"

using namespace std;

// f(x) = x^3 + x^2 - 2x - 1 = 0

float fTask1(float x) {
    return x * x * x + x * x - 2.0f * x - 1.0f;
}

float dfTask1(float x) {
    return 3.0f * x * x + 2.0f * x - 2.0f;
}

// x^3 = -x^2 + 2x + 1  =>  x = cbrt(-x^2 + 2x + 1)
float phiTask1(float x) {
    return cbrt(-x * x + 2.0f * x + 1.0f);
}

vector<IterRecord1D> solveNewton1D(float x0, float eps, int limit) {
    vector<IterRecord1D> history;
    float x = x0;
    
    IterRecord1D firstRecord = {0, x, fTask1(x), 0.0f};
    history.push_back(firstRecord);

    for (int k = 0; k < limit; ++k) {
        float fx = fTask1(x);
        float dfx = dfTask1(x);

        if (fabs(dfx) < 1e-7f) break;

        float x_next = x - fx / dfx;
        float err = fabs(x_next - x);

        history.push_back({k + 1, x_next, fTask1(x_next), err});

        if (err < eps) break;
        x = x_next;
    }
    return history;
}

vector<IterRecord1D> solveSimpleIter1D(float x0, float eps, int limit) {
    vector<IterRecord1D> history;
    float x = x0;
    
    IterRecord1D firstRecord = {0, x, fTask1(x), 0.0f};
    history.push_back(firstRecord);

    for (int k = 0; k < limit; ++k) {
        float x_next = phiTask1(x);
        float err = fabs(x_next - x);

        history.push_back({k + 1, x_next, fTask1(x_next), err});

        if (err < eps) break;
        x = x_next;
    }
    return history;
}

void task1() {
    cout << "\nx^3 + x^2 - 2x - 1 = 0\n";
    
    float eps = 1e-5f;
    int limit = 100;
    
    float x0_newton = 1.5f;
    float x0_iter = 1.5f;

    cout << "\nМетод Ньютона (начальное приближение x0 = " << x0_newton << ")\n";
    vector<IterRecord1D> histNewton = solveNewton1D(x0_newton, eps, limit);
    
    cout << left << setw(5) << "k" 
         << setw(15) << "x_k" 
         << setw(15) << "f(x_k)" 
         << setw(15) << "|x_k - x_{k-1}|" << "\n";
    cout << "\n";
    for (const auto& rec : histNewton) {
        cout << left << setw(5) << rec.k 
             << setw(15) << fixed << setprecision(6) << rec.x 
             << setw(15) << rec.fx 
             << setw(15) << rec.error << "\n";
    }

    cout << "\nМетод простой итерации (начальное приближение x0 = " << x0_iter << ")\n";
    vector<IterRecord1D> histIter = solveSimpleIter1D(x0_iter, eps, limit);
    
    cout << left << setw(5) << "k" 
         << setw(15) << "x_k" 
         << setw(15) << "f(x_k)" 
         << setw(15) << "|x_k - x_{k-1}|" << "\n";
    cout << "\n";
    for (const auto& rec : histIter) {
        cout << left << setw(5) << rec.k 
             << setw(15) << fixed << setprecision(6) << rec.x 
             << setw(15) << rec.fx 
             << setw(15) << rec.error << "\n";
    }
}


// x1^2 + x2^2 - 16 = 0
// x1 - e^(x2) + 4 = 0

Vec2 F_task2(Vec2 x) {
    return {
        x.x1 * x.x1 + x.x2 * x.x2 - 16.0f,
        x.x1 - exp(x.x2) + 4.0f
    };
}

void getJacobian_task2(Vec2 x, float J[2][2]) {
    J[0][0] = 2.0f * x.x1;
    J[0][1] = 2.0f * x.x2;
    J[1][0] = 1.0f;
    J[1][1] = -exp(x.x2);
}

// x = phi(x). 
// e^(x2) = x1 + 4 => x2 = ln(x1 + 4)
// x1^2 = 16 - x2^2 => x1 = sqrt(16 - x2^2)
Vec2 phi_task2(Vec2 x) {
    float next_x1 = sqrt(16.0f - x.x2 * x.x2);
    float next_x2 = log(x.x1 + 4.0f);
    return {next_x1, next_x2};
}

vector<IterRecord2D> solveNewton2D(Vec2 x0, float eps, int limit) {
    vector<IterRecord2D> history;
    Vec2 x = x0;
    
    Vec2 f = F_task2(x);
    IterRecord2D firstRecord = {0, x, f.x1, f.x2, 0.0f};
    history.push_back(firstRecord);

    for (int k = 0; k < limit; ++k) {
        float J[2][2];
        getJacobian_task2(x, J);
        
        float det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-7f) break;

        // J * dx = -F
        float dx1 = (-f.x1 * J[1][1] - (-f.x2) * J[0][1]) / det;
        float dx2 = (J[0][0] * (-f.x2) - J[1][0] * (-f.x1)) / det;

        Vec2 x_next = {x.x1 + dx1, x.x2 + dx2};
        float err = max(fabs(dx1), fabs(dx2));

        Vec2 f_next = F_task2(x_next);
        history.push_back({k + 1, x_next, f_next.x1, f_next.x2, err});

        if (err < eps) break;
        x = x_next;
        f = f_next;
    }
    return history;
}

vector<IterRecord2D> solveSimpleIter2D(Vec2 x0, float eps, int limit) {
    vector<IterRecord2D> history;
    Vec2 x = x0;
    
    Vec2 f = F_task2(x);
    IterRecord2D firstRecord = {0, x, f.x1, f.x2, 0.0f};
    history.push_back(firstRecord);

    for (int k = 0; k < limit; ++k) {
        Vec2 x_next = phi_task2(x);
        float err = max(fabs(x_next.x1 - x.x1), fabs(x_next.x2 - x.x2));

        Vec2 f_next = F_task2(x_next);
        history.push_back({k + 1, x_next, f_next.x1, f_next.x2, err});

        if (err < eps) break;
        x = x_next;
    }
    return history;
}

void task2() {
    cout << "\nx1^2 + x2^2 - 16 = 0\n";
    cout << "x1 - e^(x2) + 4 = 0\n";

    float eps = 1e-5f;
    int limit = 100;

    Vec2 x0 = {3.4f, 2.0f};

    cout << "\nМетод Ньютона (начальное приближение x1=" << x0.x1 << ", x2=" << x0.x2 << ")\n";
    vector<IterRecord2D> histNewton = solveNewton2D(x0, eps, limit);
    
    cout << left << setw(4) << "k" 
         << setw(12) << "x1_k" << setw(12) << "x2_k"
         << setw(12) << "f1(x)" << setw(12) << "f2(x)"
         << setw(12) << "Error" << "\n";
    cout << "\n";
    for (const auto& rec : histNewton) {
        cout << left << setw(4) << rec.k 
             << setw(12) << fixed << setprecision(6) << rec.x.x1 
             << setw(12) << rec.x.x2
             << setw(12) << rec.f1 << setw(12) << rec.f2
             << setw(12) << rec.error << "\n";
    }

    cout << "\nМетод простой итерации (начальное приближение x1=" << x0.x1 << ", x2=" << x0.x2 << ")\n";
    vector<IterRecord2D> histIter = solveSimpleIter2D(x0, eps, limit);
    
    cout << left << setw(4) << "k" 
         << setw(12) << "x1_k" << setw(12) << "x2_k"
         << setw(12) << "f1(x)" << setw(12) << "f2(x)"
         << setw(12) << "Error" << "\n";
    cout << "\n";
    for (const auto& rec : histIter) {
        cout << left << setw(4) << rec.k 
             << setw(12) << fixed << setprecision(6) << rec.x.x1 
             << setw(12) << rec.x.x2
             << setw(12) << rec.f1 << setw(12) << rec.f2
             << setw(12) << rec.error << "\n";
    }
}