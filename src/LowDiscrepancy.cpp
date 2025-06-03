#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// Constants for hQNorm
const double P0 = -0.322232431088;
const double P1 = -1.000000000000;
const double P2 = -0.342242088547;
const double P3 = -0.0204231210245;
const double P4 = -0.0000453642210148;
const double Q0 = 0.0993484626060;
const double Q1 = 0.588581570495;
const double Q2 = 0.531103462366;
const double Q3 = 0.103537752850;
const double Q4 = 0.00385607006340;

// Approximation to inverse normal CDF
double hQNorm(double p) {
    const double EPS = 1.0e-6;
    if (p >= 1.0 - EPS) p = 1.0 - EPS;
    if (p <= EPS) p = EPS;

    if (p == 0.5) return 0.0;

    double r = (p > 0.5) ? 1.0 - p : p;
    double t = sqrt(-2.0 * log(r));
    double a = ((((P4 * t + P3) * t + P2) * t + P1) * t + P0);
    double b = ((((Q4 * t + Q3) * t + Q2) * t + Q1) * t + Q0);
    double result = t + (a / b);
    return (p < 0.5) ? -result : result;
}

// Sieve of Eratosthenes for base primes
void initHalton(int dimen, vector<double>& quasi, vector<int>& base, int& offset) {
    base[0] = 2;
    if (dimen >= 2) base[1] = 3;

    int n = 3, nc = 2;
    while (nc < dimen) {
        bool isPrime = true;
        for (int i = 2; i <= sqrt(n); ++i) {
            if (n % i == 0) {
                isPrime = false;
                break;
            }
        }
        if (isPrime) base[nc++] = n;
        n += 1;
    }

    offset = 0;
    for (int nb = 0; nb < dimen; ++nb) {
        int iter = offset;
        quasi[nb] = 0.0;
        double half = 1.0 / base[nb];
        while (iter != 0) {
            int digit = iter % base[nb];
            quasi[nb] += digit * half;
            iter = (iter - digit) / base[nb];
            half /= base[nb];
        }
    }

    offset += 1;
}

void nextHalton(int dimen, vector<double>& quasi, const vector<int>& base, int& offset) {
    for (int nb = 0; nb < dimen; ++nb) {
        int iter = offset;
        quasi[nb] = 0.0;
        double half = 1.0 / base[nb];
        while (iter != 0) {
            int digit = iter % base[nb];
            quasi[nb] += digit * half;
            iter = (iter - digit) / base[nb];
            half /= base[nb];
        }
    }

    offset += 1;
}

void halton(vector<vector<double>>& qn, int n, int dimen, vector<int>& base, int& offset, int init, int transform) {
    vector<double> quasi(dimen, 0.0);

    if (init == 1) {
        initHalton(dimen, quasi, base, offset);
    }

    for (int i = 0; i < n; ++i) {
        nextHalton(dimen, quasi, base, offset);
        for (int j = 0; j < dimen; ++j) {
            qn[i][j] = (transform == 0) ? quasi[j] : hQNorm(quasi[j]);
        }
    }
}
