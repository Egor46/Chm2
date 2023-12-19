#pragma once
#include <iostream>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

using namespace std;

double FresnelInt(double x) {
	double ep = 0.0000001;
	double k = x;
	double s = k;
 	int n = 0;
	while (abs(k) >= ep) {
		k *= ((-1) * (4 * n + 1) * (x * x * x * x) * (M_PI_2 * M_PI_2)) / ((2 * n + 2) * (2 * n + 1) * (4 * n + 5));
		n++;
		s += k;
	}
		return s;
}

double Derivative(double x) {
	double ep = 0.0000001;
	double k = 1;
	double s = k;
	int n = 0;
	while (abs(k) >= ep) {
		k *= ((-1) * (M_PI_2 * x * x) * (M_PI_2 * x * x)) / ((2 * n + 2) * (2 * n + 1));
		n++;
		s += k;
	}

	return s;
}