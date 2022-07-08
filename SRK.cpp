#include "stdafx.h"
#include <iostream>
#include <iomanip>
#define USE_MATH_DEFINE
using namespace std;

double f(double, double, double);
double g(double, double, double);

int main()
{
	const int N = 1;
	int i;
	double x[N], y[N], z[N], k0[N], k1[N], k2[N], k3[N], l0[N], l1[N], l2[N], l3[N];
	double n, h;
	double a = 0.5, b = 0.6;
	x[0] = 0.5;
	y[0] = 1;
	z[0] = 1;
	h = 0.1;
	n = (b - a) / h;
	// cout << "y[0]=" << y[0] << "\n";
	for (i = 0; i < n; i++)
	{
		x[i] = a + i*h;
		k0[i + 1] = h*f(x[i], y[i], z[i]);
		l0[i + 1] = h*g(x[i], y[i], z[i]);
		k1[i + 1] = h*f((x[i] + h / 2), (y[i] + k0[i + 1] / 2), (z[i] + l0[i + 1] / 2));
		l1[i + 1] = h*g((x[i] + h / 2), (y[i] + k0[i + 1] / 2), (z[i] + l0[i + 1] / 2));
		k2[i + 1] = h*f((x[i] + h / 2), (y[i] + k1[i + 1] / 2), (z[i] + l1[i + 1] / 2));
		l2[i + 1] = h*g((x[i] + h / 2), (y[i] + k1[i + 1] / 2), (z[i] + l1[i + 1] / 2));
		k3[i + 1] = h*f((x[i] + h), (y[i] + k2[i + 1]), (z[i] + l2[i + 1]));
		l3[i + 1] = h*g((x[i] + h), (y[i] + k2[i + 1]), (z[i] + l2[i + 1]));
		y[i + 1] = y[i] + (k0[i + 1] + 2 * k1[i + 1] + 2 * k2[i + 1] + k3[i + 1]) / 6;
		z[i + 1] = z[i] + (l0[i + 1] + 2 * l1[i + 1] + 2 * l2[i + 1] + l3[i + 1]) / 6;
		// cout << "y[" << i+1 << "]=" << y[i+1] <<"\n";
		y[i] = y[i + 1];
		z[i] = z[i + 1];
	}
	cout << "y=" << y[N] << std :: setprecision(6) << "\n";
	cout << "z=" << z[N] << std :: setprecision(6) << "\n";
	
	system("pause");
	return 0;
}

double f(double x, double y, double z)
{
	double f;
	f = (2 * y - x) / z;
	return f;
}

double g(double x, double y, double z)
{
	double g;
	g = 2 * y / (z + x);
	return g;
}