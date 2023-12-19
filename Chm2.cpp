#include "Header.h"
#include <random>
#include <raylib.h>
#include <fstream>

using namespace std;

double* nodes;
double* div_diff;
int N, n;

double a = 0.4, b = 4;

void calculateDivDiff(int N) {

	for (int i = 0; i < N; i++) {
		div_diff[i] = FresnelInt(nodes[i]);
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = N - 1; j >= i; j--)
		{
			div_diff[j] = (div_diff[j] - div_diff[j - 1]) / (nodes[j] - nodes[j - i]);
		}
	}
}

double Newton(double x) {
	double res = div_diff[0];
	double temp = 1;
	for (int j = 1; j < N; j++)
	{
		temp *= (x - nodes[j - 1]);
		res += (div_diff[j] * temp);
	}
	return res;
}

void generateNodes(int N) {
	for (int i = 0; i < N; i++) {
		nodes[i] = 1.0 / 2.0 * ((b + a) + (b - a) * cos((2 * (i + 1) - 1) * M_PI / (2 * N)));
	}
}

double Lagrange(double x) {
	double result = 0;
	double prod = 1;
	for (int i = 0; i < N; i++) {
		prod = 1;
		for (int j = 0; j < N; j++) {
			if (i == j) j++;
			prod *= (x - nodes[j]) / (nodes[i] - nodes[j]);
		}
		result += FresnelInt(nodes[i]) * prod;
	}
	return result;
}

double error(double (*f)(double), double x, int N) {
	return abs(f(x) - FresnelInt(x));
}

int main() {
	cin >> N;
	cin >> n;
	nodes = new double[N];
	for (int i = 0; i < N; i++) nodes[i] = a + i * (b - a) / (N - 1);
	div_diff = new double[N];
	calculateDivDiff(N);
	double* ksi = new double[n];
	for (int i = 0; i < n; i++) ksi[i] = a+0.2 + (b - a) / (n - 1) * i;
	std::ofstream out1("out1.txt"), out3("out3.txt");
	for (int i = 0; i < n - 1; i++) {
		out1 << ksi[i] << ' ' << error(Lagrange, ksi[i], N) << endl;
		out3 << ksi[i] << ' ' << abs(FresnelInt(ksi[i]) - Newton(ksi[i])) << endl;
	}
	generateNodes(N);

	calculateDivDiff(N);
	ofstream out2("out2.txt"), out4("out4.txt");
	for (int i = 0; i < n - 1; i++) {
		out2 << ksi[i] << ' ' << error(Lagrange, ksi[i], N) << endl;
		out4 << ksi[i] << ' ' << error(Newton, ksi[i], N) << endl;
	}
	out1.close(); out2.close(); out3.close(); out4.close();
	delete[] ksi;
	delete[] nodes;
	delete[] div_diff;
}

