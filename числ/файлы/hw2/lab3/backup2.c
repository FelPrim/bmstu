#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifndef M_PI
constexpr double M_PI = acos(-1.0);
#endif

#include <lapacke.h>

constexpr double epsilon = 1e-10;

// A[N*N]*x[N] = y[N]
static void SLAE_solve(
	const int N,
	double A[static N * N],
	double x[static N],
	const double y[static N]
){
	// LAPACKE_dgesv перезаписывает A
//	double A_cpy[N * N];
//	memcpy(A_cpy, A, N*N*sizeof(double));

    int ipiv[N];
	
	memcpy(x, y, N*sizeof(double));

    int info = LAPACKE_dgesv(
        LAPACK_ROW_MAJOR,
        N,
        1,
        A,
        N,
        ipiv,
        x,
        1
    );

	if (info != 0){
		printf("LAPACKE_dgesv failed, info = %d\n", info);
        abort();
	}
}


static void print_array(const int size, const double array[static size]){
	printf("[");
	int i = 0;
	assert(size > 1);
	for (; i < size-1; ++i)
		printf("%lf, ", array[i]);
	printf("%lf]\n", array[size-1]);
}

#define PRINT_ARRAY(x, y) printf("%s: ", #y); print_array(x, y)
static void print_matrix(const int size, const double matrix[static size*size]){
	for (int i = 0; i < size; ++i)
		print_array(size, matrix + i*size);
}
#define PRINT_MATRIX(x, y) printf("%s:\n", #y); print_matrix(x, y)


constexpr int N = 13;
constexpr int n = 53;

constexpr double lambda = 1/(n - 47.0);
constexpr double a = 0;
constexpr double b = (N + 7)/((double) N);
constexpr double alpha = (N + 3) / ((double) N);
constexpr double beta = (n - 53) / 2.0;


constexpr int node_count = 5;
constexpr double h = (b-a)/node_count;

typedef double (* const math_function)(double);

static double K(double s, double tau){
	if (0 <= s && s <= tau){
		return s*(2*b - tau);
	}
	if (tau <= s && s <= (N + 7.0) / ((double) N)){
		return tau*(2*b - s);
	}
	assert(0);
}

static double y(double s){
	assert(a <= s && s <= b);
	return alpha *(s*s + beta);
}

static double Kij(int i, int j){
	auto si = a+i*h;
	auto tj = a+j*h;

	return K(si, tj);
}

static double yi(int i){
	auto si = a+i*h;
	return y(si);
}

double *x;

static double grid_fun(const double arg){
	assert(arg > a + epsilon && arg < b-h+epsilon);
	const double position = (arg - a)/h; 
	int index = (int) floor(position);
//	assert(index >= 0 && index < node_count);
	if (index == node_count-1)
		return x[node_count - 1];
	
	const double interop = position - index;
	printf("%lf %d\n", interop, index);
	assert(-epsilon < interop && interop < 1+epsilon);
	
	auto xi = x[index];
	auto xii = x[index+1];
	return xii*interop +
		   xi*(1-interop);
}


inline static double lambda_f(double t){
	return 10/39.0 * (
			2.4 * t*t
			- 4.8 * (sqrt(39/20.0) * t*sin(sqrt(20.0/39.0)*t)
					  + 39/20.0 * cos(sqrt(20.0/39.0)*t))
		);
}

inline static double lambda_f2(double t){
	return 20/6.0/13.0*(
		-sqrt(39/20.0)*t*cos(sqrt(20/39.0)*t)
		+39/20.0*sin(sqrt(20/39.0)*t)
	);
}

double C;
double analytic(){
	auto A = 4.8*(1-cos(sqrt(20/39.0)*20/13.0));
	auto B = sin(sqrt(20/39.0)*20/13.0);
	auto C = 16/13.0*400/169.0;
	auto D = lambda_f(b) - lambda_f(a);
	auto E = lambda_f2(b) - lambda_f2(a);
	// A + Bc = C + D + Ec
	auto F = C + D - A;
	auto H = B - E;
	auto C2 = F/H;
	printf("%lf, %lf, %lf, %lf, %lf\n",
			A, B, C, D, E);
	printf("Analytic: x(s) = %lf*(1-cos(sqrt(%d/%d)*s)) + %lf*sin(sqrt(%d/%d)*s)\n",
							4.8,			20,39,		C2,			20, 39);
	return C2;
}

double f(double x){
	auto omega = sqrt(20/39.0);
	return 4.8*(1-cos(omega*x)) + C*sin(omega*x);
}

int main(){
	x = malloc(node_count*sizeof(double));
	for (int i = 0; i < node_count; ++i)
		x[i] = a+i*h;
	printf("a: %lf, b: %lf, n: %d, N: %d, lambda: %lf\n",
				 a,      b,     n,     N,      lambda);
	printf("y(s) = %lf(s^2 + %lf)\n",
				 alpha,     beta);
	//PRINT_ARRAY(node_count, x);
	
	double *y = malloc(node_count*sizeof(double));
	double *F = malloc(node_count*node_count*sizeof(double));

	for (int i = 0; i < node_count; ++i)
		y[i] = yi(i);

	//PRINT_ARRAY(node_count, y);
	for (int i = 0; i < node_count; ++i){
		for (int j = 0; j < node_count; ++j){
			const double delta = (i == j)? 1: 0;
			F[i*node_count +j] = delta - lambda * Kij(i, j) * h;
		}
	}

//	PRINT_MATRIX(node_count, F);
	SLAE_solve(node_count, F, x, y);
//	PRINT_MATRIX(node_count, F);
//	PRINT_ARRAY(node_count, x);

	C = analytic();
	printf("%lf->%lf, %lf\n",
			b-h, grid_fun(b-h), f(b-h));

	double abs_error = -1;


	// Построение графиков аналитического и численного решений
	FILE *gp = popen("gnuplot -persistent", "w");
	fprintf(gp, "set encoding utf8\n");
	fprintf(gp, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
	fprintf(gp, "set output 'comparison0.png'\n");
	fprintf(gp, "set title 'Сравнение аналитического и численного решений'\n");
	fprintf(gp, "set xlabel 'x'\n");
	fprintf(gp, "set ylabel 'y'\n");
	fprintf(gp, "set grid\n");
	fprintf(gp, "plot '-' with lines linewidth 2 title 'Аналитическое решение f(x)', \\\n");
	fprintf(gp, "     '-' with linespoints linewidth 2 pointtype 7 title 'Численное решение (узлы)'\n");

	for (int i = 0; i < node_count; ++i) {
		double x_val  = a + i * h;         
		fprintf(gp, "%f %f\n", x_val, f(x_val));
		if (fabs(f(x_val) - x[i]) > abs_error)
			abs_error = fabs(f(x_val) - x[i]);
	}
	fprintf(gp, "e\n");

	for (int i = 0; i < node_count; ++i) {
		double x_node = a + i * h;          
		fprintf(gp, "%f %f\n", x_node, x[i]);
	}
	fprintf(gp, "e\n");

	fflush(gp);
	pclose(gp);

	printf("Maxerror: %lf\n", abs_error);

	free(x);
	free(F);
	free(y);
	return 0;
}
