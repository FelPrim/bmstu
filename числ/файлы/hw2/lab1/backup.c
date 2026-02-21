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

constexpr int N = 13;
constexpr int n = 53;

static double f(double tau){
	return (20+0.2*N)/
		(1+(20+0.2*N)*(1+0.05*(54-n))*(tau-1)*(tau-1));	
}

constexpr int grid_node_num = 11;

typedef double (* const math_function)(double);

static void print_array(const int size, const double array[static size]){
	printf("[");
	int i = 0;
	assert(size > 1);
	for (; i < size-1; ++i)
		printf("%lf, ", array[i]);
	printf("%lf]\n", array[size-1]);
}

#define PRINT_ARRAY(x, y) printf("%s: ", #y); print_array(x, y)

// A[N*N]*x[N] = y[N]
static void SLAE_solve(
	const int N,
	const double A[static N * N],
	double x[static N],
	const double y[static N]
){
	// LAPACKE_dgesv перезаписывает A
	double A_cpy[N * N];
	memcpy(A_cpy, A, N*N*sizeof(double));

    int ipiv[N];
	
	memcpy(x, y, N*sizeof(double));

    int info = LAPACKE_dgesv(
        LAPACK_ROW_MAJOR,
        N,
        1,
        A_cpy,
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

static void Lagrange_polinom(
	const double a,
	const double b,
	const int node_count,
	const double tau[static node_count],
	math_function f_A,
	double x[static node_count]
){
	const int k = node_count-1;
	
	double y[node_count];
	for (int i = 0; i < node_count; ++i)
		y[i] = f_A(tau[i]);
	

	double tau_matrix[node_count*node_count];
	for (int i = 0; i < node_count; ++i){
		double current = 1;
		for (int j = 0; j < node_count; ++j){
			tau_matrix[i*node_count + j] = current;
			current *= tau[i];
		}
	}

	SLAE_solve(node_count, tau_matrix, x, y);

}

static double apply_polinom(
	const int arrsz,
	const double polinom[arrsz],
	const double arg
){
	double result = 0;
	double current = 1;
	for (int i = 0; i < arrsz; ++i){
		result += polinom[i] * current;
		current *= arg;
	}
	return result;
}

static void make_grid_uniform(
	const double a,
	const double b,
	const int node_count,
	double grid[static node_count]
){
	const double diff = (b-a)/(node_count-1);
	for (int i = 0; i < node_count; ++i){
		grid[i] = a + i*diff;	
	}
}

static double apply_lagrange(
	const int node_count,
	double grid[static node_count],
	double y[static node_count],
	double x
){
	double result = 0;
	const int n = node_count;
	for (int i = 0; i < n; ++i){
		double chisl = 1;
		for (int j = 0; j < n; ++j)
			chisl *= x-grid[j];

		double znam = 1;
		for (int j = 0; j < n; ++j){
			if (j == i)
				continue;
			znam *= grid[i] - grid[j];
		}
		result += chisl/znam*y[i];
	}
	return result;
}

static void make_grid_chebyshev(
	const double a,
	const double b,
	const int node_count,
	double tau[static node_count]
){
	const int k = node_count - 1;
	for (int j = 0; j < node_count; ++j){
		tau[j] = (a+b)/2.0 -
				 (b-a)/2.0 * cos(
							(M_PI*(2.0*j+1.0)) /
							(2.0*(k + 1.0))
						 );
	}
}


int main(){
	constexpr double a = 0;
	constexpr double b = 2;

	constexpr int grid1count = 11;
	constexpr int grid2count = 21;

	double grid1[grid1count];
	make_grid_uniform(a, b, grid1count, grid1);
	
	double Y[grid1count];
	for (int i = 0; i < grid1count; ++i)
		Y[i] = f(grid1[i]);

	PRINT_ARRAY(grid1count, grid1);
	//double x[grid1count];

	//Lagrange_polinom(a, b, grid1count, grid1, f, x);
//	PRINT_ARRAY(grid1count, x);
	
	double grid2[grid2count];
	make_grid_uniform(a, b, grid2count, grid2);

	double y11[grid2count], y12[grid2count], diff[grid2count];
	for (int i = 0; i < grid2count; ++i){
		y11[i] = f(grid2[i]);
	//	y12[i] = apply_polinom(grid1count, x, grid2[i]);
		y12[i] = apply_lagrange(grid1count, grid1, Y, grid2[i]);
	//	diff[i] = fabs(y11[i] - y12[i]);
	}
	PRINT_ARRAY(grid2count, y11);
	PRINT_ARRAY(grid2count, y12);
	
	FILE *gp = popen("gnuplot -persistent", "w");
	fprintf(gp, "set encoding utf8\n");
	fprintf(gp, "set terminal pngcairo size 640,480 enhanced font 'Times New Roman, 14'\n"
				"set output '1.png'\n"
				"set title 'Интерполяционный полином Лагранжа на равномерной сетке'\n"
				"plot '-' with lines title 'f(x)', '-' with lines title 'Полином'\n"
			);
	for (int i = 0; i < grid2count; ++i){
		fprintf(gp, "%f %f\n", grid2[i], y11[i]);
	}
	fprintf(gp, "e\n");
   for (int i = 0; i < grid2count; ++i)
        fprintf(gp, "%f %f\n", grid2[i], y12[i]);
    fprintf(gp, "e\n");
	fflush(gp);
	pclose(gp);
	//PRINT_ARRAY(grid2count, diff);
	// графики f и x_k*tau^k на grid2
	
	constexpr int GCcount = 11;
	double grid_chebyshev[GCcount];
	make_grid_chebyshev(a, b, GCcount, grid_chebyshev);
	for (int i = 0; i < GCcount; ++i)
		Y[i] = f(grid_chebyshev[i]);

	PRINT_ARRAY(GCcount, grid_chebyshev);
	
//	double x_chebyshev[grid1count];

	//Lagrange_polinom(a, b, GCcount, grid_chebyshev, f, x_chebyshev);
	//PRINT_ARRAY(GCcount, x_chebyshev);

	for (int i = 0; i < grid2count; ++i){
		y11[i] = f(grid2[i]);
		y12[i] = apply_lagrange(GCcount, grid_chebyshev, Y, grid2[i]);
		//	y12[i] = apply_polinom(GCcount, x_chebyshev, grid2[i]);
	//	diff[i] = fabs(y11[i] - y12[i]);
	}
	PRINT_ARRAY(grid2count, y11);
	PRINT_ARRAY(grid2count, y12);
	//PRINT_ARRAY(grid2count, diff);

	gp = popen("gnuplot -persistent", "w");
	fprintf(gp, "set encoding utf8\n");
	fprintf(gp, "set terminal pngcairo size 640,480 enhanced font 'Times New Roman, 14'\n"
				"set output '2.png'\n"
				"set title 'Интерполяционный полином Лагранжа на чебышёвской сетке'\n"
				"plot '-' with lines title 'f(x)', '-' with lines title 'Полином'\n"
			);
	for (int i = 0; i < grid2count; ++i){
		fprintf(gp, "%f %f\n", grid2[i], y11[i]);
	}
	fprintf(gp, "e\n");
   for (int i = 0; i < grid2count; ++i)
        fprintf(gp, "%f %f\n", grid2[i], y12[i]);
    fprintf(gp, "e\n");
	fflush(gp);
	pclose(gp);

	return 0;
}
