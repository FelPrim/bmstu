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

constexpr int N = 13;
constexpr int n = 53;

static double f(double x){
	return 2*(54-n)*sin(M_PI*x)*sqrt(
			1+2*N*(x-0.5)+N*N
		);
}


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

//static void make_grid_uniform(
//	const double a,
//	const double b,
//	const int node_count,
//	double grid[static node_count]
//){
//	const double diff = (b-a)/(node_count-1);
//	for (int i = 0; i < node_count; ++i){
//		grid[i] = a + i*diff;	
//	}
//}

int find_areascount(int nodecount, int power){
	//int result = 0;
	//for (int i = 0; i < nodecount; i += power){
	//	++result;
	//}
	//return result;
	int result = (nodecount-1 + power-1)/power;
	printf("Areas: %d\n", result);
	return result;
}

double apply_s2(
	const int arr_sz,
	double splines[static arr_sz][3], 
	const double x0, const double h,
	double x){
	
	int index = trunc((x-x0)/(2*h));
	if (index >= arr_sz){
		printf("x = %f, x0 = %f, h = %f, index = %d, arr_sz = %d\n", 
				     x,      x0,      h,      index,      arr_sz);
		abort();
	}
	const double a = splines[index][0];
	const double b = splines[index][1];
	const double c = splines[index][2];
	const double x_k = x0 + h * (index + 0.5);

	return c*(x-x_k)*(x-x_k)+b*(x-x_k)+a;
}

int main(){
	constexpr double a = 0;
	constexpr double b = 2;

	constexpr double h = 0.1;
	constexpr int node_count = 21; // -> 20 областей
	

	double grid[node_count];
	double y[node_count];
	for (int i = 0; i < node_count; ++i){
		grid[i] = a + h*i;
		y[i] = f(grid[i]);
	}

	PRINT_ARRAY(node_count, grid);
	PRINT_ARRAY(node_count, y);

	// интерполяционный многочлен s2, используя квадратичные сплайны дефекта 1

	const int quadrareas = find_areascount(node_count, 2);
	double quadrsplines[quadrareas][3];
	

	// База индукции
	
	const double a1 = y[0];
	const double b1 = (y[1]-y[0])/h;
	const double c1 = 0;

	double a_k[node_count], b_k[node_count], c_k[node_count];
	a_k[0] = a1;
	b_k[0] = b1;
	c_k[0] = c1;
	// Шаг
	//
	for (int i = 0; i < node_count; ++i){
		b_k[i] = 
	}
	

	/*
	quadrsplines[0][0] = y[0];
	quadrsplines[0][1] = (y[1]-y[0])/h;
	quadrsplines[0][2] = 0;
	
	for (int i = 1; i < quadrareas; ++i){
		const double a_im = quadrsplines[i-1][0];
		const double b_im = quadrsplines[i-1][1];
		const double c_im = quadrsplines[i-1][2];

		const double a_i = y[i-1];
		const double b_i = b_im + c_im*2*h;
		const double c_i = (y[i]-y[i-1] - b_i*h) / (h * h);

		quadrsplines[i][0] = a_i;
		quadrsplines[i][1] = b_i;
		quadrsplines[i][2] = c_i;
	}
*/
	//const double a0 = (y[2] - 2*y[1] + y[0])/(2*h*h);
	//const double b0 = (y[1]-y[0] - a0*(2*h*a+h*h))/h;
	//const double c0 = y[0] - a0*a*a - b0*a;
	//
	//quadrsplines[0][0] = a0;
	//quadrsplines[0][1] = b0;
	//quadrsplines[0][2] = c0;


	for (int i = 0; i < quadrareas; ++i){
		print_array(3, quadrsplines[i]);
	}		
	
	FILE *gp = popen("gnuplot -persistnent", "w");
	fprintf(gp, "set encoding utf8\n"
				"set terminal pngcairo size 640,480 enhanced font 'Times New Roman, 14'\n"
				"set output '1.png'\n"
				"set title 'Квадратичные сплайны'\n"
				"plot '-' with lines title 'f(x)', '-' with lines title 'Сплайны'\n"
			);
	
	for (double n = a; n < b; n += 0.01)
		fprintf(gp, "%f %f\n", n, f(n));
	fprintf(gp, "e\n");

	for (double n = a; n < b; n += 0.01)
		fprintf(gp, "%f %f\n", n, apply_s2(
					quadrareas, 
					quadrsplines,
					a, h,
					n));
	fprintf(gp, "e\n");
	fflush(gp);
	pclose(gp);

	// интерполяционный многочлен s3, используя кубические сплайны дефекта 1

	// интерполяционный многочлен H3, используя кубические Эрмитовы сплайны дефекта 1
	
	// значения производной функции и сплайнов в узлах



	return 0;
}
