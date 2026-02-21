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

#define or ||

constexpr int N = 13;
constexpr int n = 53;

static double f(double x){
	return 2*(54-n)*sin(M_PI*x)*sqrt(
			1+2*N*(x-0.5)+N*N
		);
}


typedef double (* const math_function)(double);

double manual_derivative(math_function f, double x){ 
	const double epsilon = 1e-10; 
	return (f(x+epsilon) - f(x))/
				epsilon; 
}



double manual_doublederivative(math_function f, double x){ 
	const double epsilon = 1e-10; 
	return (f(x+epsilon) - 2*f(x)+f(x-epsilon))/
				(epsilon*epsilon); 
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

int find_areascount(int nodecount, int power){
	int result = (nodecount-1 + power-1)/power;
	printf("Areas: %d\n", result);
	return result;
}

double apply_s2(
	const int arr_sz,
	double s2a[static arr_sz],
	double s2b[static arr_sz],
	double s2c[static arr_sz], 
	const double x0, const double h,
	double x){
	
//	int index = trunc((x-x0)/(2*h));

	int index = (int)floor((x - x0) / h);
	if (index >= arr_sz or index < 0){
		printf("x = %f, x0 = %f, h = %f, index = %d, arr_sz = %d\n", 
				     x,      x0,      h,      index,      arr_sz);
		abort();
	}
	const double a = s2a[index];
	const double b = s2b[index];
	const double c = s2c[index];
	const double x_k = x0 + h * (index + 0.5);

	return c*(x-x_k)*(x-x_k)+b*(x-x_k)+a;
}

double derivative_s2(
	const int arr_sz,
	double s2a[static arr_sz],
	double s2b[static arr_sz],
	double s2c[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_s2(arr_sz, s2a, s2b, s2c, x0, h, x + epsilon);
	const double f_x = apply_s2(arr_sz, s2a, s2b, s2c, x0, h, x);
	return (f_xdx - f_x) / epsilon;
}


double doublederivative_s2(
	const int arr_sz,
	double s2a[static arr_sz],
	double s2b[static arr_sz],
	double s2c[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_s2(arr_sz, s2a, s2b, s2c, x0, h, x + epsilon);
	const double f_x = apply_s2(arr_sz, s2a, s2b, s2c, x0, h, x);
	const double f_xmdx = apply_s2(arr_sz, s2a, s2b, s2c, x0, h, x - epsilon);
	return (f_xdx - 2*f_x + f_xmdx) / (epsilon*epsilon);
}

double apply_s3(
	const int arr_sz,
	double s3a[static arr_sz],
	double s3b[static arr_sz],
	double s3c[static arr_sz], 
	double s3d[static arr_sz],
	const double x0, const double h,
	double x){
	
//	int index = trunc((x-x0)/(2*h));
	int index = (int)floor((x - x0) / h);
	if (index >= arr_sz or index < 0){
		printf("x = %f, x0 = %f, h = %f, index = %d, arr_sz = %d\n", 
				     x,      x0,      h,      index,      arr_sz);
		abort();
	}
	const double a = s3a[index];
	const double b = s3b[index];
	const double c = s3c[index];
	const double d = s3d[index];
	const double x_k = x0 + h * index;

	return	d*(x-x_k)*(x-x_k)*(x-x_k)+
			c*(x-x_k)*(x-x_k)+
			b*(x-x_k) + a;
}

double derivative_s3(
	const int arr_sz,
	double s3a[static arr_sz],
	double s3b[static arr_sz],
	double s3c[static arr_sz], 
	double s3d[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_s3(arr_sz, s3a, s3b, s3c, s3d, x0, h, x + epsilon);
	const double f_x = apply_s3(arr_sz, s3a, s3b, s3c, s3d, x0, h, x);
	return (f_xdx - f_x) / epsilon;
}


double doublederivative_s3(
	const int arr_sz,
	double s3a[static arr_sz],
	double s3b[static arr_sz],
	double s3c[static arr_sz], 
	double s3d[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_s3(arr_sz, s3a, s3b, s3c, s3d, x0, h, x + epsilon);
	const double f_x = apply_s3(arr_sz, s3a, s3b, s3c, s3d, x0, h, x);
	const double f_xmdx = apply_s3(arr_sz, s3a, s3b, s3c, s3d, x0, h, x - epsilon);
	return (f_xdx - 2*f_x + f_xmdx) / (epsilon*epsilon);
}

double apply_ermit(
	const int arr_sz,
	double h3a[static arr_sz],
	double h3b[static arr_sz],
	double h3c[static arr_sz], 
	double h3d[static arr_sz],
	const double x0, const double h,
	double x){
	
//	int index = trunc((x-x0)/(2*h));

	int index = (int)floor((x - x0) / h);
	if (index >= arr_sz or index < 0){
		printf("x = %f, x0 = %f, h = %f, index = %d, arr_sz = %d\n", 
				     x,      x0,      h,      index,      arr_sz);
		abort();
	}
	const double a = h3a[index];
	const double b = h3b[index];
	const double c = h3c[index];
	const double d = h3d[index];
	const double x_k = x0 + h * index;

	return	a*(x-x_k)*(x-x_k)*(x-x_k)+
			b*(x-x_k)*(x-x_k)+
			c*(x-x_k) + d;
}


double derivative_h3(
	const int arr_sz,
	double h3a[static arr_sz],
	double h3b[static arr_sz],
	double h3c[static arr_sz], 
	double h3d[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_ermit(arr_sz, h3a, h3b, h3c, h3d, x0, h, x + epsilon);
	const double f_x = apply_ermit(arr_sz, h3a, h3b, h3c, h3d, x0, h, x);
	return (f_xdx - f_x) / epsilon;
}


double doublederivative_h3(
	const int arr_sz,
	double h3a[static arr_sz],
	double h3b[static arr_sz],
	double h3c[static arr_sz], 
	double h3d[static arr_sz], 
	const double x0, const double h,
	double x){
	
	const double epsilon = 1e-10;
	const double f_xdx = apply_ermit(arr_sz, h3a, h3b, h3c, h3d, x0, h, x + epsilon);
	const double f_x = apply_ermit(arr_sz, h3a, h3b, h3c, h3d, x0, h, x);
	const double f_xmdx = apply_ermit(arr_sz, h3a, h3b, h3c, h3d, x0, h, x - epsilon);
	return (f_xdx - 2*f_x + f_xmdx) / (epsilon*epsilon);
}

int main(){

    // таблица значений функций
	// таблица значений коэффициентов
	// решение систем
	// найденные сплайны для всех промежутков
	// графики сплайнов и производных
	//

	constexpr double a = 0;
	constexpr double b = 2;

	constexpr double h = 0.1;
	constexpr int node_count = 21; // -> 20 областей
	

	double grid[node_count];
	double y[node_count];
	double y_shtrix[node_count];

	constexpr double epsilon = 1e-10;

	for (int i = 0; i < node_count; ++i){
		grid[i] = a + h*i;
		y[i] = f(grid[i]);
		y_shtrix[i] = (f(grid[i]+epsilon) - y[i])/epsilon;
	}

	PRINT_ARRAY(node_count, grid);
	PRINT_ARRAY(node_count, y);

	// интерполяционный многочлен s2, используя квадратичные сплайны дефекта 1

	const int quadrareas = 20;
	
	double s2a[quadrareas];
	double s2b[quadrareas];
	double s2c[quadrareas];

	// База индукции
	
	const double s2b0 = (y[1]-y[0])/h;
	// Доп условие: в точке (a, f(a)) наклон параболы такой же, как у производной
	const double derivative0 = (f(a+epsilon)-y[0])/epsilon;
	printf("Derivative: %lf\n", derivative0);
	const double s2c0 = (s2b0-derivative0)/h;
	const double s2a0 = y[0]+
						h/2.0*s2b0-
						h*h/4.0*s2c0;

	s2a[0] = s2a0;
	s2b[0] = s2b0;
	s2c[0] = s2c0;

	// Шаг индукции
	for (int i = 1; i < node_count-1; ++i){
		s2b[i] = (y[i+1]-y[i])/
					h;
		s2c[i] = (y[i+1]-2*y[i]+y[i-1])/
					(h*h)					- s2c[i-1];
		s2a[i] = y[i]+
				h/2.0*s2b[i]-
				h*h/4.0*s2c[i];
	}
	
	PRINT_ARRAY(node_count-1, s2a);
	PRINT_ARRAY(node_count-1, s2b);
	PRINT_ARRAY(node_count-1, s2c);


	// интерполяционный многочлен s3, используя кубические сплайны дефекта 1
	
	const int cubareas = 20;
	
	double s3a[cubareas];
	double s3b[cubareas];
	double s3c[cubareas+1]; // для метода прогонки
	double s3d[cubareas];

	// "База индукции" - не получается разделить на базу и шаг, т.к. нужно решить СЛАУ
	// Доп условия: f''(a) и f''(b) для сплайнов и оригинала совпадают
	// Эти условия задают СЛАУ. простую, но СЛАУ
	
	const double doublederivative0 = (f(a+epsilon)-2*f(a)+f(a-epsilon))/(epsilon*epsilon);

	const double doublederivative1 = (f(b+epsilon)-2*f(b)+f(b-epsilon))/(epsilon*epsilon);

	printf("Double derivatives: %lf, %lf\n", doublederivative0, doublederivative1);
	
	for (int i = 0; i < cubareas; ++i)
		s3a[i] = y[i];
	
	s3c[0] = doublederivative0/2.0;
	s3c[cubareas] = doublederivative1/2.0;
	/* 
	 * h^2 (C_{k+2} + 4 C_{k+1} + C_k) = 3 (y_{k+2} - 2y_{k+1} + y_k)
	 * C_{k+1} + 4 C_k + C_{k-1} = F_k, F_k = 3/h^2 (y_{k+1} - 2y_k + y_{k-1})
	 * C_0 = Q0, C_n = Qn
	 * */
	double alpha[cubareas], beta[cubareas];
	alpha[0] = 0;
	beta[0] = s3c[0];

	for (int i = 1; i < cubareas; ++i){
		assert(4 + alpha[i-1] != 0);
		alpha[i] = - 1 / (4 + alpha[i-1]);
		const double F_k = 3*( y[i+1] - 2*y[i] + y[i-1]) / 
			(h * h);
		beta[i] = (F_k - beta[i-1]) /
			(4 + alpha[i-1]);
	}
	for (int i = cubareas - 1; i > 0; --i){
		s3c[i] = alpha[i] * s3c[i+1] + beta[i];
		s3d[i] = (s3c[i+1] - s3c[i])/(3 * h);
		s3b[i] = (y[i+1]-y[i])/h - h*s3c[i] - h*h*s3d[i];
	}
	s3d[0] = (s3c[1] - s3c[0])/ (3 * h);
	s3b[0] = (y[1]-y[0])/h - h*s3c[0] - h*h*s3d[0];

	PRINT_ARRAY(cubareas, s3a);
	PRINT_ARRAY(cubareas, s3b);
	PRINT_ARRAY(cubareas, s3c);
	PRINT_ARRAY(cubareas, s3d);

	// интерполяционный многочлен H3, используя кубические Эрмитовы сплайны дефекта 1
	
	constexpr int ermitarea = node_count-1;
	
	double h3a[ermitarea], h3b[ermitarea], h3c[ermitarea], h3d[ermitarea];

	for (int i = 0; i < ermitarea; ++i){
		h3d[i] = y[i];
		h3c[i] = y_shtrix[i];
		
		const double A = y[i+1] - y[i] - h*y_shtrix[i];
		const double B = y_shtrix[i+1] - y_shtrix[i];

		// A = h^2 b + h^3  a
		// B = 2h  b + 3h^2 a

		h3a[i] = B/(h*h) - 2 * A / (h*h*h);
		h3b[i] = 3*A / (h * h) - B / h;
	}

	PRINT_ARRAY(ermitarea, h3a);
	PRINT_ARRAY(ermitarea, h3b);
	PRINT_ARRAY(ermitarea, h3c);
	PRINT_ARRAY(ermitarea, h3d);

	// значения производной функции и сплайнов в узлах
	
	double ys2 = apply_s2(
		quadrareas,
		s2a, s2b, s2c, a, h, 
		(a+b)/2.0 + h/3
	);
	double ys3 = apply_s3(
		cubareas,
		s3a, s3b, s3c, s3d,
		a, h,
		(a+b)/2.0 + h/3
	);
	double yh3 = apply_ermit(
		ermitarea,
		h3a, h3b, h3c, h3d, 
		a, h, 
		(a+b)/2.0 + h/3
	);
	printf("%lf\n", (a+b)/2.0+h/3);
	printf("%lf==%lf==%lf==%lf\n",
			ys2, ys3, yh3, f((a+b)/2.0+h/3));

    // таблица значений функций
	// таблица значений коэффициентов
	// решение систем
	// найденные сплайны для всех промежутков
	// графики сплайнов и производных
	
	puts("x: \t\t| f(x):\t\t| f'(x):");
	for (int i = 0; i < node_count; ++i) 
		printf("%.9f\t| %.9f\t| %.9f\n", grid[i], y[i], y_shtrix[i]);
	
	puts("\ns2(x) = {");
	for (int i = 0; i < quadrareas; ++i) {
		double x_k = (grid[i]+grid[i+1])/2;
		printf("          { %10.5f*(x-%.2f)^2 + %9.5f*(x-%.2f) + %9.5f, for x in [%.1f, %.1f)\n", 
			   s2c[i], x_k, s2b[i], x_k, s2a[i], grid[i], grid[i+1]);

	}
	puts("          }");

	puts("Coefficienti metoda podgonki:");
	puts("C_k = alpha_k*C_{k+1} + beta_k");

	puts("i:\t|alpha: \t| beta:");
	for (int i = 0; i < cubareas; ++i) 
		printf("%d\t|%.9f\t| %.9f\n", i, alpha[i], beta[i]);

	puts("\ns3(x) = {");
	for (int i = 0; i < cubareas; ++i) {
		double x_k = grid[i]; 
		printf("          { %10.5f*(x-%.1f)^3 + %10.5f*(x-%.1f)^2 + %9.5f*(x-%.1f) + %9.5f, for x in [%.1f, %.1f)\n", 
			   s3d[i], x_k, s3c[i], x_k, s3b[i], x_k, s3a[i], grid[i], grid[i+1]);
	}
	puts("    }");

	puts("");

	puts("\nh3(x) = {");
	for (int i = 0; i < ermitarea; ++i) {
		double x_k = grid[i]; 
		printf("          { %10.5f*(x-%.1f)^3 + %10.5f*(x-%.1f)^2 + %9.5f*(x-%.1f) + %9.5f, for x in [%.1f, %.1f)\n", 
			   h3a[i], x_k, h3b[i], x_k, h3c[i], x_k, h3d[i], grid[i], grid[i+1]);
	}
	puts("    }");


	// График f(x)
	FILE *gp_f = popen("gnuplot -persistent", "w");
	fprintf(gp_f, "set encoding utf8\n");
	fprintf(gp_f, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
	fprintf(gp_f, "set output 'f(x).png'\n");
	fprintf(gp_f, "set title 'f(x)'\n");
	fprintf(gp_f, "set xlabel 'x'\n");
	fprintf(gp_f, "set ylabel 'y'\n");
	fprintf(gp_f, "set grid\n");
	fprintf(gp_f, "plot '-' with lines linewidth 2 title 'f(x)'\n");

	for (double x = 0; x <= 2; x += 0.01) {
		fprintf(gp_f, "%f %f\n", x, f(x));
	}
	fprintf(gp_f, "e\n");
	fflush(gp_f);
	pclose(gp_f);

	// График s2(x)
	FILE *gp_s2 = popen("gnuplot -persistent", "w");
	fprintf(gp_s2, "set encoding utf8\n");
	fprintf(gp_s2, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
	fprintf(gp_s2, "set output 's2(x).png'\n");
	fprintf(gp_s2, "set title 's2(x) - квадратичный сплайн'\n");
	fprintf(gp_s2, "set xlabel 'x'\n");
	fprintf(gp_s2, "set ylabel 'y'\n");
	fprintf(gp_s2, "set grid\n");
	fprintf(gp_s2, "plot '-' with lines linewidth 2 title 's2(x)'\n");

	for (double x = 0; x <= 2; x += 0.01) {
		fprintf(gp_s2, "%f %f\n", x, apply_s2(quadrareas, s2a, s2b, s2c, a, h, x));
	}
	fprintf(gp_s2, "e\n");
	fflush(gp_s2);
	pclose(gp_s2);

	// График s3(x)
	FILE *gp_s3 = popen("gnuplot -persistent", "w");
	fprintf(gp_s3, "set encoding utf8\n");
	fprintf(gp_s3, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
	fprintf(gp_s3, "set output 's3(x).png'\n");
	fprintf(gp_s3, "set title 's3(x) - кубический сплайн'\n");
	fprintf(gp_s3, "set xlabel 'x'\n");
	fprintf(gp_s3, "set ylabel 'y'\n");
	fprintf(gp_s3, "set grid\n");
	fprintf(gp_s3, "plot '-' with lines linewidth 2 title 's3(x)'\n");

	for (double x = 0; x <= 2; x += 0.01) {
		fprintf(gp_s3, "%f %f\n", x, apply_s3(cubareas, s3a, s3b, s3c, s3d, a, h, x));
	}
	fprintf(gp_s3, "e\n");
	fflush(gp_s3);
	pclose(gp_s3);

	// График h3(x)
	FILE *gp_h3 = popen("gnuplot -persistent", "w");
	fprintf(gp_h3, "set encoding utf8\n");
	fprintf(gp_h3, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
	fprintf(gp_h3, "set output 'h3(x).png'\n");
	fprintf(gp_h3, "set title 'h3(x) - эрмитов сплайн'\n");
	fprintf(gp_h3, "set xlabel 'x'\n");
	fprintf(gp_h3, "set ylabel 'y'\n");
	fprintf(gp_h3, "set grid\n");
	fprintf(gp_h3, "plot '-' with lines linewidth 2 title 'h3(x)'\n");

	for (double x = 0; x <= 2; x += 0.01) {
		fprintf(gp_h3, "%f %f\n", x, apply_ermit(ermitarea, h3a, h3b, h3c, h3d, a, h, x));
	}
	fprintf(gp_h3, "e\n");
	fflush(gp_h3);
	pclose(gp_h3);

// Производная f(x)
FILE *gp_f_deriv = popen("gnuplot -persistent", "w");
fprintf(gp_f_deriv, "set encoding utf8\n");
fprintf(gp_f_deriv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_f_deriv, "set output 'f_prime(x).png'\n");
fprintf(gp_f_deriv, "set title \"f'(x)\"\n");
fprintf(gp_f_deriv, "set xlabel 'x'\n");
fprintf(gp_f_deriv, "set ylabel 'y'\n");
fprintf(gp_f_deriv, "set grid\n");
fprintf(gp_f_deriv, "plot '-' with lines linewidth 2 title \"f'(x)\"\n");

for (double x = 0; x <= 2; x += 0.01) {
    fprintf(gp_f_deriv, "%f %f\n", x, manual_derivative(f, x));
}
fprintf(gp_f_deriv, "e\n");
fflush(gp_f_deriv);
pclose(gp_f_deriv);

// Производная s2(x)
FILE *gp_s2_deriv = popen("gnuplot -persistent", "w");
fprintf(gp_s2_deriv, "set encoding utf8\n");
fprintf(gp_s2_deriv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_s2_deriv, "set output 's2_prime(x).png'\n");
fprintf(gp_s2_deriv, "set title \"s2'(x) - производная квадратичного сплайна\"\n");
fprintf(gp_s2_deriv, "set xlabel 'x'\n");
fprintf(gp_s2_deriv, "set ylabel 'y'\n");
fprintf(gp_s2_deriv, "set grid\n");
fprintf(gp_s2_deriv, "plot '-' with lines linewidth 2 title \"s2'(x)\"\n");

for (double x = 0; x <= 2; x += 0.01) {
    fprintf(gp_s2_deriv, "%f %f\n", x, derivative_s2(quadrareas, s2a, s2b, s2c, a, h, x));
}
fprintf(gp_s2_deriv, "e\n");
fflush(gp_s2_deriv);
pclose(gp_s2_deriv);

// Производная s3(x)
FILE *gp_s3_deriv = popen("gnuplot -persistent", "w");
fprintf(gp_s3_deriv, "set encoding utf8\n");
fprintf(gp_s3_deriv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_s3_deriv, "set output 's3_prime(x).png'\n");
fprintf(gp_s3_deriv, "set title \"s3'(x) - производная кубического сплайна\"\n");
fprintf(gp_s3_deriv, "set xlabel 'x'\n");
fprintf(gp_s3_deriv, "set ylabel 'y'\n");
fprintf(gp_s3_deriv, "set grid\n");
fprintf(gp_s3_deriv, "plot '-' with lines linewidth 2 title \"s3'(x)\"\n");

for (double x = 0; x <= 2; x += 0.01) {
    fprintf(gp_s3_deriv, "%f %f\n", x, derivative_s3(cubareas, s3a, s3b, s3c, s3d, a, h, x));
}
fprintf(gp_s3_deriv, "e\n");
fflush(gp_s3_deriv);
pclose(gp_s3_deriv);

// Производная h3(x)
FILE *gp_h3_deriv = popen("gnuplot -persistent", "w");
fprintf(gp_h3_deriv, "set encoding utf8\n");
fprintf(gp_h3_deriv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_h3_deriv, "set output 'h3_prime(x).png'\n");
fprintf(gp_h3_deriv, "set title \"h3'(x) - производная эрмитова сплайна\"\n");
fprintf(gp_h3_deriv, "set xlabel 'x'\n");
fprintf(gp_h3_deriv, "set ylabel 'y'\n");
fprintf(gp_h3_deriv, "set grid\n");
fprintf(gp_h3_deriv, "plot '-' with lines linewidth 2 title \"h3'(x)\"\n");

for (double x = 0; x <= 2; x += 0.01) {
    fprintf(gp_h3_deriv, "%f %f\n", x, derivative_h3(ermitarea, h3a, h3b, h3c, h3d, a, h, x));
}
fprintf(gp_h3_deriv, "e\n");
fflush(gp_h3_deriv);
pclose(gp_h3_deriv);


// Вторая производная f(x)
FILE *gp_f_dderiv = popen("gnuplot -persistent", "w");
fprintf(gp_f_dderiv, "set encoding utf8\n");
fprintf(gp_f_dderiv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_f_dderiv, "set output 'f_doubleprime(x).png'\n");
fprintf(gp_f_dderiv, "set title \"f''(x)\"\n");
fprintf(gp_f_dderiv, "set xlabel 'x'\n");
fprintf(gp_f_dderiv, "set ylabel 'y'\n");
fprintf(gp_f_dderiv, "set grid\n");
fprintf(gp_f_dderiv, "plot '-' with lines linewidth 2 title \"f''(x)\"\n");

const double eps = 1e-9;

for (double x = eps; x <= 2 - eps; x += 0.01) {
    fprintf(gp_f_dderiv, "%f %f\n", x, manual_doublederivative(f, x));
}

fprintf(gp_f_dderiv, "e\n");
fflush(gp_f_dderiv);
pclose(gp_f_dderiv);

// Вторая производная s2(x)
FILE *gp_s2_dderiv = popen("gnuplot -persistent", "w");
fprintf(gp_s2_dderiv, "set encoding utf8\n");
fprintf(gp_s2_dderiv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_s2_dderiv, "set output 's2_doubleprime(x).png'\n");
fprintf(gp_s2_dderiv, "set title \"s2''(x) - вторая производная квадратичного сплайна\"\n");
fprintf(gp_s2_dderiv, "set xlabel 'x'\n");
fprintf(gp_s2_dderiv, "set ylabel 'y'\n");
fprintf(gp_s2_dderiv, "set grid\n");
fprintf(gp_s2_dderiv, "plot '-' with lines linewidth 2 title \"s2''(x)\"\n");

for (double x = eps; x <= 2 - eps; x += 0.01) {
    fprintf(gp_s2_dderiv, "%f %f\n", x, doublederivative_s2(quadrareas, s2a, s2b, s2c, a, h, x));
}

fprintf(gp_s2_dderiv, "e\n");
fflush(gp_s2_dderiv);
pclose(gp_s2_dderiv);

// Вторая производная s3(x)
FILE *gp_s3_dderiv = popen("gnuplot -persistent", "w");
fprintf(gp_s3_dderiv, "set encoding utf8\n");
fprintf(gp_s3_dderiv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_s3_dderiv, "set output 's3_doubleprime(x).png'\n");
fprintf(gp_s3_dderiv, "set title \"s3''(x) - вторая производная кубического сплайна\"\n");
fprintf(gp_s3_dderiv, "set xlabel 'x'\n");
fprintf(gp_s3_dderiv, "set ylabel 'y'\n");
fprintf(gp_s3_dderiv, "set grid\n");
fprintf(gp_s3_dderiv, "plot '-' with lines linewidth 2 title \"s3''(x)\"\n");

for (double x = eps; x <= 2 - eps; x += 0.01) {
    fprintf(gp_s3_dderiv, "%f %f\n", x, doublederivative_s3(cubareas, s3a, s3b, s3c, s3d, a, h, x));
}

fprintf(gp_s3_dderiv, "e\n");
fflush(gp_s3_dderiv);
pclose(gp_s3_dderiv);

// Вторая производная h3(x)
FILE *gp_h3_dderiv = popen("gnuplot -persistent", "w");
fprintf(gp_h3_dderiv, "set encoding utf8\n");
fprintf(gp_h3_dderiv, "set terminal pngcairo size 1024,768 enhanced font 'Times New Roman, 14'\n");
fprintf(gp_h3_dderiv, "set output 'h3_doubleprime(x).png'\n");
fprintf(gp_h3_dderiv, "set title \"h3''(x) - вторая производная эрмитова сплайна\"\n");
fprintf(gp_h3_dderiv, "set xlabel 'x'\n");
fprintf(gp_h3_dderiv, "set ylabel 'y'\n");
fprintf(gp_h3_dderiv, "set grid\n");
fprintf(gp_h3_dderiv, "plot '-' with lines linewidth 2 title \"h3''(x)\"\n");

for (double x = eps; x <= 2 - eps; x += 0.01) {
    fprintf(gp_h3_dderiv, "%f %f\n", x, doublederivative_h3(ermitarea, h3a, h3b, h3c, h3d, a, h, x));
}

fprintf(gp_h3_dderiv, "e\n");
fflush(gp_h3_dderiv);
pclose(gp_h3_dderiv);

	return 0;
}
