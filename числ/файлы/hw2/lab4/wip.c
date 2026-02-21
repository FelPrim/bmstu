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

typedef double (* const math_function)(double);
constexpr double epsilon = 1e-10;

static void print_array(const int size, const double array[static size]){
	printf("[");
	int i = 0;
	assert(size > 1);
	for (; i < size-1; ++i)
		printf("%lf, ", array[i]);
	printf("%lf]\n", array[size-1]);
}

#define PRINT_ARRAY(x, y) printf("%s: ", #y); print_array(x, y)


// math part
constexpr int N = 13;
constexpr int n = 10;
constexpr int a = 5;
constexpr int b = 1;
constexpr int c = 3;

constexpr double A = 0;
constexpr double B = 2;
constexpr double h = (B-A)/((double) n);

static double uniform_grid(int index){ return A+h*index;}

static double f(const double x){
	const double x_2 = x*x;
	const double x_4 = x_2 * x_2;
	return ((a + 53 - N) * x_4 + (b - 53 + N) * x_2 + c) /
					((x + 1) * (x_2 + 1));
}

static double F(const double x){
	const double x_2 = x*x;
	const double ln_xp1 = log(x + 1);
	const double ln_x2p1 = log(x_2 + 1);
	const double arctg_x = atan(x);

	return 45*x_2/2.0 - 45*x - 39*ln_xp1 + 87/2.0*(
		-1/2.0*ln_x2p1 + arctg_x + ln_xp1
	);
}

double rect_approx(math_function f){
	double sum = 0;
	for (int i = 1; i <= n; ++i){
		const double x_prev = uniform_grid(i-1);
		const double x_cur  = uniform_grid(i);

		const double theta = (x_prev + x_cur)/2.0;

		sum += f(theta);
	}
	return h*sum;
}

double trap_approx(math_function f){
	double sum = 0;
	for (int i = 1; i < n; ++i){
		const double x_i = uniform_grid(i);
		sum += f(x_i);
	}
	const double result = f(uniform_grid(0)) + f(uniform_grid(n)) + 2*sum;
	return h/2*result;
}

double simp_approx(math_function f){
	constexpr int m = n/2;

	double sum_1 = 0, sum_2 = 0;
	for (int i = 1; i <= m; ++i)
		sum_1 += f(uniform_grid(2*i-1));
	for (int i = 1; i < m; ++i)
		sum_2 += f(uniform_grid(2*i));
	const double result = f(uniform_grid(0)) + f(uniform_grid(n)) + 4*sum_1 + 2*sum_2;
	return h/3 * result;
}

double df(double x){
	return (f(x+epsilon)-f(x-epsilon))/(2*epsilon);
}

double dfminus(double x){
	return -df(x);
}

double d2f(double x){
	constexpr double epsilo = sqrt(epsilon);

	return (f(x+epsilo)+f(x-epsilo) - 2*f(x))/epsilon;
}

double d2fminus(double x){
	return -d2f(x);
}

double d4f(double x){
	constexpr double epsilo = sqrt(epsilon);
	constexpr double epsil = sqrt(epsilo);

	return (  f(x-2*epsil)
			- 4*f(x-epsil)
			+ 6*f(x)
			- 4*f(x+epsil)
			+ f(x+2*epsil))/epsilon;
}

double d4fminus(double x){
	return -d4f(x);
}


double find_zero(math_function f){
	constexpr double phi = 0.6180339887498948482;

	double a = A, b = B;

	double x1 = B - phi * (B - A);
	double x2 = A + phi * (B - A);

	double f1 = f(x1);
	double f2 = f(x2);

	constexpr double tol = 1e-8;

	while (fabs(b - a) > tol){
		if (f1 < f2){
			a = x1;
			
			x1 = x2;
			f1 = f2;
			
			x2 = a + phi * (b - a);
			f2 = f(x2);
		} else {
			b = x2;

			x2 = x1;
			f2 = f1;

			x1 = b - phi * (b - a);
			f1 = f(x1);
		}
	}
	return (a+b)/2;
}

double rect_error(){
	const double x1 = find_zero(df);
	const double x2 = find_zero(dfminus);

	const double f1_ = fabs(df(x1));
	const double f2_ = fabs(df(x2));
	
	double M1;
	double x_;
	if (f1_ > f2_){
		x_ = x1;
		M1 = f1_;
	} else {
		x_ = x2;
		M1 = f2_;
	}

	printf("x: |f'(x)| = max: %lf, M1: %lf\n",
							   x_,	    M1);

	return (B - A)/2.0 * h * M1;
}


double trap_error(){
	const double x1 = find_zero(d2f);
	const double x2 = find_zero(d2fminus);

	const double f1__ = fabs(d2f(x1));
	const double f2__ = fabs(d2f(x2));
	
	double M2;
	double x_;
	if (f1__ > f2__){
		x_ = x1;
		M2 = f1__;
	} else {
		x_ = x2;
		M2 = f2__;
	}

	printf("x: |f''(x)| = max: %lf, M2: %lf\n",
							   x_,	    M2);

	return (B - A)/12.0 * h*h * M2;
}


double simp_error(){
	const double x1 = find_zero(d4f);
	const double x2 = find_zero(d4fminus);

	const double f1____ = fabs(d4f(x1));
	const double f2____ = fabs(d4f(x2));
	
	double M4;
	double x_;
	if (f1____ > f2____){
		x_ = x1;
		M4 = f1____;
	} else {
		x_ = x2;
		M4 = f2____;
	}

	printf("x: |f''''(x)| = max: %lf, M4: %lf\n",
							      x_,	   M4);

	return (B - A)/180.0 * h*h*h*h * M4;
}

int main(){
	
	printf(" f(x) = %d x^4 + %d x^2 + %d\n",
			   a+53-N,  b-53+N,       c);
	puts("______________________________");
	puts("    (x + 1)(x^2 + 1)");
	const double I_rect = rect_approx(f), 
				 I_trap = trap_approx(f), 
				 I_simp = simp_approx(f),
				 I_anal = F(B) - F(A);
	printf("Rectangles: %lf\nTrapetions: %lf\nSimpson: %lf\n",
					 I_rect,		  I_trap,		I_simp);
	printf("Analytic: %lf\n", I_anal);
	const double error_rect = rect_error(),
				 error_trap = trap_error(),
				 error_simp = simp_error();
	printf("Apriori rectangles error: %lf\n", error_rect);
	printf("Apriori trapetions error: %lf\n", error_trap);
	printf("Apriori Simpson error: %lf\n", error_simp);
	return 0;
}
