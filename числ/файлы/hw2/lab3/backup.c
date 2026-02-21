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


constexpr int node_count = 10;
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

double x[node_count];

static double grid_fun(const double arg){
	assert(arg >= a && arg <= b);
	const double position = (arg - a)/h; 
	int index = (int) floor(position);
//	assert(index >= 0 && index < node_count);
	if (index == node_count - 1)
		return x[node_count - 1];
	
	const double interop = position - index;
	assert(0 <= interop < 1);
	
	auto xi = x[index];
	auto xii = x[index+1];
	return xii*interop +
		   xi*(1-interop);
}

double analytic(){
	// x = alpha*(s^2 + beta) + lambda int_a^b K(s, tau) x(tau) dtau
	//						  + lambda( int_a^s (2b*s-tau) x(tau) dtau 
	//								  + int_s^b (2b*tau-s) x(tau) dtau)
	// x'_s = 2alpha*s + lambda( int_a^s 2b x(tau) dtau + (2b*s-s) x(s)
	//							+int_s^b -  x(tau) dtau - (2b*s-s) x(s))
	//
	// x'_s = 2alpha*s + lambda( int_a^s 2b x(tau) dtau + #############
	//							+int_s^b -  x(tau) dtau - #############)
	//
	// x'_s = 2alpha*s + lambda( 2b int_a^s x(tau) dtau - int_s^b x(tau) dtau)
	// x'' = 2alpha + lambda( 2b x(s) + x(s) )
	// x'' = 2alpha + lambda(2b + 1)x
	// x'' - lambda(2b+1) x = 2alpha
	printf("%lf, %lf\n", lambda, (double)(2*b+1));
	assert(lambda > 0 && (2*b+1) > 0);
	// x = C1 ch(omega s) + C2 sh(omega s) + X0
	auto X0 = -2*alpha/(lambda*(2*b+1));
	auto omega = sqrt(lambda/(2*b+1));
	printf("X0: %lf, omega: %lf\n", X0, omega);
	// C1 ch(omega s) + C2 sh(omega s) + X0 - lambda( int_a^s (2b*s - tau) x(tau) dtau
	//												+ int_s^b (2b*tau - s) x(tau) dtau) = y(s)
	// s = a, b:
	//
	// C1 ch(omega a) + C2 sh(omega a) + X0 - lambda(int_a^b (2b*tau - a) x(tau) dtau) = y(a)
	//			int (2b*tau - a) (C1 ch(omega tau) + C2 sh(omega tau) + X0) dtau
	//			int 2b*tau*C1 ch(omega tau) + 
	//				    -a*C1 ch(omega tau) +
	//				2b*tau*C2 sh(omega tau) +
	//					-a*C2 sh(omega tau) +
	//				2b*tau*X0 +
	//					-a*X0
	auto D0 = -a*X0;  // 1
	auto D1 = 2*b*X0; // x

	auto D2 = -a;  // C2 sh
	auto D3 = -a;  // C1 ch
	auto D4 = 2*b; // C2 x sh
	auto D5 = 2*b; // C1 x ch 
	// 1 x x^2 sh ch x*sh x*ch
	// x
	// x^2
	// sh
	// ch
	// x*sh
	// x*ch
	// ( 0 1 0 0 0 0 0) (0 )   (D0)
	// ( 0 0 2 0 0 0 0) (A0)   (D1)
	// ( 0 0 0 0 0 0 0) (A1)   (0 )
	// ( 0 0 0 0 o 1 0) (A2)   (D2 * C2)
	// ( 0 0 0 o 0 0 1) (A3) = (D3 * C1)
	// ( 0 0 0 0 0 0 o) (A4)   (D4 * C2)
	// ( 0 0 0 0 0 o 0) (A5)   (D5 * C1)
	//	
	//	A0 = D0
	//	2*A1 = D1
	//	omega*A3 + A4 = D2 * C2
	//	omega*A2 + A5 = D3 * C1
	//	omega*A5      = D4 * C2
	//	omega*A4      = D5 * C1

	auto A0 = D0;
	auto A1 = D1/2.0;
	auto A4 = D5/omega; // D5/o * C1
	auto A5 = D4/omega; // D4/o * C2
	
	// o*A3 + D5/o*C1 = D2 * C2
	// o*A2 + D4/o*C2 = D3 * C1
	//
	// o*A3 = D2 * C2 - A4 * C1
	// o*A2 = D3 * C1 - A5 * C2
	//
	// A3 = 1/o (-A4 * C1 + D2 * C2)
	// A2 = 1/o ( D3 * C1 - A5 * C2)
	auto A31 = -A4/omega;
	auto A32 = D2/omega;
	auto A210 = D3/omega;
	auto A220 = -A5/omega;

	// A0 *x + A1 * x^2 + A2 * sinh + A3 * cosh + A4 * x sinh + A5 * x cosh
	//  ok       ok       A21*C1+A22*C2 sinh
	//							+A31*C1+A32*C2 cosh
	//													C1		C2
	auto sa = sinh(omega*a);
	auto ca = cosh(omega*a);
	auto sb = sinh(omega*b);
	auto cb = cosh(omega*b);
	auto asa = a*sa;
	auto aca = a*ca;
	auto bsb = b*sb;
	auto bcb = b*cb;
	// C1 ch(omega a) + C2 sh(omega a) + X0 - lambda(int_a^b (2b*tau - a) x(tau) dtau) = y(a)
	//			int (2b*tau - a) (C1 ch(omega tau) + C2 sh(omega tau) + X0) dtau
	//
	// ca*C1 + sa*C2 + X0 - lambda*(A0*b + A1*b^2 + sb(A21*C1 + A22*C2) + cb(A31*C1 + A32*C2) + A4 * bsb * C1 + A5 * bcb * C2 -
	//						(A0*a + A1*a^2 + sa(A21*C1 + A22*C2) + ca(A31*C1 + A32*C2) + A4 * asa * C1 + A5 * aca * C2))=y(a)
	//	1*(X0 - lambda*(A0*b + A1*b^2 - A0*a - A1*a^2) )+
	//	C1*(ca - lambda*(sb*A21+cb*A31+A4*bsb - sa*A21 - ca*A31 - A4*asa))+
	//	C2*(sa - lambda*(sb*A22+cb*A32+A5*bcb - sa*A22 - ca*A32 - A5*aca))
	auto v1 = y(a) - (X0 - lambda*(A0*b + A1*b*b - A0*a - A1*a*a));
	auto M11 = ca - lambda*(sb*A210+cb*A31+A4*bsb - sa*A210 - ca*A31 - A4*asa);
	auto M12 = sa - lambda*(sb*A220+cb*A32+A5*bcb - sa*A220 - ca*A32 - A5*aca);
	// C1 ch(omega s) + C2 sh(omega s) + X0 - lambda( int_a^s (2b*s - tau) x(tau) dtau
	//												+ int_s^b (2b*tau - s) x(tau) dtau) = y(s)
	// C1 ch(omega b) + C2 sh(omega b) + X0 - lambda( int_a^b (2b*b - tau) (
	//		C1 ch(omega tau) + C2 sh(omega tau) + X0
	// ) dtau )= y(b)
	//		int (2b*b - tau) (C1 ch(omega tau) + C2 sh(omega tau) + X0)
	//	1, x, x^2, sh, ch, x*sh, x*ch
	auto D20 = 2*b*b*X0;
	auto D21 = -1*X0;
	auto D22 = 2*b*b; // C2
	auto D23 = 2*b*b; // C1
	auto D24 = -1; // C2
	auto D25 = -1; // C1

	auto A20 = D20;
	auto A21 = D21/2.0;
	auto A24 = D25/omega; // D5/o * C1
	auto A25 = D24/omega; // D4/o * C2
	
	auto A231 = -A24/omega;
	auto A232 = D22/omega;
	auto A221 = D23/omega;
	auto A222 = -A25/omega;

	auto v2 = y(b) - (X0 - lambda*(A20*b + A21*b*b - A20*a - A21*a*a));
	auto M21 = cb - lambda*(sb*A221+cb*A231+A24*bsb - sa*A221 - ca*A231 - A24*asa);
	auto M22 = sb - lambda*(sb*A222+cb*A232+A25*bcb - sa*A222 - ca*A232 - A25*aca);

	// M11 C1 + M12 C2 = v1
	// M21 C1 + M22 C2 = v2
	
	auto determinant = M11*M22 - M12*M21;
	auto d1 = v1 * M22 - M12 * v2;
	auto d2 = M11 * v2 - v1 * M21;

	auto C1 = d1/determinant;
	auto C2 = d2/determinant;

	printf("x(s) = %lf + %lf*ch(%lf*s) + %lf*sh(%lf*s)\n",
			X0, C1, omega, C2, omega);

}

int main(){
	for (int i = 0; i < node_count; ++i)
		x[i] = a+i*h;
	printf("a: %lf, b: %lf, n: %d, N: %d, lambda: %lf\n",
				 a,      b,     n,     N,      lambda);
	printf("y(s) = %lf(s^2 + %lf)\n",
				 alpha,     beta);
	PRINT_ARRAY(node_count, x);
	
	double y[node_count];
	double F[node_count*node_count];

	for (int i = 0; i < node_count; ++i)
		y[i] = yi(i);

	PRINT_ARRAY(node_count, y);
	for (int i = 0; i < node_count; ++i)
		for (int j = 0; j < node_count; ++j){
			const double delta = (i == j)? 1: 0;
			F[i*node_count +j] = delta - lambda * Kij(i, j) * h;
		}

	PRINT_MATRIX(node_count, F);
	SLAE_solve(node_count, F, x, y);
	PRINT_MATRIX(node_count, F);
	PRINT_ARRAY(node_count, x);

	analytic();

	return 0;
}
