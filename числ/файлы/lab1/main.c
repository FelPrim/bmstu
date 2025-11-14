#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

/* Теория:
 * Ax = b
 * a^i_j x^j = b^i
 * ||b|| = max{|b^i|}
 * ||A|| = max{|a^i_j|e_i}
 * Определение 1.1 cond(A) = ||A|| * ||A^-1|| - число обусловленности
 * A(x + Dx) = b + Db
 * A(Dx) = Db
 * Теорема 1.1 ||Dx||/||x|| <= cond(A)* ||Db||/||b||
 * cond(A) >= 1
 * cond(A) ~ 10 - хорошо обусловленная
 * cond(A) > 10^2 - плохо обусловленная
 * b = A x => ||b|| <= ||A||*||x|| => ||x|| >= ||b||/||A||
 * Dx = A Db => ||Dx|| <= ||A^-1||*||Db||
 * Замечание об овражности симметричной положительно определённой матрицы
 * A^T = A=>
 * ||A||=sqrt(rho(A^T*A)) = sqrt(rho(A^2)) = sqrt(lambda_max^2(A)) = lambda_max(A)
 * max(Spr(A^-1)) = 1/(lambda_min (A)), (A^-1)^T=A^-1
 * ||A^-1|| = sqrt(rho_Spr(A^T*A))=sqrt(rho_Spr(A^-2))=1/(lambda_min (A))
 * cond_e(A) = lambda_max(A)/lambda_min(A) 
 * 
 * Задание 1.1.
 * N = 13
 * alpha = (53-50)/100 = 0.03
 * { 200(1+0.5*N+alpha)*x^1 + 200*(1+0.5*N)        *x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)
 * { 200.1*(1+0.5*N)   *x^1 + 199.9*(1+0.5*N+alpha)*x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)
 * { 199.9*(1+0.5*N)   *x^1 + 200*(1+0.5*N)        *x^2 + 200.1*(1+0.5*N+alpha)*x^3=200*(3+1.5*N+alpha)
 *
 * относительная ошибка правой части 0.01
 * Приближённая СЛАУ:
 * { 200(1+0.5*N+alpha)*x^1 + 200*(1+0.5*N)        *x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)*(1+0.01)
 * { 200.1*(1+0.5*N)   *x^1 + 199.9*(1+0.5*N+alpha)*x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)*(1-0.01)
 * { 199.9*(1+0.5*N)   *x^1 + 200*(1+0.5*N)        *x^2 + 200.1*(1+0.5*N+alpha)*x^3=200*(3+1.5*N+alpha)*(1+0.01)
 *
 * Найти:
 * число обусловленности матрицы рассматриваемой СЛАУ и относительную погрешность приближённой СЛАУ
 * Прокомментриовать получившиеся результаты
 * Задание 1.2
 * N = 13, alpha = 0.02
 * lambda+alpha=0.4
 * F = arctg
 * a = 0
 * b = 1
 * На отрезке [a, b] выбрана центрально равномерная сетка с десятью узлами:
 * s1=t1=a+h/2,s2=t2=t1+h, s3=t3=t2+h, ..., s10=t10=t9+h
 * h = (b-a)/10
 * Решить приближённую СЛАУ
 * (E+lambda A)x=b+Db (5)
 * A=a^i_j, b=b^i:
 * a^i_j = F(s_i*t_j)*(b-a)/10
 * b = (E+lambda A)*x
 * x^i = 1
 * приближённая СЛАУ определяется только погрешностью Db=0.01[b^1, -b^2, b^3, -b^4, b^5, -b^6, b^7, -b^8, b^9, -b^10]
 * Найти:
 * число обусловленности матрицы рассматриваемой СЛАУ и относительную погрешность в решении приближённой СЛАУ
 * Затем, прокомментировать получившиеся результаты.
 * Кроме того, найти решение СЛАУ, которая получается из СЛАУ (5) делением каждого i-го уравнения на число b^i+Db^i
 * После этого сравнить абсолютную погрешность в решении получившейся СЛАУ с абсолютной погрешностью в решении приближённой СЛАУ
 * */
typedef struct Matrix{
	alignas(32) double * restrict data;
	const unsigned short rows;
	const unsigned short columns;
} Matrix;

typedef struct Vector{
	alignas(32) double * restrict data;
	const unsigned short length;
} Vector;

void Multiply(const double lambda, const Matrix A, Matrix result){
	for (unsigned int i = 0; i < (unsigned int)(A.rows)*A.columns; ++i)
		result.data[i] = lambda*A.data[i];
}
void Sum(const Matrix E, const Matrix la, Matrix result){
	for (unsigned int i = 0; i < (unsigned int)(E.rows)*E.columns; ++i)
		result.data[i] = la.data[i] + E.data[i];
}
void Apply(const Matrix A, const Vector x, Vector b){
	for (unsigned short i = 0; i < A.rows; ++i){
		b.data[i] = 0;
		for (unsigned short j = 0; j < A.columns; ++j){
			b.data[i] += x.data[j] * A.data[j + A.columns*i];
		}
	}
}

double Norm(const Matrix matrix){
	double result = -1;
	for (unsigned int i = 0; i < matrix.rows; ++i){
		double row_result = 0;
		for (unsigned int j = 0; j < matrix.columns; ++j){
			row_result += fabs(matrix.data[j+matrix.columns*i]);
		}
		if (row_result > result)
			result = row_result;
	}
	return result;
}

double norm(const Vector vector){
	double result = 0;
	for (unsigned short i = 0; i < vector.length; ++i){
		result += vector.data[i] * vector.data[i];
	}
	return sqrt(result);
}

typedef struct MinorStruct{
	Matrix matrix;
	const unsigned short row;
	const unsigned short column;
} MinorStruct;

inline static double xifaelsey(double x, double y, int a, int b){
	int equal = a == b;
	return equal*x+(1-equal)*y;
}

double minor(const MinorStruct input);
double minor(const MinorStruct input){

	switch (input.matrix.rows){
		case 0:
		case 1:
			return 0;
		case 2:
			return input.matrix.data[3 - input.column - 2*input.row];
		case 3:
			int i0 = input.row != 0;
			int i1 = input.row != 1;
			int i2 = input.row != 2;
			int j0 = input.column != 0;
			int j1 = input.column != 1;
			int j2 = input.column != 2;
			int a11 = 4 - j0 - 3*i0;
			int a12 = 4 + j2 - 3*i0;
			int a21 = 4 - j0 + 3*i2;
			int a22 = 4 + j2 + 3*i2;
			return input.matrix.data[a11] * input.matrix.data[a22]
			     - input.matrix.data[a12] * input.matrix.data[a21];
		default:
			// TODO
			exit(-1);
	}	
}

double determinant(const Matrix matrix){
	assert(matrix.rows < 4);
	switch (matrix.rows){
		case 1:
			return matrix.data[0];
		case 2:
			return 	  matrix.data[0]*matrix.data[3]
				- matrix.data[1]*matrix.data[2];
		case 3:
			return 	matrix.data[0]*matrix.data[4]*matrix.data[8]
				  +	matrix.data[1]*matrix.data[5]*matrix.data[6]
				  +	matrix.data[2]*matrix.data[3]*matrix.data[7]
				  -	matrix.data[2]*matrix.data[4]*matrix.data[6]
				  -	matrix.data[1]*matrix.data[3]*matrix.data[8]
				  -	matrix.data[0]*matrix.data[5]*matrix.data[7];
		default:
			//TODO
			//exit(-1);
	}
	
}

double cond(const Matrix matrix, const Matrix inversed){
	return Norm(matrix)*Norm(inversed);
}

const double EPSILON = 1e-11;
void check_inverse(const Matrix original, const Matrix inversed){
	for (int i = 0; i < original.rows; ++i){
		for (int j = 0; j < inversed.columns; ++j){
			double x = 0;
			for (int k = 0; k < original.columns; ++k){
				// original[i][k] * inversed*[k][j]
				x +=   original.data[k+original.columns*i]
					 * inversed.data[j+inversed.columns*k];
			}

			//printf("i: %d, j: %d, x: %.15lf\n", i, j, x);
			if (i == j){
				assert(fabs(x-1) < EPSILON);
			}
			else
				assert(fabs(x) < EPSILON);
		}
	}
}

void inverse(const Matrix original, Matrix result){
	assert(original.rows < 4);
	double det = determinant(original);
	printf("det: %llf\n", det);
	assert(fabs(det) > EPSILON);
	int parity = 1;
	for (unsigned short i = 0; i < result.rows; ++i){
		for (unsigned short j = 0; j < result.columns; ++j){	
			// original[j + matrix->columns*sizeof(double)*i] - original[i][j]
			// (j, i) вместо (i, j) т.к. речь идёт о транспонированной матрице

			struct MinorStruct orig = {
				.matrix = original,
				.row = j,
				.column = i
			};
			double cofactor = parity * minor(orig);
			result.data[j + result.columns*i] = cofactor/det;
			parity = -parity;
			/*
			for (unsigned short i2 = 0; i2 < original->columns; ++i2){
				if (i2 == j)
					continue;
				for (unsigned short j2 = 0; j2 < original->rows; ++j2){
					if (j2 == i)
						continue;
					det_before_sum[i2] += 
				}
			}
			*/
			//double cofactor = //-1^i+j * minor
		}
	}
	check_inverse(original, result);
}



/*
 * 
 * Задание 1.1.
 * N = 13
 * alpha = (53-50)/100 = 0.03
 * { 200(1+0.5*N+alpha)*x^1 + 200*(1+0.5*N)        *x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)
 * { 200.1*(1+0.5*N)   *x^1 + 199.9*(1+0.5*N+alpha)*x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)
 * { 199.9*(1+0.5*N)   *x^1 + 200*(1+0.5*N)        *x^2 + 200.1*(1+0.5*N+alpha)*x^3=200*(3+1.5*N+alpha)
 *
 * относительная ошибка правой части 0.01
 * Приближённая СЛАУ:
 * { 200(1+0.5*N+alpha)*x^1 + 200*(1+0.5*N)        *x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)*(1+0.01)
 * { 200.1*(1+0.5*N)   *x^1 + 199.9*(1+0.5*N+alpha)*x^2 + 200(1+0.5*N)         *x^3=200*(3+1.5*N+alpha)*(1-0.01)
 * { 199.9*(1+0.5*N)   *x^1 + 200*(1+0.5*N)        *x^2 + 200.1*(1+0.5*N+alpha)*x^3=200*(3+1.5*N+alpha)*(1+0.01)
 *
 * Найти:
 * число обусловленности матрицы рассматриваемой СЛАУ и относительную погрешность приближённой СЛАУ
 * Прокомментриовать получившиеся результаты
 * */


void Task1(){
	puts("Task1");
	const unsigned int N = 13;
	const unsigned int n = 53;
	const double alpha = (n-50)/((double) 100);
	double matrix_coeffs[9] = {
		200*(1+0.5*N+alpha), 200*(1+0.5*N)        , 200*(1+0.5*N)        ,
		200.1*(1+0.5*N)    , 199.9*(1+0.5*N+alpha), 200*(1+0.5*N)        ,
		199.9*(1+0.5*N)    , 200*(1+0.5*N)        , 200.1*(1+0.5*N+alpha)
	};
	Matrix SLE = {.data = matrix_coeffs, .rows = 3, .columns = 3};
	double vec_coeffs[3] = {
		200*(3+1.5*N+alpha),
		200*(3+1.5*N+alpha),
		200*(3+1.5*N+alpha)
	};
	Vector b = {.data = vec_coeffs, .length = 3};
	double vec2_coeffs[3] = {
		200*(3+1.5*N+alpha)*0.01,
		200*(3+1.5*N+alpha)*-0.01,
		200*(3+1.5*N+alpha)*0.01
	};
	Vector Db = {.data = vec2_coeffs, .length = 3};
	double imat_coeffs[9];
	for (unsigned int i = 0; i < 9; ++i)
		imat_coeffs[i] = matrix_coeffs[i];
	Matrix Inversed = {.data = imat_coeffs, .rows = 3, .columns = 3};
	inverse(SLE, Inversed);	
	double c = cond(SLE, Inversed);
	printf("cond: %d\n", c);
	if (c > 100)
		puts("SLE is ill conditioned");
	else
		puts("SLE is well conditioned");
	
	double norm_Db = norm(Db), norm_b = norm(b);
	double relative_error = c*norm_Db/norm_b;
	printf("Relative error: %lf, |b|: %e, |Db|: %e\n", relative_error, norm_b, norm_Db);
	puts("_________________________________________________________________________");
}

/*
 * Задание 1.2
 * N = 13, alpha = 0.02
 * lambda+alpha=0.4
 * F = arctg
 * a = 0
 * b = 1
 * На отрезке [a, b] выбрана центрально равномерная сетка с десятью узлами:
 * s1=t1=a+h/2,s2=t2=t1+h, s3=t3=t2+h, ..., s10=t10=t9+h
 * h = (b-a)/10
 * Решить приближённую СЛАУ
 * (E+lambda A)x=b+Db (5)
 * A=a^i_j, b=b^i:
 * a^i_j = F(s_i*t_j)*(b-a)/10
 * b = (E+lambda A)*x
 * x^i = 1
 * приближённая СЛАУ определяется только погрешностью Db=0.01[b^1, -b^2, b^3, -b^4, b^5, -b^6, b^7, -b^8, b^9, -b^10]
 * Найти:
 * число обусловленности матрицы рассматриваемой СЛАУ и относительную погрешность в решении приближённой СЛАУ
 * Затем, прокомментировать получившиеся результаты.
 * Кроме того, найти решение СЛАУ, которая получается из СЛАУ (5) делением каждого i-го уравнения на число b^i+Db^i
 * После этого сравнить абсолютную погрешность в решении получившейся СЛАУ с абсолютной погрешностью в решении приближённой СЛАУ
 * */

typedef double (*math_fn)(double);

void Task2(){
	puts("Task2");
	const unsigned short N = 13;
	const double alpha = 0.02;
	const double lambda = 0.4 - alpha;
	const math_fn F = atan;
	const double a = 0;
	const double b = 1;
	const unsigned int GRID = 10;
	double s[GRID];
	double h = (b-a)/GRID;
	{
	unsigned int i = 0;
	s[i] = a+h/2;
	assert(GRID > 1);
	while (i < GRID-1){
		s[i+1] = s[i]+h;
		++i;
	}
	}
	double a_ij[GRID*GRID];
	for (int i = 0; i < GRID; ++i)
		for (int j = 0; j < GRID; ++j)
			a_ij[j+GRID*i] = F(s[i]*s[j])*(b-a)/10;
	double x_i[GRID];
	for (int i = 0; i < GRID; ++i)
		x_i[i] = 1;
	double b_i[GRID];
	double E_ij[GRID*GRID] = {};
	for (int i = 0; i < GRID; ++i)
		E_ij[i + GRID*i] = 1;
	Matrix A = {
		.data = a_ij, 
		.rows = GRID,
		.columns = GRID
	};
	Vector x = {
		.data = x_i,
		.length = GRID
	};
	Vector bvec = {
		.data = b_i,
		.length = GRID
	};
	Matrix E = {
		.data = E_ij,
		.rows = GRID,
		.columns = GRID
	};
	double AplE_ij[GRID*GRID];
	Matrix Aple = {
		.data = AplE_ij,
		.rows = GRID,
		.columns = GRID
	};
	Multiply(lambda, A, Aple);
	Sum(E, Aple, Aple);
	Apply(Aple, x, bvec);
	double c = cond(Aple);
	printf("cond: %lf\n", c);
	
	double Db_i[GRID];
	{
	int sign = 1;
	for (int i = 0; i < GRID; ++i){
		Db_i[i] = sign*0.01*b_i[i];
		sign = -sign;
	}
	}
	Vector Db = {
		.data = Db_i,
		.length = GRID
	};
	double norm_Db = norm(Db), norm_b = norm(bvec);
	double relative_error = c*norm_Db/norm_b;
	printf("Relative error: %lf, |b|: %e, |Db|: %e\n", relative_error, norm_b, norm_Db);
	/*
		найти решение СЛАУ, которая получается делением каждого i-го уравнения на b[i]+Db[i]
		(E + lambda A) x = b + Db
		(E + lamdbda A)/(b+Db) x = 1
	*/
	puts("_________________________________________________________________________");
}

int main(){
	Task1();
	Task2();
	return 0;
}
