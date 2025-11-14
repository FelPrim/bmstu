#include <stdio.h>
#include <math.h>

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

double Norm(const Matrix* matrix){
	double result = -1;
	for (unsigned int i = 0; i < matrix->rows; ++i){
		double row_result = 0;
		for (unsigned int j = 0; j < matrix->columns; ++j){
			row_result += fabs(matrix[j+matrix->columns*sizeof(double)*i]);
		}
		if (row_result > result)
			result = row_result;
	}
	return result;
}

struct MinorStruct{
	Matrix matrix;
	const unsigned short row;
	const unsigned short column;
};

inline static xifaelsey(double x, double y, int a, int b){
	int equal = a == b;
	return equal*x+(1-equal)*y;
}

double minor(const MinorStruct * input);
double minor(const MinorStruct * input){

	switch (input->matrix.rows){
		case 0:
		case 1:
			return 0;
		case 2:
			return input->matrix.data[3 - input->column - 2*input->row];
		case 3:
			int i0 = input->row != 0;
			int i1 = input->row != 1;
			int i2 = input->row != 2;
			int j0 = input->column != 0;
			int j1 = input->column != 1;
			int j2 = input->column != 2;
			int a11 = 4 - j0 - 3*i0;
			int a12 = 4 + j2 - 3*i0;
			int a21 = 4 - j0 + 3*i2;
			int a22 = 4 + j2 + 3*i2;
			return input->matrix.data[a11] * input->matrix.data[a22]
			     - input->matrix.data[a12] * input->matrix.data[a21];
		case 4:
			


	}	
}

double cond(const Matrix* restrict matrix, const Matrix* inversed){
	return Norm(matrix)*Norm(inversed);
}

struct 

void adjugate(const Matrix* original, Matrix* result){
	for (unsigned short i = 0; i < result->rows; ++i){
		for (unsigned short j = 0; j < result->columns; ++j){	
			// original[j + matrix->columns*sizeof(double)*i] - original[i][j]
			double det_before_sum[result->columns];
			for (unsigned short i2 = 0; i2 < original->columns; ++i2){
				if (i2 == j)
					continue;
				for (unsigned short j2 = 0; j2 < original->rows; ++j2){
					if (j2 == i)
						continue;
					det_before_sum[i2] +=
				}
			}
			//double cofactor = //-1^i+j * minor
		}
	}
}

void inverse(const Matrix* original, Matrix* result){
		
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
		200*(3+1.5*N+alpha)*(1+0.01),
		200*(3+1.5*N+alpha)*(1-0.01),
		200*(3+1.5*N+alpha)*(1+0.01)
	};
	Vector bDb = {.data = vec2_coeffs, .length = 3};
	double imat_coeffs[9];
	for (unsigned int i = 0; i < 9; ++i)
		imat_coeffs[i] = matrix_coeffs[i];
	Matrix Inversed = {.data = imat_coeffs, .rows = 3, .columns = 3};
		
}

int main(){

	return 0;
}
