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
	double * data;
	const unsigned short rows;
	const unsigned short columns;
} Matrix;

typedef struct Vector{
	double * data;
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

inline static double xifaelsey(double x, double y, int a, int b){
	int equal = a == b;
	return equal*x+(1-equal)*y;
}



inline static double _determinant(const Matrix matrix){
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
			exit(-1);
	}
	
}

typedef struct GapeMatrix{
	Matrix matrix;
	unsigned short *rows;
	unsigned short count; // from 0 to SHORT_MAX
} GapeMatrix;

inline static unsigned short find_index(const unsigned short *sorted_array, unsigned short length, unsigned short elem){
	unsigned short left = 0, right = length - 1;
	unsigned short mid;
	while (left <= right){
		mid = left + (right - left)/2;
		if (sorted_array[mid] == elem){
			return length;
		}
		if (sorted_array[mid] < elem)
			left = mid + 1;
		else
			right = mid;
	}
	return left;
}

double recursive_determinant(const GapeMatrix);
double recursive_determinant(const GapeMatrix gape){
	puts("@");
	switch (gape.count){
		case 0:
		case 0x8000:
			exit(-1);
		case 0x8001:{
			// 0 0 0 0 0
			// 0 0 0 0 X
			// 0 0 0 0 0
			// 0 0 0 0 0
			// 0 0 0 0 0
			const int a11 =  gape.matrix.columns - 1 + gape.rows[0] * gape.matrix.columns;
			return gape.matrix.data[a11];}
		case 0x8002:{
			// 0 0 0 0 0 
			// 0 0 0 A B
			// 0 0 0 0 0
			// 0 0 0 0 0
			// 0 0 0 C D
			const int a11 = gape.matrix.columns - 2 + gape.rows[0] * gape.matrix.columns;
			const int a12 = gape.matrix.columns - 1 + gape.rows[0] * gape.matrix.columns;
			const int a21 = gape.matrix.columns - 2 + gape.rows[1] * gape.matrix.columns;
			const int a22 = gape.matrix.columns - 1 + gape.rows[1] * gape.matrix.columns;
			return 	  gape.matrix.data[a11] * gape.matrix.data[a22] 
				- gape.matrix.data[a12] * gape.matrix.data[a21];}
		case 0x8003:{
			// 0 0 A B C
			// 0 0 0 0 0
			// 0 0 D E F
			// 0 0 G H I
			// 0 0 0 0 0
			const int a11 = gape.matrix.columns - 3 + gape.rows[0] * gape.matrix.columns;
			const int a12 = gape.matrix.columns - 2 + gape.rows[0] * gape.matrix.columns;
			const int a13 = gape.matrix.columns - 1 + gape.rows[0] * gape.matrix.columns;
			const int a21 = gape.matrix.columns - 3 + gape.rows[1] * gape.matrix.columns;
			const int a22 = gape.matrix.columns - 2 + gape.rows[1] * gape.matrix.columns;
			const int a23 = gape.matrix.columns - 1 + gape.rows[1] * gape.matrix.columns;
			const int a31 = gape.matrix.columns - 3 + gape.rows[2] * gape.matrix.columns;
			const int a32 = gape.matrix.columns - 2 + gape.rows[2] * gape.matrix.columns;
			const int a33 = gape.matrix.columns - 1 + gape.rows[2] * gape.matrix.columns;
			return    gape.matrix.data[a11] * gape.matrix.data[a22] * gape.matrix.data[a33]
				+ gape.matrix.data[a12] * gape.matrix.data[a23] * gape.matrix.data[a31]
				+ gape.matrix.data[a13] * gape.matrix.data[a21] * gape.matrix.data[a32]
				- gape.matrix.data[a13] * gape.matrix.data[a22] * gape.matrix.data[a31]
				- gape.matrix.data[a12] * gape.matrix.data[a21] * gape.matrix.data[a33]
				- gape.matrix.data[a11] * gape.matrix.data[a23] * gape.matrix.data[a32];}
		default:
			const unsigned short is_remaining = 0x8000 & gape.count;
			unsigned short real_count;
			if (is_remaining)
				real_count = gape.count ^ 0x8000;
			else
				real_count = gape.count;
			double result = 0;
			int sign = 1;
			if (!is_remaining){
				const unsigned short remaining  = gape.matrix.columns - real_count;
				if (remaining > real_count + 1 && remaining > 3){
					// 0 A B C D E
					// 0 0 0 0 0 0
					// 0 F G H I J
					// 0 K L M N O
					// 0 P Q R S T
					// 0 U V W X Y
					// [1]

					unsigned short new_rows[real_count + 1];
					for (unsigned short i = 0; i < gape.matrix.columns; ++i){
						unsigned short index = find_index(gape.rows, real_count, i);
						if (index == real_count)
							// index был в массиве
							continue;
						unsigned short offset = 0;
						for (unsigned short j = 0; j < real_count; ++j){
							if (j == index){
								new_rows[j] = i; 
								++offset;
								continue;
							}
							new_rows[j] = gape.rows[j+offset];
						}
						GapeMatrix tmp = {
							.matrix = gape.matrix,
							.rows = new_rows,
							.count = real_count+1
						};
						result += sign * recursive_determinant(tmp);
						sign = -sign;
					}
				}
				else{
					// 0 0 0 A B C C2     0 A B C
					// 0 0 0 0 0 0 0      0 0 0 0
					// 0 0 0 D E F F2 или 0 D E F
					// 0 0 0 G H I I2     0 G H I
					// 0 0 0 0 0 0 0
					// 0 0 0 J K L M
					// 0 0 0 0 0 0 0 	 
					// меняем то, как храним массив: [157] -> [0346]
					unsigned short cur_rows[remaining];
					unsigned short j = 0;
					for (unsigned short i = 0; i < gape.matrix.columns; ++i){
						unsigned short index = find_index(gape.rows, real_count, i);
						if (index == real_count){
							cur_rows[j++] = i;
						}
					}
					unsigned short new_rows[remaining-1];
					for (unsigned short i = 0; i < remaining; ++i){
						unsigned short index = cur_rows[i];
						sign = index % 2;
						unsigned short offset = 0;
						for (unsigned short j = remaining; j > 0; ++j){
							if (j == index){
								++offset;
								continue;
							}
							new_rows[j+offset] = cur_rows[j];
						}
						GapeMatrix tmp = {
							.matrix = gape.matrix,
							.rows = new_rows,
							.count = (remaining-1) | 0x8000
						};
						result += sign * recursive_determinant(tmp);
					}
				}
			}
			else{
				unsigned short new_rows[real_count-1];
				for (unsigned short i = 0; i < real_count; ++i){
					unsigned short index = gape.rows[i];
					sign = index % 2;
					unsigned short offset = 0;
					for (unsigned short j = real_count; j > 0; ++j){
						if (j == index){
							++offset;
							continue;
						}
						new_rows[j+offset] = gape.rows[j];						
					}
					GapeMatrix tmp = {
						.matrix = gape.matrix,
						.rows = new_rows,
						.count = (real_count-1) | 0x8000
					};
					result += sign * recursive_determinant(tmp);

				}
			}

			return result;
		
	}

}

double determinant(const Matrix matrix){
	if (matrix.columns < 4)
		return _determinant(matrix);
	if (matrix.columns == 4){
		unsigned short arr[3] = {1, 2, 3};	
		GapeMatrix tmp = {
			.matrix = matrix,
			.rows = arr,
			.count = 0x8003
		};
		double result = recursive_determinant(tmp);
		tmp.rows[0] = 0;
		result += -recursive_determinant(tmp);
		tmp.rows[1] = 1;
		result += recursive_determinant(tmp);
		tmp.rows[2] = 2;
		result += -recursive_determinant(tmp);
		return result;
	}
	double result = 0;
	int sign = 1;
	for (unsigned short i = 0; i < matrix.columns; ++i){
		GapeMatrix tmp = {
			.matrix = matrix,
			.rows = &i,
			.count = 1
		};
		result += recursive_determinant(tmp);
		sign = -sign;
	}
	return result;
		
}

typedef struct MinorStruct{
	Matrix matrix;                  
	const unsigned short row;       
	const unsigned short column;   
} MinorStruct;

double minor(const MinorStruct input){
	unsigned short n = input.matrix.rows;
	if (n == 0) return 1.0; 
	if (n == 1) return 1.0;

	if (input.matrix.rows != input.matrix.columns){
		fprintf(stderr, "minor: matrix must be square\n");
		exit(-1);
	}

	unsigned short m = n - 1;

	double *sub = (double*)malloc(sizeof(double) * (size_t)m * (size_t)m);
	if (!sub){
		fprintf(stderr, "minor: malloc failed\n");
		exit(-1);
	}

	unsigned short si = 0;
	for (unsigned short i = 0; i < n; ++i){
		if (i == input.row) continue;
		unsigned short sj = 0;
		for (unsigned short j = 0; j < n; ++j){
			if (j == input.column) continue;
			sub[sj + si * m] = input.matrix.data[j + input.matrix.columns * i];
			++sj;
		}
		++si;
	}

	Matrix subM = { .data = sub, .rows = m, .columns = m };
	double det = determinant(subM);

	free(sub);
	return det;
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
	double det = determinant(original);
	printf("det: %.15lf\n", det);
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
	printf("cond: %.15lf\n", c);
	if (c > 100)
		puts("SLE is ill conditioned");
	else
		puts("SLE is well conditioned");
	
	double norm_Db = norm(Db), norm_b = norm(b);
	double relative_error = c*norm_Db/norm_b;
	printf("Relative error: %.15lf, |b|: %.15lf, |Db|: %.15lf\n", relative_error, norm_b, norm_Db);
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
	double I_data[GRID*GRID];
	Matrix Inversed = {
		.data = I_data,
		.rows = GRID,
		.columns = GRID
	};
	inverse(Aple, Inversed);
	double c = cond(Aple, Inversed);
	printf("cond: %.15lf\n", c);
	if (c > 100)
		puts("SLE is ill conditioned");
	else
		puts("SLE is well conditioned");
	

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
	printf("Relative error: %.15lf, |b|: %.15lf, |Db|: %.15lf\n", relative_error, norm_b, norm_Db);
	/*
		найти решение СЛАУ, которая получается делением каждого i-го уравнения на b[i]+Db[i]
		(E + lambda A) x = b + Db
		(E + lamdbda A)/(b+Db) x = 1
	*/
	for (unsigned short i = 0; i < GRID; ++i)
		for (unsigned short j = 0; j < GRID; ++j)
			Aple.data[j + i*GRID] /= bvec.data[i] + Db.data[i];
	for (unsigned short i = 0; i < GRID; ++i){
		Db.data[i] /= bvec.data[i] + Db.data[i];
		bvec.data[i] /= bvec.data[i] + Db.data[i];
	}
	inverse(Aple, Inversed);
	c = cond(Aple, Inversed);
	printf("cond: %.15lf\n", c);
	norm_Db = norm(Db);
	norm_b = norm(bvec);
	relative_error = c*norm_Db/norm_b;
	printf("Relative error: %.15lf\n", relative_error);
	puts("_________________________________________________________________________");
}

int main(){
	Task1();
	Task2();
	return 0;
}
