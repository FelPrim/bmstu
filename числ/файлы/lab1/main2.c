#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

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
			return 1;
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
	double det = _determinant(original);
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


typedef double (*math_fn)(double);

int main(){
	Task1();
	return 0;
}
