#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h> 

typedef struct Vector{
	alignas(32) double* data;
	int count;
	int capacity;
} Vector;

typedef struct UMatrix{
	alignas(32) double* data;
	int rows;
} UMatrix;

static int cap_for_align(int bytes, int alignment){
	if (bytes % alignment == 0)
		return bytes;
	else	 
		return bytes - bytes%alignment + alignment;
}

static void coords_to_local(const int i, const int j, int * const index){
	static_assert( INT_MAX == 2147483647);
	//        46340*(46341) = 2147441940
	//        46341*(46342) = 2147534622
	assert( i < 46341); // тогда нет переполнения
	assert( i >= j);
	*index = (i*(i+1))/2 + j;
}

static double get_number(
	const double * const data,
	const int i,
	const int j
){
	int local_index;
	coords_to_local(i, j, &local_index);
	return data[local_index];
}


static void Matrix_print(const int rows, const double _data[static (rows*(rows+1))/2]){
	const double * const data = (const double *) __builtin_assume_aligned(_data, 32); 
	assert( rows < 46341);
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < rows; ++j){
			if (j != 0)
				printf(", ");
			else
				printf("[ ");
			int local;
			if (i >= j)
				coords_to_local(i, j, &local);
			else
				coords_to_local(j, i, &local);
			assert(local < (rows*(rows+1))/2);
			printf("%6lf", data[local]);
		}
		puts("]");
	}
}

static void Eigenvectors_print(const int rows, const double _data[static rows*rows]){
	const double * const data = (const double *) __builtin_assume_aligned(_data, 32); 
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < rows; ++j){
			if (j != 0)
				printf(", ");
			else
				printf("[ ");
			printf("%6lf", data[j*rows + i]); // результат должен получиться транспонированным
		}
		puts("]");
	}
}

static double max_nondiag_elem(
	const int rows,
	const double _data[static (rows*(rows+1))/2],
	int *current_i, int *current_j
){
	assert( rows > 1);
	const double * const data = (const double *) __builtin_assume_aligned(_data, 32); 
	assert( rows < 46341);
	
	double current = -1;
	double elem = 0;
	for (int i = 1; i < rows; ++i)
		for (int j = 0; j < i; ++j){
			int index;
			coords_to_local(i, j, &index);
			const double value = data[index];
			const double abs_value = fabs(value);
			if (current < abs_value){
				elem = value;
				current = abs_value;
				*current_i = i;
				*current_j = j;
			}
		}
	assert(*current_i > *current_j);
	return elem;
}

static void JacobiMethod(
	const int rows,
	const double _data[static (rows*(rows+1))/2],
	double _eigenvalues[static rows],
	double _eigenvectors[static rows*rows],
	const double epsilon,
	const int maxcount,
	int *count
){
	assert( rows < 46341);
	const double * const data = (const double *) __builtin_assume_aligned(_data, 32); 
	double * const eigenvalues = (double *) __builtin_assume_aligned(_eigenvalues, 32); 
	double * const eigenvectors = (double *) __builtin_assume_aligned(_eigenvectors, 32);
	
	const int matrix_sz = (rows*(rows+1))/2;
	const int BYTE_SZ = cap_for_align(sizeof(double)*matrix_sz, 32);
	alignas(32) double * A_i;
	alignas(32) double * A_next;
	const bool dynamic_arrays_needed = BYTE_SZ > 65536;
	if (dynamic_arrays_needed){
		// большие массивы, храним в куче
		A_i = aligned_alloc(32, BYTE_SZ);
		A_next = aligned_alloc(32, BYTE_SZ);
	}
	else{
		// маленькие, можем хранить на стеке
		A_i = __builtin_alloca_with_align(BYTE_SZ, 32);
		A_next = __builtin_alloca_with_align(BYTE_SZ, 32);
	}
	assert(A_i);
	assert(A_next);
	
	memcpy(A_i, data, sizeof(double)*matrix_sz);
	memcpy(A_next, data, sizeof(double)*matrix_sz);
	
	int alpha, beta;
	for (
		double max_elem = max_nondiag_elem(rows, A_i, &beta, &alpha);
		
		fabs(max_elem) > epsilon && *count < maxcount;
		
		memcpy(A_i, A_next, sizeof(double)*matrix_sz),
		max_elem = max_nondiag_elem(rows, A_i, &beta, &alpha),
		++*count
	){
		assert(alpha < beta);

#ifndef NDEBUG
		printf("Current iteration: %d. Current fabs of max_elem: %lf.\n",
									*count+1, fabs(max_elem));
		printf("alpha: %d, beta: %d\n", alpha, beta);
		printf("A[%d]:\n", *count);
		Matrix_print(rows, A_i);
		printf("Q[%d]:\n", *count);
		Eigenvectors_print(rows, eigenvectors);
		
#endif

		const double a_aa = get_number(A_i, alpha, alpha);
		const double a_ba = max_elem;
		const double a_bb = get_number(A_i, beta, beta);
		
	//	printf("a_aa = %lf, a_bb = %lf, a_aa - a_bb = %lf\n", a_aa, a_bb, a_aa - a_bb);
	//	printf("2*a_ba / (a_aa - a_bb) = %lf\n", 2*a_ba / (a_aa - a_bb));
		const double varphi = 1/2.0 * atan(2*a_ba / (a_aa - a_bb));
		//printf("varphi = %lf\n", varphi);
		
		const double cos_phi = cos(varphi);
		const double sin_phi = sin(varphi);

	//	printf("sin(phi)=%lf, cos(phi)=%lf\n",
	//			sin_phi, cos_phi);
	
		// исправить потом
		const double b_aa = a_aa * cos_phi * cos_phi +
							a_bb * sin_phi * sin_phi +
						2 * a_ba * sin_phi * cos_phi;
		const double b_bb = a_aa * sin_phi * sin_phi +
							a_bb * cos_phi * cos_phi -
						2 * a_ba * sin_phi * cos_phi;
		const double b_ba = 0;

		
		// подсчёт собственных векторов:
		// Q[k] = Q[k-1]*Q_k

		// alpha = 0, beta = 3
		//
		// cos	0	-sin	0
		// 0	1	0		0
		// sin	0	cos		0
		// 0	0	0		1
		//
		//	0	1	2		3
		//	4	5	6		7
		//	8	9	10		11
		//	12	13	14		15
		//
		//	=>
		//
		//	0*cos +2*sin	1	0*-sin  + 2*cos		3
		//	4*cos +6*sin	5	4*-sin  + 6*cos		7
		//	8*cos +10*sin	9	8*-sin  + 10*cos	11
		//	12*cos+14*sin	13	12*-sin + 14*cos	15
		//	
		
		for (int j = 0; j < rows; ++j){
			const int aindex = alpha*rows + j;
			const int bindex = beta*rows + j;

			const double Qak = eigenvectors[aindex];
			const double Qbk = eigenvectors[bindex];

			eigenvectors[aindex] = cos_phi*Qak +
								   sin_phi*Qbk;
			eigenvectors[bindex] = -sin_phi*Qak+
								   cos_phi*Qbk;
		}

		// подсчёт собственных значений
		for (int i = 0; i < rows; ++i){
			for (int j = 0; j <= i; ++j){
				int index;
				coords_to_local(i, j, &index);

				if (i < alpha){
					continue;
				}
				if (i == alpha){
					if (j == alpha){
						A_next[index] = b_aa;
					} 
					//else if(j == beta невозможно, alpha < beta и i >= j
					else {
						// j < alpha < beta
						int bindex;
						coords_to_local( beta, j, &bindex);
						const double a_ak = A_i[index];
						const double a_bk = A_i[bindex];
						const double b_ak = a_ak*cos_phi+
											a_bk*sin_phi;
						A_next[index] = b_ak;
					}
				}
				else if (i < beta){
					if (j == alpha){
						int bindex;
						coords_to_local( beta, j, &bindex);

						const double a_ak = A_i[index];
						const double a_bk = A_i[bindex];
						const double b_ak = a_ak*cos_phi+
											a_bk*sin_phi;
						A_next[index] = b_ak;
					}
					else {
						continue;
					}
				}
				else if (i == beta){
					if (j == alpha){
						A_next[index] = b_ba;
					}
					else if (j == beta){
						A_next[index] = b_bb;
					}
					else {
						int aindex;
						coords_to_local( j, alpha, &aindex);
						const double a_ak = A_i[aindex];
						const double a_bk = A_i[index];
						const double b_bk = -a_ak*sin_phi+
											a_bk*cos_phi;

						A_next[index] = b_bk;
					}
				}
				else {
					if (j == alpha){
						int bindex;
						if (i < beta){
							coords_to_local( , , &bindex);
						}
					}
				}
			}
		}
		// b_ak = a_ak *  cos_phi + a_bk * sin_phi
		// b_kb = a_ak * -sin_phi + a_bk * cos_phi
		// b_kk = a_kk
		
		// alpha < beta
	
		for (int i = 0; i < rows; ++i){
			for (int j = 0; j < i; ++j){
				int local_index;
				coords_to_local(i, j, &local_index);
				if (i == alpha){
					if (j == beta){
						A_i[local_a] = b_ak;
						A_i[local_b] = b_bk;
					}
				}
			}
		}
		alignas(32) double a_kb[beta - alpha] = {};
		alignas(32) double a_ak[beta - alpha] = {};
		
		int local_a;
		int local_b;
		coords_to_local(alpha, 0, &local_a);
		coords_to_local(beta, 0, &local_b);
		
		// итерируем по общей части i = alpha и beta одновременно для кеш-локальности
		// т.е. i = alpha, beta; j = 0..alpha
		
		int iindex = 0;
		for (int j = 0; j <= beta; ++j, ++local_b){
			const double a_bk = A_i[local_b];
			//a_kb[j] = a_bk;
			if (j < alpha){
				const double a_ka = A_i[local_a];
				const double b_ak = a_ka * cos_phi +
									a_bk * sin_phi;
				const double b_bk = -a_ka * sin_phi +
									a_bk * cos_phi;
				
				A_i[local_a] = b_ak;
				A_i[local_b] = b_bk;

//				printf("local_a: %d, local_b: %d\n a_bk: %lf a_ka: %lf b_ak: %lf b_bk: %lf\n",
//						local_a, local_b, a_bk, a_ka, b_ak, b_bk);

				++local_a;
			} else if (j == alpha){
				const double a_ka = A_i[local_a];
				const double b_ak = b_aa;
				const double b_bk = b_ba;
									
				A_i[local_a] = b_ak;
				A_i[local_b] = b_bk;
				++local_a;
			} else{
				a_kb[iindex] = a_bk;
				++iindex;
			}
		}
		
//		printf("alpha: %d, beta: %d\n", alpha, beta);
//		printf("A[%d]:\n", *count);
//		Matrix_print(rows, A_i);
		assert(iindex == beta - alpha);
		int test_index;
		coords_to_local(alpha+1, 0, &test_index);
		assert(test_index == local_a);
		
		/*
		for (int i = 0, int local_index = 0;
			i < alpha; ++i){
			// для строк до i = alpha при таком хранении ничего не меняется
			continue;			
		}
		*/

		// для i от alpha до beta нас интересует только столбец j = alpha
		local_a += alpha;
		int jindex = 0;
		for (int i = alpha + 1; i < beta; ++i){
			const double a_ka = A_i[local_a];
			const double a_bk = a_kb[jindex];
			
			const double b_ak = a_ka * cos_phi +
								a_bk * sin_phi;

			
			a_ak[jindex] = a_ka;
			A_i[local_a] = b_ak;

			printf("local_a: %d, local_b: %d\n a_bk: %lf a_ka: %lf b_ak: %lf, jindex: %d \n",
					local_a, local_b, a_bk, a_ka, b_ak, jindex);
			local_a += i+1;
			++jindex;
		}
		printf("alpha: %d, beta: %d\n", alpha, beta);
		printf("A[%d]:\n", *count);
		Matrix_print(rows, A_i);
		
		assert(jindex == iindex-1);
		coords_to_local(beta, alpha, &test_index);
		assert(test_index == local_a);
		++local_a;
		
		int zindex = 0;
		for (int j = alpha + 1; j < beta; ++j){
			const int i = beta;
				
			const double a_ka = a_ak[zindex];
			const double a_bk = A_i[local_a];

			const double b_bk = -a_ka * sin_phi +
								a_bk * cos_phi;

			A_i[local_a] = b_bk;

			printf("local_a: %d, local_b: %d\n a_bk: %lf a_ka: %lf, ndex: %d \n",
					local_a, local_b, a_bk, a_ka, zindex);
			
			++zindex;
			++local_a;
		}

		A_i[local_a] = b_bb;

		printf("alpha: %d, beta: %d\n", alpha, beta);
		printf("A[%d]:\n", *count);
		Matrix_print(rows, A_i);

		assert(zindex == jindex);
		coords_to_local(beta+1, 0, &test_index);
		printf("local_a: %d, local_b: %d, test_index: %d\n",
				local_a, local_b, test_index);
		assert(test_index == local_a);
		local_b = local_a + beta;
		local_a += alpha;

		for (int i = beta + 1; i < rows; ++i){
			const double a_ka = A_i[local_a];
			const double a_bk = A_i[local_b];

			const double b_ak = a_ka * cos_phi +
								a_bk * sin_phi;
			const double b_bk = -a_ka * sin_phi +
								a_bk * cos_phi;
			
			A_i[local_a] = b_ak;
			A_i[local_b] = b_bk;
			
			local_a += i + 1;
			local_b += i + 1;
		}
		
		printf("alpha: %d, beta: %d\n", alpha, beta);
		printf("A[%d]:\n", *count);
		Matrix_print(rows, A_i);
		coords_to_local(rows, alpha, &test_index);
		assert(test_index == local_a);
		coords_to_local(rows, beta, &test_index);
		assert(test_index == local_b);
		
		/*
		Q(alpha, beta, varphi) :=
			Q[i, j] = 1, если i == j и (i, j) != (alpha, beta)
			Q[i, j] = 0, если i != j и (i, j) != (alpha, beta)
			Если i или j == alpha или beta:
				Q[alpha, alpha] = cos varphi
				Q[alpha,  beta] = -sin varphi
				Q[ beta, alpha] = sin varphi
				Q[ beta,  beta] = cos varphi
				Всё остальное   = 0
		*/
	}

	
	for (int i = 0; i < rows; ++i){
		// i + i(i+1)/2 = i(i+3)/2
		eigenvalues[i] = A_i[((i+3)*i)/2];
	}
	
	if (dynamic_arrays_needed)
		free(A_i);
}

int main(){
	constexpr int N = 13;
	constexpr double beta = 1+0.1*(52-53);
	constexpr double varepsilon = 0.01;
	constexpr int rows = 4;
	alignas(32) double A_data[] = {
		10*beta,//     1	   2	     3
			  1, 10*beta,//	   -3		-2
			  2, 	  -3, 10*beta,//	-1
			  3, 	  -2, 	   -1, 10*beta
	};
	alignas(32) double eigenvalues[rows] = {};
	alignas(32) double eigenvectors[rows*rows] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};

	Matrix_print(rows, A_data);
	constexpr int maxcount = 1000;
	int count = 0;
	JacobiMethod(rows, A_data, eigenvalues, eigenvectors, varepsilon, maxcount, &count);
	printf("A[%d]:\n", count);
	Matrix_print(rows, A_data);
	printf("Q[%d]:\n", count);
	Eigenvectors_print(rows, eigenvectors);
	
	return 0;
}

/*
Удалил:
typedef struct Matrix{
	alignas(64) float *data;
	int count;
	int rows;
	int columns;
} Matrix;

static inline void UMatrix_global_to_local(const UMatrix matrice, const index_t global_i, const index_t global_j, index_t *local){
	//assert(global_i > -1 && global_j > -1 && global_i <= global_j);
	*local = global_i*matrice.rows - (global_i * (global_i + 1)) / 2 + global_j;
}

static inline uint64_t int_sqrt(uint64_t n) {
    if (n <= 1) return n;
    uint64_t x = n;
    uint64_t y = (x + n / x) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

static inline void UMatrix_local_to_global(const UMatrix matrice, const index_t local, index_t *global_i, index_t *global_j) {
    const index_t N = matrice.rows;
    const index_t k = local;
    
    const uint64_t twoN_plus_one = 2 * (uint64_t)N + 1;
    const uint64_t D = twoN_plus_one * twoN_plus_one - 8 * (uint64_t)k;
    
    uint64_t sqrt_D = int_sqrt(D);
    
    index_t i = (index_t)((twoN_plus_one - sqrt_D) / 2);
    
    uint64_t k_min_i;
    do {
        k_min_i = (uint64_t)i * N - ((uint64_t)i * (i + 1)) / 2 + i;
        if (k_min_i > (uint64_t)k) {
            i--;
        } else {
            break;
        }
    } while (1);
    
    while (i + 1 < N) {
        uint64_t k_min_i1 = (uint64_t)(i + 1) * N - ((uint64_t)(i + 1) * (i + 2)) / 2 + (i + 1);
        if (k_min_i1 <= (uint64_t)k) {
            i++;
        } else {
            break;
        }
    }
    
    *global_i = i;
    *global_j = k - i * N + i * (i + 1) / 2;
    
    assert(*global_i >= 0 && *global_j >= 0 && *global_i <= *global_j && *global_j < N);
}
*/
