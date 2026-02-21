#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

constexpr bool DEBUG_EVERYTHING = true;
constexpr int DIM = 4;

static void print_array(const double array[static DIM]){
	printf("[%.16lf", array[0]);
	for (int i = 1; i < DIM; ++i)
		printf(", %.16lf", array[i]);
	puts("]");
}

static void MATRIX_PRINT(const double matrix[static DIM*DIM]){
	for (int i = 0; i < DIM; ++i)
		print_array(matrix+i*DIM);
}

static inline int get_index(const int i, const int j){
	assert(0 <= i && i < DIM);
	assert(0 <= j && j <= i);
	return i*(i+1)/2 + j;
}

void Matrix_print(const double matrix[static DIM*(DIM+1)/2]){
	int cur = 0;
	for (int i = 0; i < DIM; ++i){
		printf("[%8.5lf", matrix[cur]);
		++cur;
		for(int j = 1; j <= i; ++j){
			printf(", %8.5lf", matrix[cur]);
			++cur;
		}
		for (int j = i+1; j < DIM; ++j){
			printf(", %8.5lf", matrix[get_index(j, i)]);
		}
		printf("]\n");
	}
}


void Eigenvectors_print(const double matrix[static DIM*DIM]){
	for (int i = 0; i < DIM; ++i){
		printf("[(%8.5lf)", matrix[i]);
		for (int j = 1; j < DIM; ++j)
			printf(" (%8.5lf)", matrix[i + DIM*j]);
		printf("]\n");
		
	}
}


constexpr int N = 13;
constexpr int n = 53;
constexpr double beta = 1 + 0.1*(52-n);

constexpr double eps = 0.01;

double find_max(const double A[static DIM*(DIM+1)/2],
	int * restrict alpha, int * restrict beta){
	
	double cur = -1;
	int index = 1;
	for (int i = 1; i < DIM; ++i){
		for (int j = 0; j < i; ++j){
			const double elem = A[index];
			if (fabs(elem) > cur){
				cur = fabs(elem);
				*alpha = j;
				*beta = i;
			}
			++index;
		}
		assert(index == i*(i+3)/2);
		++index; // скипаем диагональ
	}
	assert(*alpha < *beta);
	return cur;
}

void Jacobi_method(
	double eigenvalues[static DIM*(DIM+1)/2],
	double eigenvectors[static DIM*DIM]
){
	constexpr int maxiter = 1000;
	int alpha, beta;
	double max = find_max(eigenvalues, &alpha, &beta);
	int iter = 0;
	while (max > eps && iter < maxiter){
		const double varphi = 1/2.0 * atan(
				2*eigenvalues[get_index(beta, alpha)] /
		(eigenvalues[get_index(alpha, alpha)] - eigenvalues[get_index(beta, beta)])
		);
		if (DEBUG_EVERYTHING){
			printf("Iterations: %d, max: %lf, alpha: %d, beta: %d, varphi: %lf\n", 
							  iter,		 max,	  alpha,	 beta,		varphi);

			puts("A:");
			Matrix_print(eigenvalues);
			puts("Q:");
		//	MATRIX_PRINT(eigenvectors);
			Eigenvectors_print(eigenvectors);
		}

		const double cos_phi = cos(varphi);
		const double sin_phi = sin(varphi);
		
		{
		double alpha_vec[DIM], beta_vec[DIM];
		for (int i = 0; i < DIM; ++i){
//			const double E_ai = eigenvectors[i*DIM + alpha];
//			const double E_bi = eigenvectors[i*DIM +  beta];

			const double E_ai = eigenvectors[alpha*DIM + i];
			const double E_bi = eigenvectors[beta*DIM + i];
		//	const double E_ia = eigenvectors[i + DIM*alpha];
		//	const double E_ib = eigenvectors[i + DIM*beta];
			
			alpha_vec[i] = cos_phi*E_ai + sin_phi*E_bi;
			beta_vec[i] = -sin_phi*E_ai + cos_phi*E_bi;
		}
		if (DEBUG_EVERYTHING){
			puts("alpha and beta vectors");
			print_array(alpha_vec);
			print_array(beta_vec);
		}
		for (int i = 0; i < DIM; ++i){

			//eigenvectors[alpha + DIM*i] = alpha_vec[i];
			//eigenvectors[beta + DIM*i] =  beta_vec[i];

			eigenvectors[DIM*alpha + i] = alpha_vec[i];
			eigenvectors[DIM*beta  + i] =  beta_vec[i];
		}
		}
		
		//int i = 0;
		//for (; i < alpha; ++i)
		//	continue;

		int i = alpha;
		int index_a = get_index(alpha, 0);
		int index_b = get_index(beta, 0);

		for (int k = 0; k < alpha; ++k){
			const double a_ak = eigenvalues[index_a];
			const double a_bk = eigenvalues[index_b];

			const double b_ak =  cos_phi * a_ak + sin_phi * a_bk;
			const double b_bk = -sin_phi * a_ak + cos_phi * a_bk;

			eigenvalues[index_a] = b_ak;
			eigenvalues[index_b] = b_bk;

			++index_a;
			++index_b;
		}
		assert(index_a == get_index(alpha, alpha));
		assert(index_b == get_index(beta, alpha));

		const double a_aa = eigenvalues[index_a];
		const double a_ab = eigenvalues[index_b];

		eigenvalues[index_b] = 0;
		++index_b;
		const int index_after_ba = index_b;

		// a_bk
		double storage[(beta - alpha == 1)? 1: beta - alpha - 1];
		for (int k = 0; k < beta - alpha - 1; ++k){
			storage[k] = eigenvalues[index_b];
			++index_b;
		}
		
		const double a_bb = eigenvalues[index_b];
		
		const double b_aa = cos_phi*cos_phi*a_aa +
							sin_phi*sin_phi*a_bb +
						  2*sin_phi*cos_phi*a_ab;
		const double b_bb = sin_phi*sin_phi*a_aa +
							cos_phi*cos_phi*a_bb +
						 -2*sin_phi*cos_phi*a_ab;

		assert(index_a == get_index(alpha, alpha));
		printf("index_a: %d, index_b: %d, [beta, beta]: %d\n",
				index_a, index_b, get_index(beta, beta));
		assert(index_b == get_index(beta, beta));

		eigenvalues[index_a] = b_aa;
		eigenvalues[index_b] = b_bb;

		for (int i = alpha + 1; i < beta; ++i){
			index_a += i;
			assert(index_a == get_index(i, alpha));

			const double a_bk = storage[i - alpha - 1];
			const double a_ak = eigenvalues[index_a];

			const double b_ak =  cos_phi * a_ak + sin_phi * a_bk;
			const double b_bk = -sin_phi * a_ak + cos_phi * a_bk;

			eigenvalues[index_a] = b_ak;
			eigenvalues[index_after_ba + i - alpha - 1] = b_bk;
		}
		index_a += beta;
		assert(index_a == get_index(beta, alpha));
		assert(index_after_ba + beta - alpha - 1 == get_index(beta, beta));

		for (int i = beta + 1; i < DIM; ++i){
			index_a += i;
			index_b += i;
			assert(index_a == get_index(i, alpha));
			assert(index_b == get_index(i, beta));

			const double a_ak = eigenvalues[index_a];
			const double a_bk = eigenvalues[index_b];

			const double b_ak =  cos_phi * a_ak + sin_phi * a_bk;
			const double b_bk = -sin_phi * a_ak + cos_phi * a_bk;

			eigenvalues[index_a] = b_ak;
			eigenvalues[index_b] = b_bk;
		}

		max = find_max(eigenvalues, &alpha, &beta);
		++iter;
	}

	printf("Iterations: %d, max: %lf\n", iter, max);

}


int main(){
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
	SetConsoleCP(CP_UTF8);
#endif
	
	double A[DIM*(DIM+1)/2] = {
		10*beta,
			  1, 10*beta,
			  2,	  -3, 10*beta,
			  3,	  -2,	   -1, 10*beta
	};
	double eigenvectors[DIM*DIM];
	for (int i = 0; i < DIM; ++i)
		eigenvectors[DIM*i + i] = 1;
	double eigenvalues[DIM*(DIM+1)/2];
	memcpy(eigenvalues, A, sizeof A);
	
	puts("A:");
	Matrix_print(A);
	puts("Матрица собственных векторов:");
	Eigenvectors_print(eigenvectors);

	Jacobi_method(eigenvalues, eigenvectors);

	puts("A:");
	Matrix_print(eigenvalues);
	puts("Собственные вектора:");
	Eigenvectors_print(eigenvectors);

	puts("Check");
	for (int i = 0; i < DIM; ++i){
		const double * const qi = eigenvectors+DIM*i;
		const double lambda = eigenvalues[get_index(i, i)];
		printf("i = %d, lambda = %lf\n", i, lambda);
		printf("               q_i = ");
		print_array(qi);
		double Aqi[DIM] = {};
		
		double lqi[DIM] = {};

		int index = 0;
		Aqi[0] = A[0] * qi[0];
		++index;
		
		// для векторизации операций
		double uvector[DIM];
		memset(uvector, 0, sizeof uvector);

		for (int r = 1; r < DIM; ++r){
			for (int j = 0; j < r; ++j){
				Aqi[r] += A[index] * qi[j];
				uvector[j] += A[index] * qi[r];
				++index;
			}
			assert(index == get_index(r, r));
			Aqi[r] += A[index] * qi[r];
			++index;
		}
		for (int j = 0; j < DIM; ++j)
			Aqi[j] += uvector[j];

		printf("             A*q_i = ");
		print_array(Aqi);
		
		for (int j = 0; j < DIM; ++j)
			lqi[j] += lambda*qi[j];

		printf("        lambda*q_i = ");
		print_array(lqi);

		for (int j = 0; j < DIM; ++j)
			Aqi[j] -= lqi[j];

		printf("A*q_i - lambda*q_i = ");
		print_array(Aqi);
	}
	
	return 0;
}
