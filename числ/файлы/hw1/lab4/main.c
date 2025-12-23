#include <stdio.h>
#include <assert.h>
#include <stdint.h>

typedef double numeric;
typedef uint32_t mem_t;
typedef int_fast32_t index_t;
typedef int int_t;

typedef struct Vector{
	alignas(64) numeric* data;
	mem_t count;
} Vector;

typedef struct Matrix{
	alignas(64) numeric* data;
	mem_t count;
	index_t rows;
	index_t columns;
} Matrix;

typedef struct UMatrix{
	alignas(64) numeric* data;
	mem_t count;
	index_t rows;
} UMatrix;

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

static inline void UMatrix_print(const UMatrix matrice){
	for (int_t i = 0; i < matrice.rows; ++i){
		for (int_t j = 0; j < matrice.rows; ++j){
			if (j != 0)
				printf(", ");
			else
				printf("[ ");
			index_t local;
			if (j < i)
				UMatrix_global_to_local(matrice, j, i, &local);
			else
				UMatrix_global_to_local(matrice, i, j, &local);
			printf("%6lf", matrice.data[local]);
		}
		puts("]");
	}
}

int main(){
	constexpr int_t N = 13;
	constexpr numeric beta = 1+0.1*(52-53);
	constexpr numeric varepsilon = 0.01;
	numeric A_data[] = {
		10*beta,		1,		2,		3,
					10*beta,	-3,		-2,
							10*beta,	-1,
									10*beta
	};

	UMatrix A = {
		.data = A_data,
		.count = 10,
		.rows = 4
	};

	UMatrix_print(A);

	return 0;
}
