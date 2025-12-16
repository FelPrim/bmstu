#include <stdio.h>

typedef double numeric;
typedef uint32_t mem_t;
typedef int_fast32_t index_t;
typedef int int_t;

typedef struct Vector{
	alignas(64) numeric* data;
	mem_t count;
	mem_t capacity;
} Vector;

typedef struct Matrix{
	alignas(64) numeric* data;
	mem_t count;
	mem_t capacity;
	index_t rows;
	index_t columns;
} Matrix;

typedef struct UMatrix{
	alignas(64) numeric* data;
	mem_t count;
	mem_t capacity;
	index_t rows;
} UMatrix;



int main(){
	constexpr int_t N = 13;
	constexpr numeric beta = 1+0.1*(52-53);
	constexpr numeric varepsilon = 0.01;
	UMatrix A = {
		.data = {10*beta, 1, 2, 3, 10}
	};
	return 0;
}
