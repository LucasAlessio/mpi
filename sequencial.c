#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define ITERACOES 1000

double *matriz;
double *resultado;

long double ERRO = 1e-6;

long double norma_vetor(long double *x, int n) {
	int i;
	long double soma = 0;
	
	for (i = 0; i < n; i++) soma += pow(x[i], 2);
	
	return sqrt(soma);
}


long double * jacobi(long double * A, long double * b, int n) {
	long double *x = (long double*)calloc(n, sizeof(long double));
	long double *novo_x = (long double*)calloc(n, sizeof(long double));
	
	int iter, i, j;
	long double soma;
	
	for (iter = 0; iter < ITERACOES; iter++) {		
		for (i = 0; i < n; i++) {
			soma = 0;
			for (j = 0; j < n; j++) {
				if (i != j) {
					soma += A[i * n + j] * x[j];
				}
				//printf("A[%d][%d] %Lf", i, j, A[i * n + j]);
			}
			//printf("\n");
			
			novo_x[i] = (b[i] - soma) / A[i * n + i];
		}
		
		if (fabs(norma_vetor(x, n) - norma_vetor(novo_x, n)) < ERRO) {
			free(x);
			return novo_x;
		}
		else {
			for (i = 0; i < n; i++) {
				x[i] = novo_x[i];
			}
		}
	}
	
	free(novo_x);
	return x;
}

int main() {
	int n = 4;

	long double * A = (long double*)malloc(n * n * sizeof(long double));
	long double * b = (long double*)malloc(4 * sizeof(long double));

	A[0] = 10;
	A[1] = -1;
	A[2] = 2;
	A[3] = 0;
	
	A[4] = -1;
	A[5] = 11;
	A[6] = -1;
	A[7] = 3;

	A[8] = 2;
	A[9] = -1;
	A[10] =  10;
	A[11] = -1;
	 
	A[12] = 0;
	A[13] = 3;
	A[14] = -1;
	A[15] = 8;
	
	b[0] = 6;
	b[1] = 25;
	b[2] = -11;
	b[3] = 15;

	long double * resultado = jacobi(A, b, n);	
	
	int i = 4;
	for (i = 0; i < n; i++) {
		printf("%Lf\n", resultado[i]);
	}

	free(resultado);
	
	free(A);
	free(b);

	return 0;
}
