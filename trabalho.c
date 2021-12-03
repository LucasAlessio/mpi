#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

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


void jacobi(long double * A, long double * b, int n, int id, int nproc, long double *x) {
	//long double *x = (long double*)calloc(n, sizeof(long double));
	long double *novo_x = (long double*)calloc(n, sizeof(long double));
	
	int iter, i, j, ii;
	long double soma;
	
	int ini = id * (n / nproc);
	int fim = ini + (n / nproc);
	
	for (iter = 0; iter < ITERACOES; iter++) {
		for (ii = 0, i = ini; i < fim; ii++, i++) {
			soma = 0;
			for (j = 0; j < n; j++) {
				if (i != j) {
					soma += A[i * n + j] * x[j];
				}
			}

			novo_x[ii] = (b[i] - soma) / A[i * n + i];
		}
		
		if (fabsl(norma_vetor(x, n) - norma_vetor(novo_x, n)) < ERRO) {

			//printf("[%d] ITER %d\n", id, iter);
			for (i = 0; i < n; i++) {
				//printf("[%Lf]   ", novo_x[i]);
				x[i] = novo_x[i];
			}
			//printf("\n");
			
			free(novo_x);
			return;
			//return novo_x;
		}
		else {
			for (j = 0; j < n; j++) {
				x[j] = novo_x[j];
			}
		}
	}
	
	//printf("[%d] ITER %d", id, iter);
	free(novo_x);
	return;
}

int main(int argc, char **argv) {
	int id, np, n = 4, i, j;
	
	//printf("ID %d\n", id);
	//printf("NP %d\n", np);

	long double * A = (long double*)malloc(n * n * sizeof(long double));
	long double * b = (long double*)malloc(n * sizeof(long double));

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
	A[10] = 10;
	A[11] = -1;
	 
	A[12] = 0;
	A[13] = 3;
	A[14] = -1;
	A[15] = 8;
	
	b[0] = 6;
	b[1] = 25;
	b[2] = -11;
	b[3] = 15;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	long double * resultado = NULL;
	long double * buffer = NULL;
	long double * final = NULL;
	
	buffer = (long double*)calloc(n, sizeof(long double));	
	
	//if (id == 0) {
		resultado = (long double*)calloc(n, sizeof(long double));
		final = (long double*)calloc(n, sizeof(long double));
	//}
	
	
	MPI_Bcast(A, n * n, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(buffer, n, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	
		
	//long double * buffer = (long double*)calloc(n * n, sizeof(long double));

	jacobi(A, b, n, id, np, buffer);
	
	/* for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			printf("[%d] %Lf - ", i * 4 + j, buffer[i * 4 + j]);
		}
		printf(" <------------ \n");
	} */
	
	/*printf("%d = ", id);
	for (i = 0; i < n / np; i++) {
		printf("%Lf ", buffer[i]);
	}
	printf("\n");*/
		
//	MPI_Gather(buffer, n/np, MPI_LONG_DOUBLE, resultado, n/np, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Allreduce(&buffer, &resultado, n, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (id == 0) {
		int i = n;
		for (i = 0; i < n; i++) {
			printf("%Lf\n", resultado[i]);
		}
	}
	
	MPI_Finalize();
	
	/*free(resultado);
	free(A);
	free(b);
	free(buffer);*/

	return 0;
}
