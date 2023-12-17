/*
 * CSE 351 Lab 4 (Caches and Cache-Friendly Code)
 * Part 2 - Optimizing Matrix Transpose
 *
 * Name(s): Ian Tsai 
 * NetID(s): itsai2
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[M][N], int B[N][M]);
 * and compute B = A^T.
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1 KiB direct mapped cache with a block size of 32 bytes.
 */

#include <stdio.h>
#include "support/cachelab.h"

int is_transpose(int M, int N, int A[M][N], int B[N][M]);


/*
 * transpose_submit - This is the transpose function that you will be graded
 *     on. Do not change the description string "Transpose submission", as the
 *     driver searches for that string to identify the transpose function to be
 *     graded.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[M][N], int B[N][M]) {
    int i, j, ii, jj;

    if (M == 32 || N == 32) {
	    int temp;
	    int block_size = 8;

	    for (i = 0; i < M; i+= block_size) {
		    for (j = 0; j < N; j += block_size) {
			    for (ii = 0; ii < block_size; ii++) {
				    for (jj = 0; jj < block_size; jj++) {
					    B[j + (block_size - 1) - ii][i + jj] = A[i + ii][j +jj];
                    		    }
			    }

                	    for (ii = 0; ii < block_size / 2; ii++) {
                    		    for (jj = 0; jj < block_size; jj++) {
                        		    temp = B[j + ii][i + jj];
                        		    B[j + ii][i + jj] = B[j + (block_size - 1) - ii][i + jj];
                        		    B[j + (block_size - 1) - ii][i + jj] = temp;
                    		    }
                	    }
 
                	    for (ii = 0; ii < block_size; ii++) {
                    		    for (jj = 0; jj < block_size; jj++) {
                        		    if (jj > ii) {
						    temp = B[j + ii][i + jj];
                            		    	    B[j + ii][i + jj] = B[j + jj][i + ii];
                            		    	    B[j + jj][i + ii] = temp;
                        		    }
                    		    }
                	    }
            	    }
            }
    } else {
	    int temp1[4];
	    int temp2[4];

	    for (i = 0; i < M; i += 8) {
		    for (j = 0; j < N; j += 8) {
			    for (ii = 0; ii < 4; ii++) {
				    for (jj = 0; jj < 8; jj++) {
					    B[j + 3 - ii][i + jj] = A[i + ii][j + jj];
				    }
                            }

			    for (ii = 0; ii < 2; ii++) {
				    for (jj = 0; jj < 8; jj++) {
					    temp1[0] = B[j + ii][i + jj];
                                            B[j + ii][i + jj] = B[j + 3 - ii][i + jj];
                                            B[j + 3 - ii][i + jj] = temp1[0];
                    		    }
                	    }

			    for (ii = 0; ii < 4; ii++) {
				    for (jj = 0; jj < 4; jj++) {
					    if (jj > ii) {
						    temp1[0] = B[j + ii][i + jj];
                                            	    B[j + ii][i + jj] = B[j + jj][i + ii];
                                            	    B[j + jj][i + ii] = temp1[0];
					    }
				    }
			    }

                	    for (ii = 0; ii < 4; ii++) {
                                    for (jj = 4; jj < 8; jj++) {
                                            if (jj - 4 > ii) {
                                                    temp1[0] = B[j + ii][i + jj];
                                                    B[j + ii][i + jj] = B[j + jj - 4][i + ii + 4];
                                                    B[j + jj - 4][i + ii + 4] = temp1[0];
                                            }
                                    }
                            }

                	    for (ii = 0; ii < 2; ii++) {
                                    for (jj = 4; jj < 8; jj++) {
                                            temp1[0] = B[j + ii][i + jj];
                                            B[j + ii][i + jj] = B[j + 3 - ii][i + jj];
                                            B[j + 3 - ii][i + jj] = temp1[0];
                                    }
                            }
 
                	    for (ii = 4; ii < 8; ii++) {
				    for (jj = 0; jj < 8; jj++) {
					    B[j + 11 - ii][i + jj] = A[i + ii][j + jj];
                    		    }
                	    }

               		    for (ii = 0; ii < 2; ii++) {
                    		    for (jj = 0; jj < 8; jj++) {
                        		    temp1[0] = B[j + ii + 4][i + jj];
                        		    B[j + ii + 4][i + jj] = B[j + 7 - ii][i + jj];
                        		    B[j + 7 - ii][i + jj] = temp1[0];
                    		    }
                 	    }
                
                	    for (ii = 4; ii < 8; ii++) {
                    		    for (jj = 0; jj < 4; jj++) {
                        		    if (jj > ii - 4) {
						    temp1[0] = B[j + ii][i + jj];
                            			    B[j + ii][i + jj] = B[j + jj + 4][i + ii - 4];
                            			    B[j + jj + 4][i + ii - 4] = temp1[0];
                        		    }
                    		    }
                	    }

                	    for (ii = 4; ii < 6; ii++) {
                    		    for (jj = 0; jj < 4; jj++) {
                        		    temp1[0] = B[j + ii][i + jj];
                        		    B[j + ii][i + jj] = B[j + 11 - ii][i + jj];
                        		    B[j + 11 - ii][i + jj] = temp1[0];
                    		    }
                	    }

                	    for (ii = 4; ii < 8; ii++) {
                    		    for (jj = 4; jj < 8; jj++) {
                        		    if (jj > ii) {
                            			    temp1[0] = B[j + ii][i + jj];
                            			    B[j + ii][i + jj] = B[j + jj][i + ii];
                            			    B[j + jj][i + ii] = temp1[0];
                        		    }
                    		    }
                	    }
                
                	    for (ii = 0; ii < 4; ii++) {
                    		    for (jj = 0; jj < 4; jj++) {
                        		    temp1[jj] = B[j + ii][i + jj + 4];
                        		    temp2[jj] = B[j + 7 - ii][i + jj];
                        		    B[j + 7 - ii][i + jj] = temp1[jj];
                        		    B[j + ii][i + jj + 4] = temp2[jj];
                    		    }
                	    }
		    }
	    }
    }
}


// You can define additional transpose functions below. We've defined a simple
// one below to help you get started.

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[M][N], int B[N][M]) {
    int i, j, tmp;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }

}


/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions() {
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);

}


/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[M][N], int B[N][M]) {
    int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}
