/* 
   precision.c
   @author Gerardo Veltri
   Summarize differences in outcome and expected values of matrices
   introduced by precision errors  
*/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mem.h>
#include <matrix.h>

void identityPrecision(Matrix matrix, double *stats)
{
	double curr;
	double max = 0;
	double sum = 0;
	for (int i=0;i<matrix->n;i++)
	{
		for (int j=0;j<matrix->m;j++)
		{
			if (i==j)
				curr = fabs(matrix->values[i][j] - 1.0);
			else
				curr = fabs(matrix->values[i][j]);
			max = curr > max ? curr : max;
			sum = sum + curr;
		}
	}
	sum = sum / (matrix->n * matrix->m);
	stats[0] = sum;
	stats[1] = max;
}

void matrixComparison(Matrix matrix1, Matrix matrix2, double *stats)
{
	assert(matrix1->n == matrix2->n);
	assert(matrix1->m == matrix2->m);

	double curr;
	double max = 0;
	double sum = 0;
	for (int i=0;i<matrix1->n;i++)
	{
		for (int j=0;j<matrix1->m;j++)
		{
			curr = fabs(matrix1->values[i][j] - matrix2->values[i][j]);
			max = curr > max ? curr : max;
			sum = sum + curr;
		}
	}
	sum = sum / (matrix1->n * matrix1->m);
	stats[0] = sum;
	stats[1] = max;
}
