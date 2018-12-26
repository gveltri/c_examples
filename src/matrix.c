/*
  @file matrix.c
  @author Gerardo Veltri
  Base matrix functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mem.h>
#include <matrix.h>

/* constants for rendering tables */
const int PADDING = 1;
const int PRECISION = 3;
const char *FORMATTING  = "%.3f";

void setMatrixValues(double value, char type, Matrix matrix)
{
	for (int i=0;i<matrix->n;i++)
	{
		for (int j=0;j<matrix->m;j++)
		{
			/* base matrix types */
			switch (type) {
			case 'V':
				matrix->values[i][j] = value;
				break;
			case 'I':
				if (i == j)
					matrix->values[i][j] = value;
				else
					matrix->values[i][j] = 0;
				break;
			case 'U':
				if (i <= j)
					matrix->values[i][j] = value;
				else
					matrix->values[i][j] = 0;
				break;
			case 'L':
				if (i >= j)
					matrix->values[i][j] = value;
				else
					matrix->values[i][j] = 0;
				break;
			case 'R':
				matrix->values[i][j] = (double)rand()/(double)(RAND_MAX/value);
				break;
			}
		}
	}
}

void copyMatrix(Matrix source, Matrix target)
{
	for (int i=0;i<source->n;i++)
	{
		for (int j=0;j<source->m;j++)
		{
			target->values[i][j] = source->values[i][j];
		}
	}
}

void copyColumn(Matrix source, int idx_s, Matrix target, int idx_t)
{
	for (int i=0;i<target->n;i++)
	{
		target->values[i][idx_t] = source->values[i][idx_s];
	}
}

void copyRow(Matrix source, int idx_s, Matrix target, int idx_t)
{
	for (int i=0;i<source->m;i++)
	{
		target->values[idx_t][i] = source->values[idx_s][i];
	}
}

void switchRow(Matrix matrix, int row1, int row2)
{
	double scratch[matrix->n];
	for (int i=0; i<matrix->m; i++)
	{
		scratch[i] = matrix->values[row1][i];
	}

	for (int i=0; i<matrix->m; i++)
	{
		matrix->values[row1][i] = matrix->values[row2][i];
	}

	for (int i=0; i<matrix->m; i++)
	{
		matrix->values[row2][i] = scratch[i];
	}
}

void transposeMatrix(Matrix source, Matrix target)
{
	for (int i=0;i<target->n;i++)
	{
		for (int j=0;j<target->m;j++)
		{
			target->values[i][j] = source->values[j][i];
		}
	}
}

double matrixMax(Matrix matrix, int _abs)
{
	double max, curr;

	for (int i=0; i < matrix->n; i++)
	{
		for (int j=0; j < matrix->m; j++)
		{
			if ((i==0) & (j==0))
			{
				max = matrix->values[0][0];
				if (_abs)
					max = fabs(max);
			}

			curr = matrix->values[i][j];
			if (_abs) curr = fabs(curr);

			if (max < curr) max = curr;
		}
	}
	return max;
}

double sumMatrix(Matrix matrix)
{
	double sum = 0;
	for (int i=0;i<matrix->n;i++)
	{
		for (int j=0;j<matrix->m;j++)
		{
			sum = matrix->values[i][j];
		}
	}
	return sum;
}

double meanMatrix(Matrix matrix)
{
	double sum = sumMatrix(matrix);
	return sum / (matrix->n * matrix->m);
}


int numDigits(int n)
{
	int add = 0;
	if (n < 0)
	{
		n = n * -1;
		add++;
	}

	if (n == 0)
		n = 1;
	else
		n = floor(log10(n)) + 1;

	return n+add;
}

int numDigitsDouble(double x)
{
	int n = floor(x);

	n = numDigits(n);

	return n;
}

void drawMatrix(Matrix matrix)
{
	double matrix_max = matrixMax(matrix, 1);
	int max_width = numDigitsDouble(matrix_max);
	max_width = max_width + PADDING + PRECISION;

	/* header */
	int offset = numDigits(matrix->n) + numDigits(matrix->m);
	int blank_fill;

	printf("\n");
	printf("(%d, %d) ", matrix->n, matrix->m);
	for (int i=offset; i<(max_width*matrix->m); i++)
	{
		printf("-");
	}

	for (int i=0; i<matrix->n; i++)
	{
		printf("\n");
		for (int j=0; j<max_width/2; j++)
		{
			printf("\n");
		}
		for (int j=0; j<matrix->m; j++)
		{
			double value = matrix->values[i][j];
			printf(FORMATTING, value);

			blank_fill = max_width - numDigitsDouble(value);

			for (int i=0; i<blank_fill; i++)
			{
				printf(" ");
			}
		}
	}

	printf("\n");
	for (int i=0; i<(max_width*matrix->m); i++)
	{
		printf("-");
	}
	for (int i=0; i<3; i++) printf("\n");
}

/* scaling operations */
void scaleColumn(Matrix matrix, int idx, double scalar)
{
	for (int i=0;i<matrix->n;i++)
	{
		matrix->values[i][idx] = (matrix->values[i][idx]) * scalar;
	}
}

void scaleRow(Matrix matrix, int idx, double scalar)
{
	for (int i=0;i<matrix->n;i++)
	{
		matrix->values[idx][i] = (matrix->values[idx][i]) * scalar;
	}
}

void scaleMatrix(Matrix matrix, double scalar)
{
	for (int i=0;i<matrix->m;i++)
	{
		scaleColumn(matrix, i, scalar);
	}
}

void absMatrix(Matrix matrix)
{
	for (int i=0; i<matrix->n; i++)
	{
		for (int j=0; j<matrix->m; j++)
		{
			matrix->values[i][j] = fabs(matrix->values[i][j]);
		}
	}
}

/* addition and subtraction */
void addColumn(Matrix target, int idx1, Matrix source, int idx2)
{
	assert(source->n == target->n);
	for (int i=0; i<target->n; i++)
	{
		target->values[i][idx1] = target->values[i][idx1] + source->values[i][idx2];
	}
}

void addMatrix(Matrix target, Matrix source)
{
	assert(source->m == target->m);
	for (int i=0; i<target->m; i++)
	{
		addColumn(target, i, source, i);
	}
}

void subtractRow(Matrix target, int idx1, Matrix source, int idx2)
{
	assert(source->m == target->m);
	for (int i=0; i<target->m; i++)
	{
		target->values[idx1][i] = target->values[idx1][i] - source->values[idx2][i];
	}
}


void subtractColumn(Matrix target, int idx1, Matrix source, int idx2)
{
	assert(source->n == target->n);
	for (int i=0; i<target->n; i++)
	{
		target->values[i][idx1] = target->values[i][idx1] - source->values[i][idx2];
	}
}

void subtractMatrix(Matrix target, Matrix source) {
	assert(source->m == target->m);
	for (int i=0; i<target->m; i++)
	{
		subtractColumn(target, i, source, i);
	}
}

/* scale and add */
void addRowScalarMultiple(Matrix target, int idx_t, double scalar, Matrix source, int idx_s)
{
	assert(source->m == target->m);
	for (int i=0; i<target->m; i++)
	{
		target->values[idx_t][i] = (source->values[idx_s][i] * scalar) + \
			target->values[idx_t][i];
	}
}

/* dot product */
double dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2)
{
	double x = 0.0;

	switch (orient)
	{
	case 'R':
		assert(matrix1->m == matrix2->m);
		for (int i=0; i<matrix1->m; i++)
		{
			x = x + (matrix1->values[idx1][i] * matrix2->values[idx2][i]);
		}
		break;
	case 'C':
		assert(matrix1->n == matrix2->n);
		for (int i=0; i<matrix1->n; i++)
		{
			x = x + (matrix1->values[i][idx1] * matrix2->values[i][idx2]);
		}
		break;
	}

	return x;
}

double dotProductV(Matrix matrix1, Matrix matrix2)
{
	double x = dotProduct('C', matrix1, 0, matrix2, 0);
	return x;
}

/* norm and normalization */
double norm(char orient, Matrix matrix, int idx)
{
	double x = sqrt(dotProduct(orient, matrix, idx, matrix, idx));
	return x;
}

double normV(Matrix matrix) {
	double x = norm('C', matrix, 0);
	return x;
}

void normalizeColumn(Matrix matrix, int idx) {
	double _norm = norm('C', matrix, idx);

	/* if norm is zero, zero vector */
	if (_norm != 0) 
		scaleColumn(matrix, idx, 1/_norm);
}

/* 
    project

    project column vector of source 1 onto source 2 and save to column target[idx_t]

    target[idx_t] <- proj_{source2[idx2]}source1[idx1]
*/
void simpleProject(Matrix source1, int idx1, Matrix source2, int idx2,
             Matrix target, int idx_t)
{
	assert(source1->n == source2->n);
	assert(source2->n == target->n);

	copyColumn(source2, idx2, target, idx_t);
	double st_dot = dotProduct('C', source1, idx1, source2, idx2);
	double norm_squared = dotProduct('C', source2, idx2, source2, idx2);

	/* projection of 0 vector onto any vector is the 0 vector */
	if (norm_squared == 0)
		return;

	st_dot = st_dot / norm_squared;
	scaleColumn(target, idx_t, st_dot);
}

/* 

   add scaled projection to matrix
   a faster, more efficient implementation of projection with scalar addition

   target[idx_t] <- (proj_scalar * proj_{source2[idx2]}source1[idx1]) + 
                    tscalar * target[idx_t]

*/
void project(Matrix source1, int idx1, Matrix source2, int idx2, double proj_scalar,
             Matrix target, int idx_t, double tscalar)
{
	assert(source1->n == source2->n);
	assert(source2->n == target->n);
	
	double st_dot = dotProduct('C', source1, idx1, source2, idx2);
	double norm_squared = dotProduct('C', source2, idx2, source2, idx2);

	if (norm_squared == 0)
		st_dot = 0;
	else
		st_dot = st_dot / norm_squared;
	
	for (int i=0; i<source1->n; i++)
	{
		target->values[i][idx_t] = (proj_scalar * source2->values[i][idx2] * st_dot) + \
			(tscalar * target->values[i][idx_t]);
	}
}

/*  
    outerMatrix

    outer product of self
*/
void outerMatrix(Matrix source, int idx_s, Matrix target)
{
	assert(source->n == target->n);
	assert(target->n == target->m);

	for (int i=0; i<target->n; i++)
	{
		for (int j=0; j<target->m; j++)
		{
			target->values[i][j] = source->values[i][idx_s] * source->values[j][idx_s];
		}
	}
}


/* 
    simpleMultiplyMatrices

    simplified version of multiplyMatrices
*/
void simpleMultiplyMatrices(Matrix source1, Matrix source2, Matrix target)
{
	assert(source1->m == source2->n);
	assert(source1->n == target->n);
	assert(source2->m == target->m);
	for (int i=0; i<source1->n; i++)
	{
		for (int j=0; j<source2->m; j++)
		{
			double value = 0;
			for (int k=0; k<source1->m; k++)
			{
				value = value + (source1->values[i][k] * source2->values[k][j]);
			}
			target->values[i][j] = value;
		}
	}
}


/* 
   multiplyMatrices

   target <- transpose1(source1)transpose2(source2) + (tscalar * target)

   multiply two matrices and write to target matrix
   parameters for transposing either source matrix and
   adding a scalar multiple of the target to the result
   inspired by GEMM of BLAS

 */
void multiplyMatrices(Matrix source1, int transpose1, Matrix source2, int transpose2,
		      Matrix target, double tscalar)
{
	int iterations;
	if (transpose1 & transpose2)
	{
		assert(source1->n == source2->m);
		assert(source1->m == target->n);
		assert(source2->n == target->m);

		iterations = source1->n;
	}
	else if (transpose1)
	{
		assert(source1->n == source2->n);
		assert(source1->m == target->n);
		assert(source2->m == target->m);

		iterations = source1->n;
	}
	else if (transpose2)
	{
		assert(source1->m == source2->m);
		assert(source1->n == target->n);
		assert(source2->n == target->m);

		iterations = source1->m;
	}
	else
	{
		assert(source1->m == source2->n);
		assert(source1->n == target->n);
		assert(source2->m == target->m);

		iterations = source1->m;
	}

	for (int i=0; i<target->n; i++)
	{
		for (int j=0; j<target->m; j++)
		{
			double value = 0;
			for (int k=0; k<iterations; k++)
			{
				if (transpose1 & transpose2)
					value = value			\
						+ (source1->values[k][i] * source2->values[j][k]);
				else if (transpose1)
					value = value			\
						+ (source1->values[k][i] * source2->values[k][j]);
				else if (transpose2)
					value = value			\
						+ (source1->values[i][k] * source2->values[j][k]);
				else
					value = value			\
						+ (source1->values[i][k] * source2->values[k][j]);
			}
			target->values[i][j] = value + (tscalar * target->values[i][j]);
		}
	}
}
