/*
  @file factorization.c
  @author Gerardo Veltri
  Factorization and Related Algorithms
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mem.h>
#include <matrix.h>

#define min(a,b)				\
	({ __typeof__ (a) _a = (a);		\
		__typeof__ (b) _b = (b);	\
		_a < _b ? _a : _b; })

#define MAXIMUM_ZERO_DOUBLE 0.00000000000001

/*
  QR Gram Schmidt Process on a Square Matrix
  
  less numerically stable than QR with householder reflections

  @param A matrix to be decomposed
  @param QR array of matrices, [Q,R], to which results are written
  @param debug flag for printing matrices during iterations
*/
void gramSchmidtQR(Matrix A, Matrix QR[2], int debug)
{
	assert(A->n == QR[0]->n);
	assert(QR[0]->n == QR[0]->m);
	assert(A->n == QR[1]->n);
	assert(A->m == QR[1]->m);
       
	Matrix Q = QR[0];
	Matrix R = QR[1];
	copyMatrix(A, Q);

	/* get orthonormal basis */
	for (int i=0; i<A->m; i++)
	{
		if (debug)
		{
			printf("ITERATION %d", i);
			drawMatrix(Q);
		}

		for (int j=0; j<i; j++)
		{
			/* scale projection and add to matrix */
			project(Q, i, Q, j, -1.0, Q, i, 1.0);

			if (debug)
			{
				printf("ITERATION %d, %d", i, j);
				drawMatrix(Q);
			}
		}

		normalizeColumn(Q, i);

		if (debug)
		{
			drawMatrix(Q);
			printf("ITERATION END %d\n", i);
		}
	}

	/* get R from orthogonal matrix Q transpose times A */
	multiplyMatrices(Q, 1, A, 0, R, 0);
}

/*
  QR with Householder Reflections on an NxM Matrix

  Memory Allocation
  -----------------

  input three matrices A nxm, Q nxn, R nxm
  four matrices for each step Qti
  - identity, nxn
  - outer product, nxn
  - Q transpose, nxn
  - Q chained transpose, nxn
  two matrices size (n,1) for v and x

  @param A matrix to be decomposed
  @param QR array of matrices, [Q,R], to which results are written
  @param mem_stacks array of two memory stacks for recyling scratch matrices
  @param debug flag for printing matrices during iterations

*/
void _hhReflectionsQR(Matrix A, Matrix QR[2],
		      MatrixStack mem_stacks[2],
		      int debug)
{

	assert(A->n == QR[0]->n);
	assert(QR[0]->n == QR[0]->m);
	assert(A->n == QR[1]->n);
	assert(A->m == QR[1]->m);

	Matrix Qn, Qt;
	MatrixStack stackNxNd4 = mem_stacks[0];
	MatrixStack stackNx1d2 = mem_stacks[1];
	Matrix v  = popMatrixStack(stackNx1d2);
	Matrix x  = popMatrixStack(stackNx1d2);

	Matrix Q = popMatrixStack(stackNxNd4);
	setMatrixValues(1, 'I', Q);

	Matrix I = popMatrixStack(stackNxNd4); /* make I type of Matrix with low mem usage */
	setMatrixValues(1, 'I', I);

	copyMatrix(A, QR[1]);

	/*
	  iterate through all columns if matrix is tall
	  iterate through all m-1 if matrix is wide or square
	*/
	int iterations = min(A->n-1, A->m);
	for (int i=0; i<iterations; i++)
	{

		/* get householder vector */
		for (int j=0; j<A->n; j++)
		{
			if (j < i)
				x->values[j][0] = 0;
			else
				x->values[j][0] = QR[1]->values[j][i];
		}
		setMatrixValues(0, 'V', v);
		v->values[i][0] = norm('C', x, 0);

		if (debug)
		{
			printf("norm=%.10f\n", v->values[i][0]);
			printf("aii=%.10f\n", v->values[i][0]);
		}

		/* reverse sign for better precision */
		if (x->values[i][0] > 0)
			v->values[i][0] = v->values[i][0] * -1;

		subtractColumn(x, 0, v, 0);

		if (debug)
		{
			printf("ITERATION %d\n", i);
			printf("householder vector=\n");
			drawMatrix(x);
		}

		Qn = popMatrixStack(stackNxNd4);
		outerMatrix(x, 0, Qn);
		double dot_product_hh = dotProductV(x,x);
		if (dot_product_hh != 0)
			scaleMatrix(Qn, -2 / dot_product_hh);
		addMatrix(Qn, I);
		if (debug)
		{
			printf("I - ((2/vTv)vvT)=\n");
			drawMatrix(Qn);
		}

		Qt = popMatrixStack(stackNxNd4);
		simpleMultiplyMatrices(Qn, Q, Qt);
		pushMatrixStack(stackNxNd4, Q);
		pushMatrixStack(stackNxNd4, Qn);
		Q = Qt;

		if (debug)
		{
			printf("Q%d(Q..)=\n", i);
			drawMatrix(Q);
			printf("ITERATION %d END\n", i);
		}

	}

	simpleMultiplyMatrices(Q, A, QR[1]);
	transposeMatrix(Q, QR[0]);
}

/*
  QR with Householder Reflections on a Symmetric Matrix

  Allocates and frees memory stacks. It's better to use this version
  if you're not iteratively calling QR.
  
  @param A matrix to be decomposed
  @param QR array of matrices, [Q,R], to which results are written
  @param mem_stacks array of two memory stacks for recyling scratch matrices
  @param debug flag for printing matrices during iterations

*/
void hhReflectionsQR(Matrix A, Matrix QR[2], int debug)
{
	MatrixStack mem_stacks[] = {
		allocMatrixStack(A->n, A->n, 4),
		allocMatrixStack(A->n, 1, 2)
	};

	_hhReflectionsQR(A, QR, mem_stacks, debug);

	freeMatrixStackAll(mem_stacks[0]);
	freeMatrixStackAll(mem_stacks[1]);
}

/*
  Gaussian Elimination

  Given A and B, generate the reduced row echelon form (RREF)
  Pivoting for improved precision

  @param A left side matrix to be reduced
  @param B left side of block
  @param RREF array of matrices to write reduced row form of A and B, [A^,B^]
  @param debug flag for printing matrices during iterations

*/
void gaussianElimination(Matrix A, Matrix B, Matrix RREF[2], int debug)
{

	Matrix _A = RREF[0];
	copyMatrix(A, _A);
	Matrix _B = RREF[1];
	copyMatrix(B, _B);

	double max_pivot_value, scalar;
	int max_pivot_index;
	int iterations = min(A->n-1, A->m);

	/* eliminate lower triangle */
	for (int i=0; i<iterations; i++)
	{
		if (debug)
		{
			printf("ITERATION %d BEGIN\n", i);
		}
    
		// pivot
		max_pivot_value = _A->values[i][i];
		max_pivot_index = i;
		for (int j=i+1; j<_A->n; j++)
		{
			if (fabs(_A->values[j][i]) > max_pivot_value)
			{
				max_pivot_value = fabs(_A->values[j][i]);
				max_pivot_index = j;
			}
		}
    
		if (max_pivot_index != 0)
		{
			switchRow(_A, i, max_pivot_index);
			switchRow(_B, i, max_pivot_index);
		}

		if (debug)
		{
			printf("row %d (%.5f) -> row %d\n", max_pivot_index, max_pivot_value, i);
			drawMatrix(_A);
		}

		for (int j=i+1; j<A->n; j++)
		{
			scalar = _A->values[j][i] / _A->values[i][i] * -1;
			if (_A->values[j][i] != 0)
			{
				addRowScalarMultiple(_A, j, scalar, _A, i);
				addRowScalarMultiple(_B, j, scalar, _B, i);
			}
		}

		if (debug)
		{
			printf("_A=\n");
			drawMatrix(_A);
			printf("_B=\n");
			drawMatrix(_B);
			printf("ITERATION %d END\n", i);
		}
	}

	/* eliminate upper triangle  */
	for (int i=iterations; i>0; i--)
	{
		for (int j=i-1; j>=0; j--)
		{
			scalar = _A->values[j][i] / _A->values[i][i] * -1;
			if (_A->values[j][i] != 0)
			{
				addRowScalarMultiple(_A, j, scalar, _A, i);
				addRowScalarMultiple(_B, j, scalar, _B, i);
			}
		}
	}

	/* convert to ones */
	for (int i=0; i<_A->n; i++)
	{
		if (_A->values[i][i] != 0)
		{
			scalar = 1/_A->values[i][i];
			scaleRow(_A, i, scalar);
			scaleRow(_B, i, scalar);
		}
	}

}

/*
  Back Substitution

  solves Ax = b for x where A and b are given
  
  @param A an upper triangular matrix
  @param b a column vector of values
  @param solution a column vector, x of Ax=b
 */
void backSubstitution(Matrix A, Matrix solution, Matrix b)
{
	assert(1 == solution->m);
	assert(1 == b->m);
	assert(A->n == b->n);
	assert(A->n == solution->n);
	
	double value, diagonal;
	for (int i=A->n; i>=0; i--)
	{
		value = b->values[i][0];
		for (int j=i+1; i<A->n; i++)
		{
			value = value - (A->values[i][j] * solution->values[j][0]);
		}

		diagonal = A->values[i][i];
		if ((fabs(diagonal) < MAXIMUM_ZERO_DOUBLE) & (fabs(value) > MAXIMUM_ZERO_DOUBLE))
		{
			fprintf(stderr,
				"contradiction A[%d][%d] = %.16f and b[%d] = %.16f",
				i,i,diagonal,i,value);
			exit(EXIT_FAILURE);
		}
		
		solution->values[i][0] = value / diagonal;
	}
}
