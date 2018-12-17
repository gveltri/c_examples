/*
  @file factorization.c
  @author Gerardo Veltri
  Factorization and Related Algorithms
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mem.h>
#include <matrix.h>

/*
  QR Gram Schmidt Process on a Square Matrix

  Memory Allocation
  -----------------

  input three matrices A, Q, R
  one scratch matrix of size n,1
  one scratch matrix of size n,n to store transpose

*/
void _gramSchmidtQR(Matrix A, Matrix QR[2], MatrixStack mem_stacks[2], int debug)
{
  assert(A->m == A->n);
  assert(A->m == QR[0]->m);
  assert(A->m == QR[0]->n);

  MatrixStack stackMxNd1 = mem_stacks[0];
  MatrixStack stackNx1d1 = mem_stacks[1];

  Matrix Q_t = popMatrixStack(stackMxNd1);
  Matrix cur_proj = popMatrixStack(stackNx1d1);

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

      project(Q, i, Q, j, cur_proj, 0);
      subtractColumn(Q, i, cur_proj, 0);

      if (debug)
      {
        printf("ITERATION %d, %d", i, j);
        drawMatrix(cur_proj);
      }
    }

    normalizeColumn(Q, i);

    if (debug)
    {
      printf("ITERATION END %d\n", i);
      drawMatrix(Q);
    }

  }

  transposeMatrix(Q, Q_t); /* combine into one step with transposeMultiply */
  multiplyMatrices(Q_t, A, R);

}

void gramSchmidtQR(Matrix A, Matrix QR[2], int debug)
{
  MatrixStack mem_stacks[] = {
    allocMatrixStack(A->m, A->n, 1),
    allocMatrixStack(A->n, 1, 1),
  };

  _gramSchmidtQR(A, QR, mem_stacks, debug);

  freeMatrixStackAll(mem_stacks[0]);
  freeMatrixStackAll(mem_stacks[1]);
}


/*
  QR with Householder Reflections on a Symmetric Matrix

  Memory Allocation
  -----------------

  input three matrices A, Q, R
  four matrices of size (n,n) for each step Qi
   - identity
   - outer product
   - Q transpose
   - Q result
  two matrices size (n,1) for v and x

 */
void _hhReflectionsQR(Matrix A, Matrix QR[2],
		     MatrixStack mem_stacks[2],
		     int debug) {

  assert(A->m == A->n);
  assert(A->m == QR[0]->m);
  assert(A->m == QR[0]->n);

  Matrix Qn, Qt;
  MatrixStack stackNxNd4 = mem_stacks[0];
  MatrixStack stackNx1d2 = mem_stacks[1];
  Matrix v  = popMatrixStack(stackNx1d2);
  Matrix x  = popMatrixStack(stackNx1d2);

  Matrix Q = popMatrixStack(stackNxNd4);
  setMatrixValues(1, 'I', Q);

  Matrix I = popMatrixStack(stackNxNd4); /* make I type of Matrix with low mem usage */
  setMatrixValues(1, 'I', I);

  copyMatrix(A,QR[1]);

  /*
    iterate through all but last column
    for symmetric matrix
   */
  for (int i=0; i<A->m-1; i++)
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
    multiplyMatrices(Qn, Q, Qt);
    pushMatrixStack(stackNxNd4, Q);
    pushMatrixStack(stackNxNd4, Qn);
    Q = Qt;

    multiplyMatrices(Q, A, QR[1]);

    if (debug)
    {
      printf("Q%d(Q..)=\n", i);
      drawMatrix(Q);
      printf("Q%d(Q..) * A=\n", i);
      drawMatrix(QR[1]);
      printf("ITERATION %d END\n", i);
    }

  }

  transposeMatrix(Q, QR[0]);
}

void hhReflectionsQR(Matrix A, Matrix QR[2],
		     int debug)
{
  MatrixStack mem_stacks[] = {
    allocMatrixStack(A->n, A->m, 4),
    allocMatrixStack(A->n, 1, 2)
  };

  _hhReflectionsQR(A, QR, mem_stacks, debug);

  freeMatrixStackAll(mem_stacks[0]);
  freeMatrixStackAll(mem_stacks[1]);
}


/*
  Gaussian Elimination

*/
void gaussianElimination(Matrix A, Matrix B, Matrix RREF[2], int debug)
{

  Matrix _A = RREF[0];
  copyMatrix(A, _A);
  Matrix _B = RREF[1];
  copyMatrix(B, _B);

  Matrix a_scratch = allocMatrix(1, _A->m);
  Matrix b_scratch = allocMatrix(1, _B->m);

  double max_pivot_value, scalar;
  int max_pivot_index;
  for (int i=0; i<_A->m; i++)
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

    if (debug)
    {
      printf("row %d (%.5f) -> row %d\n", max_pivot_index, max_pivot_value, i);
    }
    if (max_pivot_index != 0)
    {
      switchRow(_A, 0, max_pivot_index);
      switchRow(_A, 0, max_pivot_index);
    }

    for (int j=i+1; j<_A->n-i; j++)
    {
      if (_A->values[j][i] != 0)
      {
        scalar = _A->values[j][i] / _A->values[i][i];
        copyRow(_A, i, a_scratch, 0);
        scaleMatrix(a_scratch, scalar);
        subtractRow(_A, j, a_scratch, 0);

        copyRow(_B, i, b_scratch, 0);
        scaleMatrix(b_scratch, scalar);
        subtractRow(_B, j, b_scratch, 0);
      }
    }
    if (debug)
    {
      printf("ITERATION %d END\n", i);
      drawMatrix(_A);
      drawMatrix(_B);
    }
  }
  freeMatrix(b_scratch);
  freeMatrix(a_scratch);
}
