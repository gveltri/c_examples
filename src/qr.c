/*
  @file qr.c
  @author Gerardo Veltri
  QR Decomposition
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
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
  double _norm;
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

  pushMatrixStack(stackNx1d1, cur_proj);
  pushMatrixStack(stackMxNd1, Q_t);
}

void gramSchmidtQR(Matrix A, Matrix QR[2], int debug)
{
  MatrixStack mem_stacks[] = {
    allocMatrixStack(A->m, A->n, 1),
    allocMatrixStack(A->n, 1, 1),
  };

  _gramSchmidtQR(A, QR, mem_stacks, debug);

  freeMatrixStack(mem_stacks[0]);
  freeMatrixStack(mem_stacks[1]);
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

  pushMatrixStack(stackNxNd4, Q);
  pushMatrixStack(stackNxNd4, I);

  pushMatrixStack(stackNx1d2,v);
  pushMatrixStack(stackNx1d2,x);
}

void hhReflectionsQR(Matrix A, Matrix QR[2],
		     int debug)
{
  MatrixStack mem_stacks[] = {
    allocMatrixStack(A->n, A->m, 4),
    allocMatrixStack(A->n, 1, 2)
  };
  
  _hhReflectionsQR(A, QR, mem_stacks, debug);

  freeMatrixStack(mem_stacks[0]);
  freeMatrixStack(mem_stacks[1]);
}





int main()
{
  MatrixStack stack = allocMatrixStack(5,5,4);
  
  Matrix A = popMatrixStack(stack);
  Matrix QR[2] = {popMatrixStack(stack),
                  popMatrixStack(stack)};
  Matrix _A = popMatrixStack(stack);
  
  setMatrixValues(10, 'R', A);

  printf("A=\n");
  drawMatrix(A);
  
  hhReflectionsQR(A, QR, 0);
  
  printf("Q=\n");
  drawMatrix(QR[0]);
  printf("R=\n");
  drawMatrix(QR[1]);

  multiplyMatrices(QR[0], QR[1], _A);
  printf("QR=\n");
  drawMatrix(_A);

  subtractMatrix(A, _A);
  absMatrix(A);
  double mean = meanMatrix(A);
  printf("Q - QR=\n");
  drawMatrix(A);

  printf("Mean Error=\n");
  printf("%.6f\n", mean);

  pushMatrixStack(stack, QR[0]);
  pushMatrixStack(stack, QR[1]);
  pushMatrixStack(stack, A);
  pushMatrixStack(stack, _A);
  freeMatrixStack(stack);

  return 0;
}
