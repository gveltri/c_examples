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

  Memory Complexity
  -----------------

  input three matrices A, Q, R
  one scratch matrix of size n,1
  one scratch matrix of size n,n to store transpose

*/
void gramSchmidtQR(Matrix A, Matrix QR[2], int debug)
{
  assert(A->m == A->n);
  assert(A->m == QR[0]->m);
  assert(A->m == QR[0]->n);

  Matrix cur_proj = allocMatrix(A->n, 1);
  Matrix Q = QR[0];
  Matrix R = QR[1];
  copyMatrix(A,Q);

  /* get orthonormal basis */
  float _norm;
  for (int i=0; i<A->m; i++)
  {

    if (debug)
    {
      printf("ITERATION %d", i);
      draw2DMatrix(Q);
    }

    for (int j=0; j<i; j++)
    {

      project(Q, i, Q, j, cur_proj, 0);
      subtractColumn(Q, i, cur_proj, 0);

      if (debug)
      {
        printf("ITERATION %d, %d", i, j);
        draw2DMatrix(cur_proj);
      }
    }

    normalizeColumn(Q, i);

    if (debug)
    {
      printf("ITERATION END %d\n", i);
      draw2DMatrix(Q);
    }

  }

  freeMatrix(cur_proj);
  Matrix Q_t = allocMatrix(A->m, A->n);
  transposeMatrix(Q, Q_t); /* combine into one step with transposeMultiply */
  multiplyMatrices(Q_t, A, R);
  freeMatrix(Q_t);
}


/*
  QR with Householder Reflections on a Symmetric Matrix

  Memory Complexity
  -----------------

  input three matrices A, Q, R
  four matrices of size (n,n) for each step Qi
   - identity
   - outer product
   - Q transpose
   - Q result
  one matrix size (n,1) for v and vT

 */
void hhReflectionsQR(Matrix A, Matrix QR[2], int debug) {

  assert(A->m == A->n);
  assert(A->m == QR[0]->m);
  assert(A->m == QR[0]->n);

  Matrix Qn, Qt;
  MatrixStack stack = allocMatrixStack(A->n, A->n, 4);
  Matrix v  = allocMatrix(A->n, 1);

  Matrix Q = popMatrixStack(stack);
  setMatrixValues(1, 'I', Q);

  Matrix I = popMatrixStack(stack); /* make I type of Matrix with low mem usage */
  setMatrixValues(1, 'I', I);

  for (int i=0; i<A->m; i++)
  {
    setMatrixValues(0, 'V', v);
    v->values[i][0] = norm('C', A, i) * -1;

    if (debug)
    {
      printf("ITERATION %d\n", i);
      printf("norm = %.10f\n", v->values[i][0]);
      draw2DMatrix(v);
    }

    addColumn(v, 0, A, i);

    if (debug)
    {
      draw2DMatrix(v);
    }

    normalizeColumn(v, 0);

    if (debug)
    {
      draw2DMatrix(v);
    }

    Qn = popMatrixStack(stack);
    outerMatrix(v, 0, Qn);
    if (debug)
    {
        draw2DMatrix(Qn);
    }

    scaleMatrix(Qn, -2); /* this step earlier to reduce complexity */
    addMatrix(Qn, I);
    if (debug)
    {
      draw2DMatrix(Qn);
    }

    Qt = popMatrixStack(stack);
    transposeMatrix(Qn, Qt);

    multiplyMatrices(Q, Qt, Qn);
    pushMatrixStack(stack, Q);
    Q = Qn;

    pushMatrixStack(stack, Qt);

    if (debug)
    {
      draw2DMatrix(Q);
      printf("ITERATION %d END\n", i);
    }

  }

  copyMatrix(Q, QR[0]);

  transposeMatrix(QR[0], Q);
  multiplyMatrices(Q, A, QR[1]);

  freeMatrixStack(stack);
  freeMatrix(v);
}

int main()
{
  MatrixStack stack = allocMatrixStack(4,4,4);
  Matrix A = popMatrixStack(stack);
  setMatrixValues(5, 'R', A);

  printf("A=\n");
  draw2DMatrix(A);

  Matrix QR[2] = {popMatrixStack(stack),
                  popMatrixStack(stack)};
  Matrix _A = popMatrixStack(stack);

  hhReflectionsQR(A, QR, 1);

  printf("Q=\n");
  draw2DMatrix(QR[0]);
  printf("R=\n");
  draw2DMatrix(QR[1]);

  multiplyMatrices(QR[0], QR[1], _A);
  printf("QR=\n");
  draw2DMatrix(_A);

  subtractMatrix(A, _A);
  absMatrix(A);
  float mean = meanMatrix(A);
  printf("Q - QR=\n");
  draw2DMatrix(A);

  printf("Mean Error=\n");
  printf("%.20f\n", mean);

  freeMatrixStack(stack);

  return 0;
}
