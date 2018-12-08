/*
  @file qr.c
  @author Gerardo Veltri
  QR Decomposition
*/
#include <stdio.h>
#include <math.h>
#include <mem.h>
#include <matrix.h>

/* Gram Schmidt Process on a Square Matrix */
void gramSchmidt(Matrix A, Matrix QR[2], int debug)

{
  Matrix cur_proj = allocMatrix(A->n, 1);
  Matrix Q = QR[0];
  Matrix R = QR[1];
  copyMatrix(A,Q);

  draw2DMatrix(A);

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

      project(A, i, Q, j, cur_proj, 0);
      subtractColumn(Q, i, cur_proj, 0);

      if (debug)
      {
        printf("ITERATION %d, %d", i, j);
        draw2DMatrix(cur_proj);
      }
    }

    _norm = norm('C', Q, i);

    if (_norm != 0) /* if norm is zero, zero vector? */
      scaleColumn(Q, i, 1/_norm);

    if (debug)
      {
        printf("%.10f", _norm);
        printf("ITERATION END %d\n", i);
        draw2DMatrix(Q);
      }

  }

  freeMatrix(cur_proj);
  Matrix Q_t = allocMatrix(A->m, A->n);
  transposeMatrix(Q, Q_t);
  draw2DMatrix(Q_t);
  multiplyMatrices(Q_t, A, R);
  freeMatrix(Q_t);
}

int main()
{
  MatrixStack stack = allocMatrixStack(4,4,4);
  Matrix A = popMatrixStack(stack);
  setMatrixValues(1, 'V', A);

  printf("A=\n");
  draw2DMatrix(A);

  Matrix QR[2] = {popMatrixStack(stack),
                  popMatrixStack(stack)};
  Matrix _A = popMatrixStack(stack);

  gramSchmidt(A, QR, 1);

  printf("Q=\n");
  draw2DMatrix(QR[0]);
  printf("R=\n");
  draw2DMatrix(QR[1]);

  multiplyMatrices(QR[0], QR[1], _A);
  printf("QR=\n");
  draw2DMatrix(_A);

  freeMatrixStack(stack);

  return 0;
}
