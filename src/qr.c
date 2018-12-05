/*
  @file qr.c
  @author Gerardo Veltri
  QR Decomposition
*/
#include <stdio.h>
#include <matrix.h>
#include <stdlib.h>
#include <math.h>

/* Gram Schmidt Process on a Square Matrix */
Matrix *gramSchmidt(Matrix A, Matrix *QR)
{
  Matrix Q = copyMatrix(A);

  /* get orthonormal basis */
  float _norm;
  for (int i=0; i<A->m; i++)
  {
    printf("ITERATION %d", i);
    draw2DMatrix(Q);
    for (int j=0; j<i; j++)      
    {
      Matrix cur_proj;
      cur_proj = project(A, i, Q, j);
      subtractColumn(Q, i, cur_proj, 0);
      printf("ITERATION %d, %d", i, j);
      draw2DMatrix(cur_proj);
      freeMatrix(cur_proj);
    }
    _norm = norm('C', Q, i);
    printf("ITERATION END %d\n", i);
    draw2DMatrix(Q);
    if (_norm != 0) /* if norm is zero, zero vector? */
      scaleColumn(Q, i, 1/_norm);
  }

  Matrix Q_t = transposeMatrix(Q);
  Matrix R = multiplyMatrices(Q_t, A);
  
  freeMatrix(Q_t);
  QR[0] = Q;
  QR[1] = R;
  return QR;
}

int main()
{
  Matrix matrix1, matrix2;

  matrix1 = makeMatrix(4, 4, 'V', 1);

  draw2DMatrix(matrix1);

  Matrix *QR = malloc(sizeof(Matrix) * 2);

  gramSchmidt(matrix1, QR);
  
  printf("Q:\n");
  draw2DMatrix(QR[0]);
  printf("R:\n");
  draw2DMatrix(QR[1]);

  matrix2 = multiplyMatrices(QR[0], QR[1]);
  printf("QR:\n");
  draw2DMatrix(matrix2);

  freeMatrix(matrix1);
  freeMatrix(QR[0]);
  freeMatrix(QR[1]);
  freeMatrix(matrix2);
  free(QR);

  return 0;
}

