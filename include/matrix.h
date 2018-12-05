/*
  @file matrix.h
  @author Gerardo Veltri
  Base matrix struct and functions
*/
#ifndef HEADER_FILE
#define HEADER_FILE

typedef struct _Matrix_ {
  int n, m;
  float **values;
} *Matrix;

Matrix makeMatrix(int n, int m, char type, float scalar);
Matrix copyMatrix(Matrix matrix);
Matrix transposeMatrix(Matrix source);
void freeMatrix(Matrix matrix);
float matrixMax(Matrix matrix);
void draw2DMatrix(Matrix matrix);

void scaleColumn(Matrix matrix, int idx, float scalar);
void subtractColumn(Matrix matrix1, int idx1, Matrix Matrix2, int idx2);

float dotProduct(char orient, Matrix matrix1, int idx1, Matrix matrix2, int idx2);
float dotProductV(Matrix matrix1, Matrix matrix2);
float norm(char orient, Matrix matrix, int idx);
float normV(Matrix matrix);
Matrix project(Matrix source, int idx1, Matrix target, int idx2);

Matrix multiplyMatrices(Matrix matrix1, Matrix matrix2);

#endif
