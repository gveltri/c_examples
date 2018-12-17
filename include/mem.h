/*
  @file mem.h
  @author Gerardo Veltri
  Memory management of matrix
*/
#ifndef MEM_HEADER
#define MEM_HEADER

typedef struct _Matrix_ {

  int n; /* columns */
  int m; /* rows */

  double **values;

} *Matrix;

typedef struct _MatrixStack_ {

  int n;
  int m;
  int depth;
  int cur_depth;

  Matrix *top;
  Matrix *matrices;

} *MatrixStack;


Matrix allocMatrix(int n, int m);
void freeMatrix(Matrix matrix);

MatrixStack allocMatrixStack(int n, int m, int depth);
Matrix popMatrixStack(MatrixStack stack);
void pushMatrixStack(MatrixStack stack, Matrix matrix);
void freeMatrixStackAll(MatrixStack stack);
void freeMatrixStack(MatrixStack stack);

#endif
