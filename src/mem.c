/*
  @file mem.c
  @author Gerardo Veltri
  Memory functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mem.h>

void allocValues(Matrix matrix)
{
  for (int i=0;i<matrix->n;i++)
  {
    matrix->values[i] = malloc(matrix->m*sizeof(double));
  }
}

Matrix allocMatrix(int n, int m)
{
  Matrix matrix = malloc(sizeof(struct _Matrix_));

  matrix->n = n;
  matrix->m = m;
  matrix->values = malloc(n*sizeof(double*));
  allocValues(matrix);

  return matrix;
}

void freeMatrix(Matrix matrix)
{
  for (int i=0; i<matrix->n; i++)
  {
    free(matrix->values[i]);
  }

  free(matrix->values);
  free(matrix);
}

MatrixStack allocMatrixStack(int n, int m, int depth)
{
  MatrixStack stack = malloc(sizeof(struct _MatrixStack_));

  stack->n = n;
  stack->m = m;
  stack->depth = depth;
  stack->cur_depth = depth;
  stack->matrices = malloc((stack->depth)*sizeof(Matrix));
  for (int i=0; i<stack->depth; i++)
  {
    stack->matrices[i] = allocMatrix(n,m);
  }

  stack->top = (stack->matrices)+(stack->depth)-1;

  return stack;
}

Matrix popMatrixStack(MatrixStack stack)
{
  Matrix popped = *(stack->top);
  stack->top--;
  stack->cur_depth--;

  return popped;
}

void pushMatrixStack(MatrixStack stack, Matrix matrix)
{
  stack->top++;
  *(stack->top) = matrix;
  stack->cur_depth++;
}

void freeMatrixStack(MatrixStack stack)
{
  if (stack->cur_depth == 0)
  {
    free(stack);
    return;
  }

  do
  {
    freeMatrix(*(stack->top--));
    stack->cur_depth--;
  } while (stack->cur_depth > 0);

  free(stack);
}
