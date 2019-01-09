/*
  @file mem.c
  @author Gerardo Veltri
  Memory functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mem.h>

Matrix allocMatrix(int n, int m)
{
        Matrix matrix = malloc(sizeof(struct _Matrix_));

        matrix->n = n;
        matrix->m = m;
        matrix->values = malloc(n*m*sizeof(double));

        return matrix;
}

void freeMatrix(Matrix matrix)
{
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
    assert(stack->cur_depth > 0);
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

/*
  freeMatrixStackAll

  frees all initialized matrices from stack
  frees stack
*/
void freeMatrixStackAll(MatrixStack stack)
{
    while (stack->depth > 0)
        {
            freeMatrix(*(stack->matrices++));
            stack->depth--;
        }

    free(stack);
}

/*
  freeMatrixStack

  frees remaining matrices from stack
  frees stack
*/
void freeMatrixStack(MatrixStack stack)
{
    while (stack->cur_depth > 0)
        {
            freeMatrix(*(stack->top--));
            stack->cur_depth--;
        }

    free(stack);
}
