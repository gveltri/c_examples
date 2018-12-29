/*
  @file precision.h
  @author Gerardo Veltri
  Summarize differences in outcome and expected values of matrices
  introduced by precision errors  
*/
#ifndef PRECISION_HEADER
#define PRECISION_HEADER

void identityPrecision(Matrix matrix, double *stats);
void matrixComparison(Matrix matrix1, Matrix matrix2, double *stats);

#endif
