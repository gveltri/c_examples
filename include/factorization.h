/*
  @file factorization.h
  @author Gerardo Veltri
  Factorization and Related Algorithms
*/
#ifndef FACTORIZATION_HEADER
#define FACTORIZATION_HEADER

void _gramSchmidtQR(Matrix A, Matrix QR[2], MatrixStack mem_stacks[2], int debug);
void gramSchmidtQR(Matrix A, Matrix QR[2], int debug);

void _hhReflectionsQR(Matrix A, Matrix QR[2],
                      MatrixStack mem_stacks[2],
                      int debug);
void hhReflectionsQR(Matrix A, Matrix QR[2],
                     int debug);

void gaussianElimination(Matrix A, Matrix B, Matrix RREF[2], int debug);

#endif
