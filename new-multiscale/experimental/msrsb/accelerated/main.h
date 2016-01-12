#pragma once
#ifdef USEMEX
#include <mex.h>
#endif

int main(int argc, char* argv[]);
int flatFileBasis(std::string fn);
void printMatrix(ConnMatrix * mat);
void printGrid(Grid * grid);

#ifdef USEMEX
void mexFunction(int nlhs, mxArray *[], int nrhs, const mxArray *prhs[]);
void buildMatrixFromMxSparse(ConnMatrix * mat, const mxArray * A);
#endif
void printBasis(Grid * grid, double * basis);
