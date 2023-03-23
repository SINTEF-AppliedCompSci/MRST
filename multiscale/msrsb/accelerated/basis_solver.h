/*
 * Copyright 2016 (c) SINTEF ICT, Applied Mathematics.
 * Maintained by Olav Moyner <olav.moyner@sintef.no>
 */
#pragma once
#include "data_defs.h"
#include <string>
void getBasis(Grid * grid, ConnMatrix * mat, double * basis, double tol, int N, double relax);
int computeBasis(Grid * grid, ConnMatrix * mat, double * basis, double tol, int maxiter, double omega);
int readInfo(Grid * grid, std::string file_path);
int readConnMatrix(Grid * grid, ConnMatrix * mat, std::string file_path);
int readBasisOperator(Grid * grid, double * basis, std::string file_path);
int openFlatFile(std::string file_path, std::ifstream * file);
int writeBasisOperator(Grid * grid, std::string file_path, double * basis);
void setupCoarseMapping(Grid * grid, ConnMatrix * mat);
void renormalize(Grid * grid, ConnMatrix * mat, double * basis);


int testfn(Grid * grid);
