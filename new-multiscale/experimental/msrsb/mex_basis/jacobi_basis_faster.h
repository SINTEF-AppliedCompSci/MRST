int compute_jacobi_basis(double* I, int nf, int nc, int nel, int ngb,
        int *cells, int *bndIndicator, int *cellPos,
        int *GB_globalCells, int *GBLocal, int *GBLocalMap, int* cellIsActive, 
        struct submatrix *subsys, double tol, int maxiter, double omega);