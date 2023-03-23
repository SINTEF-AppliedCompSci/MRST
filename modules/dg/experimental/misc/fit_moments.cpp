#include <iostream>
#include <mex.h>
#include <matrix.h>

void fit_moments(double *moments, double *points){
    /* Fit moments                                                       */
    int m; // Number of cubature points
    int n; // Number of moments 
    
    while (m > 0 && reduced){
        /* Buil matrix                                                   */
        P       = build_basis_matix(points, m, n);
        weights = least_squares_nonneg(P, moments);
    }
}

double build_basis_matrix(int points, int m, int n, const mxArray *psi){
    /* Build basis matrix                                                */
    double P[m][n];
    for (m_no = 0; m_no < m; m_no++){
        for n_no = 0; n_no < n; n_no++){
            P[m_no][n_no] = psi[m_no](points[n_no]);
        }
    }
    return P;
}

double significance(double P, double weights){
    /* Get significance of cubature points                               */    
}

void least_squares_nonneg(){
    /* Solve least squares problem with nonnegative constraints          */
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    const mxArray *moment_ptr;
    double *moments;
    const mxArray *points;

    int n_cells;
    n_cells = mxGetNumberOfElements(prhs[0]);
    for (int c_no = 0; c_no < n_cells; c_no++){
        
        moment_ptr = mxGetCell(prhs[0],c_no);
        
        int n_moments;
        n_moments = mxGetM(moment_ptr);
        
        moments = mxGetDoubles(moment_ptr);
        const mwSize *dims;
        
        for (int m_no = 0; m_no < n_moments; m_no++){
        std::cout << moments[m_no] << " ";
        }
        std::cout << "\n";
    }
    
}