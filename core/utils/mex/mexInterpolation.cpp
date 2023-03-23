// Base script for various interpolation routines

#include <cmath>
#include <mex.h>
#include <array>

/* Binary search follows */
size_t binary_search(const double * X, const double x, size_t n){
    size_t upper = n-1;
    size_t lower = 0;
    size_t mid;
    if(x >= X[upper]){
        // Start of last interval
        return n - 2;
    }else if(x <= X[lower]){
        return 0;
    }
    while(true){
        // We are inside the interval, just return lower bound
        if(upper - lower == 1){
            return lower;
        }
        mid = (upper + lower)/2;
        if(X[mid] > x){
            // Midpoint is larger, set as upper bound
            upper = mid;
        }else{
            // Upper bound is larger, set as lower bound
            lower = mid;
        }
    }
}

void interp1_binary_search(const double * X, const double * F, const double * x, double * Fx, double * dFdx, size_t n, size_t m){
    /* Interpolate table of length n with m samples, using binary search */
    #pragma omp parallel for schedule(static, 512)
    for(int ix = 0; ix < m; ix++){
        size_t pos = binary_search(X, x[ix], n);
        double der = (F[pos+1] - F[pos])/(X[pos+1] - X[pos]);
        Fx[ix] = F[pos] + der*(x[ix] - X[pos]);
        dFdx[ix] = der;
    }
}

/* Constant size search */
void interp1_equal_width_search(const double X0, const double dX, const double * F, const double * x, double * Fx, double * dFdx, size_t n, size_t m){
    /* Interpolate table of length n with m samples, using binary search */
    #pragma omp parallel for schedule(static, 512)
    for(int ix = 0; ix < m; ix++){
        double xi = x[ix] - X0;
        int pos;
        if(xi >= 0 && xi <= (n-1)*dX ){
            pos = xi/dX;
        }else if(xi < 0){
            // Out of bounds, lower
            pos = 0;
        }else{
            // Out of bounds, upper
            pos = n - 2;
        }
        double F_i = F[pos];
        double der = (F[pos+1] - F_i)/dX;
        // mexPrintf("%d: Value %f, position %d (%f -> %f)\n", ix, x[ix], pos, X0 + pos*dX, X0 + pos*dX + dX);
        Fx[ix] = F_i + der*(xi - pos*dX);
        dFdx[ix] = der;
    }
}

/* Binned search follows */
void interp1_binning(const double * X, const double * x, size_t * bins, size_t n, size_t m){
    /* Assign each interpolation value to a bin */
    #pragma omp parallel for
    for(int i = 0; i < m; i++){
        if(x[i] <= X[0]){
            bins[i] = 0;
        }else if(x[i] >= X[n-1]){
            bins[i] = n-2;
        }
    }
    #ifdef _WIN32
        #pragma omp parallel for
    #else
        #pragma omp parallel for collapse(2)
    #endif
    for(int bin = 1; bin < n; bin++){
        for(int i = 0; i < m; i++){
            if(x[i] > X[bin-1] & x[i] <= X[bin]){
                bins[i] = bin-1;
            }
        }
    }
}

void interp1_binned_search(const double * X, const double * F, const double * x, double * Fx, double * dFdx, size_t n, size_t m){
    /* Interpolate table of length n with m samples, using binned search */
    size_t * bins = new size_t[m];
    interp1_binning(X, x, bins, n, m);
    #pragma omp parallel for
    for(int ix = 0; ix < m; ix++){
        // mexPrintf("%d: Value %f, bin %d\n", ix, x[ix], bins[ix]);
        size_t pos = bins[ix];
        double der = (F[pos+1] - F[pos])/(X[pos+1] - X[pos]);
        Fx[ix] = F[pos] + der*(x[ix] - X[pos]);
        dFdx[ix] = der;
    }
    delete[] bins;
}


