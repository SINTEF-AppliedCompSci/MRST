/*
  mex version of inhull for 3D
*/

#include "mex.h"
#include <math.h>
#include <stdio.h>

void printvec(const double a[3])
{
    printf("%f %f %f\n", a[0], a[1], a[2]);
}

void printmtx(const double *a, size_t n)
{
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < 3; ++j)
      printf("%f ", a[i+j*n]);
    printf("\n");
  }
}
void printmtxi(const double *a, size_t n)
{
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < 3; ++j)
      printf("%i ", (int)a[i+j*n]);
    printf("\n");
  }
}

double** dalloc(int rows, int cols, double init)
{
    int i, j;
    double** mat = (double **)malloc(rows * sizeof(double*));
    for (i = 0; i < rows; ++i) {
        mat[i] = (double *)malloc(cols * sizeof(double));
        for (j = 0; j < cols; ++j)
            mat[i][j] = init;
    }
    return mat;
}

void dfree(int rows, double **mat)
{
    int i;
    for (i = 0; i < rows; ++i)
        free((void *)mat[i]);
    free((void *)mat);
}

void dreshape(const double* a, int r, int c, double** b)
{
    int i, j;
    for (i = 0; i < r; ++i) {
        for (j = 0; j < c; ++j) {
            b[i][j] = a[j*r+i];
        }
    }
}

void cross(const double u[3], const double v[3], double w[3])
{
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}

double norm(const double w[3])
{
    return sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
}

double dot(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void mean(const double* const* a, int r, double* m)
{
    int i, d;
    m[0] = m[1] = m[2] = 0.0;

    for (i = 0; i < r; ++i)
        for (d = 0; d < 3; ++d)
            m[d] += a[i][d];

    for (d = 0; d < 3; ++d)
        m[d] /= r;
}

void mean3(const double a[3], const double b[3], const double c[3], double* m)
{
    int d;
    m[0] = m[1] = m[2] = 0.0;
    for (d = 0; d < 3; ++d) {
        m[d] += a[d];
        m[d] += b[d];
        m[d] += c[d];
        m[d] /= 3.0;
    }
}

void inhull(const double* const* testpts,
            mwSize ntestpts,
            const double* const* xyz,
            mwSize nxyz,
            const double* const* tess,
            mwSize ntess,
            double tol,
            int *is_inside,
            int verbose)
{
    /* This function can be optimized */
    const int dim = 3;
    double** nrmls = dalloc(ntess, dim, 0.0);
    double** xyzmeans = dalloc(ntess, dim, 0.0);
    double length;
    double center[3];
    int d;
    mwSize i, j;

    /* Compute center of tesselation */
    mean(xyz, nxyz, center);

    /* Compute normals for each face */
    for (i = 0; i < ntess; ++i)
    {
        int a, b, c;
        double ab[3], ac[3];
        a = (int)tess[i][0];
        b = (int)tess[i][1];
        c = (int)tess[i][2];
        for (d = 0; d < 3; ++d) {
            ab[d] = xyz[a][d] - xyz[b][d];
            ac[d] = xyz[a][d] - xyz[c][d];
        }
        cross(ab, ac, nrmls[i]);
        length = norm(nrmls[i]);
        mean3(xyz[a], xyz[b], xyz[c], xyzmeans[i]);

        if (length > tol)
        {
            /* Normalize and compute vector used for normal orientation */
            double s[3];
            for (d = 0; d < 3; ++d) {
                nrmls[i][d] /= length;
                s[d] = center[d] - xyzmeans[i][d];
            }
            /* Flip normals */
            if (dot(s, nrmls[i]) < 0.0)
                for (d = 0; d < 3; ++d)
                    nrmls[i][d] *= -1;
        }
        else {
            if (verbose)
            {
                printf("mex_inhull found degenerate triangle %d:\n", i);
                printvec(ab);
                printvec(ac);
                printvec(nrmls[i]);
            }
        }
    }

    /* Check if the test points j are inside */
    for (j = 0; j < ntestpts; ++j) {
        /* Assume point inside */
        is_inside[j] = 1;
        for (i = 0; i < ntess; ++i)
        {
            double s[3];
            double dotp;
            for (d = 0; d < 3; ++d)
                s[d] = xyzmeans[i][d] - testpts[j][d];
            dotp = dot(s, nrmls[i]);
            if (verbose)
                printf("%d dotp=%f\n", i, dotp);
            if (dotp > tol) {
                is_inside[j] = 0;
                break;
            }
        }
    }

    dfree(ntess, nrmls);
    dfree(ntess, xyzmeans);
}

/* Gateway fcn */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *testpts;
    size_t ntestpts;
    size_t dim;
    double *xyz;
    size_t nxyz;
    double *tess;
    size_t ntess;
    double tol;
    int *is_inside;
    size_t i,j;
    double **testpts_;
    double **xyz_;
    double ** tess_;
    int verbose;

    if (nrhs < 4) {
        mexErrMsgTxt("Four or five inputs are required.");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("One output is required.");
    }

    testpts = mxGetPr(prhs[0]);
    ntestpts = mxGetM(prhs[0]); /* Rows */
    dim = mxGetN(prhs[0]); /* Columns */
    if (dim != 3) {
        mexErrMsgTxt("dim (number of columns) must be 3");
    }

    xyz = mxGetPr(prhs[1]);
    if (mxGetN(prhs[1]) != dim) {
        mexErrMsgTxt("Test points dimension and points defining the hull are not the same");
    }
    nxyz = mxGetM(prhs[1]);

    tess = mxGetPr(prhs[2]);
    if (mxGetN(prhs[2]) != dim) {
        mexErrMsgTxt("tesselation must be simplices in dimension d-1");
    }
    ntess = mxGetM(prhs[2]);

    tol = mxGetScalar(prhs[3]);
    if (nrhs == 5)
        verbose = mxGetScalar(prhs[4]);
    else
        verbose = 0;

    if (verbose) {
        printf("ntestpts %i\n", ntestpts);
        printmtx(testpts, ntestpts);
        printf("nxyz %i\n", nxyz);
        printmtx(xyz, nxyz);
        printf("ntess %i\n", ntess);
        printmtxi(tess, ntess);
        printf("tol=%f\n", tol);
    }


    plhs[0] = mxCreateNumericMatrix((mwSize)ntestpts, 1, mxINT32_CLASS, mxREAL);
    is_inside = (int*)mxGetPr(plhs[0]);

    /* Reshape from array to matrix */
    testpts_ = dalloc(ntestpts, dim, 0.0);
    dreshape(testpts, ntestpts, dim, testpts_);
    xyz_ = dalloc(nxyz, dim, 0.0);
    dreshape(xyz, nxyz, dim, xyz_);
    tess_ = dalloc(ntess, dim, 0.0);
    dreshape(tess, ntess, dim, tess_);

    /* Adjust for zero indexing */
    for (i = 0; i < ntess; ++i)
        for (j = 0; j < dim; ++j)
            tess_[i][j]--;

    /* Call main fcn */
    inhull((const double* const*)testpts_, (mwSize)ntestpts,
           (const double* const*)xyz_, (mwSize)nxyz,
           (const double* const*)tess_, (mwSize)ntess,
           tol,
           is_inside,
           verbose);

    dfree(ntestpts, testpts_);
    dfree(nxyz, xyz_);
    dfree(ntess, tess_);
}
