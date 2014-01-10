#include <mex.h>
#include <umfpack.h>
#include "call_umfpack.h"


#if 0
static void
sortRows(int n, mwSignedIndex* ia, mwSignedIndex* ja,
         double *sa)
{
    int i,j,k,row;

    for(row=0; row<n; ++row)
    {
        for(k=ia[row]; k<ia[row+1]; ++k)
        {
            for(i=k+1; i<ia[row+1]; ++i)
            {
                if(ja[k] > ja[i])
                {
                    mwSignedIndex itmp = ja[i];
                    ja[i] = ja[k];
                    ja[k] = itmp;

                    double dtmp = sa[i];
                    sa[i] = sa[k];
                    sa[k] = dtmp;
                }
            }
        }
    }
}
#endif


struct CSCMatrix {
    mwSignedIndex  n;
    mwSignedIndex  nnz;

    mwSignedIndex *p;
    mwSignedIndex *i;
    double        *x;
};


/* ---------------------------------------------------------------------- */
static void
csc_deallocate(struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    if (csc != NULL) {
        if (csc->x != NULL) { mxFree(csc->x); }
        if (csc->i != NULL) { mxFree(csc->i); }
        if (csc->p != NULL) { mxFree(csc->p); }

        mxFree(csc);
    }
}


/* ---------------------------------------------------------------------- */
static struct CSCMatrix *
csc_allocate(mwSignedIndex n, mwSignedIndex nnz)
/* ---------------------------------------------------------------------- */
{
    struct CSCMatrix *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        new->p = mxMalloc((n + 1) * sizeof *new->p);
        new->i = mxMalloc(nnz     * sizeof *new->i);
        new->x = mxMalloc(nnz     * sizeof *new->x);

        if ((new->p == NULL) || (new->i == NULL) || (new->x == NULL)) {
            csc_deallocate(new);
            new = NULL;
        } else {
            new->n   = n;
            new->nnz = nnz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
csr_to_csc(const int        *ia,
           const int        *ja,
           const double     *sa,
           struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    mwSignedIndex i, nz;

    /* Clear garbage, prepare for counting */
    for (i = 0; i <= csc->n; i++) { csc->p[i] = 0; }

    /* Count column connections */
    for (nz = 0; nz < csc->nnz; nz++) {
        csc->p[ ja[nz] + 1 ] += 1;
    }

    /* Define column start pointers */
    for (i = 1; i <= csc->n; i++) {
        csc->p[0] += csc->p[i];
        csc->p[i]  = csc->p[0] - csc->p[i];
    }

    mxAssert (csc->p[0] == csc->nnz,
              "Internal error defining CSC start pointers.");

    /* Fill matrix whilst defining column end pointers */
    for (i = nz = 0; i < csc->n; i++) {
        for (; nz < ia[i + 1]; nz++) {
            csc->i[ csc->p[ ja[nz] + 1 ] ] = i;      /* Insertion sort */
            csc->x[ csc->p[ ja[nz] + 1 ] ] = sa[nz]; /* Insert mat elem */

            csc->p        [ ja[nz] + 1 ]  += 1;      /* Advance col ptr */
        }
    }

    mxAssert (csc->p[csc->n] == csc->nnz,
              "Internal error defining CSC end pointers.");

    csc->p[0] = 0;
}


/* ---------------------------------------------------------------------- */
static void
solve_umfpack(struct CSCMatrix *csc, const double *b, double *x)
/* ---------------------------------------------------------------------- */
{
    void *Symbolic, *Numeric;
    double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

    umfpack_dl_defaults(Control);

    umfpack_dl_symbolic(csc->n, csc->n, csc->p, csc->i, csc->x,
                        &Symbolic, Control, Info);
    umfpack_dl_numeric (csc->p, csc->i, csc->x,
                        Symbolic, &Numeric, Control, Info);

    umfpack_dl_free_symbolic(&Symbolic);

    umfpack_dl_solve(UMFPACK_A, csc->p, csc->i, csc->x, x, b,
                     Numeric, Control, Info);

    umfpack_dl_free_numeric(&Numeric);
}


/*---------------------------------------------------------------------------*/
void
callMWUMFPACK(int n, int* ia, int* ja, double *sa, double *b, double *x)
/*---------------------------------------------------------------------------*/
{
    struct CSCMatrix *csc;

    csc = csc_allocate(n, ia[n]);

    if (csc != NULL) {
        csr_to_csc(ia, ja, sa, csc);

        solve_umfpack(csc, b, x);
    }

    csc_deallocate(csc);
}
