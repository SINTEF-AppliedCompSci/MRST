#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <mex.h>

#include "coarse_conn.h"
#include "coarse_sys.h"
#include "partition.h"
#include "sparse_sys.h"

#include "mrst_objects.h"
#include "mrst_api.h"
#include "mrst_objects.h"
#include "call_umfpack.h"


/* ------------------------------------------------------------------ */
static int
verify_faces_structure(mxArray *faces)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(faces, "neighbors") >= 0);
    ok = ok && (mxGetFieldNumber(faces, "areas"    ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "normals"  ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "centroids") >= 0);

    return ok;
}


/* ------------------------------------------------------------------ */
static int
verify_cells_structure(mxArray *cells)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(cells, "facePos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "faces"    ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "volumes"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "centroids") >= 0);

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_grid(const mxArray *G)
/* ---------------------------------------------------------------------- */
{
    int nodes_ok = 0, faces_ok = 0, cells_ok = 0, field_no;

    mxArray *pm;

    if (mxIsStruct(G)) {
        nodes_ok = mxGetFieldNumber(G, "nodes") >= 0;

        field_no = mxGetFieldNumber(G, "faces");
        faces_ok = field_no >= 0;
        if (faces_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            faces_ok = verify_faces_structure(pm);
        }

        field_no = mxGetFieldNumber(G, "cells");
        cells_ok = field_no >= 0;
        if (cells_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            cells_ok = verify_cells_structure(pm);
        }
    }

    return nodes_ok && faces_ok && cells_ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_rock(const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return mxIsStruct(rock) &&
        (mxGetFieldNumber(rock, "perm") >= 0);
}


/* ---------------------------------------------------------------------- */
static int
verify_cg(const mxArray *cg)
/* ---------------------------------------------------------------------- */
{
    int      ok;
    mxArray *sub, *fld;

    ok = !mxIsEmpty(cg) && mxIsStruct(cg);

    if (ok) {
        ok = (sub = mxGetField(cg, 0, "cells")) != NULL;

        ok = ok && !mxIsEmpty(sub) && mxIsStruct(sub);

        ok = ok &&
            ((fld = mxGetField(sub, 0, "facePos")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));

        ok = ok &&
            ((fld = mxGetField(sub, 0, "faces")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));
    }

    if (ok) {
        ok = (sub = mxGetField(cg, 0, "faces")) != NULL;

        ok = ok && !mxIsEmpty(sub) && mxIsStruct(sub);
        ok = ok                                               &&
            ((fld = mxGetField(sub, 0, "neighbors")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));

        ok = ok                                                &&
            ((fld = mxGetField(sub, 0, "subfacePos")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));

        ok = ok                                              &&
            ((fld = mxGetField(sub, 0, "subfaces")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));
    }

    if (ok) {
        ok = ok                                              &&
            ((fld = mxGetField(cg, 0, "partition")) != NULL) &&
            (mxIsDouble(fld) || mxIsInt32(fld));
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_src(const mxArray *src)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = mxIsEmpty(src);        /* Support empty src */

    if (!ok) {
        ok =        mxIsStruct(src);
        ok = ok && (mxGetFieldNumber(src, "cell") >= 0);
        ok = ok && (mxGetFieldNumber(src, "rate") >= 0) &&
            mxIsDouble(mxGetField(src, 0, "rate"));
    }

    return ok;
}


/* cs = mex_generate_coarsesystem(G, rock, CG, Lt, src) */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 5);

    ok = ok && verify_grid(prhs[0]);                            /* G */
    ok = ok && verify_rock(prhs[1]);                            /* rock */
    ok = ok && verify_cg  (prhs[2]);                            /* CG */
    ok = ok && !mxIsEmpty (prhs[3]) && mxIsDouble(prhs[3]);     /* Lt */
    ok = ok && verify_src (prhs[4]);                            /* src */

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
assign_int_vector(const mxArray *M_v, int *v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (i = 0; i < n; i++) { v[i] = pd[i]; }
    } else {
        pi = mxGetData(M_v);

        memcpy(v, pi, n * sizeof *v);
    }

    for (i = 0; i < n; i++) { v[i] -= 1; }
}


/* ---------------------------------------------------------------------- */
static void
assign_int_vector_c2m(const int *v, mxArray *M_v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (i = 0; i < n; i++) { pd[i] = v[i] + 0; }
    } else {
        pi = mxGetData(M_v);

        for (i = 0; i < n; i++) { pi[i] = v[i] + 0; };
    }
}


/* ---------------------------------------------------------------------- */
static void
assign_double_vector_c2m(const double *v, mxArray *M_v)
/* ---------------------------------------------------------------------- */
{
    size_t n;
    double *pd;

    mxAssert (mxIsDouble(M_v), "Attempt to assign DOUBLE to non-DOUBLE");

    n  = mxGetNumberOfElements(M_v);
    pd = mxGetPr(M_v);

    memcpy(pd, v, n * sizeof *pd);
}


/* ---------------------------------------------------------------------- */
static mwSize
count_dofs(size_t nblocks, struct coarse_sys *cs)
/* ---------------------------------------------------------------------- */
{
    size_t p, e;
    mwSize d, ndof = 0;

    for (p = 0, e = cs->blkdof_pos[ nblocks ]; p < e; p++) {
        d    = cs->blkdof[ p ];
        ndof = (d > ndof) ? d : ndof;
    }

    return ndof + 1;
}


/* ---------------------------------------------------------------------- */
static mxArray *
create_return_value(grid_t                 *g,
                    int                    *p,
                    struct coarse_topology *ct,
                    struct coarse_sys      *cs)
/* ---------------------------------------------------------------------- */
{
    int b, ndof, sum_ndof2, max_npair, npair, max_bcells, bcells;

    const char *fields[] = {"Dof2Conn", "blkDofPos", "basisPos",
                            "cellIPPos", "blkDof", "basis",
                            "cellIP", "BI", "maxBlkCells", "workSize",
                            "sumNdof2", "b2cPos", "b2c"};
    int nfields = sizeof(fields) / sizeof(fields[0]);

    mxArray *ret, *fld, *fld2;

    mwSize dims[] = { 1, 1 };
    mwSize ndims = sizeof(dims) / sizeof(dims[0]);

    if ((g == NULL) || (ct == NULL) || (cs == NULL)) {
        ret = mxCreateDoubleMatrix(0, 0, mxREAL); /* Empty */
    } else {
        ret = mxCreateStructArray(ndims, dims, nfields, fields);

        if (ret != NULL) {
            fld = mxCreateNumericMatrix(count_dofs(ct->nblocks, cs), 1,
                                        mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                assign_int_vector_c2m(cs->dof2conn, fld);
                mxSetField(ret, 0, "Dof2Conn", fld);
            }

            fld = mxCreateNumericMatrix(ct->nblocks + 1, 1,
                                        mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                assign_int_vector_c2m(cs->blkdof_pos, fld);
                mxSetField(ret, 0, "blkDofPos", fld);
            }

            fld = mxCreateNumericMatrix(ct->nblocks + 1, 1,
                                        mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                assign_int_vector_c2m(cs->basis_pos, fld);
                mxSetField(ret, 0, "basisPos", fld);
            }

            fld = mxCreateNumericMatrix(ct->nblocks + 1, 1,
                                        mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                assign_int_vector_c2m(cs->cell_ip_pos, fld);
                mxSetField(ret, 0, "cellIPPos", fld);
            }

            fld = mxCreateNumericMatrix(cs->blkdof_pos[ct->nblocks], 1,
                                        mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                assign_int_vector_c2m(cs->blkdof, fld);
                mxSetField(ret, 0, "blkDof", fld);
            }

            fld = mxCreateDoubleMatrix(cs->basis_pos[ct->nblocks], 1,
                                       mxREAL);
            if (fld != NULL) {
                assign_double_vector_c2m(cs->basis, fld);
                mxSetField(ret, 0, "basis", fld);
            }

            fld = mxCreateDoubleMatrix(cs->cell_ip_pos[ct->nblocks], 1,
                                       mxREAL);
            if (fld != NULL) {
                assign_double_vector_c2m(cs->cell_ip, fld);
                mxSetField(ret, 0, "cellIP", fld);
            }

            sum_ndof2 = max_npair = 0;
            for (b = 0; b < ct->nblocks; b++) {
                ndof = cs->blkdof_pos[b + 1] - cs->blkdof_pos[b];
                npair = ndof * (ndof + 1) / 2;

                if (npair > max_npair) max_npair = npair;

                sum_ndof2 += ndof * ndof;
            }

            fld = mxCreateDoubleMatrix(sum_ndof2, 1, mxREAL);
            if (fld != NULL) {
                assign_double_vector_c2m(cs->Binv, fld);
                mxSetField(ret, 0, "BI", fld);
            }

            fld = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                *(int *)mxGetData(fld) = sum_ndof2;
                mxSetField(ret, 0, "sumNdof2", fld);
            }

            max_bcells = 0;

            fld  = mxCreateNumericMatrix(ct->nblocks + 1, 1,
                                         mxINT32_CLASS, mxREAL);
            fld2 = mxCreateNumericMatrix(g->number_of_cells, 1,
                                         mxINT32_CLASS, mxREAL);
            if ((fld != NULL) && (fld2 != NULL)) {
                int *b2c_pos, *b2c;

                b2c_pos = mxGetData(fld );
                b2c     = mxGetData(fld2);

                partition_invert(g->number_of_cells, p, b2c_pos, b2c);

                for (b = 0; b < ct->nblocks; b++) {
                    bcells = b2c_pos[b + 1] - b2c_pos[b];

                    if (bcells > max_bcells) max_bcells = bcells;
                }

                mxSetField(ret, 0, "b2cPos", fld );
                mxSetField(ret, 0, "b2c"   , fld2);
            }

            fld = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                *(int *)mxGetData(fld) = max_bcells;
                mxSetField(ret, 0, "maxBlkCells", fld);
            }

            fld = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            if (fld != NULL) {
                *(int *)mxGetData(fld) = max_bcells + max_npair;
                mxSetField(ret, 0, "workSize", fld);
            }
        }
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int *
mrst_partition(const mxArray *cg)
/* ---------------------------------------------------------------------- */
{
    size_t   nc;
    int     *p;
    mxArray *M_p;

    M_p = mxGetField(cg, 0, "partition");

    nc  = mxGetNumberOfElements(M_p);
    p   = mxMalloc(nc * sizeof *p);

    if (p != NULL) {
        assign_int_vector(M_p, p);
    }

    return p;
}


/* ---------------------------------------------------------------------- */
static void
destroy_coarse_topology(struct coarse_topology *ct)
/* ---------------------------------------------------------------------- */
{
    if (ct != NULL) {
        if (ct->subfaces   != NULL) { mxFree(ct->subfaces  ); }
        if (ct->subfacepos != NULL) { mxFree(ct->subfacepos); }
        if (ct->neighbours != NULL) { mxFree(ct->neighbours); }
        if (ct->blkfaces   != NULL) { mxFree(ct->blkfaces  ); }
        if (ct->blkfacepos != NULL) { mxFree(ct->blkfacepos); }

        mxFree(ct);
    }
}


/* ---------------------------------------------------------------------- */
static struct coarse_topology *
build_coarse_topology(const mxArray *cg)
/* ---------------------------------------------------------------------- */
{
    int    n_alloc;
    size_t nel;

    struct coarse_topology *ct;

    const mxArray *sub;
    const mxArray *fld;

    ct = mxMalloc(1 * sizeof *ct);

    if (ct != NULL) {
        ct->nblocks    =           ct->nfaces   = -1  ;

        ct->neighbours = NULL;
        ct->blkfacepos = NULL;     ct->blkfaces = NULL;
        ct->subfacepos = NULL;     ct->subfaces = NULL;
        
        n_alloc = 0;
        ct->blkfacepos = getCellFacePos(cg);
        if (ct->blkfacepos != NULL) {
            n_alloc     += 1;
            ct->nblocks  = getNumberOfCells(cg);
        }

        ct->blkfaces = getCellFaces(cg);
        if (ct->blkfaces != NULL) {
            n_alloc += 1;
        }

        ct->neighbours = getFaceCellNeighbors(cg);
        if (ct->neighbours != NULL) {
            n_alloc += 1;
        }

        sub = mxGetField(cg , 0, "faces");
        fld = mxGetField(sub, 0, "subfacePos");
        nel = mxGetNumberOfElements(fld);
        ct->subfacepos = mxMalloc(nel * sizeof *ct->subfacepos);
        if (ct->subfacepos != NULL) {
            assign_int_vector(fld, ct->subfacepos);

            n_alloc    += 1;
            ct->nfaces  = nel - 1;
        }

        fld = mxGetField(sub, 0, "subfaces");
        nel = mxGetNumberOfElements(fld);
        ct->subfaces = mxMalloc(nel * sizeof *ct->subfaces);
        if (ct->subfaces != NULL) {
            assign_int_vector(fld, ct->subfaces);

            n_alloc += 1;
        }

        if (n_alloc != 5) {
            destroy_coarse_topology(ct);
            ct = NULL;
        }
    }

    return ct;
}


/* ---------------------------------------------------------------------- */
static void
local_solver(struct CSRMatrix *A, double *b, double *x)
/* ---------------------------------------------------------------------- */
{
    static unsigned long NSYS = 0;
    static double TOTAL_T = 0.0;

    double solve_t;
    struct rusage t0, t1;

    mexPrintf("Solving BF system %5lu (m=%5lu, nnz=%6lu) ... ",
              NSYS + 1, (unsigned long)A->m,
              (unsigned long)(A->ia[A->m]));

    getrusage(RUSAGE_SELF, &t0);

    callMWUMFPACK(A->m, A->ia, A->ja, A->sa, b, x);

    getrusage(RUSAGE_SELF, &t1);

    solve_t  =  t1.ru_utime.tv_sec  - t0.ru_utime.tv_sec;
    solve_t += (t1.ru_utime.tv_usec - t0.ru_utime.tv_usec) / 1.0e6;

    mexPrintf("done (solve=%.2e, total=%.3e)\n",
              solve_t, TOTAL_T + solve_t);

    NSYS++;
    TOTAL_T += solve_t;
}


/*
 * cs = mex_generate_coarsesystem(G, rock, CG, Lt, src)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int    *p;
    double *perm;
    grid_t *g;

    double *src, *totmob;
    char errmsg[1023 + 1];

    struct coarse_topology *ct;
    struct coarse_sys      *sys;

    if (args_ok(nlhs, nrhs, prhs)) {
        g = mrst_grid(prhs[0]);

        if (g != NULL) {
            perm = mrst_perm     (g->dimensions, prhs[1]);
            p    = mrst_partition(               prhs[2]);

            src    = mrst_src    (g->number_of_cells, prhs[4]);
            totmob = mxGetPr     (prhs[3]);

            ct  = NULL;
            sys = NULL;

            if ((perm != NULL) && (p != NULL)) {
                ct = build_coarse_topology(prhs[2]);
            }

            if (ct != NULL) {
                sys = coarse_sys_construct(g, p, ct, perm,
                                           src, totmob,
                                           local_solver);
            }

            if (sys != NULL) {
                plhs[0] = create_return_value(g, p, ct, sys);

                coarse_sys_destroy(sys);
                destroy_coarse_topology(ct);

                mrst_src_deallocate (src);    mxFree(p);
                mrst_perm_deallocate(perm);
            } else {
                plhs[0] = create_return_value(NULL, NULL, NULL, NULL);
            }
        }

        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tCS = %s(G, rock, CG, Lt, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
