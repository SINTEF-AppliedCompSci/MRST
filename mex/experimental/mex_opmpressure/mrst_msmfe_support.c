#include <mex.h>

#include "coarse_sys.h"

#include "mrst_msmfe_support.h"

/* ---------------------------------------------------------------------- */
int
verify_mex_cs(const mxArray *cs)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =        !mxIsEmpty(cs) && mxIsStruct(cs);

    /* set_coarse_sys_meta() */
    ok = ok && (mxGetFieldNumber(cs, "blkDofPos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "maxBlkCells") >= 0);
    ok = ok && (mxGetFieldNumber(cs, "workSize"   ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "sumNdof2"   ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "b2cPos"     ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "b2c"        ) >= 0);

    /* set_coarse_sys() */
    ok = ok && (mxGetFieldNumber(cs, "Dof2Conn"   ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "cellIPPos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "blkDof"     ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "basis"      ) >= 0);
    ok = ok && (mxGetFieldNumber(cs, "cellIP"     ) >= 0);

    return ok;
}


/* ---------------------------------------------------------------------- */
void
set_coarse_sys_meta(const mxArray *cs, struct coarse_sys_meta *csys_meta)
/* ---------------------------------------------------------------------- */
{
    mxArray *fld;

    fld = mxGetField(cs, 0, "blkDofPos");
    csys_meta->nblocks = (int)(mxGetM(fld) - 1);

    fld = mxGetField(cs, 0, "maxBlkCells");
    csys_meta->max_bcells = *(const int *)mxGetData(fld);

    fld = mxGetField(cs, 0, "workSize");
    csys_meta->work_size = *(const int *)mxGetData(fld);

    fld = mxGetField(cs, 0, "sumNdof2");
    csys_meta->sum_ndof2 = *(const int *)mxGetData(fld);

    fld = mxGetField(cs, 0, "b2cPos");
    csys_meta->b2c_pos = mxGetData(fld);

    fld = mxGetField(cs, 0, "b2c");
    csys_meta->b2c = mxGetData(fld);
}


/* ---------------------------------------------------------------------- */
void
set_coarse_sys(const mxArray *cs, struct coarse_sys *csys)
/* ---------------------------------------------------------------------- */
{
    mxArray *fld;

    fld = mxGetField(cs, 0, "Dof2Conn");
    csys->dof2conn = mxGetData(fld);

    fld = mxGetField(cs, 0, "blkDofPos");
    csys->blkdof_pos = mxGetData(fld);

    fld = mxGetField(cs, 0, "basisPos");
    csys->basis_pos = mxGetData(fld);

    fld = mxGetField(cs, 0, "cellIPPos");
    csys->cell_ip_pos = mxGetData(fld);

    fld = mxGetField(cs, 0, "blkDof");
    csys->blkdof = mxGetData(fld);

    fld = mxGetField(cs, 0, "basis");
    csys->basis = mxGetPr(fld);

    fld = mxGetField(cs, 0, "cellIP");
    csys->cell_ip = mxGetPr(fld);
}
