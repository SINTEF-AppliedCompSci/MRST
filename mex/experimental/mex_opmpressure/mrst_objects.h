#ifndef MRST_OBJECTS_H_INCLUDED
#define MRST_OBJECTS_H_INCLUDED

#include <stddef.h>

#include <mex.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "flow_bc.h"
#include "well.h"


struct well_data {
    double *WI, *wdp;
};


struct mrst_well {
    well_t           *wdesc;
    well_control_t   *wctrl;
    struct well_data *wdata;
};


flowbc_t *
mrst_flowbc(size_t nf, const mxArray *BC);

void
mrst_flowbc_deallocate(flowbc_t *bc);

double *
mrst_src(size_t nc, const mxArray *SRC);

void
mrst_src_deallocate(double *src);

double *
mrst_perm(int d, const mxArray *rock);

void
mrst_perm_deallocate(double *perm);

size_t
mrst_well_count_totperf(const mxArray *W);

struct mrst_well *
mrst_well(const mxArray *W);

void
mrst_well_deallocate(struct mrst_well *W);
#ifdef __cplusplus
}
#endif

#endif  /* MRST_OBJECTS_H_INCLUDED */
