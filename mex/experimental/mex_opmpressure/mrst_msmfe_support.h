#ifndef MRST_MSMFE_SUPPORT_H_INCLUDED
#define MRST_MSMFE_SUPPORT_H_INCLUDED

struct coarse_sys;

struct coarse_sys_meta {
    int     nblocks;
    int     max_bcells;

    size_t  work_size;
    size_t  sum_ndof2;

    int    *b2c_pos;
    int    *b2c;
};


int verify_mex_cs(const mxArray *cs);

void set_coarse_sys_meta(const mxArray *cs,
                         struct coarse_sys_meta *csys_meta);

void set_coarse_sys(const mxArray *cs, struct coarse_sys *csys);

#endif  /* MRST_MSMFE_SUPPORT_H_INCLUDED */
