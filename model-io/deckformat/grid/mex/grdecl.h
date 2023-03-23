#ifndef GRDECL_H_INCLUDED
#define GRDECL_H_INCLUDED

struct grdecl{
    int           dims[3];
    const double *coord;
    const double *zcorn;
    const int    *actnum;
};


#endif /* GRDECL_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
