/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TWOPHASE_H_INCLUDED
#define TWOPHASE_H_INCLUDED



struct UnstructuredGrid;
struct SolverData  {
    struct UnstructuredGrid *grid;
    const double            *darcyflux;   /* one flux per face  in cdata::grid*/
    const double            *porevolume;  /* one volume per cell */
    const double            *source;      /* one source per cell */
    const int               *satnum;      /* saturation region of each cell */
    double                   dt;
    double                  *saturation;      /* one per cell */
    double                  *fractionalflow;  /* one per cell */
};

struct NonlinearSolverCtrl;


void
solvecell (void *data, struct NonlinearSolverCtrl *ctrl, int cell);

void
destroy_solverdata(struct SolverData *d);

struct SolverData *
init_solverdata(struct UnstructuredGrid *grid, const double *darcyflux,
                const double *porevolume, const double *source,
                const int *satnum, double dt, double *saturation);

#endif /* TWOPHASE_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
