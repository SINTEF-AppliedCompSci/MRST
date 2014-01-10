/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TWOPHASETRANSPORT_H_INCLUDED
#define TWOPHASETRANSPORT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

struct UnstructuredGrid;
void twophasetransport(
    const double *porevolume,
    const double *source,
    double dt,
    struct UnstructuredGrid *grid,
    const double *darcyflux,
    const int *satnum,
    double *saturation);

#ifdef __cplusplus
}
#endif

#endif /* TWOPHASETRANSPORT_H_INCLUDED */
