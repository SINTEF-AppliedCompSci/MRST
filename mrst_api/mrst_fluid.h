/*===========================================================================
//
// File: mrst_fluid.h
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011-2021 SINTEF Digital, Mathematics & Cybernetics.
*/

#ifndef MRST_FLUID_H_HEADER
#define MRST_FLUID_H_HEADER

#include <stddef.h>
#include <mex.h>

struct mrst_fluid {
    int            (*nphase)  (struct mrst_fluid *fluid);
    const double * (*visc)    (struct mrst_fluid *fluid);
    const double * (*dens)    (struct mrst_fluid *fluid);
    void           (*relperm) (struct mrst_fluid *fluid,
                               size_t nsat, const double *sat,
                               double *kr, double *dkr);

    void           (*destroy) (struct mrst_fluid *fluid);

    void            *fluid_data;
};

struct mrst_fluid *
mrst_fluid_create(const mxArray *FLUID);

#endif  /* MRST_FLUID_H_HEADER */
