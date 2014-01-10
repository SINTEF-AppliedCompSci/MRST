/*===========================================================================

 File: system.h

 Created: Tue Nov 15 12:59:17 CET 2011

 Author: Knut-Andreas Lie <Knut-Andreas.Lie@sintef.no>

 Revision: $Id$

 Description:
   Implementation of the Corey fluid object for a single rock type

===========================================================================*/

#ifndef FLUID_H
#define FLUID_H

#ifdef MATLAB_MEX_FILE
void   init_fluid(const mxArray *);
#else
void   init_fluid(void);
#endif

double fluxfun   (double sw, const int);
double dfluxfun  (double sw, const int);
int    getNumSatRegions(void);

#endif /* FLUID_H */
