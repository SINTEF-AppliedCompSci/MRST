/*===========================================================================

 File: system.h

 Created: Tue Nov 15 12:59:17 CET 2011

 Author: Knut-Andreas Lie <Knut-Andreas.Lie@sintef.no>

 Revision: $Id$

 Description:
   Implementation of the Corey fluid object for a single stone type

===========================================================================*/

#ifndef FLUID_H
#define FLUID_H

void   init_fluid(const mxArray *);
void   end_fluid();
double fluxfun (double sw, const int);
double dfluxfun(double sw, const int);
int    getNoSatRegions(void);

#endif /* FLUID_H */
