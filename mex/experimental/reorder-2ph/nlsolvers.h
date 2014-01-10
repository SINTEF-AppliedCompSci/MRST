/*------------------------------------------------------------
 File       : nlsolvers.h
 Created    : Thu Oct 16 09:38:47 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id$


 Description


 Changes

------------------------------------------------------------*/

#ifndef NLSOLVERS_H
#define NLSOLVERS_H

double bisection   (double (*)(double, void*), void*, double, int, int*);
double ridder      (double (*)(double, void*), void*, double, int, int*);
double regulafalsi (double (*)(double, void*), void*, double, int, int*);

#endif /* NLSOLVERS_H */
