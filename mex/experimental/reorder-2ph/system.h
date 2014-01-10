/*===========================================================================

 File: system.h

 Created: Mon Jun  2 17:23:16 2008

 Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>

 Revision: $Id$

 Description:

===========================================================================*/

#ifndef SYSTEM_H
#define SYSTEM_H


/*                                                      *
 *   System of type                                     *
 *                                                      *
 *     S - So + dtBS + dtAf(S) + dtCf(S) =  dt*Q        *
 *   or                                                 *
 *     (1 + dtB)S + dt(A + C)f(S) =  So + dt*Q          *
 *                                                      *
 *                                                      *
 * Sources are actually positive sources.  Negative     *
 * ones are added to the diagonal of A.  Need not       *
 * distinguich between A and C either. We would         *
 * then get                                             *
 *                                                      *
 *   (1+dtB)*S + dt*AC*f(S) = So dt*Q                   *
 *                                                      *
 * in the paper we use                                  *
 *                                                      *
 *   (1+dtA)u^n + dt/phi(B+V)f(u^n) = dt/phiQ +u^{n-1}  *
 *                                                      *
 *                                                      *
 *                                                      */

typedef struct {
  int       n;         /* dimension                     */
  double   *u;         /* new solution                  */
  double   *u_old;     /* old solution                  */
  int      *satnum;    /* saturation regions            */

  /*
   *  Make changes here to avoid copy of flux matrix.
   */
  sparse_t *V;         /* flux matrix(change to straigh input flux matrix. */
  double   *diagonal;  /* sum flux out + negative source. */

  double   *Qp;        /* Source terms.                 */
  double    dt;        /* time step                     */
  double   *f;         /* flux function                 */

  /* Settings */
  int       maxit;
  double    tolerance;
} system_t;

typedef struct
{
  /* Logging */
  int      *iterations; /* Num. iterations in each cell. */
  double   *res;        /* Final residual in each cell.  */
} report_t;


void init_system(system_t *sys, double dt, sparse_t *flux, double *sources,
		 double *volumes, double *u_old, int *satnum,
		 char *solver, int maxit, double tolerance,double *work);

int solveScalar(int cell, void *data, void* report);

#endif /* SYSTEM_H */
