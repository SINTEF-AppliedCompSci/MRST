/*
Copyright (C) 2012 (c) Jostein R. Natvig <jostein natvig at gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef NLSOLVERS_H
#define NLSOLVERS_H

struct NonlinearSolverCtrl
{
   enum     Method {RIDDERS, REGULAFALSI, BISECTION} method;
   double  nltolerance;
   int     maxiterations;
   double  initialguess;
   int     iterations;    /* set by solver */
   double  residual;      /* set by solver */
};

double find_zero  (double (*G)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl);
double bisection   (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);
double ridders     (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);
double regulafalsi (double (*)(double, void*), void*, struct NonlinearSolverCtrl *ctrl);

#endif /* NLSOLVERS_H */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
