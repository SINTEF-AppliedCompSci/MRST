/*=========================================================================
//
// File: mxgrdecl.c
//
// Created: Fri Jun 19 08:48:21 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2012-05-15 22:24:02 +0200 (Tue, 15 May 2012) $
//
// $Revision: 967 $
//
//=======================================================================*/

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <mex.h>

#include "grdecl.h"
#include "mxgrdecl.h"

/* Get COORD, ZCORN, ACTNUM and DIMS from mxArray.       */
/*-------------------------------------------------------*/
void mx_init_grdecl(struct grdecl *g, const mxArray *s)
{
    int i,n;
    size_t numel;
    mxArray *cartdims=NULL, *actnum=NULL, *coord=NULL, *zcorn=NULL;

    if (!mxIsStruct(s)
        || !(cartdims = mxGetField(s, 0, "cartDims"))
        || !(coord    = mxGetField(s, 0, "COORD"))
        || !(zcorn    = mxGetField(s, 0, "ZCORN"))
        )
    {
        char str[]="Input must be a single MATLAB struct with fields\n"
            "'cartDims', 'COORD' and 'ZCORN'. ACTNUM may be included.\n";
        mexErrMsgTxt(str);
    }


    numel = mxGetNumberOfElements(cartdims);
    if (!mxIsNumeric(cartdims) || numel != 3){
        mexErrMsgTxt("cartDims field must be 3 numbers");
    }

    if (mxIsDouble(cartdims)) {
        double *tmp = mxGetPr(cartdims);
        for (i = 0; i < 3; ++i) {
            g->dims[i] = (int) tmp[i];
        }
    }
    else if (mxIsInt32(cartdims)) {
        int *tmp = mxGetData(cartdims);
        memcpy(g->dims, tmp, 3 * sizeof *g->dims);
    }

    n = g->dims[0];
    for (i = 1; i < 3; i++) { n *= g->dims[ i ]; }


    if ((actnum = mxGetField(s, 0, "ACTNUM")) != NULL) {
        numel = mxGetNumberOfElements(actnum);
        if ((! mxIsInt32(actnum)) || (numel != (size_t)(n))) {
            mexErrMsgTxt("ACTNUM field must be nx*ny*nz numbers int32");
        }
        g->actnum = mxGetData(actnum);
    }
    else {
        g->actnum = NULL;
    }


    numel = mxGetNumberOfElements(coord);
    if ((! mxIsDouble(coord)) ||
        numel != (size_t)(6*(g->dims[0]+1)*(g->dims[1]+1))) {
        mexErrMsgTxt("COORD field must have 6*(nx+1)*(ny+1) doubles.");
    }
    g->coord = mxGetPr(coord);


    numel = mxGetNumberOfElements(zcorn);
    if ((! mxIsDouble(zcorn)) ||
        numel != (size_t)(8*g->dims[0]*g->dims[1]*g->dims[2])) {
        mexErrMsgTxt("ZCORN field must have 8*nx*ny*nz doubles.");
    }
    g->zcorn = mxGetPr(zcorn);
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
