/*****************************************************************************/
/*                                                                           */
/*  tricall.c                                                                */
/*                                                                           */
/*                                                                           */
/*  Gateway routine for the Triangle mesh generator made by                  */
/*  Jonathan R. Shewchuck.  Basic building block in 2D triangular mesh       */
/*  generating tools for Matlab.                                             */
/*                                                                           */
/*                                                                           */
/*  Copyright - All rights reserved.                                         */
/*  Jostein Roald Natvig.                                                    */
/*  SINTEF ICT, Applied Mathematics.                                         */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
  #define REAL float
#else
  #define REAL double
#endif

#include <stdio.h>
#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "triangle.h"


#ifndef NULL
#define NULL 0
#endif

#define DELETED -2
#define le 0 /* Left edge */
#define re 1 /* Right edge */

#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree

void import_struct (struct triangulateio *out, const mxArray  *in );
void export_struct (struct triangulateio *in,        mxArray **out);




void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

  int     n, i;
  int     dims[2];
  int 	 *ptriangles;
  double *ppoints;
  int     status, buflen;
  char   *inbuf;

  struct triangulateio in, mid, out, vorout;

  memset(&in,     0, sizeof(struct triangulateio));
  memset(&out,    0, sizeof(struct triangulateio));
  memset(&vorout, 0, sizeof(struct triangulateio));
  import_struct(&in, prhs[0]);


  buflen   = mxGetNumberOfElements (prhs [1])  + 1;
  inbuf    = (char *) mxCalloc(buflen, sizeof(char));
  status   = mxGetString(prhs[1], inbuf, buflen);
  if (status != 0) mexWarnMsgTxt ("Not enough memory.  String is truncated.");

  triangulate (inbuf, &in, &out, &vorout);

  export_struct (&out,    &plhs[0]);
  export_struct (&vorout, &plhs[1]);


  /* free(out.pointlist); */

} /* mexFunction() */


void import_struct (struct triangulateio *out, const mxArray *in)
{
  int        nfields, i;
  const char **fnames;
  const int  *dims;
  mxArray    *field;


  if(!mxIsStruct (in)) mexErrMsgTxt ("Input to getiostruct must be a struct");

  nfields = mxGetNumberOfFields (in);
  fnames  = mxMalloc (nfields * sizeof(const char *));

  for (i=0; i<nfields; ++i)
  {
    field      = mxGetFieldByNumber     (in, 0, i);
    fnames [i] = mxGetFieldNameByNumber (in,    i);
    if (field == NULL)
    {
      mexErrMsgTxt ("Above field is empty!");
    }

    else if (!strncmp(fnames [i], "pointlist", 9))
    {
      out->pointlist      = (REAL *) mxGetPr(field);
      out->numberofpoints = mxGetN (field);
    }

    else if (!strncmp(fnames [i], "pointattributelist", 18))
    {
      out->pointattributelist  = (REAL *) mxGetPr(field);
    }

    else if (!strncmp(fnames [i], "pointmarkerlist", 15))
    {
      out->pointmarkerlist = (int *) mxGetData(field);
    }

    else if (!strncmp(fnames [i], "trianglelist", 12))
    {
      out->trianglelist      = (int *) mxGetData(field);
      out->numberoftriangles = mxGetN (field);
    }

    else if (!strncmp(fnames [i], "triangleattributelist", 21))
    {
      out->triangleattributelist  = (REAL *) mxGetPr(field);
    }

    else if (!strncmp(fnames [i], "trianglearealist", 16))
    {
      out->trianglearealist  = (REAL *) mxGetPr(field);
    }

    else if (!strncmp(fnames [i], "segmentlist", 11))
    {
      out->segmentlist      = (int *) mxGetData(field);
      out->numberofsegments = mxGetN (field);
      /*int i;
      fprintf(stderr, "number of segments = %d\n",
	      out->numberofsegments);
      for (i=0; i< out->numberofsegments; ++i)
   	fprintf(stderr, "%d %d\n", out->segmentlist[2*i],
		out->segmentlist[2*i+1]);*/
    }

    else if (!strncmp(fnames [i], "segmentmarkerlist", 17))
    {
      out->segmentmarkerlist  = (int *) mxGetData(field);
    }

    else if (!strncmp(fnames [i], "edgelist", 8))
    {
      out->edgelist      = (int *) mxGetData(field);
      out->numberofedges = mxGetN (field);
    }

    else if (!strncmp(fnames [i], "edgemarkerlist", 14))
    {
      out->edgemarkerlist  = (int *) mxGetData(field);
    }

    else if (!strncmp(fnames [i], "neighborlist", 12))
    {
      out->neighborlist  = (int *) mxGetData(field);
    }

    else if (!strncmp(fnames [i], "normlist", 8))
    {
      out->normlist  = (REAL *) mxGetPr(field);
    }
    else if (!strncmp(fnames [i], "holelist", 8))
    {
      out->holelist      = (REAL *) mxGetPr(field);
      out->numberofholes =  mxGetN (field);
    }

    else if (!strncmp(fnames [i], "regionlist", 10))
    {
      out->regionlist      = (REAL *) mxGetPr(field);
      out->numberofregions = mxGetM (field);
    }
  }

  mxFree(fnames);
}

void export_struct(struct triangulateio *in, mxArray **output)
{
  int          nfields, ifield, ndims;
  mwSize       *dims;
  int          i;
  REAL         *pr_dbl;
  int          *pr_int;
  mxClassID    id;
  mxComplexity complexflag;
  mxArray      *tmp, *field, *out;
  const char   *fnames[] = {"pointlist",
			    "pointattributelist",
			    "pointmarkerlist",
			    "trianglelist",
			    "triangleattributelist",
			    "trianglearealist",
			    "segmentlist",
			    "segmentmarkerlist",
			    "edgelist",
			    "edgemarkerlist",
			    "neighborlist",
			    "normlist",
			    "holelist",
			    "regionlist"};

/*   if(out != NULL) */
/*     mexWarnMsgTxt ("Anything out points to in export_struct will be lost"); */
  ndims = 2;
  dims =  mxCalloc (ndims, sizeof(*dims));
  dims [0] = 1;
  dims [1] = 1;

  nfields = 14;

  out = mxCreateStructArray (2, dims, nfields, fnames);

  for (ifield = 0; ifield < nfields; ifield++)
  {

    /* POINTLIST */
    if (!strncmp(fnames [ifield], "pointlist", 9) &&
	in->pointlist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberofpoints;
      dims [0] = 2;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->pointlist[i];
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* POINTATTRIBUTELIST */
    else if (!strncmp(fnames [ifield], "pointattributelist", 18)
	     && in->pointattributelist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [0] = in->numberofpoints;
      dims [1] = in->numberofpointattributes;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < in->numberofpointattributes*in->numberofpoints; i++)
	pr_dbl [i] = in->pointattributelist[i];
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* POINTMARKERLIST */
    else if (!strncmp(fnames [ifield], "pointmarkerlist", 15)
	     && in->pointmarkerlist !=NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [0] = in->numberofpoints;
      dims [1] = 1;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *) mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->pointmarkerlist[i];
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* TRIANGLELIST */
    else if (!strncmp(fnames [ifield], "trianglelist", 12)
	     && in->trianglelist !=NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberoftriangles;
      dims [0] = in->numberofcorners;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->trianglelist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* TRIANGLEATTRIBUTELIST */
    else if (!strncmp(fnames [ifield], "triangleattributelist", 21)
	     && in->triangleattributelist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [0] = in->numberoftriangleattributes;
      dims [1] = in->numberoftriangles;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->triangleattributelist[i];
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* TRIANGLEAREALIST */
    else if (!strncmp(fnames [ifield], "trianglearealist", 16)
	     && in->trianglearealist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberoftriangles;
      dims [0] = 1;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->trianglearealist[i];
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* SEGMENTLIST */
    else if (!strncmp(fnames [ifield], "segmentlist", 11)
	     && in->segmentlist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [0] = 2;
      dims [1] = in->numberofsegments;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->segmentlist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* SEGMENTMARKERLIST */
    else if (!strncmp(fnames [ifield], "segmentmarkerlist", 17)
	     && in->segmentmarkerlist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberofsegments;
      dims [0] = 1;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->segmentmarkerlist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* EDGELIST */
    else if (!strncmp(fnames [ifield], "edgelist", 8)
	     && in->edgelist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberofedges;
      dims [0] = 2;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->edgelist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* EDGEMERKERLIST */
    else if (!strncmp(fnames [ifield], "edgemarkerlist", 14)
	     && in->edgemarkerlist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberofedges;
      dims [0] = 1;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = in->edgemarkerlist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* NEIGHBORLIST */
    else if (!strncmp(fnames [ifield], "neighborlist", 12)
	     && in->neighborlist != NULL)
    {
      /* mexPrintf("%s\n",fnames[ifield]); */
      dims [1] = in->numberoftriangles;
      dims [0] = 3;
      field = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      pr_int = (int *)mxGetData(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_int [i] = (int) in->neighborlist[i] + 1;
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* NORMLIST */
    else if (!strncmp(fnames [ifield], "normlist", 8)
	     && in->normlist != NULL && 0)
    {
      mexPrintf("%s\n",fnames[ifield]);
      dims [1] = in->numberofpoints;
      dims [0] = 2;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->normlist[i];
      mxSetFieldByNumber (out, 0, ifield, field);

    }
    /* HOLELIST */
    else if (!strncmp(fnames [ifield], "holelist", 8)
	     && in->holelist != NULL)
    {
      mexPrintf("%s\n",fnames[ifield]);
      dims [1] = in->numberofholes;
      dims [0] = 2;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->holelist[i];
      mxSetFieldByNumber (out, 0, ifield, field);
    }
    /* REGIONLIST */
    else if (!strncmp(fnames [ifield], "regionlist", 10)
	     && in->regionlist != NULL)
    {
      mexPrintf("%s\n",fnames[ifield]);
      dims [1] = in->numberofregions;
      dims [0] = 4;
      field = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      mexPrintf("regions: %d\n",in->numberofregions);
      pr_dbl = mxGetPr(field);
      for (i = 0; i < dims [0] * dims [1]; i++)
	pr_dbl [i] = in->regionlist[i];
      mxSetFieldByNumber (out, 0, ifield, field);
    }
  }
  *output=out;
}
