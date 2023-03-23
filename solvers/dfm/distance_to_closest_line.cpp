/*********************************************************************
 * distance_to_closest_line.cpp
 * 
 * For a set of points [px,py] find the closest line
 * [x0,y0,x1,y1] and return its index as well as its
 * distance. 
 *
 *
 * Copyright 2011-2012 University of Bergen
 * This file is licenced under the GNU General Public Licence v3.0
 ********************************************************************/
#include <matrix.h>
#include <mex.h> 
#include <math.h>

double getd(double px, double py, double x0,
        double y0,double x1,double y1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *a_in_m, *b_in_m, *c_in_m, *d_in_m, *e_in_m, *f_in_m,*g_out_m,*h_out_m;
    const mwSize *dim_points_m;
    const mwSize *dim_edges_m;
    double *a, *b, *c, *d, *e, *f, *g, *h;
    int num_points, num_edges;
    int i,j;
  

//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);
    c_in_m = mxDuplicateArray(prhs[2]);
    d_in_m = mxDuplicateArray(prhs[3]);
    e_in_m = mxDuplicateArray(prhs[4]);
    f_in_m = mxDuplicateArray(prhs[5]);

//figure out dimensions
    dim_points_m = mxGetDimensions(prhs[0]);
    num_points = (int)dim_points_m[0];
    
    dim_edges_m = mxGetDimensions(prhs[2]);
    num_edges = (int)dim_edges_m[0];

//associate outputs
    g_out_m = plhs[0] = mxCreateDoubleMatrix(num_points,1,mxREAL);
    h_out_m = plhs[1] = mxCreateDoubleMatrix(num_points,1,mxREAL);

//associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_in_m);
    d = mxGetPr(d_in_m);
    e = mxGetPr(e_in_m);
    f = mxGetPr(f_in_m);
    g = mxGetPr(g_out_m);
    h = mxGetPr(h_out_m);

//do something
    
    double dist_new;
    
    
    for(i=0;i<num_points;i++)
    {
        g[i] = 99999;
        for(j=0;j<num_edges;j++)
        {
            dist_new = getd(a[i],b[i],c[j],d[j],e[j],f[j]); 
            if (dist_new < g[i])
            {
                g[i] = dist_new;
                h[i] = j+1;
            }
            
        }
    }

    return;

}

double getd(double px, double py, double x0,
        double y0,double x1,double y1)
{
    
    double vx = x1 - x0;
    double vy = y1 - y0;
    double ax = px - x0;
    double ay = py - y0;
    double bx = px - x1;
    double by = py - y1;
    double normv = vx * vx + vy*vy;
    double t = (ax * vx + ay*vy)/normv;
    
    
    
    double projx = t * vx + x0;
    double projy = t * vy + y0;
    
    double cx = px - projx;
    double cy = py - projy;
    
    
    double a = ax * ax + ay * ay;
    double b = bx * bx + by * by;
    double c = cx * cx + cy * cy;
    
    a = sqrt (a);
    b = sqrt (b);
    c = sqrt (c);
    
    double d = 0;
    if (t < 0 )
    {
        d = d +a;
    }
    if (t > 1 )
    {
        d = d +b;
    }
    if (t >= 0 && t <=1 )
    {
        d = d +c;
    }
    return d;
}

