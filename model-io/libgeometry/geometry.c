/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#include <math.h>
#include <stdio.h>
#include "geometry.h"

/* ------------------------------------------------------------------ */
static void
cross(const double u[3], const double v[3], double w[3])
/* ------------------------------------------------------------------ */
{
   w[0] = u[1]*v[2]-u[2]*v[1];
   w[1] = u[2]*v[0]-u[0]*v[2];
   w[2] = u[0]*v[1]-u[1]*v[0];
}

/* ------------------------------------------------------------------ */
static double
norm(const double w[3])
/* ------------------------------------------------------------------ */
{
   return sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
}



/* ------------------------------------------------------------------ */
void
compute_face_geometry(int ndims, double *coords, int nfaces,
                      int *nodepos, int *facenodes, double *fnormals,
                      double *fcentroids, double *fareas)
/* ------------------------------------------------------------------ */
{

   /* Assume 3D for now */
   int f;
   double x[3];
   double u[3];
   double v[3];
   double w[3];

   int i,k;
   int node;

   double cface[3]  = {0};
   double n[3]  = {0};
   const double twothirds = 0.666666666666666666666666666667;
   for (f=0; f<nfaces; ++f)
   {
      int    num_face_nodes;
      double area = 0.0;
      for(i=0; i<ndims; ++i) x[i] = 0.0;
      for(i=0; i<ndims; ++i) n[i] = 0.0;
      for(i=0; i<ndims; ++i) cface[i] = 0.0;

      /* average node */
      for(k=nodepos[f]; k<nodepos[f+1]; ++k)
      {
         node = facenodes[k];
         for (i=0; i<ndims; ++i) x[i] += coords[3*node+i];
      }
      num_face_nodes = nodepos[f+1] - nodepos[f];
      for(i=0; i<ndims; ++i) x[i] /= num_face_nodes;



      /* compute first vector u (to the last node in the face) */
      node = facenodes[nodepos[f+1]-1];
      for(i=0; i<ndims; ++i) u[i] = coords[3*node+i] - x[i];


      /* Compute triangular contrib. to face normal and face centroid*/
      for(k=nodepos[f]; k<nodepos[f+1]; ++k)
      {
         double a;

         node = facenodes[k];
         for (i=0; i<ndims; ++i) v[i] = coords[3*node+i] - x[i];

         cross(u,v,w);
         a = 0.5*norm(w);
         area += a;
         if(!(a>0))
         {
            fprintf(stderr, "Internal error in compute_face_geometry.");
         }

         /* face normal */
         for (i=0; i<ndims; ++i) n[i] += w[i];

         /* face centroid */
         for (i=0; i<ndims; ++i)
            cface[i] += a*(x[i]+twothirds*0.5*(u[i]+v[i]));

         /* Store v in u for next iteration */
         for (i=0; i<ndims; ++i) u[i] = v[i];
      }

      /* Store face normal and face centroid */
      for (i=0; i<ndims; ++i)
      {
         /* normal is scaled with face area */
         fnormals  [3*f+i] = 0.5*n[i];
         fcentroids[3*f+i] = cface[i]/area;
      }
      fareas[f] = area;
   }
}


/* ------------------------------------------------------------------ */
void
compute_cell_geometry(int ndims, double *coords,
                      int *nodepos, int *facenodes, int *neighbors,
                      double *fnormals,
                      double *fcentroids,
                      int ncells, int *facepos, int *cellfaces,
                      double *ccentroids, double *cvolumes)
/* ------------------------------------------------------------------ */
{

   int i,k, f,c;
   int face,node;
   double x[3];
   double u[3];
   double v[3];
   double w[3];
   double xcell[3];
   double ccell[3];
   double cface[3]  = {0};
   const double twothirds = 0.666666666666666666666666666667;

   int ndigits;

   ndigits = ((int) (log(ncells) / log(10.0))) + 1;


   for (c=0; c<ncells; ++c)
   {
      int num_faces;
      double volume = 0.0;

      for(i=0; i<ndims; ++i) xcell[i] = 0.0;
      for(i=0; i<ndims; ++i) ccell[i] = 0.0;


      /*
       * Approximate cell center as average of face centroids
       */
      for(f=facepos[c]; f<facepos[c+1]; ++f)
      {
         face = cellfaces[f];
         for (i=0; i<ndims; ++i) xcell[i] += fcentroids[3*face+i];
      }
      num_faces = facepos[c+1] - facepos[c];

      for(i=0; i<ndims; ++i) xcell[i] /= num_faces;



      /*
       * For all faces, add tetrahedron's volume and centroid to
       * 'cvolume' and 'ccentroid'.
       */
      for(f=facepos[c]; f<facepos[c+1]; ++f)
      {
         int num_face_nodes;

         for(i=0; i<ndims; ++i) x[i] = 0.0;
         for(i=0; i<ndims; ++i) cface[i] = 0.0;

         face = cellfaces[f];

         /* average face node x */
         for(k=nodepos[face]; k<nodepos[face+1]; ++k)
         {
            node = facenodes[k];
            for (i=0; i<ndims; ++i) x[i] += coords[3*node+i];
         }
         num_face_nodes = nodepos[face+1] - nodepos[face];
         for(i=0; i<ndims; ++i) x[i] /= num_face_nodes;



         /* compute first vector u (to the last node in the face) */
         node = facenodes[nodepos[face+1]-1];
         for(i=0; i<ndims; ++i) u[i] = coords[3*node+i] - x[i];


         /* Compute triangular contributions to face normal and face centroid */
         for(k=nodepos[face]; k<nodepos[face+1]; ++k)
         {
            double tet_volume, subnormal_sign;
            node = facenodes[k];
            for (i=0; i<ndims; ++i) v[i] = coords[3*node+i] - x[i];

            cross(u,v,w);

            tet_volume = 0;
            for(i=0; i<ndims; ++i){
                tet_volume += w[i] * (x[i] - xcell[i]);
            }
            tet_volume *= 0.5 / 3;

            subnormal_sign=0.0;
            for(i=0; i<ndims; ++i){
                subnormal_sign += w[i]*fnormals[3*face+i];
            }

            if (subnormal_sign < 0.0) {
                tet_volume = - tet_volume;
            }
            if (neighbors[2*face + 0] != c) {
                tet_volume = - tet_volume;
            }

            volume += tet_volume;

            /* face centroid of triangle  */
            for (i=0; i<ndims; ++i) cface[i] = (x[i]+twothirds*0.5*(u[i]+v[i]));

            /* Cell centroid */
            for (i=0; i<ndims; ++i) ccell[i] += tet_volume * 3/4.0*(cface[i] - xcell[i]);


            /* Store v in u for next iteration */
            for (i=0; i<ndims; ++i) u[i] = v[i];
         }
      }

      if (! (volume > 0.0)) {
          fprintf(stderr,
                  "Internal error in mex_compute_geometry(%*d): "
                  "negative volume\n", ndigits, c);
      }

      for (i=0; i<ndims; ++i) ccentroids[3*c+i] = xcell[i] + ccell[i]/volume;
      cvolumes[c] = volume;
   }
}
