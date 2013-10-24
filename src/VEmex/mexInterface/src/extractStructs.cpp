/*
-------------------------------------------------------------------------------
File Name :

Purpose :

Creation Date : 2010-11-01

Last Modified : ma. 13. des. 2010 kl. 14.57 +0100

Created By :  Martin Ertsaas (martiert@student.matnat.uio.no)
-------------------------------------------------------------------------------
*/

#include "extractStructs.h"

Well<double> *extractWell(const mxArray *well_struct)
{
    if(mxGetField(well_struct, 0, "cells")) {
        int     numcells = mxGetM(mxGetField(well_struct, 0, "cells"));

        if(numcells > 1) {
            int     *cells      = (int*) mxGetPr(mxGetField(well_struct, 0,
                        "cells"));
            double  *val        = mxGetPr(mxGetField(well_struct, 0, "val"));
            double  *h          = mxGetPr(mxGetField(well_struct, 0, "h"));

            Well<double> *well = new Well<double>(numcells, cells, val, h);

            return well;
        }

        int     cell        = (int)    mxGetScalar(mxGetField(well_struct, 0,
                    "cells"));
        double  val         = mxGetScalar(mxGetField(well_struct, 0, "val"));
        double  h           = mxGetScalar(mxGetField(well_struct, 0, "h"));


        Well<double> *well = new Well<double>(numcells, &cell, &val, &h);

        return well;
    }

    Well<double> *well = new Well<double>();
    return well;
}

BC<double> *extractBC(const mxArray *bc_struct)
{
    if(mxGetNumberOfFields(bc_struct) == 0)
        return new BC<double>();

    int num   = mxGetM(mxGetField(bc_struct, 0, "face"));

    double *face_arr    = mxGetPr(mxGetField(bc_struct, 0, "face"));
    double *value       = mxGetPr(mxGetField(bc_struct, 0, "value"));
    double *h           = mxGetPr(mxGetField(bc_struct, 0, "h"));

    int *face = (int*) mxMalloc(num*sizeof(int));
    char **type = (char**) mxMalloc(num*sizeof(char*));

#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        face[i] = (int) face_arr[i] - 1;
        type[i] = (char*) mxMalloc(20*sizeof(char));
        mxArray *cell = mxGetCell(mxGetField(bc_struct, 0, "type"), i);
        mxGetString(cell, type[i], 20*sizeof(char));
    }

    BC<double> *bc = new BC<double>(num, face, value, h, type);

    for(int i = 0; i < num; i++)
        mxFree(type[i]);

    mxFree(type);
    mxFree(face);

    return bc;
}

opt_t *extractOptions(const mxArray *opts[], int length)
{
    opt_t *options = (opt_t*) malloc(sizeof(opt_t));

    options->verbose        = false;
    options->intVert_poro   = false;
    options->semi_implicit  = false;
    options->no_dif         = false;
    options->central        = false;
    options->grav_upwind    = false;
    options->computedt      = true;
    options->intVert        = true;
    options->wells          = NULL;
    options->bc             = NULL;
    options->heightWarn     = 1.0e-16;
    options->dt             = 0.0;
    options->flag           = 1;

    char *buf;
    mwSize buflen;
    for(int i = 5; i < length; i += 2) {
        buflen = mxGetN(opts[i])*sizeof(mxChar)+1;
        buf = (char*) mxMalloc(buflen);

        int status = mxGetString(opts[i], buf, buflen);

        if(status)
            mexErrMsgTxt("Expected a string, got something else. Aborting\n");


        if(!strcmp(buf, "dt"))
            options->dt             = (int) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "verbose"))
            options->verbose        = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "computedt"))
            options->computedt      = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "intVert"))
            options->intVert        = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "intVert_poro"))
            options->intVert_poro   = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "semi_implicit"))
            options->semi_implicit  = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "wells"))
            options->wells          = extractWell(opts[i+1]);
        else if(!strcmp(buf, "src"))
            printf("Src not implemented yet. \n");
        else if(!strcmp(buf, "heightWarn"))
            options->heightWarn     = (double) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "no_dif"))
            options->no_dif         = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "central"))
            options->central        = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "gravity_upwind"))
            options->grav_upwind    = (bool) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "gravity"))
            options->gravity        = (double) mxGetScalar(opts[i+1]);
        else if(!strcmp(buf, "flag"))
            options->flag           = (int) mxGetScalar(opts[i+1]);

        mxFree(buf);
    }

    return options;
}

cells_t<double> *extractCellStruct(const mxArray *cell_struct)
{
    cells_t<double> *cells = (cells_t<double>*) mxMalloc(
            sizeof(cells_t<double>));

    cells->num              = (int) mxGetScalar(mxGetField(cell_struct, 0,
                "num"));
    cells->H                = mxGetPr(mxGetField(cell_struct, 0, "H"));
    cells->volumes          = mxGetPr(mxGetField(cell_struct, 0, "volumes"));
    cells->z                = mxGetPr(mxGetField(cell_struct, 0, "z"));
    cells->normals          = mxGetPr(mxGetField(cell_struct, 0, "normals"));
    cells->centroids        = mxGetPr(mxGetField(cell_struct, 0, "centroids"));

    if(mxIsDouble(mxGetField(cell_struct, 0, "faces"))) {
        double *faces = mxGetPr(mxGetField(cell_struct, 0, "faces"));
        cells->faces = (int*) mxMalloc(4*cells->num*sizeof(int));

        for(int i = 0; i < 4*cells->num; i++) {
            cells->faces[i] = (int) faces[i];
        }
    } else {
        cells->faces           = (int*) mxGetPr(mxGetField(cell_struct, 0,
                    "faces"));
    }

    mxArray* numFaces       = mxGetField(cell_struct, 0, "numFaces");

    if(numFaces) {
        cells->numFaces         = (int*)    mxGetPr(numFaces);
    } else {
        cells->numFaces = (int*) mxMalloc(cells->num*sizeof(int));
        for(int i = 0; i < cells->num; i++)
            cells->numFaces[i] = 4;
    }

    double *columnPos       = mxGetPr(mxGetField(cell_struct, 0, "columnPos"));
    cells->columnPos        = (int*)    mxMalloc((cells->num+1)*sizeof(int));

    cells->columnPos[0] = (int) columnPos[0];

#pragma omp parallel for
    for(int i = 0; i < cells->num; i++) {
        cells->columnPos[i+1]       = (int) columnPos[i+1];
    }

    return cells;
}

void freeCellStruct(cells_t<double> *cells)
{
    mxFree(cells->columnPos);
    mxFree(cells);
}

faces_t *extractFaceStruct(const mxArray *face_struct)
{
    faces_t *faces = (faces_t*) mxMalloc(sizeof(faces_t));

    faces->num          = (int)     mxGetScalar(mxGetField(face_struct, 0,
                "num"));
    faces->nodes        = (int*)    mxGetPr(mxGetField(face_struct, 0,
                "nodes"));
    faces->neighbors    = (int*)    mxGetPr(mxGetField(face_struct, 0,
                "neighbors"));
    faces->nodes        = (int*)    mxGetPr(mxGetField(face_struct, 0,
                "nodes"));

    return faces;
}

void freeFaceStruct(faces_t *faces)
{
    mxFree(faces);
}

nodes_t<double> *extractNodeStruct(const mxArray *node_struct)
{
    nodes_t<double> *nodes = (nodes_t<double>*) mxMalloc(
            sizeof(nodes_t<double>));

    nodes->num      = (int) mxGetScalar(mxGetField(node_struct, 0, "num"));
    nodes->coords   = mxGetPr(mxGetField(node_struct, 0, "coords"));

    return nodes;
}

void freeNodeStruct(nodes_t<double> *nodes)
{
    mxFree(nodes);
}

columns_t<double> *extractColumnStruct(const mxArray *column_struct)
{
    columns_t<double> *columns = (columns_t<double>*) mxMalloc(
            sizeof(columns_t<double>));

    columns->num    = mxGetM(mxGetField(column_struct, 0, "cells"));
    double *cells   = mxGetPr(mxGetField(column_struct, 0, "cells"));
    columns->dz     = mxGetPr(mxGetField(column_struct, 0, "dz"));
    columns->z      = mxGetPr(mxGetField(column_struct, 0, "z"));

    columns->cells  = (int*) mxMalloc(columns->num*sizeof(int));

#pragma omp parallel for
    for(int i = 0; i < columns->num; i++)
        columns->cells[i] = (int) cells[i];

    return columns;
}

void freeColumnStruct(columns_t<double> *columns)
{
    mxFree(columns->cells);
    mxFree(columns);
}
