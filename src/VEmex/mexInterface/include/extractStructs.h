#ifndef EXTRACTSTRUCTS_H
#define EXTRACTSTRUCTS_H

// Hack to compile on windows with Visual Studio
#if (_MSC_VER >= 1600)
#include <yvals.h>
#define __STDC_UTF_16__
#endif

#include "mex.h"

#include <omp.h>

#include "State.h"
#include "Grid2D.h"
#include "Rock.h"
#include "Fluid.h"
#include "Well.h"
#include "BC.h"
#include "Options.h"

Well<double>    *extractWell(const mxArray *well_struct);
BC<double>      *extractBC(const mxArray *bc_struct);
opt_t   *extractOptions(const mxArray *opts[], int length);

cells_t<double> *extractCellStruct(const mxArray *cell_struct);
void freeCellStruct(cells_t<double> *cells);

faces_t *extractFaceStruct(const mxArray *face_struct);
void freeFaceStruct(faces_t *faces);

nodes_t<double> *extractNodeStruct(const mxArray *node_struct);
void freeNodeStruct(nodes_t<double> *nodes);

columns_t<double> *extractColumnStruct(const mxArray *column_struct);
void freeColumnStruct(columns_t<double> *columns);

#include "tmpl/extractStructs_tmpl.h"
#endif
