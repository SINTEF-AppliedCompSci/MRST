template<class Real>
State<Real> *extractStateStruct(const mxArray *state_struct)
{
    if(mxGetNumberOfFields(state_struct) < 5)
        mexErrMsgTxt("The state struct must have at least 5 fields (pressure, \
            flux, s, h and h_max)");

    mxArray *flux_array     = mxGetField(state_struct, 0, "flux");
    mxArray *h_array        = mxGetField(state_struct, 0, "h");
    mxArray *h_max_array    = mxGetField(state_struct, 0, "h_max");
    
    if(!(flux_array && h_array && h_max_array))
        mexErrMsgTxt("State struct must contain flux, h and h_max fields");

    double *flux        = mxGetPr(flux_array);
    double *h           = mxGetPr(h_array);
    double *h_max       = mxGetPr(h_max_array);

    int numCells = mxGetM(h_array);
    int numFaces = mxGetM(flux_array);

    return new State<Real>(numCells, numFaces, flux, h, h_max);
}

template<class Real>
Grid2D<Real>  *extractGrid2DStruct(const mxArray *grid_struct)
{
    if(mxGetNumberOfFields(grid_struct) < 6)
        mexErrMsgTxt("The grid struct must contain at least 7 fields (cells, \
            faces, nodes, cartDims, type, columns)");

    mxArray *cells_struct   = mxGetField(grid_struct, 0, "cells");
    mxArray *faces_struct   = mxGetField(grid_struct, 0, "faces");
    mxArray *nodes_struct   = mxGetField(grid_struct, 0, "nodes");
    mxArray *columns_struct = mxGetField(grid_struct, 0, "columns");
    mxArray *dims_array     = mxGetField(grid_struct, 0, "cartDims");

    if(!(cells_struct && faces_struct && nodes_struct && dims_array &&
                columns_struct))
        mexErrMsgTxt("Grid struct must containt cells, faces, nodes, dims,\
                columns and top faces");

    cells_t<double>     *cells      = extractCellStruct(cells_struct);
    faces_t             *faces      = extractFaceStruct(faces_struct);
    nodes_t<double>     *nodes      = extractNodeStruct(nodes_struct);
    columns_t<double>   *columns    = extractColumnStruct(columns_struct);

    double *dim_tmp = mxGetPr(dims_array);
    int dims[2];
    dims[0] = (int) dim_tmp[0];
    dims[1] = (int) dim_tmp[1];

    Grid2D<Real> *grid = new Grid2D<Real>(cells, faces, nodes, columns, dims);

    freeCellStruct(cells);
    freeFaceStruct(faces);
    freeNodeStruct(nodes);
    freeColumnStruct(columns);

    return grid;
}

template<class Real>
Rock<Real>    *extractRockStruct(const mxArray *rock_struct)
{
    int cells       = mxGetM(mxGetField(rock_struct, 0, "poro"));
    double *perm    = mxGetPr(mxGetField(rock_struct, 0, "perm"));
    double *poro    = mxGetPr(mxGetField(rock_struct, 0, "poro"));

    Rock<Real> *rock = new Rock<Real>(cells, perm, poro);

    return rock;
}

template<class Real>
Fluid<Real>     *extractFluidStruct(const mxArray *fluid_struct)
{
    double *mu_tmp  = mxGetPr(mxGetField(fluid_struct, 0, "mu"));
    double *rho_tmp = mxGetPr(mxGetField(fluid_struct, 0, "rho"));
    double *kwm_tmp = mxGetPr(mxGetField(fluid_struct, 0, "kwm"));

    Real mu[2]    = {(Real) mu_tmp[0], (Real) mu_tmp[1]};
    Real rho[2]   = {(Real) rho_tmp[0], (Real) rho_tmp[1]};
    Real kwm[2]   = {(Real) kwm_tmp[0], (Real) kwm_tmp[1]};
    Real sr       = (Real) mxGetScalar(mxGetField(fluid_struct, 0, "res_gas"));
    Real sw       = (Real) mxGetScalar(mxGetField(fluid_struct, 0, "res_water"));
    Real flux     = (Real) mxGetScalar(mxGetField(fluid_struct, 0,
                "fluxInterp"));

    return new Fluid<Real>(mu, rho, sr, sw, flux, kwm);
}
