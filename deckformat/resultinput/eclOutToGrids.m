function [G, Gs, rock, info] = eclOutToGrids(init, grid, varargin)
% Reads eclipse output files and converts to MRST-compatible datastructures
%
% SYNOPSIS:
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
% RETURNS:
% G -    MRST compatible grid representing the eclipse grid, but which is 
%        suitable for plotting ("topology and geometry").
%
% Gs -   Grid directly representing the eclipse output ("topolgoy"), which 
%        is suitable for computing TOF etc.
%
% rock - Rock parameters (permeability / porosity / etc.)
%
% info - 
%
% EXAMPLE:
%
% SEE ALSO:

opt = struct('outputRock',     false, ...
             'outputSimGrid',  true);
             
opt     = merge_options(opt, varargin{:});


% units
units = {'metric', 'field', 'lab'};
unit  = units{init.INTEHEAD.values(3)};
u     = getUnitFacs(unit);

% get grid stuff
cartDims = grid.GRIDHEAD.values(2:4)';
anGrid   = grid.ACTNUM.values;
coord    = grid.COORD.values;
zcorn    = grid.ZCORN.values;

% ACTNUM or PORV > 0 
actNum   = and(anGrid, init.PORV.values > 0);

grdecl = struct('cartDims', cartDims, ...
                'COORD'   , convertFrom(coord, u.length), ...
                'ZCORN'   , convertFrom(zcorn, u.length), ...
                'ACTNUM'  , actNum );

dispif(mrstVerbose, 'Creating MRST-grid ...');
try 
    mrstModule add libgeometry mex opm_gridprocessing
    G = mprocessGRDECL(grdecl, 'SplitDisconnected', false);
    G = mcomputeGeometry(G);
catch
    dispif(mrstVerbose, '(calling efficent mex-routines failed ...)');
    G = processGRDECL(grdecl, 'SplitDisconnected', false);
    G = computeGeometry(G);
end
dispif(mrstVerbose, 'done.\n' )

% check for consistency
if G.cells.num ~= nnz(actNum)
    tmp = zeros(size(actNum));
    tmp(actNum) = (1:nnz(actNum))';
    eMap    = tmp(G.cells.indexMap);
    eMapInv = zeros(nnz(actNum),1);
    eMapInv(eMap) = (1:G.cells.num)';
else
    eMap = ':';
    eMapInv = ':';
end   
G.cells.eMap = eMap; 
G.cells.eMapInv = eMapInv; 

%G.cells.centroids(:,3) = convertFrom(init.DEPTH.values(eMap), u.length);
%G.cells.PORV           = convertFrom(init.PORV.values(G.cells.indexMap), u.resvol);
% probably don't need DX/DY/DZ
%G.DX    = convertFrom(init.DX.values(eMap), u.length);
%G.DY    = convertFrom(init.DY.values(eMap), u.length);
%G.DZ    = convertFrom(init.DZ.values(eMap), u.length);

% rock should be of size equal to Gs
if or(opt.outputRock, opt.outputSimGrid) 
    rock.poro = init.PORO.values;
    rock.ntg  = init.NTG.values;
    rock.perm = convertFrom([init.PERMX.values, ...
                             init.PERMY.values, ...
                             init.PERMZ.values], u.perm);
else
    rock = [];
end

if opt.outputSimGrid
    dispif(mrstVerbose, 'Creating simulation grid ...');
    [N, ix, T] = getActiveNeighbors(init, grid, actNum, cartDims);
    %[N, T] = getTrans(init, grid, actNum, cartDims);
    Gs = simGridTPFA(G, rock, 'neighbors', N, ...
            'porv',   convertFrom(init.PORV.values(actNum), u.resvol), ...
            'depth',  convertFrom(init.DEPTH.values, u.length), ...
            'actnum', actNum);
    % set index field with descriptive name ix ...
    Gs.faces.ix = ix;
    % add aquifer-cells if present
    if isfield(init, 'AQUIFERN')
        Gs.cells.aquifer = init.AQUIFERN.values;
    end
    dispif(mrstVerbose, 'done\n')
else
    Gs = [];
end

if nargout > 3
    info.T = convertFrom(T, u.trans);
end
end

function u = getUnitFacs(unit)
switch unit
    case 'metric'
        u.length  = meter;
        u.resvol  = meter^3;
        u.perm    = milli*darcy;
        u.trans   = centi*poise * meter^3 / (day * barsa);
    case 'field'
        u.length  = ft;
        u.resvol  = stb;
        u.perm    = milli*darcy;
        u.trans   = centi*poise * stb / (day * psia);
    case 'lab'
        u.length  = centi*meter;
        u.resvol  = (centi*meter)^3;
        u.perm    = milli*darcy;
        u.trans   = centi*poise * (centi*meter)^3 / (hour * atm);
    otherwise
        error(['Unknown unit:', unit])
end
end

function [N, ix, T] = getActiveNeighbors(init, grid, actNum, cartDims)
nCells = nnz(actNum);

%full grid to active grid mapping
actInx = zeros(prod(cartDims), 1);
actInx(actNum) = (1:nCells)';

% temporary matrix M to find indices of active neighbors in X, Y and Z-direction,
M = zeros(cartDims+1);
M(1:end-1, 1:end-1, 1:end-1) = reshape(actInx, cartDims);

NX = M(2:end, 1:end-1, 1:end-1); NX = NX(actNum);
NY = M(1:end-1, 2:end, 1:end-1); NY = NY(actNum);
NZ = M(1:end-1, 1:end-1, 2:end); NZ = NZ(actNum);
clear M

%non-neighbouring connections (indices given in full grid -> map to active):
if isfield(init, 'TRANNNC')
    if isfield(grid, 'NNC1')
        NNC1    = actInx(grid.NNC1.values);
        NNC2 	= actInx(grid.NNC2.values);
    else
        NNC1    = actInx(init.NNC1.values);
        NNC2 	= actInx(init.NNC2.values);
    end
    TRANNNC = init.TRANNNC.values;
else
    NNC1    = [];
    NNC2 	= [];
    TRANNNC = [];
end

% produce neighbors and corresponding maps for I,J,K and NNC
% I - direction
NI  = [(1:nCells)', NX(:)];
TI  = init.TRANX.values; 
inx = and(TI>0, prod(NI,2)>0);
ix.iIe = inx;
ix.iI  = (1:nnz(inx))';
ix.sI  = 1;
% J - direction
NJ  = [(1:nCells)', NY(:)];
TJ  = init.TRANY.values; 
inx = and(TJ>0, prod(NJ,2)>0);
ix.iJe = inx;
ix.iJ  = (1:nnz(inx))' + numel(ix.iI);
ix.sJ  = 1;
% K - direction
NK  = [(1:nCells)', NZ(:)];
TK  = init.TRANZ.values; 
inx = and(TK>0, prod(NK,2)>0);
ix.iKe = inx;
ix.iK  = (1:nnz(inx))' + numel(ix.iI) + numel(ix.iJ);
ix.sK  = 1;
% NNCs
NN = [NNC1, NNC2];
TN = TRANNNC;
isPos = NNC1 < NNC2;
if ~isempty(NN)
    NN(~isPos,:) = NN(~isPos, [2, 1]);
    %inx = and(TN>0, prod(NN,2)>0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inx = TN>0;
else
    inx = [];
end
ix.iNe = inx;
ix.iN  = (1:nnz(inx))' + numel(ix.iI) + numel(ix.iJ) + numel(ix.iK);
ix.sN  = 2*isPos(inx) -1;

% produce neighbour list N
N = [NI(ix.iIe,:); ...
     NJ(ix.iJe,:); ...
     NK(ix.iKe,:); ...
     NN(ix.iNe,:)];

% also transmissibilities if required
if nargout == 3
    T = [TI(ix.iIe,:); ...
         TJ(ix.iJe,:); ...
         TK(ix.iKe,:); ...
         TN(ix.iNe,:)];
end
end


function [N, T] = getTrans(init, grid, actNum, cartDims)
nCells = nnz(actNum);

%full grid to active grid mapping
actInx = zeros(prod(cartDims), 1);
actInx(actNum) = (1:nCells)';

% temporary matrix M to find indices of active neighbors in X, Y and Z-direction,
M = zeros(cartDims+1);
M(1:end-1, 1:end-1, 1:end-1) = reshape(actInx, cartDims);

NX = M(2:end, 1:end-1, 1:end-1); NX = NX(actNum);
NY = M(1:end-1, 2:end, 1:end-1); NY = NY(actNum);
NZ = M(1:end-1, 1:end-1, 2:end); NZ = NZ(actNum);

%non-neighbouring connections (indices given in full grid -> map to active):
if isfield(init, 'TRANNNC')
    if isfield(grid, 'NNC1')
        NNC1    = actInx(grid.NNC1.values);
        NNC2 	= actInx(grid.NNC2.values);
    else
        NNC1    = actInx(init.NNC1.values);
        NNC2 	= actInx(init.NNC2.values);
    end
    TRANNNC = init.TRANNNC.values;
else
    NNC1    = [];
    NNC2 	= [];
    TRANNNC = [];
end


% produce neighbour list N, and transmissibilities T
N = [(1:nCells)' NX; ...
     (1:nCells)' NY; ...
     (1:nCells)' NZ; ...
     NNC1 NNC2];
T = [init.TRANX.values; ...
     init.TRANY.values; ...
     init.TRANZ.values; ...
     TRANNNC];

% remove 0-transmissibilities and non-active neighbors (appearing as zeros)
inx = and(T>0, prod(N,2)>0);
N = N(inx,:);
T = T(inx);
end