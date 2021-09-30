function [G, rock, N, T] = initGridFromEclipseOutput(init, grid, varargin)
% Reads eclipse output files and converts to MRST-compatible datastructures
%
% SYNOPSIS:
%   [G, rock, N, T] = initGridFromEclipseOutput(init, grid, varargin)
%
% REQUIRED PARAMETERS:
%   init    - structure obtained by reading INIT-file
%   grid    - structure obtined by reading EGRID (or GRID) - file
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
% outputSimGrid - if true, output minimal grid representing the eclipse
%                 output ("topolgoy") suitable for computing TOF etc.,
%                 but not for plotting. Output as second component of G.
%                 Default is false
% RETURNS:
% G     -   if outputSimGrid = false, standard MRST-grid
%           if outputSimGrid = true, G(1) is standard MRST-grid, G(2) is
%           simulation grid
% rock  -   structure compatible with G if nargout < 5, compatible with Gs
%           otherwise, but compatible with both unless special situations
%           where e.g., an aquifer is represented in zero-volume cells
% N, T  -   Neighbour + transmissiblity list typically not compatible with G,
%           but with Gs
% Gs    -   Minimal grid directly

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('outputSimGrid', false);
opt     = merge_options(opt, varargin{:});
outputRock  = nargout >= 2;
outputTrans = nargout >= 3;

% units
units = {'metric', 'field', 'lab'};
unit  = units{init.INTEHEAD.values(3)};
u     = getUnitSystem(unit);

% get grid stuff
gridFromDeck = isfield(grid, 'cartDims');
if gridFromDeck
    cartDims = grid.cartDims;
    anGrid = grid.ACTNUM;
    coord  = grid.COORD;
    zcorn  = grid.ZCORN;
else
    cartDims = grid.GRIDHEAD.values(2:4)';
    anGrid   = grid.ACTNUM.values;
    coord    = grid.COORD.values;
    zcorn    = grid.ZCORN.values;
end

% ACTNUM or PORV > 0
actNum   = and(anGrid, init.PORV.values > 0);

grdecl = struct('cartDims', cartDims, ...
                'COORD'   , convertFrom(coord, u.length), ...
                'ZCORN'   , convertFrom(zcorn, u.length), ...
                'ACTNUM'  , actNum );

dispif(mrstVerbose, 'Creating MRST-grid ...');
try
    mrstModule add libgeometry
    G = mprocessGRDECL(grdecl, 'SplitDisconnected', false);
    G = mcomputeGeometry(G);
catch
    dispif(mrstVerbose, '(calling efficent mex-routines failed ...)');
    G = processGRDECL(grdecl, 'SplitDisconnected', false);
    G = computeGeometry(G);
end
dispif(mrstVerbose, 'done.\n' )

% check for consistency (zero volume cells might represent aquifers etc.)
if G.cells.num ~= nnz(actNum)
    consistentGrids = false;
    tmp = zeros(size(actNum));
    tmp(actNum) = (1:nnz(actNum))';
    eMap    = tmp(G.cells.indexMap);
    eMapInv = zeros(nnz(actNum),1);
    eMapInv(eMap) = (1:G.cells.num)';
else
    consistentGrids = true;
    eMap = ':';
    eMapInv = ':';
end
G.cells.eMap = eMap;
G.cells.eMapInv = eMapInv;
if isfield(init, 'DEPTH')
    G.cells.centroids(:,3) = convertFrom(init.DEPTH.values(eMap), u.length);
end
if isfield(init, 'PORV')
    G.cells.PORV = convertFrom(init.PORV.values(G.cells.indexMap), u.resvolume);
end
if all(isfield(init, {'DX', 'DY', 'DZ'}))
    G.cells.DX    = convertFrom(init.DX.values(eMap), u.length);
    G.cells.DY    = convertFrom(init.DY.values(eMap), u.length);
    G.cells.DZ    = convertFrom(init.DZ.values(eMap), u.length);
end

% if we are outputting simgrid or we have consistency between G and
% ECLIPSE-grid, rock is directly given by init
if opt.outputSimGrid || (consistentGrids && outputRock)
    rock.poro = init.PORO.values;
    if isfield(init, 'NTG')
        rock.ntg  = init.NTG.values;
    end
    rock.perm = convertFrom([init.PERMX.values, ...
        init.PERMY.values, ...
        init.PERMZ.values], u.perm);
elseif outputRock % we want rock to match G
    rock.poro = init.PORO.values(eMap);
    rock.ntg  = init.NTG.values(eMap);
    rock.perm = convertFrom([init.PERMX.values(eMap), ...
        init.PERMY.values(eMap), ...
        init.PERMZ.values(eMap)], u.perm);
end


% transmissibilities and neighbor-list
[N, T] = deal([]);
if outputTrans && ~opt.outputSimGrid
    % use simple short routine
    [N, T] = getActiveNeighbors(init, grid, actNum, cartDims);
elseif opt.outputSimGrid
    % use more detailed routine to produce neccessary indices
    [N, T, ix] = getActiveNeighbors(init, grid, actNum, cartDims);
end
T = convertFrom(T, u.trans);

if opt.outputSimGrid
    dispif(mrstVerbose, 'Creating simulation grid ...');
    Gs = simGridTPFA(G, rock, 'neighbors', N, ...
            'porv',   convertFrom(init.PORV.values(actNum), u.resvolume), ...
            'depth',  convertFrom(init.DEPTH.values, u.length), ...
            'actnum', actNum);
    % set index field with descriptive name ix ...
    Gs.faces.ix = ix;
    % add aquifer-cells if present
    if isfield(init, 'AQUIFERN')
        Gs.cells.aquifer = init.AQUIFERN.values;
    end
    Gs.cells.DX    = convertFrom(init.DX.values(eMap), u.length);
    Gs.cells.DY    = convertFrom(init.DY.values(eMap), u.length);
    Gs.cells.DZ    = convertFrom(init.DZ.values(eMap), u.length);
    Gs.cells.PORV   = convertFrom(init.PORV.values(actNum), u.resvolume);
    dispif(mrstVerbose, 'done\n')
    G = {G, Gs};
end

end

function [N, T, ix] = getActiveNeighbors(init, grid, actNum, cartDims)
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

% if transmissibilities are not given in init-file, but NNC's are given,
% then set TRANNNC to inf and issue warning
if ~isfield(init, 'TRANNNC') && (isfield(grid, 'NNC1') || isfield(init, 'NNC1'))
    warning('Transmissibilities for NNCs not given, values are set to inf!');
    if isfield(grid, 'NNC1')
        n = numel(grid.NNC1.values);
    else
        n = numel(init.NNC1.values);
    end
    init.TRANNNC.values = inf(n,1);
end

%non-neighbouring connections (indices given in full grid -> map to active):
if isfield(init, 'TRANNNC')
    if isfield(grid, 'NNC1')
        NNC1    = actInx(grid.NNC1.values);
        NNC2    = actInx(grid.NNC2.values);
    else
        NNC1    = actInx(init.NNC1.values);
        NNC2    = actInx(init.NNC2.values);
    end
    TRANNNC = init.TRANNNC.values;
else
    NNC1    = [];
    NNC2        = [];
    TRANNNC = [];
end
if nargout < 3
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
else % do more thorough job and produce index maps
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
end
