function [G, rock, N, T] = eclOut2mrst(init, grid)
%Undocumented Utility Function

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

units = {'metric', 'field', 'lab'};
unit  = units{init.INTEHEAD.values(3)};
u     = getUnitFacs(unit);

isInput = isfield(grid, 'cartDims');
if isInput
    cartDims = grid.cartDims;
    anGrid = grid.ACTNUM;
    coord  = grid.COORD;
    zcorn  = grid.ZCORN;
else
    cartDims = grid.GRIDHEAD.values(2:4)';
    anGrid = grid.ACTNUM.values;
    coord  = grid.COORD.values;
    zcorn  = grid.ZCORN.values;
end

actNum = and(anGrid, init.PORV.values > 0);
grdecl = struct('cartDims', cartDims, ...
                'COORD'   , convertFrom(coord, u.length), ...
                'ZCORN'   , convertFrom(zcorn, u.length), ...
                'ACTNUM'  , int32(actNum));
try 
    mrstModule add libgeometry mex opm_gridprocessing
    G = mprocessGRDECL(grdecl, 'SplitDisconnected', false);
    G = mcomputeGeometry(G);
catch
    G = processGRDECL(grdecl, 'SplitDisconnected', false);
    G = computeGeometry(G);
end
G.cells.centroids(:,3) = convertFrom(init.DEPTH.values, u.length);
G.PORV  = convertFrom(init.PORV.values(actNum), u.resvol);
G.DX    = convertFrom(init.DX.values, u.length);
G.DY    = convertFrom(init.DY.values, u.length);
G.DZ    = convertFrom(init.DZ.values, u.length);

rock.poro = init.PORO.values;
rock.ntg  = init.NTG.values;
rock.perm = convertFrom([init.PERMX.values, init.PERMX.values, init.PERMX.values], u.perm);

[N, T] = getTrans(init, grid, actNum, cartDims);
T      = convertFrom(T, u.trans);

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
