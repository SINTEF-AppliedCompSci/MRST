function [T, A, q] = computeTimeOfFlight_DFM(state, G, rock,  varargin)
%Compute time of flight using finite-volume scheme.
%
% MODIFIED FROM computeTimeOfFlight.m to allow for cell 2 cell connections
% 	Copyright (C) in portion 2013 IRIS AS
%
%
% SYNOPSIS:
%    T        = computeTimeOfFlight(state, G, rock)
%    T        = computeTimeOfFlight(state, G, rock, 'pn1', pv1, ...)
%   [T, A]    = computeTimeOfFlight(...)
%   [T, A, q] = computeTimeOfFlight(...)
%
% DESCRIPTION:
%   Compute time of flight by solving
%
%       \nabla·(vT) = \phi
%
%   using a first-order finite-volume method with upwind flux.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.  Must contain valid cell interface fluxes,
%           'state.flux'.
%
% OPTIONAL PARAMETERS:
%   wells - Well structure as defined by function 'addWell'.  May be empty
%           (i.e., wells = []) which is interpreted as a model without any
%           wells.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = []) which is
%           interpreted as a reservoir model without explicit sources.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions
%           to the reservoir flow.  May be empty (i.e., bc = []) which is
%           interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   'reverse' - Reverse the fluxes and rates.
%
%   'tracer'  - Cell-array of cell-index vectors for which to solve tracer
%               equation. One equation is solved for each vector with
%               tracer injected in cells given indices. Each vector adds
%               one additional RHS to the original tof-system. Output given
%               as additional columns in T.
% RETURNS:
%   T - Cell values of a piecewise constant approximation to time-of-flight
%       computed as the solution of the boundary-value problem
%
%           (*)    \nabla·(vT) = \phi
%
%       using a finite-volume scheme with single-point upwind approximation
%       to the flux.
%
%   A - Discrete left-hand side of (*), a G.cells.num-by-G.cells.num matrix
%       whose entries are
%
%           A_ij = min(F_ij, 0), and
%           A_ii = sum_j max(F_ij, 0) + max(q_i, 0),
%
%       where F_ij = -F_ji is the flux from cell i to cell j
%
%           F_ij = A_ij·n_ij·v_ij.
%
%       and n_ij is the outward-pointing normal of cell i for grid face ij.
%
%       OPTIONAL.  Only returned if specifically requested.
%
%   q - Aggregate source term contributions (per grid cell) from wells,
%       explicit sources and boundary conditions.  These are the
%       contributions referred to as 'q_i' in the definition of the matrix
%       elements, 'A_ii'.  Measured in units of m^3/s.
%
%       OPTIONAL.  Only returned if specifically requested.
%
% SEE ALSO:
%   `simpleTimeOfFlight`, `solveIncompFlow`.

%{
Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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


opt = struct('bc', [], 'src', [], 'wells', [], 'reverse', false, 'tracer', {{}});
opt = merge_options(opt, varargin{:});

assert (isempty(opt.src) || ~isempty(opt.src.sat), ...
   'Source terms must have a valid ''sat'' field.');
assert (isempty(opt.bc) || ~isempty(opt.bc.sat), ...
   'Boundary conditions must have a valid ''sat'' field.');
assert (isfield(rock, 'poro')         && ...
        numel(rock.poro)==G.cells.num,   ...
        ['The rock input must have a field poro with porosity ',...
         'for each cell in the grid.']);
assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');

tr = opt.tracer;
if ~iscell(tr), tr = {tr}; end

% Find external sources of inflow
q  = computeTransportSourceTerm(state, G, opt.wells, opt.src, opt.bc);

if opt.reverse,
   q = -q;
   state.flux = -state.flux;
end

% Build upwind flux matrix with v_ii = max(flux_ij, 0) + max(qi, 0) and
% v_ij = min(flux_ij, 0).
i  = ~any(G.faces.neighbors==0, 2);
n  = double(G.faces.neighbors(i,:));
flux = state.flux(i);

if isfield(G.cells,'neighbors')
    n = [n ; G.cells.neighbors];
    flux = [flux; state.fluxc2c];

end
nc = G.cells.num;
qp = max(q, 0);
A  = sparse(n(:,1), n(:,2),  max(flux, 0), nc, nc)...
   + sparse(n(:,2), n(:,1), -min(flux, 0), nc, nc);
A  = -A' + spdiags(sum(A)' + qp, 0, nc, nc);

pv = poreVolume(G, rock);

% Subtract divergence of non-source cells from rhs
cellNo = rldecode(1:G.cells.num,diff(G.cells.facePos),2)';
cflux = faceFlux2cellFlux(G, state.flux);

if isfield(G.cells,'neighbors')
    cellNo = [cellNo; G.cells.neighbors(:)];
    cflux = [cflux; state.fluxc2c; -state.fluxc2c];
end

div = accumarray(cellNo, cflux);
% The q below is added to the line below to make it a no-op for source
% cells, div will be equal to q for such cells.
pv = pv - div + q;

% build RHSs for tracer equations
numTrRHS = numel(tr);
numCells = cellfun(@numel, tr);
TrRHS    = full( sparse(vertcat(tr{:}), ...
                 rldecode((1:numTrRHS)', numCells(:)), ...
                 1, nc, numTrRHS ));
TrRHS    = bsxfun(@times, TrRHS, full(qp));

% Time of flight for a divergence-free velocity field.
T  = A \ [pv TrRHS];
end
