function [T, A, q] = computeTimeOfFlight(state, G, rock,  varargin)
%Compute time of flight using finite-volume scheme.
%
% SYNOPSIS:
%    T        = computeTimeOfFlight(state, G, rock)
%    T        = computeTimeOfFlight(state, G, rock, 'pn1', pv1, ...)
%   [T, A]    = computeTimeOfFlight(...)
%   [T, A, q] = computeTimeOfFlight(...)
%
% DESCRIPTION:
%   Compute time-of-flight by solving
%
%       v·\nabla T = \phi
%
%   using a first-order finite-volume method with upwind flux. Here, 'v' is
%   the Darcy velocity, '\phi' is the porosity, and T is time-of-flight.
%   The time-of-flight T(x) is the time it takes an inert particles that is
%   passively advected with the fluid to travel from the nearest inflow
%   point to the point 'x' inside the reservoir. For the computation to
%   make sense, inflow points must be specified in terms of inflow
%   boundaries, source terms, wells, or combinations of these.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir state structure, must contain valid cell interface
%           fluxes, 'state.flux'. Typically, 'state' will contain the
%           solution produced by a flow solver like 'incompTPFA'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
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
%   reverse - Boolean variable. If true, the reverse time-of-flight will be
%           computed, i.e., the travel time from a cell inside the
%           reservoir to the nearest outflow point (well, source, or
%           outflow boundary). Default: FALSE, i.e, compute forward tof.
%
%   tracer - Cell-array of cell-index vectors for which to solve a
%           stationary tracer equation:
%
%             \nabla·(v C) = 0,      C(inflow,t)=1
%
%           One equation is solved for each vector with tracer injected in
%           the cells with the given indices. Each vector adds one
%           additional right-hand side to the original tof-system. Output
%           given as additional columns in T.
%
%
% RETURNS:
%   T - Cell values of a piecewise constant approximation to time-of-flight
%       computed as the solution of the boundary-value problem
%
%           (*)    v · \nabla T = \phi
%
%       using a finite-volume scheme with single-point upwind approximation
%       to the flux.
%
%   A - Discrete left-hand side of (*), a G.cells.num-by-G.cells.num matrix
%       whose entries are
%
%           A_ij = min(F_ij, 0), and
%           A_ii = sum_j max(F_ij, 0) + max(2*q_i, 0),
%
%       where F_ij = -F_ji is the flux from cell i to cell j
%
%           F_ij = A_ij·n_ij·v_ij.
%
%       and n_ij is the outward-pointing normal of cell i for grid face ij.
%       The discretization uses a simple model for cells containing inflow.
%       If q_i denotes the rate and V_i the volume of cell i, then T_i is
%       set to half the time it takes to fill the cell, T_i = V_i/(2q_i). 
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
%   simpleTimeOfFlight, solveIncompFlow.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


opt = struct('bc',          [], ...
             'src',         [], ...
             'wells',       [], ...
             'reverse',     false,...
             'allowInf',    false, ...
             'tracer',      {{}});
opt = merge_options(opt, varargin{:});

assert (~all([isempty(opt.src), isempty(opt.bc), isempty(opt.wells)]), ...
    'Must have inflow described as boundary conditions, sources, or wells');
assert (isfield(rock, 'poro')         && ...
        numel(rock.poro)==G.cells.num,   ...
        ['The rock input must have a field poro with porosity ',...
         'for each cell in the grid.']);
assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');

tr = opt.tracer;
if ~iscell(tr), tr = {tr}; end

% Find external sources of inflow: q contains the contribution from from
% sources (src) and wells (W), while qb contains the cell-wise sum of
% fluxes from boundary conditions.
[q,qb] = computeSourceTerm(state, G, opt.wells, opt.src, opt.bc);

pv = poreVolume(G, rock);

if opt.reverse,
   q  = -q;
   qb = -qb;
   state.flux = -state.flux;
end

% Build upwind flux matrix in which we define v_ji = max(flux_ij, 0) and
% v_ij = -min(flux_ij, 0). Then the diagonal of the discretization matrix
% is obtained by summing rows in the upwind flux matrix. This will give the
% correct diagonal in all cells except for those with a positive fluid
% source. In these cells, the average time-of-flight is set to half the
% time it takes to fill the cell, which means that the diagonal entry
% should be equal twice the fluid rate inside these cells.
i  = ~any(G.faces.neighbors==0, 2);
n  = double(G.faces.neighbors(i,:));
nc = G.cells.num;
qp = max(q+qb, 0);

out = min(state.flux(i), 0);
in  = max(state.flux(i), 0);

% Cell wise total out/inflow
outflow = accumarray([n(:, 2); n(:, 1)], [out; -in]);
inflow  = accumarray([n(:, 2); n(:, 1)], [in; -out]);

% The diagonal entries are equal to the sum of outfluxes, plus the absolute
% value of any source terms to ensure that the cell has the value after it
% is half-filled
d = -(outflow - abs(q + qb));

% Any cells that have no inflow and no source terms are problematic for the
% formulation. These cells can occur from a non-converged pressure solution
% or if the problem is not incompressible.
isSingular = full(inflow == 0 & q == 0 & qb == 0);
% Set inflow/outflow for those cells to zero
in(isSingular(n(:, 1)))  = 0;
out(isSingular(n(:, 2))) = 0;

% Set diagonal to 1 and pv to one so that the cell gets infinite tof
% without producing an ill-posed linear system.
d(isSingular) = 1;
pv(isSingular) = inf;

% Inflow flux matrix
A  = sparse(n(:,2), n(:,1),  in, nc, nc)...
   + sparse(n(:,1), n(:,2), -out, nc, nc);
A = -A + spdiags(d, 0, nc, nc);

% Subtract the divergence of the velocity minus any source terms from the
% diagonal to account for compressibility effects. Inflow/outflow from
% boundary conditions are accounted for in the divergence, and hence we
% only need to subtract q (and not qb).
div = accumarray(gridCellNo(G), faceFlux2cellFlux(G, state.flux));
A   = A - spdiags(div-q, 0, nc, nc);

% Build right-hand sides for tracer equations. Since we have doubled the
% rate in any cells with a positive source, we need to also double the rate
% on the right-hand side here.
numTrRHS = numel(tr);
TrRHS = zeros(nc,numTrRHS);
for i=1:numTrRHS,
   TrRHS(tr{i},i) = 2*qp(tr{i});
end

% Time of flight for a divergence-free velocity field.
T  = A \ [pv TrRHS];

if ~opt.allowInf
    % User does not want inf values in output (could be problematic for
    % some automatic post processing). Set those values to the maximum
    % non-infinite value found in the remaining cells.
    bad = isinf(T(:, 1));
    T(bad, 1) = max(T(~bad, 1));
end

end


function [q, qb] = computeSourceTerm(state, G, W, src, bc)
   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source term (in m^3/s).

   % Contribution from wells
   if ~isempty(W),
      qi = [qi; vertcat(W.cells)];
      qs = [qs; vertcat(state.wellSol.flux)];
   end

   % Contribution from sources
   if ~isempty(src),
      qi = [qi; src.cell];
      qs = [qs; src.rate];
   end

   % Assemble all source and sink contributions to each affected cell.
   q = sparse(qi, 1, qs, G.cells.num, 1);

   % Contribution from boundary conditions
   if ~isempty(bc),
      ff    = zeros(G.faces.num, 1);

      isDir = strcmp('pressure', bc.type);
      i     = bc.face(isDir);
      if ~isempty(i)
         ff(i) = state.flux(i) .* (2*(G.faces.neighbors(i,1)==0) - 1);
      end

      isNeu = strcmp('flux', bc.type);
      ff(bc.face(isNeu)) = bc.value(isNeu);

      is_outer = ~all(double(G.faces.neighbors) > 0, 2);
      qb = sparse(sum(G.faces.neighbors(is_outer,:), 2), 1, ...
         ff(is_outer), G.cells.num, 1);
   else
      qb = sparse(G.cells.num,1);
   end
end
