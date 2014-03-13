function q = computeTransportSourceTerm(state, G, wells, src, bc)
%Compute source terms for transport
%
% SYNOPSIS:
%   q = computeTransportSourceTerm(state, G, wells, src, bc)
%
% PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.
%
%   wells - Well structure as defined by function 'addWell'.  May be empty
%           (i.e., W = []) which is interpreted as a model without any
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
% RETURNS:
%   q     - Aggregate source term contributions (per grid cell) suitable
%           for passing to a transport solver.  Measured in units of m^3/s.
%
% SEE ALSO:
%   twophaseJacobian, computeTimeOfFlight.

% TODO:
%   - implement gravity effects for pressure boundary and wells

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


   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source term (in m^3/s).

   if ~isempty(wells),
      [i, s] = contrib_wells(wells, state.wellSol);
      qi = [qi; i];
      qs = [qs; s];
   end

   if ~isempty(src), assert (~isempty(src.sat))
      [i, s] = contrib_src(src);
      qi = [qi; i];
      qs = [qs; s];
   end

   if ~isempty(bc), assert (~isempty(bc.sat))
      [i, s] = contrib_bc(G, state, bc);
      qi = [qi; i];
      qs = [qs; s];
   end

   %-----------------------------------------------------------------------
   % Assemble all source and sink contributions to each affected cell. ----
   %
   q = sparse(qi, 1, qs, G.cells.num, 1);
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_wells(W, wellSol)
   % Wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';

   qi = vertcat(W.cells);
   qs = vertcat(wellSol.flux);

   % Injection perforations have positive flux (into reservoir).
   %
   comp      = rldecode(vertcat(W.compi), nperf);
   inj_p     = qs > 0;
   if any(inj_p)
      qs(inj_p) = qs(inj_p) .* comp(inj_p,1);
   end
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_src(src)
   % Explicit sources defined by (e.g.) 'addSource'.

   qi = src.cell;
   qs = src.rate;

   % Injection sources have positive rate into reservoir.
   %
   in = find(src.rate > 0);
   if ~isempty(in),
      qs(in) = qs(in) .* src.sat(in,1);
   end
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_bc(G, state, bc)
   % Contributions from boundary conditions as defined by 'addBC'.

   N = getNeighbourship(G, 'Geometrical', true);
   N = double(N(bc.face, :));

   assert (all(sum(N == 0, 2) == 1), ...
           'Boundary condition supplied on internal face?');

   bdryCell   = @(i) sum(N(i, :), 2);
   bcFluxSign = @(i) 2*(N(i, 1) == 0) - 1;

   % 1) Dirichlet BCs.  Retrieve rate from state.flux.
   isDir = strcmp('pressure', bc.type);
   qi    = bdryCell  (isDir);
   qs    = bcFluxSign(isDir) .* state.flux(bc.face(isDir));

   % 2) Neumann BCs.  Retrieve rate from BC specification itself.
   isNeu = strcmp('flux', bc.type);
   qi    = [ qi ; bdryCell(isNeu) ];
   qs    = [ qs ; bc.value(isNeu) ];

   % Injection BCs have positive rate (flux) into reservoir.
   %
   is_inj = qs > 0;
   if any(is_inj),
      i = [ reshape(find(isDir), [], 1) ; ...
            reshape(find(isNeu), [], 1) ];

      assert (numel(is_inj) == numel(i), 'Internal error in contrib_bc');

      qs(is_inj) = qs(is_inj) .* bc.sat(i(is_inj), 1);
   end
end
