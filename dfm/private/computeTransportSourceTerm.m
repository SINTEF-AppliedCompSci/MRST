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
%   `twophaseJacobian`, `computeTimeOfFlight`.

% TODO:
%   - implement gravity effects for pressure boundary and wells

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
      is_int = all(double(G.faces.neighbors) > 0, 2);
      [i, s] = contrib_bc(G, state, bc, is_int);
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
   qs(inj_p) = qs(inj_p) .* comp(inj_p,1);
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

function [qi, qs] = contrib_bc(G, state, bc, is_int)
   % Contributions from boundary conditions as defined by 'addBC'.

   qs = zeros([G.faces.num, 1]);
   dF = false([G.faces.num, 1]);

   isDir = strcmp('pressure', bc.type);
   isNeu = strcmp('flux',     bc.type);

   dF(bc.face(isDir))      = true;
   cfIx                    = dF(G.cells.faces(:,1));

   cflux = faceFlux2cellFlux(G, state.flux);
   qs(G.cells.faces(cfIx,1)) = -cflux(cfIx);
   qs(bc.face(isNeu))      = bc.value(isNeu);

   % Injection BC's have positive rate (flux) into reservoir.
   %
   is_inj = qs > 0;
   if any(is_inj),
      qs(is_inj) = qs(is_inj) .* bc.sat(is_inj(bc.face), 1);
   end

   is_outer = ~is_int;

   qi = sum(G.faces.neighbors(is_outer,:), 2);
   qs = qs(is_outer);
end
