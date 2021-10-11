function q = computeTransportSourceTerm(state, G, wells, src, bc, varargin)
%Compute source term contributions for transport
%
% SYNOPSIS:
%   q = computeTransportSourceTerm(state, G, wells, src, bc)
%   q = computeTransportSourceTerm(state, G, wells, src, bc, 'pn1', pv1, ...)
%
% PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from function 'initState', or the results from a
%           call to function 'incompTPFA'.
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
%   'pn'/pv -
%           List of 'key'/value pairs defining optional parameters.  The
%           supported options are:
%             - use_compi --
%                 Whether or not to include the composition of injected
%                 fluids into the source term.  This is appropriate for
%                 two-phase flow only, and simply scales the rate of fluid
%                 sources with component one of 'q.compi'.  Sinks, i.e,
%                 those terms for which the source rate is negative, remain
%                 unaffected.
%
%                 LOGICAL.  Default value: use_compi=TRUE (include
%                 composition of injected fluid into fluid sources).
%
% RETURNS:
%   q - Source term structure containing details of each individual,
%       derived fluid source and sink term in the model--be they from
%       wells, boundary conditions or synthetic fluid sources (i.e., the
%       'src' parameter).
%
%       Structure containing the following fields
%         - cell --
%             m-by-1 (numeric) array of cell indices, each specifying which
%             cell is affected by a particular contribution.
%
%         - flux --
%             m-by-1 numeric array of source/sink rates.  In particular,
%             'flux(i)' specifies the rate of the 'i'th contributions
%             (affecting 'cell(i)').
%
%         - compi --
%             Empty or m-by-np numeric array ('np' is number of fluid
%             phases) specifying the compositions of the injected fluids.
%
%             This field is empty if use_compi==FALSE.
%
%             The underlying information is derived (extracted) from
%             appropriate structure fields ('.comp_i' for wells, '.sat' for
%             'src' and 'bc).  Note that it is an error to specify
%
%                use_compi=TRUE
%
%             unless all pertinent source term objects (non-empty 'wells',
%             non-empty 'src', or non-empty 'bc') carry sufficient
%             information to specify the composition of the injected fluid.
%
%             The value of compi(p,:) is undefined (typically all zero) for
%             production terms (i.e., fluid sinks) characterised by
%
%                flux(p) < 0
%
% NOTE:
%   - This function is typically appropriate only for single- and two-phase
%     flow and transport problems.
%
%   - MRST does not support specifying both injection and production (i.e.,
%     fluid sources and fluid sinks) in the same cell.  Enforcing that
%     restriction is deferred to function assembleTransportSource.
%
% SEE ALSO:
%   `computeTimeOfFlight`, `explicitTransport`, `private/assembleTransportSource`

% TODO:
%   - implement gravity effects for pressure boundary and wells

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

   opt = struct('use_compi', true);
   opt = merge_options(opt, varargin{:});

   check_input(wells, bc, src, opt)

   c = [];  % Cells to which sources are connected
   r = [];  % Actual strength of source term (in m^3/s).
   i = [];  % Injection composition (size NUMEL(qs)-by-SIZE(sat,2))

   if ~isempty(wells)
      [c, r, i] = ...
         concat(c, r, i, opt, @() contrib_wells(wells, state.wellSol));
   end

   if ~isempty(src)
      [c, r, i] = ...
         concat(c, r, i, opt, @() contrib_src(src));
   end

   if ~isempty(bc), assert (~isempty(bc.sat))
      [c, r, i] = ...
         concat(c, r, i, opt, @() contrib_bc(G, state, bc));
   end

   %-----------------------------------------------------------------------
   % Assemble all source and sink contributions to each affected cell. ----
   %
   q = struct('cell', c, 'flux', r, 'compi', i);
end

%--------------------------------------------------------------------------

function check_input(wells, bc, src, opt)
   if opt.use_compi

      if ~isempty(wells) && any(cellfun('isempty', { wells.compi }))
         error('ComputeTransportSourceTerm:EmptyWellCompi', ...
              ['Wells must have valid ''.comp_i'' fields to ', ...
               'specify (''use_compi'',TRUE).']);
      end

      if ~isempty(bc) && isempty(bc.sat)
         error('ComputeTransportSourceTerm:EmptyBCSat', ...
              ['Boundary conditions must have valid ''.sat'' ', ...
               'field to specify (''use_compi'',TRUE).']);
      end

      if ~isempty(src) && isempty(src.sat)
         error('ComputeTransportSourceTerm:EmptySrcSat', ...
              ['Explicit source terms must have valid ''.sat'' ', ...
               'field to specify (''use_compi'',TRUE).']);
      end

   end
end

%--------------------------------------------------------------------------

function [qi, qs, compi] = concat(qi, qs, compi, opt, contrib)
   [i, s, c] = contrib();

   qi = [ qi ; i ];
   qs = [ qs ; s ];

   if opt.use_compi
      assert (isempty(compi) || size(c, 2) == size(compi, 2), ...
             ['Injection composition incompatible with existing ', ...
              'source terms']);

      compi = [ compi ; c ];
   end
end

%--------------------------------------------------------------------------

function [i, s, c] = contrib_wells(W, wellSol)
   % Wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';

   i = vertcat(W.cells);
   s = vertcat(wellSol.flux);
   c = rldecode(vertcat(W.compi), nperf);
end

%--------------------------------------------------------------------------

function [i, s, c] = contrib_src(src)
   % Explicit sources defined by (e.g.) 'addSource'.

   i = src.cell;
   s = src.rate;
   c = src.sat;
end

%--------------------------------------------------------------------------

function [qi, qs, c] = contrib_bc(G, state, bc)
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

   if ~isempty(bc.sat)
      % Reorder composition to match cells ('qi').

      i = [ reshape(find(isDir), [], 1) ; ...
            reshape(find(isNeu), [], 1) ];

      bc.sat = bc.sat(i, :);
   end

   % Empty if ISEMPTY(bc.sat)
   c = bc.sat;
end
