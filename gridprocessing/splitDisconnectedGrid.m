function G = splitDisconnectedGrid(G, varargin)
%Split grid into disconnected components
%
% SYNOPSIS:
%   G = splitDisconnectedGrid(G)
%   G = splitDisconnectedGrid(G, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G - Grid structure.
%
% KEYWORD ARGUMENTS:
%
%   Verbose - Whether or not to display progress information
%             Logical.  Default value: `Verbose = mrstVerbose()`.
%
% RETURNS:
%   G - Array of grid structures, one element for each disconnected grid
%       component, sorted in order of decreasing number of cells.  If the
%       grid consists of a single component (i.e., if there are no
%       disconnected grid components), then the return value `G` is equal
%       to the input parameter `G`.
%
%       Function `splitDisconnectedGrid` will also consider any explicit
%       non-neighbouring connections represented in a field `nnc` stored in
%       the top-level structure of `G` when determining whether or not any
%       disconnected components exist.
%
% SEE ALSO:
%   `processGRDECL`, `processPINCH`, `extractSubgrid`

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


   opt = struct('verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   verbose = opt.verbose;

   % Check if grid is connected
   [a, c, c] = dmperm(adjacency(G));                                   %#ok

   ncomp = numel(c) - 1;
   if ncomp > 1,
      dispif(verbose, '\nGrid has %d disconnected components\n', ncomp);

      % Partition grid into connected subgrids
      g = arrayfun(@(i) extractSubgrid(G, a(c(i) : c(i + 1) - 1)), ...
                   1 : ncomp, 'UniformOutput', false);

      % Return grids in order of decreasing number of cells.
      cartDims = G.cartDims;
      [i, i]   = sort(- cellfun(@(g) g.cells.num, g));                 %#ok
      G        = [ g{i} ];

      [ G.cartDims ] = deal(cartDims);
   end
end

%--------------------------------------------------------------------------

function A = adjacency(G)
   N = double(G.faces.neighbors(~ any(G.faces.neighbors == 0, 2), :));

   if isfield(G, 'nnc'),
      N = [ N ; G.nnc.cells ];
   end

   I = [ N(:,1) ; N(:,2) ; (1 : G.cells.num) .' ];
   J = [ N(:,2) ; N(:,1) ; (1 : G.cells.num) .' ];
   A = sparse(I, J, 1, G.cells.num, G.cells.num);
end
