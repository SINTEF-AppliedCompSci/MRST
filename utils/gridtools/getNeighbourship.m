function [N, varargout] = getNeighbourship(G, varargin)
%Retrieve neighbourship relation ("graph") from grid
%
% SYNOPSIS:
%    N         = getNeighbourship(G)
%    N         = getNeighbourship(G, kind)
%    N         = getNeighbourship(G, kind, incBdry)
%   [N, isnnc] = getNeighbourship(...)
%
% PARAMETERS:
%   G       - MRST grid as defined by `grid_structure`.
%
%   kind    - What kind of neighbourship relation to extract.  String.  The
%             following options are supported:
%               - 'Geometrical':
%                 Extract geometrical neighbourship relations.  The
%                 geometric connections correspond to physical,
%                 geometric interfaces and are the ones listed in
%                 `G.faces.neighbors`.
%
%               - 'Topological':
%                 Extract topological neighbourship relations.  In
%                 addition to the geometrical relations of `Geometrical`
%                 these possibly include non-neighbouring connections
%                 resulting from pinch-out processing or explicit NNC
%                 lists in an ECLIPSE input deck.
%
%                    Additional connections will only be defined if the
%                    grid `G` contains an `nnc` sub-structure.
%
%   incBdry - Flag to indicate whether or not to include boundary
%             connections.  A boundary connection is a connection in which
%             one of the connecting cells is the outside (i.e., cell zero).
%             LOGICAL scalar.  Default value: incBdry = FALSE (do NOT
%             include boundary connections).
%
% RETURNS:
%   N     - Neighbourship relation.  An m-by-2 array of cell indices that
%           form the connections, geometrical or otherwise.  This array has
%           similar interpretation to the field `G.faces.neighbors`, but
%           may contain additional connections if kind='Topological'.
%
%   isnnc - An m-by-1 LOGICAL array indicating whether or not the
%           corresponding connection (row) of N is a geometrical connection
%           (i.e., a geometric interface from the grid `G`).
%
%           Specifically, isnnc(i) is TRUE if N(i,:) comes from a
%           non-neighbouring (i.e., non-geometrical) connection.
%
% NOTE:
%   If the neighbourship relation is later to be used to compute the graph
%   adjacency matrix using function `getConnectivityMatrix`, then `incBdry`
%   must be `false`.
%
% SEE ALSO:
%   `processGRDECL`, `processPINCH`, `getConnectivityMatrix`.

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


   get_topo = false;  % Geometric neighbourship by default.
   if nargin > 1 && ~isempty(varargin{1}) && ischar(varargin{1}),
      get_topo = strncmpi(varargin{1}, 'Topological', 1);
   end

   incBdry = false;
   if nargin > 2 && (numel(varargin{2}) == 1),
      incBdry = logical(varargin{2});
   end

   % Geometric neighbourship (default)
   N = G.faces.neighbors;

   if ~incBdry,
      % Exclude boundary connections.
      N = N(all(N ~= 0, 2), :);
   end

   if nargout > 1,
      % Caller requested NNC flag.
      % Geometric neighbourships aren't NNCs in MRST.
      varargout{1} = false([size(N, 1), 1]);
   end

   if get_topo && isfield(G, 'nnc') && isstruct(G.nnc) && ...
         isfield(G.nnc, 'cells') && isnumeric(G.nnc.cells) && ...
         size(G.nnc.cells, 2) == 2,

      nnc = [ G.nnc ];
      c   = vertcat(nnc.cells);

      if ~strcmp(class(N), class(c)),
         N = cast(N, class(c));
      end

      N = [ N ; c ];

      if nargout > 1,
         varargout{1} = [ varargout{1} ; true([size(c, 1), 1]) ];
      end
   end
end
