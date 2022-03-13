function [conn, ind, cpos, fconn] = coarseConnections(G, p, varargin)
%Derive connection/topology structure on coarse grid
%
% SYNOPSIS:
%    conn                    = coarseConnections(G, p)
%    conn                    = coarseConnections(G, p, Ic)
%   [conn, ind, cpos, fconn] = coarseConnections(...)
%
% PARAMETERS:
%   G  - Fine-scale grid structure over which to derive coarse topology.
%
%   p  - Partition vector for fine-scale entities (cells).
%
%   Ic - Discrete indicator function.  An m-by-n numeric array of function
%        values.  There must be one row in 'Ic' for each connection (row)
%        in the neighbourship definition. In other words, if we set
%        N=G.faces.neighbors, then m==SIZE(N,1). 
%        OPTIONAL.  Default value: Ic = ONES([m, 1]).
%
% RETURNS:
%   conn - 
%        Coarse-scale connection structure represented as a neighbourship
%        definition (m-by-2 array of coarse-scale block pairs).  If we set
%        q = [0; p], then the coarse-scale connection is defined as the
%        intersection of the unique block pairs defined by q(N+1) with the
%        corresponding unique applicable indicator tuples defined by Ic.
%
%   ind - 
%        Actual indicator value (including half-face tag if present)
%        associated with each coarse-scale connection.
%
%   cpos, fconn -
%        Packed data-array representation of coarse->fine connection
%        mapping.  Specifically, the elements
%
%              fconn(cpos(i) : cpos(i + 1) - 1)
%
%        are the fine-scale connections (i.e., rows of the neighbourship
%        definition N) that constitute coarse-scale connection i.
%
% NOTE:
%   This function uses SORTROWS.
%
% SEE ALSO:
%   `generateCoarseGrid`, `rlencode`, `sortrows`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


   N = G.faces.neighbors;
   if (nargin > 2) && isnumeric(varargin{1}),
      Ic = varargin{1};
   else
      Ic = ones([size(N,1), 1]);
   end

   assert (size(Ic,1) == size(N,1), ...
           'Connection indicators must be defined for each connection');

   % Enumerate connections.
   p = [0; p];
   c = (1 : size(N,1)) .';

   pN = p(N + 1);
   i  = pN(:,1) ~= pN(:,2);

   pN = sort(double(pN(i,:)), 2);
   c  = c(i);
   Ic = double(Ic(i,:));

   % Sort and save connections
   A      = sortrows([pN, Ic, c]);
   [B, n] = rlencode(A(:, 1 : end-1));
   conn  = B(:,1:2);
   ind   = B(:,3:end);
   cpos  = cumsum([1; reshape(n, [], 1)]);
   fconn = A(:, end);

end
