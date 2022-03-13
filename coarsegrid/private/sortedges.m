function [map, sign, success] = sortedges(edges, pos)
% SORTEDGES(E, POS) sorts the edges given by rows in E.
%
% SYNOPSIS:
%   [map, sign, success] = sortedges(e, pos)
%
% PARAMETERS:
%   edges   - An n-by-2 array of edges given by pairs of node numbers; one edge
%             per row.
%
%   pos     - Row positions in e, marking the beginning of each section in e to
%             be sorted separately.
%
%
% RETURNS:
%   map    - Row permutation of edges such that edge edge(j(i),:) is
%            followed by edge edges(j(i+1),:).
%
%   sign   - Sign of each edge, i.e., if E=edge(j,:), then
%            "E(s<0,:)=E(s<0,[2,1]);" renders E fully sorted such that
%            E(k, 2) == E(k+1, 1) within each section.
%
%   success - #sections-by-1 logical array indicating which sections have
%             been sorted successfully.
%
% NOTE:
%    1) This is the Pure Matlab Version.
%
%    2) The sequence of edges in each section must form a one-dimensional
%       connencted graph, that may form a closed loop or not.   Branches
%       are considered an input error.
%
%       If branches, i.e., nodes with more than two edges or disconnected
%       edge sets are encountered, a warning message is issued. Such edges
%       do not have a well-defined order.
%
% EXAMPLE:
%
%   p     = [1;4;6];
%   e     = [2,3; 2,1; 3,1; 1,3; 2,3];
%   [m,s] = sortedges(e, p);
%   E     = e(sub2ind(size(e), [m,m], [1-min(s,0), 1+max(0,s)]));
%
% SEE ALSO:
%

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

% Written by Jostein Natvig, SINTEF Applied Mathematics.

   if nargin == 1,
      pos = [1; numel(edges)/2 + 1];
   end
   edges = double(edges);
     n   = diff(pos);
     N   = max(n);
     act = true(numel(n), 1);


     % Allocate space for 2n-1 indices and signs, for each section.
     m    = 2*n-1;
     mpos = cumsum([1;m]);

     map  = -ones (mpos(end), 1, 'int32');
     sign =  zeros(mpos(end), 1, 'int32');

     % p is the current-write-position. Set current-write-position to the
     % middle of each section of size 2n-1 s.t. it is possible to shift p n
     % positions forward or n positions backwards.
     p    = mpos(1:end-1)+n;

     % Write first edge index and sign.
     map (p) = pos(1:end-1);
     sign(p) = 1;

     % Store next node number to look for in forward and backward directions
     next       = edges(pos(1:end-1), 2);
     prev       = edges(pos(1:end-1), 1);

     % Wipe edge
     edges(pos(1:end-1),:) = nan;

     % Increment forward and backward write positions
     pforward  = p + 1;
     pbackward = p - 1;

     for i = 2:N,
        for j = 1:N-1,
           % Update query edges
           q = min(pos(1:end-1) + j, pos(2:end)-1);

           % Does edge q have any node equal to 'next'?
           ix = any(next(:, [1,1]) == edges(q,:), 2) & act;
           if any(ix),
              % Write
              map (pforward(ix)) = q(ix);
              sign(pforward(ix)) = 2*(next(ix) == edges(q(ix), 1)) - 1;

              % Next node to look for
              next(ix) = (next(ix) == edges(q(ix), 1)).*edges(q(ix), 2) + ...
                         (next(ix) == edges(q(ix), 2)).*edges(q(ix), 1) ;

              % Update write position
              pforward(ix) = pforward(ix) + 1;
           end
           % Delete the edges
           edges(q(ix),:)=nan;


           % Does edge q have any node equal to 'prev'?
           ix = any(prev(:, [1,1]) == edges(q,:), 2) & act;
           if any(ix),
              % Write
              map (pbackward(ix)) = q(ix);
              sign(pbackward(ix)) = 2*(prev(ix) == edges(q(ix), 2)) - 1;

              % Next node to look for
              prev(ix) = (prev(ix) == edges(q(ix), 1)).*edges(q(ix), 2) + ...
                         (prev(ix) == edges(q(ix), 2)).*edges(q(ix), 1) ;
              % Update write position
              pbackward(ix) = pbackward(ix) - 1;
           end
           % Delete the edges
           assert(~any(any(isnan(edges(q(ix),:)))))
           edges(q(ix),:)=nan;

           % Activate sections:
           %   section where not all edges have been assiged a position yet.
           act = pforward-pbackward-2 < n;
        end
     end

     edgeno  = rldecode(1:numel(pos)-1, diff(pos), 2)';
     success = accumarray(edgeno, all(isnan(edges), 2), ...
                          [numel(pos)-1, 1], @all);

     map  = map (map  >  0);
     sign = sign(sign ~= 0);
end
