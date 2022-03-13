function mrg = mergeSingleNeighbour(bN, has_src, is_ext)
%Undocumented Internal Utility Function

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

   assert (~any(any(bN < 1)), ...
          ['Neighbourship definition must include ', ...
           'internal connections only.']);

   conn   = blockConnectivity(bN);
   nneigh = cellfun('prodofsize', conn);

   mrg = (1 : numel(has_src)) .';

   has_one_neigh = nneigh == 1;
   has_force     = has_src | is_ext;

   need_merge = find(has_one_neigh & ~has_force);

   while ~isempty(need_merge)
      next_merge = [];

      for b = reshape(need_merge, 1, [])
         [mrg, next_merge, conn] = merge_block(mrg, next_merge, b, conn);
      end

      need_merge = next_merge(~has_force(next_merge));
   end
end

%--------------------------------------------------------------------------

function [mrg, next, conn] = merge_block(mrg, next, b, conn)
   n = blockNeighbours(conn, b);
   n = unique(mrg(n));

   assert (numel(n) == 1, 'Internal error defining block neighbours.');

   mrg(mrg == b) = n;       % Merge 'b' into 'n'
   conn{b} = [];            % Quash the rebellion
   conn{n} = mrg(conn{n});  % Remove (now) undefined 'b' from n's neighs.

   u = unique(blockNeighbours(conn, n));
   if numel(u) == 1
      next = [next, n];
   end
end
