function mrg = mergeSingleNeighbour(bN, has_src, is_ext)
   assert (~any(any(bN < 1)), ...
          ['Neighbourship definition must include ', ...
           'internal connections only.']);

   conn   = blockConnectivity(bN);
   nneigh = cellfun('prodofsize', conn);

   mrg = (1 : numel(has_src)) .';

   has_one_neigh = nneigh == 1;
   has_force     = has_src | is_ext;

   need_merge = find(has_one_neigh & ~has_force);

   while ~isempty(need_merge),
      next_merge = [];

      for b = reshape(need_merge, 1, []),
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
   if numel(u) == 1,
      next = [next, n];
   end
end
