function [vals, ixs] = tsparsemul(avals, aind, bvals, bind)
% will multiply/contract in last index of aind and bind.
% Prerequisite: aind and bind are sorted wrt. last index   
      
   % remove entries that do not contribute
   common_ix = intersect(aind(:,end), bind(:,end));
   num_unique_ix = numel(common_ix);
   
   a_keep_ix = ismember(aind(:, end), common_ix);
   b_keep_ix = ismember(bind(:, end), common_ix);
   
   avals = avals(a_keep_ix); aind = aind(a_keep_ix, :);
   bvals = bvals(b_keep_ix); bind = bind(b_keep_ix, :);
   
   % identify location of each relevant index
   astart = [1; find(diff(aind(:,end)))+1];
   aend   = [astart(2:end)-1; numel(aind(:,end))];

   bstart = [1; find(diff(bind(:,end)))+1];
   bend   = [bstart(2:end)-1; numel(bind(:, end))];

   % split up matrices
   asplit = arrayfun(@(x) aind(astart(x):aend(x), :), 1:num_unique_ix, ...
                     'uniformoutput', false);
   bsplit = arrayfun(@(x) bind(bstart(x):bend(x), :), 1:num_unique_ix, ...
                     'uniformoutput', false);
   
   avalsplit = arrayfun(@(x) avals(astart(x):aend(x), :), 1:num_unique_ix, ...
                     'uniformoutput', false);
   bvalsplit = arrayfun(@(x) bvals(bstart(x):bend(x), :), 1:num_unique_ix, ...
                     'uniformoutput', false);
   
   % function producing all permutation of rows in matrices u and v
   loc_kron = @(u, v) [repmat(u, size(v, 1), 1), repelem(v, size(u, 1), 1)];
   
   % make array with all indices that are to be multiplied together
   ixarr = arrayfun(@(ix) loc_kron(asplit{ix}, bsplit{ix}), 1:num_unique_ix, ...
               'uniformoutput', false);
   ixarr = vertcat(ixarr{:});

   % produce corresponding permutations also for the value vectors (which may
   % be ADI, so we must be a bit more careful)
   
   avalarr = arrayfun(@(ix) repmat(avalsplit{ix}, numelem(bvalsplit{ix}), 1), ...
                      1:num_unique_ix, 'uniformoutput', false);
   avalarr = vertcat(avalarr{:});
   
   bvalarr = ...
       arrayfun(@(ix) bvalsplit{ix}(reshape(repmat(1:numelem(bvalsplit{ix}),...
                                                   numelem(avalsplit{ix}), 1), ...
                                            [], 1)), ...
                1:num_unique_ix, 'uniformoutput', false);
   bvalarr = vertcat(bvalarr{:});
   
   vals = avalarr .* bvalarr;
   ixs = ixarr;
   ixs(:, end) = [];
   ixs(:, size(aind,2)) = [];
   
   if isempty(ixs)
      % the result is an intrinsic scalar
      vals = sum(vals);
      ixs = [];
      return
   end
   
   % here we must do the assembly
   tmp = mat2cell(ixs, size(ixs, 1), ones(size(ixs,2), 1));
   ix_expanded = sub2ind(max(ixs), tmp{:});
   
   map = accumarray([ix_expanded(:), (1:numel(ix_expanded))'], 1, [], [], [], ...
                    true);
   kk = find(sum(map, 2));
   vals = map * vals;
   vals = vals(kk);
   [tmp{:}] = ind2sub(max(ixs), kk);
   ixs = [tmp{:}];
   
end


function ne = numelem(x)
   if isa(x, 'ADI')
      ne = numel(x.val);
   else
      ne = numel(x);
   end
end
