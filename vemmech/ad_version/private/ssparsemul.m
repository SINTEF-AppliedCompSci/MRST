function [rows, cols, vals] = ssparsemul(ixs1, v1, ixs2, v2)
%Undocumented Utility Function

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

   im = ixs1(:,1); ik1 = ixs1(:,2);
   in = ixs2(:,2); ik2 = ixs2(:,1);

   M = max(im);
   N = max(in);

   % [im, ik1, v1] = find(m1);
   % [ik2, in, v2] = find(m2);
   % M = size(m1, 1);
   % N = size(m2, 2);

   tmp = diff(in);
   tmp_ix = find(tmp);
   
   k_for_n = [ones(in(1), 1); 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(in)+1];

   k_for_n = [k_for_n; repmat(k_for_n(end), N+1-numel(k_for_n), 1)];
   
   tmp = diff(ik1);
   tmp_ix = find(tmp);
   m_for_k = [ones(ik1(1), 1); 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(ik1)+1];

   m_for_k = [m_for_k; repmat(m_for_k(end), max(ik2)+1 - numel(m_for_k), 1)];
   
   if isa(v1, 'ADI') || isa(v2, 'ADI')
      % temporal storage arrays must be ADI too
      ALLOC_MULT = 5;
      model = v1;
      if ~isa(model, 'ADI')
         model = v2;
      elseif (numel(in) > numel(im)) && isa(v2, 'ADI')
         model = v2;
      end
      final_m = zeros(ALLOC_MULT * numel(model), 1);
      final_n = zeros(ALLOC_MULT * numel(model), 1);
      %final_v = repmat(model, ALLOC_MULT, 1) * 0; % zero, in ADI sense
      
      tmp_v = model * 0;
      tmp_m = value(tmp_v);
   else
      ALLOC_SIZE = 100 * max(numel(im), numel(in));
      final_m = zeros(ALLOC_SIZE, 1);
      final_n = zeros(ALLOC_SIZE, 1);
      %final_v = zeros(ALLOC_SIZE, 1);
      
      %ALLOC_SIZE_LOC = 1000;
      ALLOC_SIZE_LOC = numel(im); % likely overkill, but ok
      tmp_m = zeros(ALLOC_SIZE_LOC, 1);
      tmp_v = zeros(ALLOC_SIZE_LOC, 1);
   end
      
   count = 1;
   next_entry = 1;

   %for n = 1:N
   n_unique = unique(in(:))';
   accum_v = cell(numel(n_unique), 1);
   accum_counter = 1;
   for n = n_unique
      next_entry_loc = 1;
      for k = ik2(k_for_n(n):k_for_n(n+1)-1)'
         val = v2(count); count = count + 1;
         mix0 = m_for_k(k);
         mix1 = m_for_k(k+1);
         num_new = mix1-mix0;
         if num_new > 0
            tmp_m(next_entry_loc:next_entry_loc + num_new - 1) = im(mix0:mix1-1);
            tmp_v(next_entry_loc:next_entry_loc + num_new - 1) = val * v1(mix0:mix1-1);
         end
         next_entry_loc = next_entry_loc + num_new;
      end
      
      % ensure uniqueness of entries
      tmp_m_loc = tmp_m(1:next_entry_loc-1);

      % @@ the following method seems to work efficiently, and works with ADI
      ltable = find(accumarray(tmp_m_loc, 1, [], [], [], true));
      ltable_inv = sparse(ltable, 1, 1:numel(ltable));
      
      tmp_m_loc = ltable_inv(tmp_m_loc);
      spmat = sparse(tmp_m_loc, 1:numel(tmp_m_loc), 1);
      
      um = ltable(:);
      tot = spmat * tmp_v(1:next_entry_loc-1);
      nu = numel(um);
            
      % check if reallocation is needed (if so, double array size)
      if next_entry + nu > size(final_m, 1)
         final_m = [final_m; 0 * final_m]; %#ok
         final_n = [final_n; 0 * final_n]; %#ok
         %final_v = [final_v; 0 * final_v]; %#ok
      end
      
      final_m(next_entry:next_entry+nu-1) = um;
      final_n(next_entry:next_entry+nu-1) = n;
      %final_v(next_entry:next_entry+nu-1) = tot;
      accum_v{accum_counter} = tot;
      accum_counter = accum_counter + 1;
      
      next_entry = next_entry + nu;
   end
   
   %keyboard;
   rows = final_m(1:next_entry-1);
   cols = final_n(1:next_entry-1);
   %vals = final_v(1:next_entry-1);
   vals = vertcat(accum_v{:});
   
   
   % res = sparse(final_m(1:next_entry-1), ...
   %              final_n(1:next_entry-1), ...
   %              final_v(1:next_entry-1), M, N);
end
