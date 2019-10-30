function res = ssparsemul(m1, m2)

   [im, ik1, v1] = find(m1);
   [ik2, in, v2] = find(m2);
   
   tmp = diff(in);
   tmp_ix = find(tmp);
   
   k_for_n = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(in)+1];
   
   tmp = diff(ik1);
   tmp_ix = find(tmp);
   m_for_k = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(ik1)+1];

   M = size(m1, 1);
   N = size(m2, 2);
   
   ALLOC_SIZE = 100 * max(nnz(m1), nnz(m2));
   next_entry = 1;
   final_m = zeros(ALLOC_SIZE, 1);
   final_n = zeros(ALLOC_SIZE, 1);
   final_v = zeros(ALLOC_SIZE, 1);
   
   ALLOC_SIZE_LOC = 1000;
   tmp_m = zeros(ALLOC_SIZE_LOC, 1);
   tmp_v = zeros(ALLOC_SIZE_LOC, 1);
   
   count = 1;
   for n = 1:N
      next_entry_loc = 1;
      for k = ik2(k_for_n(n):k_for_n(n+1)-1)'
         val = v2(count); count = count + 1;
         mix0 = m_for_k(k);
         mix1 = m_for_k(k+1);
         num_new = mix1-mix0;
         
         tmp_m(next_entry_loc:next_entry_loc + num_new - 1) = im(mix0:mix1-1);
         tmp_v(next_entry_loc:next_entry_loc + num_new - 1) = val * v1(mix0:mix1-1);
      
         next_entry_loc = next_entry_loc + num_new;
      end
      
      % ensure uniqueness of entries
      tmp_m_loc = tmp_m(1:next_entry_loc-1);

      % @@ the following method seems to work efficiently, and should be
      % applicable for ADI variables too, since tmp_v is used for nothing but
      % matrix multiplication
      ltable = find(accumarray(tmp_m_loc, 1, [], [], [], true));
      ltable_inv = sparse(ltable, 1, 1:numel(ltable));
      
      tmp_m_loc = ltable_inv(tmp_m_loc);
      
      spmat = sparse(tmp_m_loc, 1:numel(tmp_m_loc), 1);
      
      um = ltable(:);
      tot = spmat * tmp_v(1:next_entry_loc-1);
      nu = numel(um);
      
      
      % @@ method works and is efficient, but sparse might not work on adi
      % tmp_m_loc = tmp_m(1:next_entry_loc-1);
      % [um, ~, tot] = find(sparse(tmp_m_loc, 1:numel(tmp_m_loc), 1) * ... ...
      %                     sparse(tmp_v(1:next_entry_loc-1)));
      % nu = numel(um);
      
      % @@ method works and is efficient, but accumarray might not work for ADI
      % aarr = accumarray(tmp_m(1:next_entry_loc-1), tmp_v(1:next_entry_loc-1), [], [], ...
      %                   [], true);
      % [um, ~, tot] = find(aarr);
      % nu = numel(um);

      
      
      % tmp_m_loc = tmp_m(1:next_entry_loc-1);
      % [um, ~, idx] = unique(tmp_m_loc);
      % nu = numel(um);
      % tot = zeros(nu, 1);
      % for ii = 1:nu
      %    tot(ii) = sum(tmp_v(idx==ii));
      % end
      
      % check if reallocation is needed (if so, double array size)
      if next_entry + nu > size(final_m, 1)
         cur_size = size(final_m);
         final_m = [final_m; zeros(cur_size)]; %#ok
         final_n = [final_n; zeros(cur_size)]; %#ok
         final_v = [final_v; zeros(cur_size)]; %#ok
      end
      
      final_m(next_entry:next_entry+nu-1) = um;
      final_n(next_entry:next_entry+nu-1) = n;
      final_v(next_entry:next_entry+nu-1) = tot;
      next_entry = next_entry + nu;
   end
   
   res = sparse(final_m(1:next_entry-1), ...
                final_n(1:next_entry-1), ...
                final_v(1:next_entry-1));
end

% function res = ssparsemul(m1, m2)

%    [im, ik1, v1] = find(m1);
%    [ik2, in, v2] = find(m2);
   
%    tmp = diff(ik1);
%    tmp_ix = find(tmp);
%    m_for_k = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(ik1)+1];

%    M = size(m1, 1);
%    N = size(m2, 2);
   
%    ALLOC_SIZE = 100 * max(nnz(m1), nnz(m2));
%    next_entry = 1;
%    final_m = zeros(ALLOC_SIZE, 1);
%    final_n = zeros(ALLOC_SIZE, 1);
%    final_v = zeros(ALLOC_SIZE, 1);
   
%    for el = 1:numel(in)
   
%       val = v2(el);
%       n = in(el);
%       k = ik2(el);
      
%       mix0 = m_for_k(k);
%       mix1 = m_for_k(k+1);
%       num_new = mix1-mix0;
      
%       % check if reallocation is needed (if so, double array size)
%       if next_entry + num_new > size(final_m, 1)
%          cur_size = size(final_m);
%          final_m = [final_m; zeros(cur_size)]; %#ok
%          final_n = [final_n; zeros(cur_size)]; %#ok
%          final_v = [final_v; zeros(cur_size)]; %#ok
%       end

%       final_m(next_entry:next_entry+num_new-1) = im(mix0:mix1-1);
%       final_n(next_entry:next_entry+num_new-1) = n;
%       final_v(next_entry:next_entry+num_new-1) = val * v1(mix0:mix1-1);
      
%       next_entry = next_entry + num_new;
      
%    end
   
%    res = accumarray([final_m(1:next_entry-1),final_n(1:next_entry-1)], ...
%                     final_v(1:next_entry-1), [M, N], [], [], true);
% end
