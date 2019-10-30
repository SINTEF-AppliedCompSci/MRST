function res = ssparsemul(m1, m2)

   [im, ik1, v1] = find(m1);
   [ik2, in, v2] = find(m2);

   % tmp = diff(in);
   % tmp_ix = find(tmp);
   
   %k_for_n = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(in)+1];
   %k_for_n = [1; 1+repmat(tmp_ix', tmp(tmp_ix), 1); numel(in)+1];
   %k_for_n = [1; 1+find(diff(in)); numel(in)+1];
   
   tmp = diff(ik1);
   tmp_ix = find(tmp);
   m_for_k = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(ik1)+1];
   %m_for_k = [1; 1+repmat(tmp_ix', tmp(tmp_ix), 1); numel(ik1)+1];
   %m_for_k = [1; 1+find(diff(ik1)); numel(ik1)+1];

   M = size(m1, 1);
   N = size(m2, 2);
   
   ALLOC_SIZE = 100000;
   next_entry = 1;
   final_m = zeros(ALLOC_SIZE, 1);
   final_n = zeros(ALLOC_SIZE, 1);
   final_v = zeros(ALLOC_SIZE, 1);
   %entries = zeros(ALLOC_SIZE, 3);
   
   for el = 1:numel(in)
   
      val = v2(el);
      n = in(el);
      k = ik2(el);
      
      %m_ix = m_for_k(k):m_for_k(k+1)-1;
      m = im(m_for_k(k):m_for_k(k+1)-1);
      mval = v1(m_for_k(k):m_for_k(k+1)-1);
      
      % check if reallocation is needed (if so, double array size)
      if next_entry + m > size(final_m, 1)
         cur_size = size(final_m);
         final_m = [final_m; zeros(cur_size)]; %#ok
         final_n = [final_n; zeros(cur_size)]; %#ok
         final_v = [final_v; zeros(cur_size)]; %#ok
         %entries = [entries; zeros(size(entries))]; %#ok
      end

      final_m(next_entry:next_entry+numel(m)-1) = m;
      final_n(next_entry:next_entry+numel(m)-1) = n;
      final_v(next_entry:next_entry+numel(m)-1) = mval * val;
      
      % entries(next_entry:next_entry+numel(m)-1,:) = ...
      %     [m, n * ones(size(m)), mval * val];
      
      next_entry = next_entry + numel(m);
      
   end
   
   % for n = 1:N
   
   %    kset = ik2(k_for_n(n):k_for_n(n+1)-1);
      
   %    for k = kset(:)'

   %       m = im(m_for_k(k):m_for_k(k+1)-1);

   %       % check if reallocation is needed (if so, double array size)
   %       if next_entry + m > size(entries, 1)
   %          entries = [entries; zeros(size(entries))]; %#ok
   %       end
         
   %       entries(next_entry:next_entry+numel(m)-1, :) = ...
   %           [m, n * ones(size(m)), v1(m) * v2(n)];
         
   %       next_entry = next_entry + numel(m);
         
   %    end
      
   % end
   
   %entries = entries(1:next_entry-1, :);
   
   res = accumarray([final_m(1:next_entry-1),final_n(1:next_entry-1)], ...
                    final_v(1:next_entry-1), [M, N], [], [], true);
end
