function res = ssparsemul(m1, m2)

   [im, ik1, v1] = find(m1);
   [ik2, in, v2] = find(m2);
   
   tmp = diff(ik1);
   tmp_ix = find(tmp);
   m_for_k = [1; 1+rldecode(tmp_ix, tmp(tmp_ix)); numel(ik1)+1];

   M = size(m1, 1);
   N = size(m2, 2);
   
   ALLOC_SIZE = max(nnz(m1), nnz(m2));
   next_entry = 1;
   final_m = zeros(ALLOC_SIZE, 1);
   final_n = zeros(ALLOC_SIZE, 1);
   final_v = zeros(ALLOC_SIZE, 1);
   
   for el = 1:numel(in)
   
      val = v2(el);
      n = in(el);
      k = ik2(el);
      
      mix0 = m_for_k(k);
      mix1 = m_for_k(k+1);
      num_new = mix1-mix0;
      
      % check if reallocation is needed (if so, double array size)
      if next_entry + num_new > size(final_m, 1)
         cur_size = size(final_m);
         final_m = [final_m; zeros(cur_size)]; %#ok
         final_n = [final_n; zeros(cur_size)]; %#ok
         final_v = [final_v; zeros(cur_size)]; %#ok
      end

      final_m(next_entry:next_entry+num_new-1) = im(mix0:mix1-1);
      final_n(next_entry:next_entry+num_new-1) = n;
      final_v(next_entry:next_entry+num_new-1) = val * v1(mix0:mix1-1);
      
      next_entry = next_entry + num_new;
      
   end
   
   res = accumarray([final_m(1:next_entry-1),final_n(1:next_entry-1)], ...
                    final_v(1:next_entry-1), [M, N], [], [], true);
end
