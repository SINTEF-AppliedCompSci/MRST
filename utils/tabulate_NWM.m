function tab = tabulate_NWM(u)
% Use 'accumarray' with val = 1 to count the number of identical subscripts 
% in u. Equivalent to the matlab function tabulate(u).
% Written by Knut-Andreas Lie
    v = accumarray(u, 1);
    tab = [(1:numel(v))', v];
end