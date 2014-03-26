function J = lMultDiag(d, J1)
  % J = lMultDiag(d, J1)  
    n = numel(d);
    D = sparse((1:n)', (1:n)', d, n, n);
    J = cell(1, numel(J1));
    for k = 1:numel(J)
        J{k} = D*J1{k};
    end
end