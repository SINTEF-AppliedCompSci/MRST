function [iC, kC, C] = multiSparseProduct(iA, jA, A, jB, kB, B)

    [m_iA, ~, b_iA] = unique(iA, 'rows');
    [m_kB, ~, b_kB] = unique(kB, 'rows');

    [m_j, ~, b_j] = unique([jA; jB], 'rows');
    
    b_jA = b_j(1 : size(jA, 1));
    b_jB = b_j(size(jA, 1) + (1 : size(jB, 1)));
    
    sA = sparse(b_iA, b_jA, A);
    sB = sparse(b_jB, b_kB, B);
    sC = sA * sB;

    [b_iC, b_kC, C] = find(sC);

    iC = m_iA(b_iC, :);
    kC = m_kB(b_kC, :);
    
end
