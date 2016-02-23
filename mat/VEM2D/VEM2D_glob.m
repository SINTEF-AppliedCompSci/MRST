function [S, b, hG] = VEM2D_glob(G, f)

    Nc = G.cells.num;
    Ne = G.faces.num;
    Nn = G.nodes.num;
    Ndof = Nn + Ne + Nc;
    
    S = sparse(Ndof, Ndof);
    b = zeros(Ndof, 1);
    
    hG = 0;
    hK = 0;
    for K = 1:Nc; 
        [Sl, bl, dofVec] = VEM2D_loc(G, K, f);
        S(dofVec, dofVec) = S(dofVec, dofVec) + Sl;
        b(dofVec) = b(dofVec) + bl;
        hG = max(hG,hK);
    end

end