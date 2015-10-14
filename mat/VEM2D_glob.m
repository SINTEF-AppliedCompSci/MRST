function [S, b] = VEM2D_glob(G, f, g)

    Nc = G.cells.num;
    Ne = G.faces.num;
    Nn = G.nodes.num;
    Ndof = Nn + Ne + Nc;
    
    S = sparse(Ndof, Ndof);
    b = zeros(Ndof, 1);
    
    for K = 1:Nc;
        
        [Sl, bl, dofVec] = VEM2D_loc(G, K, f);
        S(dofVec, dofVec) = S(dofVec, dofVec) + Sl;
        b(dofVec) = b(dofVec) + bl;

    end

end