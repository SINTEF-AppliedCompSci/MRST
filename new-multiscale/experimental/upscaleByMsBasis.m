function Tc = upscaleByMsBasis(CG, T, varargin)
    G = CG.parent;
    A = getIncomp1PhMatrix(G, T);
    basis = getMultiscaleBasis(CG, A, varargin{:});
    
    A_c = basis.R*A*basis.B;
    
    [ic, jc, vc] = find(A_c);
    Tc = ones(CG.faces.num, 1)*sqrt(eps)*max(abs(vc));

    for i = 1:CG.cells.num
        fa = gridCellFaces(CG, i);
        for j = 1:numel(fa)
            f = fa(j);
            
            c = CG.faces.neighbors(f, :);
            c = c(c~=i);
            if c == 0
                continue
            end
            Tc(f) = abs(vc(ic == i & jc == c));
        end
    end
end