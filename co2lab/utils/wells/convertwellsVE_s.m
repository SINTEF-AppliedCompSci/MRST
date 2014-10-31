function W = convertwellsVE_s(W, G, Gt, rock2D,inner_product)
% Convert wells in 3D grid G to wells suitable for VE simulations in 2D
% 
% SYNOPSIS:
% W = convertwellsVE(W, G, Gt, rock2D)
%
%
% DESCRIPTION:
%
% Converts wells into a suitable format for vertical equilibrium
% simulations on a top grid. Wells must be definable by one and only
% one logical index in ij-plane.

    [i, j, k] = ind2sub(G.cartDims, G.cells.indexMap);          %#ok<NASGU>
    [i2, j2] = ind2sub(Gt.cartDims, Gt.cells.indexMap);
    W_2D = [];
    for nw = 1:numel(W)
        w = W(nw);
        % Find i/j positions of well in 3D grid
        wi = i(w.cells); wj = j(w.cells);
        assert(numel(unique(wi)) == 1 || numel(unique(wj)) == 1)
        ind2d = find(i2 == wi(1) & j2 == wj(1));
        W_2D = addWell(W_2D, Gt, rock2D, ind2d,'Type', w.type, 'Val', w.val, 'Radius', w.r,...
                      'Comp_i', w.compi, 'InnerProduct', inner_product, 'name', w.name,'Sign',w.sign); 
                 
    end
    for nw = 1:numel(W)
       % Set VE specific parameters       
        W_2D(nw).WI= W_2D(nw).WI.*Gt.cells.H(W_2D(nw).cells);              %#ok
        W_2D(nw).dZ = Gt.cells.H(W_2D(nw).cells)*0.0;%#ok
        if( isfield(W(nw),'bhpLimit') )
          W_2D(nw).bhpLimit=W(nw).bhpLimit; %#ok<AGROW>
        end
    end
    W = W_2D;
end
