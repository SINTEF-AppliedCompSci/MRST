function W = convertwellsVE(W, G, Gt, rock2D, varargin)
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
%
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if nargin > 4
        inner_product = varargin{1};
    else
       % @@ Should the default rather be 'ip_tpf' here, since this is the new
       % default for 'addWell'?
        inner_product = 'ip_simple'; 
    end
    [i, j, k] = ind2sub(G.cartDims, G.cells.indexMap); %#ok<NASGU>
    [i2, j2] = ind2sub(Gt.cartDims, Gt.cells.indexMap);
    W_2D = [];
    for nw = 1:numel(W)
        w = W(nw);
        % Find i/j positions of well in 3D grid
        wi = i(w.cells); wj = j(w.cells);
        assert(numel(unique(wi)) == 1 || numel(unique(wj)) == 1)
        ind2d = find(i2 == wi(1) & j2 == wj(1));
        W_2D = addWell(W_2D, Gt, rock2D, ind2d, ...
                       'Type'         , w.type        , ...
                       'Val'          , w.val         , ...
                       'Radius'       , w.r(1)           , ...
                       'Comp_i'       , w.compi       , ...
                       'InnerProduct' , inner_product , ...
                       'name'         , w.name        , ...
                       'refDepth'     , w.refDepth    , ...
                       'Sign'         , w.sign);           
        %'WI'           , sum(w.WI)     , ... 
    end
    
    for nw = 1:numel(W)
       % Set VE specific parameters
       W_2D(nw).WI= W_2D(nw).WI.*Gt.cells.H(W_2D(nw).cells);              %#ok
       W_2D(nw).h = Gt.cells.H(W_2D(nw).cells);                           %#ok
       W_2D(nw).dZ = Gt.cells.H(W_2D(nw).cells)*0.0;                      %#ok
       if( isfield(W(nw),'bhpLimit') )
          W_2D(nw).bhpLimit = W(nw).bhpLimit; %#ok<AGROW>
       end
    end
    W = W_2D;
end
