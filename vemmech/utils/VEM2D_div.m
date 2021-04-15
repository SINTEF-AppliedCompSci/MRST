function [div] = VEM2D_div(G, varargin)
% Discrete divergence operator for the virtual element method in 2D
%
% SYNOPSIS:
%   function [div] = VEM2D_div(G)
%
% DESCRIPTION: Computes a discrete divergence operator in 2D. This discrete div operator is a mapping from node-valued
% displacement vector to cell-valued 2D vector. Node-valued displacement vectors correspond to the degrees of freedom that
% determine for each cell a displacement function over the cell via the virtual basis function. The discrete divergence
% operator that is assembled here computes the divergence of the L^2 projection, cell-wise, of this displacement
% function. For more detail, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
%
% PARAMETERS:
%   G - Grid structure
%
% OPTIONAL PARAMETERS:
%
%   'extraFaceDof' - This option has to be set to get a stable divergence
%                    operator when extra degrees of freedom have been introduced
%                    on the edges to avoid numerical locking.
% RETURNS:
%   div - matrix corresponding to the discrete divergence operator.
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('extraFaceDof', false); 
    opt = merge_options(opt, varargin{:}); 

    qf     = calculateQF(G);
    qf     = reshape(qf', [], 1);
    dofs   = mcolon(G.griddim * (G.cells.nodes - 1) + 1, ...
                    G.griddim * (G.cells.nodes - 1) + G.griddim);
    ndofs  = G.nodes.num * G.griddim;
    cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos) * G.griddim);
    
    if (opt.extraFaceDof)
        ncf        = diff(G.cells.facePos); 
        fcellno    = rldecode([1:G.cells.num]', ncf); 
        faces      = G.cells.faces(:, 1); 
        fasgn      = (2 * (G.faces.neighbors(faces, 2) == fcellno) - 1) .* ...
            G.faces.areas(faces); 
        cellno     = [cellno; fcellno]; dofs = [dofs'; ndofs + faces]; 
        qf         = [qf; fasgn]; 
        ndofs      = ndofs + G.faces.num; 
    end    
    
    div = sparse(cellno, dofs, qf, G.cells.num, ndofs);

end
