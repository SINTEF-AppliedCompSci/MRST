function [div] = VEM2D_div(G)
%
%
% SYNOPSIS:
%   function [div] = VEM2D_div(G)
%
% DESCRIPTION:
% Computes a discrete div operator in 2D. This discrete div operator is a
% mapping from node-valued displacement vector to cell-valued 2D
% vector. Node-valued displacement vectors correspond to the degrees of freedom
% that determine for each cell a displacement function over the cell via the
% virtual basis function. The discrete div operator computes the L^2
% projection, cell-wise,  of this displacement function. For more detail, see
% ecmor paper.
%
% PARAMETERS:
%   G        - Grid structure
%
% RETURNS:
%   div - matrix corresponding to the discrete div operator.
%
% EXAMPLE:
%
% SEE ALSO:
%

    qf     = calculateQF(G);
    qf     = reshape(qf', [], 1);
    dofs   = mcolon(G.griddim * (G.cells.nodes - 1) + 1, ...
                    G.griddim * (G.cells.nodes - 1) + G.griddim);
    ndofs  = G.nodes.num * G.griddim;
    cellno = rldecode([1:G.cells.num]', diff(G.cells.nodePos) * G.griddim);
    div = sparse(cellno, dofs, qf, G.cells.num, ndofs);

end
