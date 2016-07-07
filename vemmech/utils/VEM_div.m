function [div] = VEM_div(G)
%
%
% SYNOPSIS:
%   function [div] = VEM_div(G)
%
% DESCRIPTION: Computes a discrete div operator. This discrete div operator is a
% mapping from node-valued displacement vector to cell-valued
% vector. Node-valued displacement vectors correspond to the degrees of freedom
% that determine for each cell a displacement function over the cell via the
% virtual basis function. The discrete div operator that is assembled here
% computes the L^2 projection, cell-wise, of this displacement function. For
% more detail, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
%
% PARAMETERS:
%   G - Grid struture
%
% RETURNS:
%   div - matrix corresponding to the discrete div operator.
%
% EXAMPLE:
%
% SEE ALSO:
%

    if (G.griddim == 3)
        div = VEM3D_div(G);
    else
        assert(G.griddim == 2)
        div = VEM2D_div(G);
    end
end
