function [div] = VEM_div(G)
% Discrete divergence operator for the virtual element method
%
% SYNOPSIS:
%   function [div] = VEM_div(G)
%
% DESCRIPTION: Computes a discrete divergence operator. This discrete divergence operator is a
% mapping from node-valued displacement vector to cell-valued
% vector. Node-valued displacement vectors correspond to the degrees of freedom
% that determine for each cell a displacement function over the cell via the
% virtual basis function. The discrete divergence operator that is assembled here
% computes the L^2 projection, cell-wise, of this displacement function. For
% more detail, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1].
%
% PARAMETERS:
%   G - Grid struture
%
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

    if (G.griddim == 3)
        div = VEM3D_div(G);
    else
        assert(G.griddim == 2)
        div = VEM2D_div(G);
    end
end
