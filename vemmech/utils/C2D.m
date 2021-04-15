function D = C2D(C, G,varargin)
% Compute the matrix of normalized strain energies (D) from the elasticity
% tensor (C) associated with grid (G).
%
% convert Stiffness tensor from Voigt notation to plain notations (no
% factors are applied to the off-diagonal coefficients of the strain
% tensor). The latter is used in our formulation of VEM. The resulting tensor
% D is such that, for S in plain notation, we have: S'DS = energy. 

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

    opt = struct('inv', false);
    opt = merge_options(opt, varargin{:});
    
    if(G.griddim == 3)
        nlin = 6; 
    else
        assert(G.griddim == 2)
        nlin = 3; 
    end

    factor = [ones(G.griddim, G.griddim),            ones(G.griddim, nlin - G.griddim) * 2; ...
              ones(nlin - G.griddim, G.griddim) * 2, ones(nlin - G.griddim, nlin - G.griddim) * ...
              4]; 
    factor = reshape(factor, 1, []);
    if(opt.inv)
        factor = 1./factor;
    end
    D = bsxfun(@times, C, factor);
end

