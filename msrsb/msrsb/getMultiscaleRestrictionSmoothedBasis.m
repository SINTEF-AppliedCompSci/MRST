function [I, interactionMap, res] = getMultiscaleRestrictionSmoothedBasis(A, CG, varargin)
%Set up basis function for the MsRSB (multiscale restricted smoothed basis) method
%
% SYNOPSIS:
%   I = getMultiscaleRestrictionSmoothedBasis(A, CG)
%
% DESCRIPTION:
%   This function constructs basis functions iteratively by Jacobi-like
%   iterations to a initial interpolator. The method is formulated as
%   matrix-matrix products to get reasonable speed in Matlab.
%
% REQUIRED PARAMETERS:
%   A   - System matrix.
%   CG  - Coarse grid for which to construct basis functions. Should have
%   interaction regions stored in CG.cells.interaction beforehand.
%
%
% RETURNS:
%   I   - Basis function interpolator
%
% SEE ALSO:
%   `getMultiscaleBasis`

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
    opt = struct('enforceCenter', false, ...
                 'iterations',    round(10*CG.parent.cells.num/CG.cells.num), ...
                 'incrementTol',  1e-3, ...
                 'autostop',      true, ...
                 'verbose',       mrstVerbose(), ...
                 'useConstant',   true, ...
                 'omega',         2/3, ...
                 'interpolator',  [], ...
                 'limitSupport',  true ...
                 );
    opt = merge_options(opt, varargin{:});
    centers = mapCenters(CG);
    CG.cells.centers = centers;
    A = A - diag(sum(A, 2));
    
    if ~isfield(CG.cells, 'interaction')
        CG = storeInteractionRegion(CG);
    end
    
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);

    G = CG.parent;
    
    isNode = false(G.cells.num, 1);
    if isempty(opt.interpolator)
        if opt.useConstant
            I_0 = controlVolumeRestriction(CG.partition)';
        else
            I_0 = linearInterpolator(CG);
        end
    else
        I_0 = opt.interpolator;
    end
    
    if opt.enforceCenter
        isNode(centers) = true;

        rhs = -A(~isNode, isNode)*speye(numel(centers));
        A = A(~isNode, ~isNode);
        interactionMap = interactionMap(~isNode, :);
    else
        rhs = 0*I_0;
    end
    I = I_0(~isNode, :);
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    
    res = nan(min(opt.iterations, 1000), 1);
    w = opt.omega;
    i = 1;
    while i < opt.iterations
        def = D_inv*(A*I - rhs);
        if opt.limitSupport
             def = removeOutOfBounds(def, interactionMap);
        end
        nval = full(sum(abs(def), 1));
        if size(res, 1) < i
            res = [res; nan*res];
        end
        
        res(i) = norm(nval, 2);
        if i > 1 && (opt.autostop && res(i) > res(i-1) || abs(res(i)-res(i-1)) < opt.incrementTol)
            dispif(opt.verbose, 'Finished after %d iterations\n', i);
            break
        end
        dispif(opt.verbose, '%d of %d iterations in basis\n', i, opt.iterations)
        
        I = I - w*def;
        I = bsxfun(@rdivide, I, sum(I, 2));
        i = i + 1;
    end
    I_0(~isNode, :) = I;
    I = I_0;
end

function increment = removeOutOfBounds(increment, interactionMap)
    increment = increment.*interactionMap;
end
