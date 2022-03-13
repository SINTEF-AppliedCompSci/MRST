function [dof, tol] = ensureMinimumFractionDG(disc, dof, state, tol)
% Set a minimum value on a composition matrix

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

    if nargin == 3
        tol = 1e-8;
    end
    
    [~, xc, cc] = disc.getCubature((1:disc.G.cells.num)', 'volume');
    [~, xf, cf] = disc.getCubature((1:disc.G.cells.num)', 'surface');
    
    nDof = disc.basis.nDof;
    for cNo = 1:size(dof,2)
        
        zc = disc.evaluateDGVariable(xc, cc, state, dof);
        [~, nc] = rlencode(cc);
        [zMinc, zMaxc] = getMinMax(zc, nc);
        
        zf = disc.evaluateDGVariable(xf, cf, state, dof);
        [~, nf] = rlencode(cf);
        [zMinf, zMaxf] = getMinMax(zf, nf);
        
        zMin = min(zMinc, zMinf);
        
        if any(zMin < tol)
            z = disc.getCellMean(dof, state);
            theta = [(z - tol)./(z - zMin), ones(disc.G.cells.num,1)];
            theta(~isfinite(theta) | abs(theta) > 1) = 1;
            %             theta(abs(theta) > 1) = 1;
            theta = min(theta, [], 2);

            for dofNo = 2:nDof
                ix = disc.getDofIx(state, dofNo, (1:G.cells.num)', true);
                dof(ix(ix>0),:) = dof(ix(ix>0),:).*theta(ix>0);
            end

            ind = state.degree > 0;
            dofIx = disc.getDofIx(state, 1, ind);
            dof(dofIx,:) = (dof(dofIx,:) - z(ind,:)).*theta(ind) + z(ind,:);
        end
        
%         z = disc.getCellMean(dof(:,cNo), state);
%         ix = disc.getDofIx(state, 1, z < tol);
%         dof(ix,cNo) = tol;
%         ix = disc.getDofIx(state, 2:nDof, z < tol);
%         dof(ix,cNo) = 0;
    end

%     dof = bsxfun(@rdivide, dof, sum(dof, 2));
end
