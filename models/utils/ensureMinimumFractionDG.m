function [dof, tol] = ensureMinimumFractionDG(disc, dof, state, tol)
% Set a minimum value on a composition matrix
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