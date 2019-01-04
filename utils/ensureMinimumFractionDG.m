function [dof, tol] = ensureMinimumFractionDG(disc, dof, state, tol)
% Set a minimum value on a composition matrix
    if nargin == 3
        tol = 1e-8;
    end
    nDof = disc.basis.nDof;
    for cNo = 1:size(dof,2)
        z = disc.getCellMean(dof(:,cNo), state);
        ix = disc.getDofIx(state, 1, z < tol);
        dof(ix,cNo) = tol;
        ix = disc.getDofIx(state, 2:nDof, z < tol);
        dof(ix,cNo) = 0;
    end

%     dof = bsxfun(@rdivide, dof, sum(dof, 2));
end