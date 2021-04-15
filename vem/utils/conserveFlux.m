function [state, r] = conserveFlux(state, G, rock, varargin)
%Postprocess nonconservative flux field.
%
% SYNOPSIS:
%   [state, r] = conserveFlux(state, G, rock)
%   [state, r] = conserveFlux(state, G, rock, 'pn1', vn1, ...)
%
% DESCRIPTION:
%   Postprocesses nonconsrevative flux field using a TPFA-like scheme.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure, result from
%            a previous call to f.ex. function 'incompVEM' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%
%   G       - Grid structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.
%
% OPTIONAL PARAMETERS:
%   bc      - Boundary condition structure as defined by function 'addBC'.
%             This structure accounts for all external boundary conditions
%             to the reservoir flow.  May be empty (i.e., bc = []) which is
%             interpreted as all external no-flow (homogeneous Neumann)
%             conditions.
%
%   src     - Explicit source contributions as defined by function
%             'addSource'.  May be empty (i.e., src = []) which is
%             interpreted as a reservoir model without explicit sources.
%
%   faceWeights - The choice of face weights for the L2 norm. String.
%             Default vale = 'permWeighted'.
%             Supported value are:
%               - permWeighted : Each face is weighted by the inverse of
%                 the sum of the inverses of the permeability in its
%                 neighbor cells.
%               - tpf          : Face weights equals the TPFA
%                                tranmissibilities.
%               - ones         : All face weights equals one.
%
%   tol     - Tolerance for the residuals
%             r_i = \int \Omega_i q dx - \int_\Omega_i v \cdot n ds,
%             If the function fails to construc a flux field with residuals
%             norm(r)/norm(rhs) < tol, a warning is displayed. Scalar.
%             Default value = 1e-14.
%
% RETURNS:
%   state - Update reservoir solution structure with locally conservative
%           fluxes.
%   r     - Residuals.
%
% SEE ALSO:
%   incompVEM.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.

%% Main function

%   Merge input parameters
opt = struct('bc'              , []            , ...
             'src'             , []            , ...
             'faceWeights'     , 'permWeighted', ...
             'tol'             , 1e-14              );

opt = merge_options(opt, varargin{:});

%   Identify neumann boundary faces.
neu = false(G.faces.num, 1);
bf = boundaryFaces(G);
neu(bf) = true;
if ~isempty(opt.bc)
    neu(opt.bc.face(strcmp(opt.bc.type, 'pressure'))) = false;
end

%   Face orientations
f = G.cells.faces(:,1);
ncf = diff(G.cells.facePos);
fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));

%   Set rhs.
rhs = zeros(G.cells.num,1);
if ~isempty(opt.src)
    rhs(opt.src.cell) = opt.src.rate;
end

%   Matrix P sums faces per cell.
[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
P = sparse(ii, jj, 1);

%   Calculate residulas.
r = rhs-P*(state.flux(f).*fSgn);

%   If above tolerance, apply postprocessing.
den = rhs;
if nnz(rhs) == 0
    den = state.flux;
end
if norm(r)/norm(den) < opt.tol
    warning('Flux already conservative. No need for postprocessing.');
else

    %   Calculate weights omega
    omega = computeFaceWeights(G, rock, opt);

    %   Build matrix systems as if all half-face normals have
    %   orientation out of the cell.

    %   Off-diagonal element (i,j) equals minus the face area of common
    %   face for cell i and j.
    c = G.faces.neighbors(f,:);
    c(fSgn == -1,:) = c(fSgn == -1, 2:-1:1);
    nz = c(:,2) ~= 0;
    B = sparse(c(nz,1), c(nz,2), -G.faces.areas(f(nz))./omega(nz), ...
                                             G.cells.num, G.cells.num);

    %   Diagonal element (i,i) equals sum of face areas for cell i.
    I = ones(numel(f),1);
    I(neu(f)) = 0;
    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
    ca = sparse(ii,jj,I)*(G.faces.areas(f)./omega);
    B = B + spdiags(ca, 0, G.cells.num, G.cells.num);

    if nnz(neu) == numel(bf)
        B(1,1) = 2*B(1,1);
    end

    beta = B\r;
    beta = rldecode(beta, ncf,1);
    beta(neu(f)) = 0;

    beta = sparse(f, 1:numel(f), I)*(beta./omega.*fSgn);
    flux = state.flux + G.faces.areas.*beta;

    r = rhs-P*(flux(f).*fSgn);
    if norm(r)/norm(den) > opt.tol;
        warning('Could not construct conservative flux field');
    end

    state.flux = flux;

end

end

%% Helper funciton.

%--------------------------------------------------------------------------

function omega = computeFaceWeights(G, rock, opt)
%   Computes the face weights.

    if strcmp(opt.faceWeights, 'permWeighted')
        
        f = G.cells.faces(:,1);
        ncf = diff(G.cells.facePos);
        fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));
        
        K = permTensor(rock, G.griddim)';
        [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1));
        K = sparse(ii, jj, K(:));
        
        fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f))';
        [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1), ncf);
        fnMat = sparse(ii,jj, fn(:));
        delta = fnMat'*K*fnMat;
        delta = spdiags(delta, 0);

        ii = f;
        jj = (1:numel(f))';
        omega = (sparse(ii, jj,1)*delta);

        for i = 1:G.faces.num
            d = delta(f == i);
            omega(i) = omega(i)/(numel(d)*prod(d));
        end
        
        omega = omega(f);
        
    elseif strcmp(opt.faceWeights, 'tpf')
        
        omega = computeTrans(G, rock)./G.faces.areas(G.cells.faces(:,1));
        
    elseif strcmp(opt.faceWeights, 'ones')
        
        omega = ones(numel(G.cells.faces(:,1)),1);
        
    end
    
end

%--------------------------------------------------------------------------