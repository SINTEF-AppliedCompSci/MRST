function [S, operators] = VEM_assemble(G, C, varargin)
% Assemble the virtual element stiffness matrix and intermediate operators
%
% SYNOPSIS:
%   function [S, operators] = VEM_assemble(G, C, varargin)
%
% DESCRIPTION:
%   Compute the full VEM stiffness matrix; optionally returns intermediate
%   operators
% 
%   Assembly based on Paulino's paper in 3D. The notations follow Paulino's
%   paper.
%
%   Vectorized version. 
% 
%   S = |E| W_c D W_c^T + (I - P_P)^T s (I - P_P)
%
% PARAMETERS:
%   G        - Grid (2D or 3D) C        - matrix where each row represents
%              the elasticity tensor for the grid cell corresponding to
%              that row
%   varargin - options are:
%              'blocksize' - size of blocks (# of cells) for vectorized
%              calc.
%
% RETURNS:
%   S         - full VEM system matrix, including rows/columns for
%   Dirichlet nodes. operators - will contain the fields:
%               - D      - matrix of normalized strain energies for each
%                          cell
%               - WC     - matrix for projection of basis functions onto
%                          space of constant strain modes 'c' 
%               - assemb -
%
% EXAMPLE:
%
% SEE ALSO: `VEM_linElast`

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

    opt = struct('blocksize', 1000, ...
                 'alpha_scaling', 1, ...
                 'experimental_scaling', false, ...
                 'S', []);
    opt = merge_options(opt, varargin{:});

    if (G.griddim == 3)
        qc_all = calculateQC(G);
    else
        qf_all = calculateQF(G);
    end

    D_all = C2D(C, G);

    % Make block sizes
    
    if nargout > 1
        nblocks = 1;
    else
        nblocks = ceil(G.cells.num/opt.blocksize);
    end
    bsizes = floor(G.cells.num/nblocks)*(ones(nblocks, 1));
    rest = G.cells.num-sum(bsizes);
    bsizes(1:rest) = bsizes(1:rest)+1;
    assert(sum(bsizes) == G.cells.num);
    bcellpos = cumsum([1;bsizes]);
    blocks = cell(numel(bsizes), 1);
    for i = 1 : (numel(bcellpos)-1)
        blocks{i} = [bcellpos(i) : bcellpos(i + 1) - 1]';     
    end

    % Initialize a zero matrix to be added in case of block assembly
    ndofs = G.nodes.num*G.griddim;
    num_loc_nodes = sum(diff(G.cells.nodePos));
    S = sparse(ndofs, ndofs);

    for iblocks = 1:numel(blocks)
        cells = blocks{iblocks};
        % Here every thing should be done only for the cells in block
        inodes = mcolon(G.cells.nodePos(cells), (G.cells.nodePos(cells+1)-1))';
        linodes = [1:numel(inodes)]';
        nodes = G.cells.nodes(inodes)';
        nlc = G.cells.nodePos(cells+1)-G.cells.nodePos(cells);
        cellnum = rldecode(cells, nlc);
        lcellnum = rldecode([1:numel(cells)]', nlc);
        
        BB = nan(numel(cells), 1);
        for i = 1:G.griddim
            BB(:, i) = accumarray(lcellnum, G.nodes.coords(nodes, i), [numel(cells), 1]);
        end
        fac = accumarray(lcellnum, 1, [numel(cells), 1]);
        BB = bsxfun(@rdivide, BB, fac);

        % Calculate coordinates relative to center
        XB = G.nodes.coords(nodes, :)-BB(lcellnum, :);
        
        % Dimension of Voigt matrix and therefore also dimension of of the local linear
        % space
        if(G.griddim == 2)
            nlin = 3;        
        else
            nlin = 6; % when G.griddim == 3
        end
        
        nn = numel(nodes);
        zz = zeros(nn, 2);
        z = zeros(nn, 1);

        % Construction of NC and NR, which give the nodal representations of the local
        % linear transformations.
        
        NC = zeros(G.griddim*nn, nlin);
        NR = zeros(G.griddim*nn, nlin);
        if(G.griddim == 3)
            NC(:, 1:3) = [reshape([XB(:, 1), zz]', [], 1), reshape([z, XB(:, 2), z]', [], 1), reshape([zz, XB(:, 3)]', [], 1)];
        else
            NC(:, 1:2) = [reshape([XB(:, 1), z]', [], 1), reshape([z, XB(:, 2)]', [], 1)];
        end
        NR(:, 1:G.griddim) = repmat(eye(G.griddim), nn, 1);
        
        if(G.griddim == 3)
            perm = {[2, 1, 3], [1, 3, 2], [3, 2, 1]};
            nrv = {[1, -1, 0], [0, 1, -1], [-1, 0, 1]};
            nrc = {[1, 1, 0], [0, 1, 1], [1, 0, 1]};
        else
            perm = {[2, 1]};
            nrv = {[1, -1]};
            nrc = {[1, 1]};
        end
        for i = 1:numel(perm);
            tmp = XB(:, perm{i});
            NR(:, G.griddim+i) = reshape(bsxfun(@times, tmp, nrv{i})', [], 1);
            NC(:, G.griddim+i) = reshape(bsxfun(@times, tmp, nrc{i})', [], 1);
        end
        
        WC = zeros(G.griddim*nn, nlin);
        WR = zeros(G.griddim*nn, nlin);
        nlcl = rldecode(nlc, nlc);
        if (G.griddim == 3)
            % XX is 2*q_i in [Gain et al : doi:10.1016/j.cma.2014.05.005] eqs 77
            XX = bsxfun(@rdivide, qc_all(inodes, :), G.cells.volumes(cellnum));
            WC(:, 1:3) = [reshape([XX(:, 1), zz]', [], 1), reshape([z, XX(:, 2), z]', [], 1), reshape([zz, XX(:, 3)]', [], 1)];
            XX = XX/2; % now XX is q_i
            WR(:, 1:3) = [reshape([1./nlcl, zz]', [], 1), reshape([z, 1./nlcl, z]', [], 1), reshape([zz, 1./nlcl]', [], 1)];
        else
            % XX is 2*q_i in [Gain et al : doi:10.1016/j.cma.2014.05.005] eqs 77
            XX = bsxfun(@rdivide, qf_all(inodes, :), G.cells.volumes(cellnum));
            WC(:, 1:2) = [reshape([XX(:, 1), z]', [], 1), reshape([z, XX(:, 2)]', [], 1)];
            XX = XX/2; % now XX is q_i
            WR(:, 1:2) = [reshape([1./nlcl, z]', [], 1), reshape([z, 1./nlcl]', [], 1)];
        end
        
        for i = 1:numel(perm)
            tmp = XX(:, perm{i});
            WR(:, G.griddim+i) = reshape(bsxfun(@times, tmp, nrv{i})', [], 1);
            WC(:, G.griddim+i) = reshape(bsxfun(@times, tmp, nrc{i})', [], 1);
        end
        
        % Make block matrices for the operators efficently do local cell operations
        [ib, jb] = blockDiagIndex(nlin*ones(size(nlc)), G.griddim*nlc);
        mat_sz = numel(nodes)*G.griddim;
        WR = sparse(ib, jb, reshape(WR', [], 1), nlin*numel(cells), mat_sz)';
        WC = sparse(ib, jb, reshape(WC', [], 1), nlin*numel(cells), mat_sz)';
        NR = sparse(ib, jb, reshape(NR', [], 1), nlin*numel(cells), mat_sz)';
        NC = sparse(ib, jb, reshape(NC', [], 1), nlin*numel(cells), mat_sz)';
        
        % Define the projection operators in nodal basis
        PR = NR*WR'; % Rigid body rotation
        PC = NC*WC'; % Other linear displacements
        PP = PR+PC;  % Projection to linear displacement
        
        inda = 1 : nlin;
        indd = nlin*(inda-1) + inda;
        trD = sum(D_all(cells, indd), 2); % trace of D_all

        D = reshape(D_all(cells, :)', nlin, [])';
        [i, j] = blockDiagIndex(nlin*ones(numel(cells), 1), nlin*ones(numel(cells), 1));
        D = sparse(i, j, reshape(D', [], 1), nlin*numel(cells), nlin*numel(cells));

        vol = G.cells.volumes(cells);
        volfac = rldecode(vol, nlc*G.griddim);
        nldofs = G.griddim*numel(inodes);
        volmap = sparse(1:nldofs, 1:nldofs, volfac, nldofs, nldofs); % not always necessary 
        assert(size(PR, 1) == nldofs);
        I = sparse(1:nldofs, 1:nldofs, 1, nldofs, nldofs);

        if isempty(opt.S)
            % Calculate the inner product by adding the consistent part on linear
            % displacement and the stabilization part for the orthogonal of linear
            % displacements.
            DNC = diag(NC'*NC);
            if(~opt.experimental_scaling)
                trDNC = sum(reshape(DNC, nlin, []), 1)';  
                alpha = rldecode((trD.*vol)./(trDNC), nlc*G.griddim); % Not properly scaled
                if(any(alpha<0))
                    warning('Natural scaling variable is negative');
                end
                alpha = abs(alpha);               
            else
                %assert(G.griddim == 2);
                trDNCn = reshape(DNC, nlin, []);
                trDNCn = sum(1./trDNCn, 1)';% use all
                alpha = rldecode((1/nlin^2)*trD.*vol.*trDNCn, nlc*G.griddim); 
            end
            alpha = alpha .* opt.alpha_scaling;
            SE = sparse(1:nldofs, 1:nldofs, alpha, nldofs, nldofs);
        else
            if iblocks == 1
                SE = opt.S;
            else
                SE = 0 * opt.S;
            end
        end
        
        KH = volmap*WC*D*WC' + (I - PP)'*SE*(I - PP);

        % Proceed with the global assembly
        dofsg = mcolon(G.griddim*(nodes - 1) + 1, G.griddim*(nodes - 1) + G.griddim);
        dofsl = mcolon(G.griddim*(linodes - 1) + 1, G.griddim*(linodes - 1) + G.griddim);
        assemb = sparse(dofsg, dofsl, 1, G.griddim*G.nodes.num, numel(linodes)*G.griddim);
        SL = assemb*KH*assemb';
        S = S + SL;
        if nargout > 1
            operators.I      = I;
            operators.volmap = volmap; 
            operators.D      = D;
            operators.WC     = WC;
            operators.KH     = KH;
            operators.SE     = SE;
            operators.PP     = PP;
            operators.assemb = assemb;
        end
    end
end
