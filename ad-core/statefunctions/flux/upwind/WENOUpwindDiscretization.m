classdef WENOUpwindDiscretization < UpwindDiscretization
    % Weighted essentially non-oscillatory scheme for upwinding
    properties
        dim;                        % Dimension of discretization
        includeBoundary = true;     % Include boundary values in interpolation
        interpolateReference = true;% Interpolate in a scaled reference space
        useMatrixProducts = true;   % Faster implementation for many cases
        compactStencil  = true;     % Take only a single connection for multiple neighbors
        cap = true;                 % Limit values to positive
    end
    
    properties (Access = protected)
        interp_setup = struct();
    end
    
    methods
        function weno = WENOUpwindDiscretization(model, dim, varargin)
            weno@UpwindDiscretization();
            if nargin == 1
                weno.dim = model.G.griddim;
            else
                weno.dim = dim;
            end
            weno = merge_options(weno, varargin{:});
            weno = weno.setupInterpolators(model);
        end
        
        function fv = faceUpstream(weno, model, state, flag, cellvalue)
            switch weno.dim
                case 1
                    fv = weno.getFaceValue1D(model, state, flag, cellvalue);
                case {2, 3}
                    fv = weno.getFaceValueTri(model, state, flag, cellvalue);
                otherwise
                    error('Dimension %d too high', weno.dim);
            end
        end
        
        
        function disc = setupInterpolators(disc, model)
            switch disc.dim
                case 1
                    disc = disc.setupInterpolators1D(model);
                case 2
                    disc = disc.setupInterpolatorsTri(model);
                case 3
                    disc = disc.setupInterpolatorsTri(model);
                otherwise
                    error('Dimension too high');
            end
        end
        % --- Helpers for evaluation --- %
        function [vq, faces, w, state] = getFaceValueTri(disc, model, state, flag, cellvalue)
            % WENO using triangulation - 2D or 3D, general grids
            N = model.operators.N;
            cells = flag.*N(:, 1) + ~flag.*N(:, 2);
            
            % Face centroids expanded to triangles
            fc = model.G.faces.centroids(model.operators.internalConn, 1:disc.dim);
            [vq, faces, w] = disc.interpolateTriWENO(model, state, cells, fc, cellvalue);
        end
        
        function [q_interp, pointNo, w] = interpolateTriWENO(disc, model, state, cells, points, q)
            % Interpolate using WENO reconstruction for a set of points,
            % with reconstruction taken from weighted polynomial for each
            % cell.
            % Find the triangles corresponding to each cell together with
            % each set of linear weights and the number of cells in the
            % support
            nf = size(points, 1);
            tri_ix = vertcat(disc.interp_setup.cell_support{cells});
            linear_weights = vertcat(disc.interp_setup.linear_weights{cells});
            supCount = disc.interp_setup.cell_support_count(cells);
            
            
            % For each triangle, list the corresponding cell
            centers = rldecode(cells, supCount);
            % For each triangle, the corresponding point index
            pointNo = rldecode((1:nf)', supCount);
            
            % Each cell has a scaling factor for the local coordinates
            scaling = disc.interp_setup.scaling;
            scale = scaling.scale(centers, 1:disc.dim);
            points = points(pointNo, :);
            
            % Distance from point to center centroid where values are
            % defined
            dx = points - model.G.cells.centroids(centers, 1:disc.dim);
            
            % Apply shift to local coordinate system
            if ~isempty(scaling.mapping)
                dx0 = dx;
                for d = 1:disc.dim
                    dx(:, d) = sum(scaling.mapping(centers, :, d).*dx0, 2);
                end
            end
            dx = dx./scale;
            
            % Compute gradient on each triangle. Sigma{i}(j) is the gradient
            % in direction i, on triangle j.
            sigma = cell(1, disc.dim);
            [sigma{:}] = deal(0);
            
            if disc.useMatrixProducts
                Q_tri = sparse((1:numel(tri_ix))', tri_ix, 1, ...
                    numel(tri_ix), size(disc.interp_setup.C, 1));
                Qc = sparse((1:numel(cells))', cells, 1, numel(cells), model.G.cells.num);
            end
            nw = numel(pointNo);
            M = sparse(pointNo, (1:nw)', 1);
            
            if disc.useMatrixProducts
                basis_mat = disc.interp_setup.tri_basis_matrices;
                for d = 1:disc.dim
                    sigma{d} = basis_mat{d}*q;
                end
            else
                for d = 1:disc.dim
                    b = disc.interp_setup.tri_basis{d};
                    for l = 1:disc.dim+1
                        loc_cells = disc.interp_setup.C(:, l);
                        ccl = disc.interp_setup.tri_cells(loc_cells);
                        ds = b(:, l).*q(ccl);
                        
                        sigma{d} = sigma{d} + ds;
                    end
                end
            end
            
            % Get the cell values and interpolate
            if disc.useMatrixProducts
                q_c = Qc*q;
            else
                q_c = q(cells);
            end
            % We obtain the interpolated value for each cell. This is done
            % by computing q(x) = grad q(x0) dot (x - x0)
            if disc.useMatrixProducts
                ns = size(value(sigma{1}), 1);
                ss = vertcat(sigma{:});
                dd = size(value(ss), 1);
                I = repmat((1:nw)', 1, disc.dim);
                J = tri_ix + (0:disc.dim-1)*ns;
                V = dx;
                TF = sparse(I, J, V, nw, dd);
                ws = TF*ss;
            else
                % This branch never executes, old equivialent code
                ws = 0;
                
                for d = 1:disc.dim
                    if disc.useMatrixProducts
                        dq = Q_tri*sigma{d}.*dx(:, d);
                    else
                        dq = sigma{d}(tri_ix).*dx(:, d);
                    end
                    ws = ws + dq;
                end
            end
            
            % Get the gradient norm for each triangle,
            % sigma_norm = (dq/dx)^2 + (dq/dy)^2 + (dq/dz)^2
            sigma_norm = 0;
            for i = 1:disc.dim
                ds = sigma{i}.^2;
                sigma_norm = sigma_norm + ds;
            end
            if disc.useMatrixProducts
                sigma_norm_tri = Q_tri*sigma_norm;
            else
                sigma_norm_tri = sigma_norm(tri_ix);
            end
            % Compute smoothness indicator beta
            L = 2;
            E = 1e-12;
            beta = (sigma_norm_tri + E).^(-L);
            beta = beta.*linear_weights;
            
            % For each cell, the weight of each interpolated value is the
            % smoothness of each gradient relative to the sum of gradients.
            if 0
                w = beta./((M'*M)*beta);
            else
                beta_tot = M*beta;
                w = beta./beta_tot(pointNo);
            end
            w(~isfinite(value(w))) = 0;
            
            dvq = M*(ws.*w);
            q_interp = q_c + dvq;
            if disc.cap
                bad = q_interp <= 0;
                q_interp(bad) = q_c(bad);
            end
        end
        
        function qv = getFaceValue1D(disc, model, state, flag, cellvalue)
            G = model.G;
            x = G.cells.centroids(:, 1);
            xq = G.faces.centroids(model.operators.internalConn, 1);
            N = model.operators.N;
            mid = flag.*N(:, 1) + ~flag.*N(:, 2);
            
            left = disc.interp_setup.cellNeighbors(mid, 1);
            right = disc.interp_setup.cellNeighbors(mid, 2);

            [q_L, B_L] = interpolate1D(disc, left, mid, x, cellvalue, xq);
            [q_R, B_R] = interpolate1D(disc, mid, right, x, cellvalue, xq);
            
            B_T = B_L + B_R;
            
            w_L = B_L./B_T;
            w_R = B_R./B_T;
            
            qv = w_L.*q_L + w_R.*q_R;
        end

        % --- Helpers for setup --- %
        function disc = setupInterpolatorsTri(disc, model)
            G = model.G;
            nc = G.cells.num;
            disp('Getting supports.')
            tic()
            [C, pts, cells, basis, supports, linear_weights, scaling] = disc.getTriangulation(model);
            toc();
            disp('Inverting systems...');
            tic()
            
            disc.interp_setup.tri_cells = cells;
            disc.interp_setup.tri_basis = basis;
            disc.interp_setup.tri_points = pts;
            disc.interp_setup.linear_weights = linear_weights;
            disc.interp_setup.cell_support = supports;
            disc.interp_setup.scaling = scaling;
            disc.interp_setup.C = C;
            
            disc.interp_setup.cell_support_count = cellfun(@numel, disc.interp_setup.cell_support);
            
            
            [n_tri, d_tri] = size(C);
            basis_mat = cell(1, disc.dim);
            I = repmat((1:n_tri)', 1, d_tri);
            J = disc.interp_setup.tri_cells(disc.interp_setup.C);
            for d = 1:disc.dim
                M = sparse(I, J, disc.interp_setup.tri_basis{d});
                basis_mat{d} = M;
            end
            disc.interp_setup.tri_basis_matrices = basis_mat;
        end
        
        function [C, pts, cells, grad_basis, supports, linear_weights, scaling] = getTriangulation(disc, model)
            % C - Connectivity matrix for triangulation
            % pts - Points used in the triangulation
            % supports - Cell array for each grid cell with supported
            % triplets
            % linear_weights cell array with linear weights
            % scaling - scaling factors for tri
            
            
            % Set up triangulation
            G = model.G;
            nc = G.cells.num;
            pts = G.cells.centroids(:, 1:disc.dim);
            if disc.includeBoundary
                [bfaces, bcells, bf_cent] = getBoundaryFacesMerged(G, disc.dim);
                pts = [pts; bf_cent];
            else
                [bfaces, bcells] = deal([]);
            end
            
            is_boundary_face = false(G.faces.num, 1);
            is_boundary_face(bfaces) = true;
            
            N_global = G.faces.neighbors;
            l = N_global(:, 1);
            r = N_global(:, 2);
            N_global(l == 0, 1) = N_global(l == 0, 2);
            N_global(r == 0, 2) = N_global(r == 0, 1);
            
            
            bnd_map = zeros(G.faces.num, 1);
            bnd_map(bfaces) = 1:numel(bfaces);
            
            cells = [(1:nc)'; bcells];
            supports = cell(nc, 1);
            linear_weights = cell(nc, 1);
            
            doSVD = true;
            if doSVD
                coord_mapping = zeros(nc, disc.dim, disc.dim);
            else
                coord_mapping = [];
                M = 1;
            end
            
            min_values = zeros(nc, disc.dim);
            ranges = ones(nc, disc.dim);
            pts_tri = pts;
            
            % Local approach
            if 1
                N = model.operators.N;
            else
                N = neighboursByNodes(G);
            end
            [tri, coord_inverse] = deal(cell(nc, 1));
            gbases = cell(nc, disc.dim);
            offset = 0;
            
            for i = 1:nc
                if mod(i, 100) == 0
                    fprintf('%d of %d...\n', i, nc);
                end
                fpos = G.cells.facePos(i):G.cells.facePos(i+1)-1;
                local_faces = G.cells.faces(fpos, 1);
                if size(G.cells.faces, 2) > 1 && ...
                        any(ismember(G.type, {'tensorGrid', 'processGRDECL'})) && ...
                        disc.compactStencil
                    % Pick only one neighbor in each cardinal direction,
                    % for corner-point type grids
                    tags = G.cells.faces(fpos, 2);
                    
                    areas = G.faces.areas(local_faces);
                    act = false(numel(local_faces), 1);
                    for j = 1:max(tags)
                        tmp = areas.*(tags == j);
                        [~, ix] = max(tmp);
                        act(ix) = true;
                    end
                    local_faces = local_faces(act);
                end
                L = N_global(local_faces, 1);
                R = N_global(local_faces, 2);
                
                L_ok = L ~= i;
                neighbors = L_ok.*L + ~L_ok.*R;
                local_cells = [i; neighbors];
                
                
                global_cells = local_cells;
                is_bf_local = find(is_boundary_face(local_faces))+1;
                global_cells(is_bf_local) = bnd_map(local_faces(is_bf_local-1)) + nc;
                
                lpts = pts_tri(global_cells, :);
                
                if doSVD
                    if 0
                        % nodes
                        nn = gridCellNodes(G, i);
                        fc = G.nodes.coords(nn, :);
                    elseif 1
                        % faces
                        fc = G.faces.centroids(local_faces, :);
                    elseif 1
                        fc = lpts;
                    else
                        nn = gridCellFaces(G, i);
                        fc = [lpts; G.faces.centroids(nn, :)];
                    end
                    min_pts = G.cells.centroids(i, :);
                    [U,S,V] = svd(fc - min_pts);
                    range = diag(S)';
                    M = inv(V);
                    lpts_scaled = projectPoints(lpts, min_pts, M, range);
                    
                    coord_inverse{i} = M;
                    for d = 1:disc.dim
                        coord_mapping(i, :, d) = M(d, :);
                    end
                elseif 1
                    [lpts_scaled, range, min_pts] = disc.scalePoints(lpts);
                else
                    if 1
                        % nodes
                        nn = gridCellNodes(G, i);
                        fc = G.nodes.coords(nn, :);
                    else
                        % faces
                        nn = gridCellFaces(G, i);
                        fc = G.faces.centroids(nn, :);
                    end
                    min_pts = min(fc);
                    range = max(fc) - min_pts;
                    lpts_scaled = disc.scalePoints(lpts, range, min_pts);
                end
                
                if numel(global_cells) >= disc.dim + 1
                    if 1
                        gpts = pts(global_cells(2:end), :);
                        centers = projectPoints(gpts, min_pts, M, range);
                        pcenter = zeros(1, disc.dim);
                        fpoints = [pcenter; centers];
                        
                        t = delaunayn(fpoints);
                    else
                        t = delaunayn(lpts_scaled);
                    end
                    T = global_cells(t);
                    if size(t, 1) == 1
                        T = T';
                    end
                    
                    % Remove triangles not touching current cell
                    keep = any(T == i, 2);
                    T = T(keep, :);
                    t = t(keep, :);
                    
                    % Get basis
                    if disc.interpolateReference
                        ranges(i, :) = range;
                        min_values(i, :) = min_pts;
                        pts_basis = lpts_scaled;
                    else
                        pts_basis = lpts;
                    end
                    b = disc.getGradientBasis(t, pts_basis);
                    [gbases{i, :}] = deal(b{:});
                    [vol, cent] = disc.getTriVolume(t, pts_basis);
                    lw = vol;
                    linear_weights{i} = lw;
                    
                else
                    assert(0);
                end
                nt = size(T, 1);
                tri{i} = T;
                supports{i} = reshape((offset+1):(offset+nt), [], 1);
                
                offset = offset + nt;
            end
            C = vertcat(tri{:});
            grad_basis = cell(disc.dim, 1);
            for d = 1:disc.dim
                grad_basis{d} = vertcat(gbases{:, d});
            end
            
            scaling.scale = ranges;
            scaling.shift = -min_values;
            scaling.mapping = coord_mapping;
            scaling.inverses = coord_inverse;
        end
        
        function disc = setupInterpolators1D(disc, model)
            % Take x-coordinate since MRST does not have 1D grids
            G = model.G;
            x = G.cells.centroids(:, 1);
            
            N = model.operators.N;
            
            % Each cell has two neighbors, ordered from left to right
            disc.interp_setup.cellNeighbors = zeros(G.cells.num, 2);
            xf = G.faces.centroids(model.operators.internalConn, 1);
            
            L = N(:, 1);
            R = N(:, 2);
            
            smaller = xf < x(R);
            disc.interp_setup.cellNeighbors(R(smaller), 1) = L(smaller);
            disc.interp_setup.cellNeighbors(R(~smaller), 2) = L(~smaller);
            smaller = xf < x(L);
            disc.interp_setup.cellNeighbors(L(smaller), 1) = R(smaller);
            disc.interp_setup.cellNeighbors(L(~smaller), 2) = R(~smaller);
        end
        
        function [vq, smoothness] = interpolate1D(disc, l, r, x, v, xq)
            % Constant interpolation for boundaries
            l(l == 0) = r(l == 0);
            r(r == 0) = l(r == 0);
            
            x_r = x(r);
            x_l = x(l);
            v_l = v(l);
            v_r = v(r);
            
            dx = x_r - x_l;
            bad = dx == 0;
            dx(bad) = 1;
            
            sigma = (v_r - v_l)./dx;
            sigma(bad) = 0;
            
            vq = v_l + sigma.*(xq - x_l);
            
            E = 1e-7;
            L = 2;
            smoothness = (dx.*sigma.^2 + E).^(-L);
        end
        
        function [w, cells] = getSingleWeights(disc, tri_num, center, point)
            s = disc.interp_setup;
            assert(numel(tri_num) == 1);
            tri_pt_num = s.C(tri_num, :);
            cells = s.tri_cells(tri_pt_num);
            dx = disc.G.cells.centroids(center, :) - point;
            
            nw = numel(cells);
            
            w = zeros(nw, 1);
            for i = 1:nw
                for d = 1:disc.dim
                    b = disc.interp_setup.tri_basis{d}(tri_num, :);
                    w(i) = w(i) + b(i).*dx(d);
                end
            end
        end
        
        function [basis, C_inv] = getGradientBasis(disc, C, pts)
            [n_tri, d_tri] = size(C);
            basis = cell(1, disc.dim);
            [basis{:}] = deal(zeros(n_tri, d_tri));
            
            C_inv = cell(1, disc.dim);
            ws = warning('query', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:nearlySingularMatrix');
            for i = 1:n_tri
                if all(C(i, :) == C(i, 1))
                    % This is an upwind "triangle"
                    bi = zeros(disc.dim, disc.dim + 1);
                else
                    [bi, C_inv{i}] = getBasis(pts(C(i, :)', :), d_tri);
                end
                
                for j = 1:disc.dim
                    basis{j}(i, :) = bi(j, :);
                end
            end
            warning(ws.state, ws.identifier);
        end
        
        function [vol, cent] = getTriVolume(disc, C, pts)
            p1 = pts(C(:, 1), :);
            p2 = pts(C(:, 2), :);
            p3 = pts(C(:, 3), :);
            v1 = p2 - p1;
            v2 = p3 - p1;
            cent = p1 + p2 + p3;
            if disc.dim == 2
                z = zeros(size(v1, 1), 1);
                vol = 0.5*sqrt(sum(cross([v1, z], [v2, z]).^2, 2));
            else
                p4 = pts(C(:, 4), :);
                cent = cent + p4;
                v3 = p4 - p1;
                vol = (1/6)*sqrt(sum(sum(cross(v1, v2).*v3, 2).^2, 2));
            end
            cent = cent/(disc.dim+1);
        end
        
        function disc = getSubsetDiscretization(disc0, cells)
            G = disc0.G;
            if islogical(cells)
                cells = find(cells);
            end
            dof_cells = cells;
            
            
            tri = vertcat(disc0.interp_setup.cell_support{cells});
            supp = disc0.interp_setup.tri_cells(disc0.interp_setup.C(tri, :));
            cells = unique(supp(:));
            
            
            disc = disc0;
            
            disc.G.cells.num = numel(cells);
            disc.G.cells.centroids = disc.G.cells.centroids(cells, :);
            
            keep = false(G.cells.num, 1);
            keep(cells) = true;
            
            renum = nan(G.cells.num, 1);
            renum(cells) = 1:numel(cells);
            
            
            keep_all = false(G.cells.num, 1);
            keep_all(dof_cells) = true;
            keep_face = all(keep_all(disc.N), 2);
            
            disc.N = renum(disc.N(keep_face, :));
            if sum(keep_face) == 1
                disc.N = reshape(disc.N, 1, []);
            end
            
            is = disc.interp_setup;
            is.cell_support = is.cell_support(keep);
            is.linear_weights = is.linear_weights(keep);
            is.cell_support_count = is.cell_support_count(keep);
            
            is.scaling.scale = is.scaling.scale(keep, :);
            is.scaling.shift = is.scaling.shift(keep, :);
            
            offsets = cumsum([1; disc0.interp_setup.cell_support_count]);
            
            tri_keep = mcolon(offsets(cells), offsets(cells+1)-1)';
            
            tri_renum = nan(size(tri_keep));
            tri_renum(tri_keep) = 1:numel(tri_keep);
            
            is.cell_support = cellfun(@(x) tri_renum(x), is.cell_support, 'UniformOutput', false);
            
            
            is.C = is.C(tri_keep, :);
            
            for i = 1:numel(is.tri_basis_matrices)
                is.tri_basis_matrices{i} = is.tri_basis_matrices{i}(tri_keep, keep);
                is.tri_basis{i} = is.tri_basis{i}(tri_keep, :);
            end
            
            is.tri_cells = renum(is.tri_cells);
            disc.interp_setup = is;
        end
    end
    methods (Static)
        function [pts_scaled, range, mv] = scalePoints(pts, range, mv)
            if nargin < 3
                mv = min(pts);
                if nargin < 2
                    range = max(pts) - mv;
                end
            end
            pts_scaled = (pts - mv)./range;
        end
    end
end

function [basis, C_inv] = getBasis(tri_pts, Tdim)
C = [tri_pts'; ones(1, Tdim)];
C_inv = inv(C);
C_inv(~isfinite(C_inv)) = 0;
basis = C_inv(:, 1:Tdim-1)';
end

function pts = projectPoints(pts, center, M, range)
pts = pts - center;
for d = 1:size(pts, 1)
    pts(d, :) = M*pts(d, :)';
end
pts = pts./range;
end

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
