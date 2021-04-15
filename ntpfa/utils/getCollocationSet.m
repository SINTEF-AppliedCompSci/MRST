function coSet = getCollocationSet(G, rock)
%Undocumented Utility Function

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

    if nargin == 1
        rock = makeRock(G, 1, 1);
    end
    
    n_hf = size(G.cells.faces, 1);
    dim = G.griddim;
    
    jumpFlag = getJumpFlag(G, rock);
    faceSign = zeros(n_hf, 1);
    N = G.faces.neighbors;
    % This builds a collocation set for the nonlinear two-point flux
    % approximation. Generally, we construct the necessary coefficients and
    % vectors as a preprocessing step. 
    % 
    % Most of these sets are constructed twice: Once for the cell to face
    % case, and once for the face to cell case.
    
    % Vectors from d.o.f. to interpolation points
    T_all = cell(dim, 2);
    [T_all{:}] = deal(zeros(n_hf, dim));
    
    % Coefficient matrices for each set of vectors
    [c_cell, c_face]       = deal(zeros(n_hf, dim));
    % Active indicators for each set of vectors/coefficients. A point is
    % said to be active if it is defined in the cell where we want a
    % degree of freedom for a two point method (plus the face between two
    % cells).
    [act_cell, act_face]   = deal(false(n_hf, dim));
    % Global indices for the pressure values we need to sample for each
    % point. Note that the interpretation of the index depends on the type
    % array.
    [ix_cell, ix_face]     = deal(zeros(n_hf, dim));
    % Type array indicating what the indices actually point to:
    %         1: Cell
    %         2: Face
    %         3: Edge ??
    [type_cell, type_face] = deal(zeros(n_hf, dim));
    % L is the face normal vectors, scaled by cell permeability and face
    % area. One for each half face.
    L = zeros(n_hf, dim);
    L_face = zeros(G.faces.num, dim);
    
    cellNo = rldecode(1 : G.cells.num, ...
                     diff(G.cells.facePos), 2) .';
    faceNo = G.cells.faces(:, 1);
    isBF = false(G.faces.num, 1);
    isBF(boundaryFaces(G)) = true;
    for i = 1:n_hf
        c = cellNo(i);
        f = faceNo(i);
        % Self faces are faces belonging to current cell
        selfFaces = gridCellFaces(G, c);
        adjFaces = adjacentFacesForFace(G, f, dim - 1);
        
        % Compute permeability scaled normal vector
        sgn = 1 - 2*(G.faces.neighbors(f, 2) == c);
        faceSign(i) = sgn;
        normal = G.faces.normals(f, :);
        K = getK(rock, c, dim);
        l = normal*sgn*K;
        L(i, :) = l;
        L_face(f, :) = normal*K;
        
        % Cell to half-face collaction
        cpt = G.cells.centroids(c, :);
        
        isJump = jumpFlag(selfFaces);
        fa = selfFaces(isJump);
        cells = N(selfFaces(~isJump), :);
        cells(cells == c) = 0;
        cells = sum(cells, 2);

        [T, type, coeff, ix, successCell] = collocateSet(G, cpt, cells, fa, [], l);
        assert(successCell)
        
        for j = 1:dim
            T_all{j, 1}(i, :) = T(j, :);
        end
        type_cell(i, :) = type;
        ix_cell(i, :) = ix;
        c_cell(i, :) = coeff;%./factor;
        act_cell(i, :) = ix == f & type == 2;
        
        % Half-face to cell collaction
        fpt = G.faces.centroids(f, :);
        hfFaces =  intersect(adjFaces, selfFaces);
        
        isJump = jumpFlag(hfFaces);
        fa = hfFaces(isJump);
        cells = N(hfFaces(~isJump), :);
        cells(cells == c) = 0;
        cells = sum(cells, 2);
        cells = [cells; c];
        
        [T, type, coeff, ix, successFace] = collocateSet(G, fpt, cells, fa, [], -l);

        
        if ~successFace
            selfNodes = unique(gridFaceNodes(G, selfFaces));
            dispif(mrstVerbose > 1, 'Half-face %d of %d...\n', i, n_hf);
            dispif(mrstVerbose > 1, 'Extending definition to nodes, unable to find face basis');
            % Should probably be more careful about these nodes
            [T, type, coeff, ix, successFace] = collocateSet(G, fpt, c, selfFaces, selfNodes, -l);
            assert(successFace)
        end
        
        for j = 1:dim
            T_all{j, 2}(i, :) = T(j, :);
        end
        type_face(i, :) = type;
        ix_face(i, :) = ix;
        c_face(i, :) = coeff;%./factor;
        act_face(i, :) = ix == c & type == 1;
    end
    coSet.jumpFace = jumpFlag;
    coSet.T = T_all;
%     coSet.T_norm = cellfun(@(x) sqrt(sum(x.^2, 2)), T_all, 'UniformOutput', false);
    coSet.C = {c_cell, c_face};
    coSet.types = {type_cell, type_face};
    coSet.globalIndices = {ix_cell, ix_face};
    coSet.l = L;
    coSet.active = {act_cell, act_face};
    
%     exclude_cell = struct('ix', faceNo, 'type', 2);
%     exclude_face = struct('ix', cellNo, 'type', 1);
    % Construct mapping operators
    ops = cell(dim, 2);
    Pc = speye(G.cells.num);
    Pf = getFaceFromCellInterpolator(G, rock);
    Pn = getNodeFromCellInterpolator(G);
    
    for i = 1:dim
        for j = 1:2
            if j == 1
                exclude = G.faces.neighbors(faceNo, 2);
            else
                exclude = G.faces.neighbors(faceNo, 1);
            end
            gi = coSet.globalIndices{j}(:, i);
            ti = coSet.types{j}(:, i);
            ops{i, j} = getInterpolationOperator(Pc, Pf, Pn, gi, ti, exclude);
        end
    end
    coSet.pressureOperators = ops;
    
    
    coSet = storeFaceSet(G, rock, coSet, faceSign, L_face);
    
    
    intx = all(G.faces.neighbors > 0, 2);
    ops = cell(dim, 2);
    for i = 1:dim
        for j = 1:2
            if j == 1
                exclude = G.faces.neighbors(intx, 2);
            else
                exclude = G.faces.neighbors(intx, 1);
            end
            gi = coSet.faceSet.globalIndices{j}(:, i);
            ti =  coSet.faceSet.types{j}(:, i);
            ops{i, j} = getInterpolationOperator(Pc, Pf, Pn, gi, ti, exclude);
        end
    end
    coSet.faceSet.pressureOperators = ops;
end

function K = getK(rock, cell, dim)
    k = rock.perm(cell, :);
    switch numel(k)
        case 1
            % Scalar perm
            K = k;
        case dim
            % Diagonal tensor
            K = diag(k);
        case 3*(dim - 1)
            % Full symmetric tensor
            if dim == 2
                K = [k(1), k(2); ...
                     k(2), k(3)];
            else
                K = [k(1), k(2), k(3); ...
                     k(2), k(4), k(5); ...
                     k(3), k(5), k(6)];
            end
        otherwise
            error('What sorcery is this?!');
    end
end

function [T, types, coefficients, globIx, ok] = collocateSet(G, centerpt, cells, faces, nodes, l)
    [pts, typ, glob] = getCandidatePoints(G, cells, faces, nodes);
    [T, ix, coefficients, ok] = getSetForHF(centerpt, l, pts);
    globIx = glob(ix);
    types = typ(ix);
end

function [T, ixSet, coeff, done] = getSetForHF(center, l, pts)
    dim = size(pts, 2);
%     assert(dim == 3);
    N = size(pts, 1);
    
    t = bsxfun(@minus, pts, center);
    t_v = sqrt(sum(t.^2, 2));
    t_u = bsxfun(@rdivide, t, t_v);
    
    l_v = norm(l, 2);
    l_u = l./l_v;
    
    xl = center + l_u;
    xi = bsxfun(@plus, center, t_u);
    
    
    dist = bsxfun(@minus, xl, xi);
%     dist = cross(xi, repmat(xl, size(xi, 1), 1));
    dist = sum(dist.^2, 2);
    dist = sqrt(sum(dist.^2, 2));
    
    [v, ix] = sort(dist);
    
    
    done = false;
    ptr = 1;
    
    if dim == 3
    	nComb = (N-1)*(N-2)*N;
    else
        nComb = (N-1)*N;
    end
    [coefficents, allsets] = deal(nan(nComb, dim));
    if dim == 3
        for i = 1:(N-2)
            if done
                break;
            end
            for j = (i+1):(N-1)
                if done
                    break;
                end
                for k = (j+1):N
                    I = ix(i);
                    J = ix(j);
                    K = ix(k);

                    t1 = t_u(I, :);
                    t2 = t_u(J, :);
                    t3 = t_u(K, :);

                    [C, D] = getCoefficients3D(t1, t2, t3, l_u);
                    coefficents(ptr, :) = C;

                    allsets(ptr, :) = [I, J, K];

                    if all(C >= 0) %&& abs(D) > 1e-8
                        if all(C <= 1)
                            done = true;
                            break
                        end
                    else
                        coefficents(ptr, :) = inf;
                    end
                    ptr = ptr + 1;
                end
            end
        end
    else
        for i = 1:(N-1)
            if done
                break;
            end
            for j = (i+1):N
                if done
                    break;
                end
                I = ix(i);
                J = ix(j);
                t1 = t_u(I, :);
                t2 = t_u(J, :);
                [C, D] = getCoefficients2D(t1, t2, l_u);
                coefficents(ptr, :) = C;
                allsets(ptr, :) = [I, J];
                if all(C >= 0)% && abs(D) > 1e-8
                    if all(C <= 1)
                        done = true;
                        break
                    end
                else
                    coefficents(ptr, :) = inf;
                end
                ptr = ptr + 1;
            end
        end
    end


    
    if done
        ixSet = allsets(ptr, :);
        coeff = coefficents(ptr, :);
    else
%         assert(ptr == nComb) 
        % Use fallback 
        mx = max(coefficents, [], 2);
        [v, sb] = min(mx);
        coeff = coefficents(sb, :);
        ixSet = allsets(sb, :);
        done = all(isfinite(coeff));
    end

    T = t(ixSet', :);
    for i = 1:numel(coeff)
        coeff(i) = coeff(i)./norm(T(i, :), 2);
    end
%     coeff = bsxfun(@rdivide, coeff, sqrt(sum(T.^2, 2)));;
%     T = bsxfun(@rdivide, T, sqrt(sum(T.^2, 2)));
    assert(all(coeff >= 0));
end

function D = computeDCoefficient3D(a, b, c)
    D = (dot(cross(a, b), c))/(norm(a, 2)*norm(b, 2)*norm(c, 2));
end

function D = computeDCoefficient2D(a, b)
    D = norm(a - b);
%     D = cross(a, b)
%     D = norm(a - b, 2)/(norm(a,2)*norm(b, 2));
end

function [C, D] = getCoefficients2D(t1, t2, l)
    D = nan;
    C = [0, 0];
    
    t1 = [t1, 0];%./norm(t1, 2);
    t2 = [t2, 0];%./norm(t2, 2);
    l = [l, 0];%./norm(l, 2);
    n = [0, 0, 1];
    C(1) = dot(n, cross(t2, l))/dot(n, cross(t2, t1));
    C(2) = dot(n, cross(t1, l))/dot(n, cross(t1, t2));
%     C = round(C, 8);
%     
%     l_new = t1*C(1) + t2*C(2);
%     C = C./norm(l_new, 2);
    l_new = t1*C(1) + t2*C(2);
%     assert(norm(l - l_new)/norm(l) < 1e-12);
%     C
%     D1 = computeDCoefficient2D(l, t2);
%     D2 = computeDCoefficient2D(t1, l);
% 
%     D = computeDCoefficient2D(t1, t2);
%     
%     C = [D1, D2]./D;
end

function flag = getJumpFlag(G, rock)
    flag = true(G.faces.num, 1);
    N = G.faces.neighbors;
    intx = all(N > 0, 2);
    flag(intx) = ~all(rock.perm(N(intx, 1), :) == rock.perm(N(intx, 2), :), 2);
end

function [C, D] = getCoefficients3D(t1, t2, t3, l)
    D1 = computeDCoefficient3D(l, t2, t3);
    D2 = computeDCoefficient3D(t1, l, t3);
    D3 = computeDCoefficient3D(t1, t2, l);

    D = computeDCoefficient3D(t1, t2, t3);
    
    C = [D1, D2, D3]./D;
end


function [pts, types, ix] = getCandidatePoints(G, cells, faces, nodes)
    ix = [cells; faces; nodes];
    ix = ix(:);
    pts = [G.cells.centroids(cells, :); ...
           G.faces.centroids(faces, :); ...
           G.nodes.coords(nodes, :)];
    types = rldecode((1:3)', [numel(cells); numel(faces); numel(nodes)]);
end

function coSet = storeFaceSet(G, rock, coSet, faceSign, L_face)
    dim = G.griddim;
    [C_l, C_r, types_l, types_r, glob_l, glob_r, active_l, active_r] = ...
                                            deal(zeros(G.faces.num, dim));
    left = faceSign == 1;
    right = ~left;
    
    lf = G.cells.faces(left);
    rf = G.cells.faces(right);
    
    C_l(lf, :) = coSet.C{1}(left, :);
    C_r(rf, :) = coSet.C{1}(right, :);

    types_l(lf, :) = coSet.types{1}(left, :);
    types_r(rf, :) = coSet.types{1}(right, :);
    
    glob_l(lf, :) = coSet.globalIndices{1}(left, :);
    glob_r(rf, :) = coSet.globalIndices{1}(right, :);
    
    active_l(lf, :) = coSet.active{1}(left, :);
    active_r(rf, :) = coSet.active{1}(right, :);
    
    
    
    faceSet = struct();
    faceSet.C = {C_l, C_r};
    faceSet.types = {types_l, types_r};
    faceSet.globalIndices = {glob_l, glob_r};
    faceSet.active = {active_l, active_r};
    
    % Remove external
    intx = all(G.faces.neighbors > 0, 2);
    for i = 1:2
        faceSet.C{i} = faceSet.C{i}(intx, :);
        faceSet.types{i} = faceSet.types{i}(intx, :);
        faceSet.globalIndices{i} = faceSet.globalIndices{i}(intx, :);
        faceSet.active{i} = faceSet.active{i}(intx, :);
    end
    faceSet.l = L_face(intx, :);

    exclude_left = struct('ix', G.faces.neighbors(intx, 2), 'type', 1);
    exclude_right = struct('ix', G.faces.neighbors(intx, 2), 'type', 1);

    coSet.faceSet = faceSet;
end
