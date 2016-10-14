function coSet = getCollactionSet(G, rock)
    if nargin == 1
        rock = makeRock(G, 1, 1);
    end
    
    n_hf = size(G.cells.faces, 1);
    dim = G.griddim;
    
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

    cellNo = rldecode(1 : G.cells.num, ...
                     diff(G.cells.facePos), 2) .';
    faceNo = G.cells.faces(:, 1);
    
    for i = 1:n_hf
        c = cellNo(i);
        f = faceNo(i);
        % Self faces are faces belonging to current cell
        selfFaces = gridCellFaces(G, c);
        adjFaces = adjacentFacesForFace(G, f, dim - 1);
        
        % Compute permeability scaled normal vector
        sgn = 1 - 2*(G.faces.neighbors(f, 2) == c);
        normal = G.faces.normals(f, :);
        K = getK(rock, c, dim);
        l = normal*sgn*K;
        L(i, :) = l;
        
        % Cell to half-face collaction
        cpt = G.cells.centroids(c, :);
        [T, type, coeff, ix, successCell] = collocateSet(G, cpt, [], selfFaces, [], l);
        assert(successCell)
        
        for j = 1:dim
            T_all{j, 1}(i, :) = T(j, :);
        end
        type_cell(i, :) = type;
        ix_cell(i, :) = ix;
        c_cell(i, :) = coeff;
        act_cell(i, :) = ix == f & type == 2;
        
        % Half-face to cell collaction
        fpt = G.faces.centroids(f, :);
        hfFaces =  intersect(adjFaces, selfFaces);
%         hfFaces = adjFaces;
        [T, type, coeff, ix, successFace] = collocateSet(G, fpt, c, hfFaces, [], -l);

        
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
        c_face(i, :) = coeff;
        act_face(i, :) = ix == c & type == 1;
    end
    coSet.T = T_all;
    coSet.T_norm = cellfun(@(x) sqrt(sum(x.^2, 2)), T_all, 'UniformOutput', false);
    coSet.C = {c_cell, c_face};
    coSet.types = {type_cell, type_face};
    coSet.globalIndices = {ix_cell, ix_face};
    coSet.l = L;
    coSet.active = {act_cell, act_face};
    
    % Construct mapping operators
    ops = cell(dim, 2);
    Pc = speye(G.cells.num);
    Pf = getFaceFromCellInterpolator(G);
    Pn = getNodeFromCellInterpolator(G);
    
    for i = 1:dim
        for j = 1:2
            gi = coSet.globalIndices{j}(:, i);
            ti =  coSet.types{j}(:, i);
            ops{i, j} = getInterpolationOperator(Pc, Pf, Pn, gi, ti);
        end
    end
    coSet.pressureOperators = ops;
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
    
    t1 = [t1, 0]./norm(t1, 2);
    t2 = [t2, 0]./norm(t2, 2);
    l = [l, 0]./norm(l, 2);
    n = [0, 0, 1];
    C(1) = dot(n, cross(t2, l))/dot(n, cross(t2, t1));
    C(2) = dot(n, cross(t1, l))/dot(n, cross(t1, t2));
    C = round(C, 8);
%     
%     l_new = t1*C(1) + t2*C(2);
%     C = C./norm(l_new, 2);
%     l_new = t1*C(1) + t2*C(2);
%     norm(l - l_new)/norm(l)
%     C
%     D1 = computeDCoefficient2D(l, t2);
%     D2 = computeDCoefficient2D(t1, l);
% 
%     D = computeDCoefficient2D(t1, t2);
%     
%     C = [D1, D2]./D;
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