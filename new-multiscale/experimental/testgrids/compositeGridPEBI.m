function varargout = compositeGridPEBI(dims, pdims, varargin)
    nx = dims(1);
    ny = dims(2);
    
    dx = pdims(1)/(nx - 1);
    dy = pdims(2)/(ny - 1);

    
    opt = struct('makePEBI', true, ...
                 'randomPerturb', 0.5,...
                 'padding', 1,  ...
                 'uniform', true, ...
                 'wells',   [], ...
                 'radius',  3*sqrt(dx*dy), ...
                 'wellseg', 20, ...
                 'radnum', 10, ...
                 'growthfactor', 1.5, ...
                 'extraPts',    [], ...
                 'lines',       {{}}, ...
                 'randomizePoints', false);
         
    opt = merge_options(opt, varargin{:});

    
    vx = 0:dx:pdims(1);
    vy = 0:dy:pdims(2);
    [X, Y] = meshgrid(vx, vy);


    [ii, jj] = meshgrid(1:nx, 1:ny);

    nedge = opt.padding;
    exterior = (ii <= nedge | ii > nx - nedge) | ...
               (jj <= nedge | jj > ny - nedge);


    if opt.randomizePoints
        rx = opt.randomPerturb*rand(size(X)).*dx;
        ry = opt.randomPerturb*rand(size(Y)).*dy;
    else
        if opt.uniform
            rx = 0.5*dx*(mod(jj, 2) == 0);
            ry = zeros(size(rx));
        else
            rx = 0.25*dx*(mod(jj, 2) == 0);
            ry = 0.25*dy*(mod(ii, 2) == 0);
        end
    end

    interior = ~exterior;

    X(interior) = X(interior) + rx(interior);
    Y(interior) = Y(interior) + ry(interior);

    if opt.randomizePoints
        X(interior) = min(X(interior), pdims(1) - nedge*dx);
        X(interior) = max(X(interior), nedge*dx);

        Y(interior) = min(Y(interior), pdims(2) - nedge*dy);
        Y(interior) = max(Y(interior), nedge*dy);
    end

    Pts = [X(:), Y(:)];
    Pts0 = Pts;
    if ~isempty(opt.wells)
        R = expandToSize(opt.radius, opt.wells);
        WS = expandToSize(opt.wellseg, opt.wells);
        NR = expandToSize(opt.radnum, opt.wells);
        GF = expandToSize(opt.growthfactor, opt.wells);
        for i = 1:size(opt.wells, 1)
            Pts = insertWellRefinement(opt.wells(i, :), Pts, R(i), WS(i), NR(i), GF(i));
        end
    end
    
    if ~isempty(opt.extraPts)
        Pts = [Pts; opt.extraPts];
        Pts = unique(Pts, 'rows');
    end
    
    C = [];
    faultType = zeros(size(Pts, 1), 1);
    for i = 1:numel(opt.lines)
        l = opt.lines{i};
        
        assert(all(size(l) == [2 2]));
        p1 = l(1, :);
        p2 = l(2, :);
        
        v = p2 - p1;
        
        dists = norm(v, 2)/norm([dx, dy], 2);
        dists = 2.5*max(ceil(dists), 2);
        
        % Expand into a vector
        l = p1;
        for j = 1:dists
            l = [l; p1 + v*j/dists];
        end
        
        if opt.makePEBI
            n1 = v./norm(v, 2);
            n2 = [-n1(2), n1(1)];
            
            factor = 1.5*norm(n2.*[dx, dy], 2)/2;
            
            left = bsxfun(@plus, l, n2.*factor);
            right = bsxfun(@minus, l, n2.*factor);
            
            nl = size(l, 1);
            
            [Pts, new, removed] = replacePointsByHull(Pts, [left; right]);
            faultType = faultType(~removed);
            faultType = [faultType; rldecode([1; 2], [nl; nl])];
        else
            N = size(Pts, 1);

            Pts = [Pts; l];
            v = (1:size(l, 1) - 1)';
            C = [C; [v, v + 1] + N];
        end
    end
    
    if ~isempty(C)
        C0 = C;
        [Pts, ij, ik] = unique(Pts, 'rows', 'stable');
        C = ik(C);
        Tri = delaunayTriangulation(Pts, C);
    else
        Tri = delaunayTriangulation(Pts);
    end
    
    G = triangleGrid(Pts, Tri.ConnectivityList);
    if opt.makePEBI
        G = pebi(G);
        
        N = G.faces.neighbors + 1;
        faultType = [0; faultType];
        
        ft1 = faultType(N(:, 1));
        ft2 = faultType(N(:, 2));
        G.faces.tag = double(ft1 ~= ft2 & ft1 > 0 & ft2 > 0);
        indicator = ~ismember(Pts, Pts0, 'rows');
    else
        % Should tags for triangular grids here at some point.
        indicator = false(G.cells.num, 1);
    end
    
    varargout{1} = G;
    if nargout > 1
        varargout{2} = indicator;
    end
end

function [Pts, isWell] = insertWellRefinement(pt, Pts, r, N, nr, growthfactor)
    Pw = [];
    rad = growthfactor.^-(1:nr);

    rad = rad - min(rad);
    rad = rad./max(rad);
    
    rad = r.*rad;
    for r = rad
        [x,y,z] = cylinder(r,N); 
        Pw = [Pw; [x(1,:); y(1, :)]'];
    end
    Pw = [Pw; [0, 0]];
    
    Pw = bsxfun(@plus, Pw, pt);
    
    Pw = unique(Pw, 'rows');
    [Pts, isWell] = replacePointsByHull(Pts, Pw);
end

function x = expandToSize(x, target)
    if numel(x) == 1
        x = repmat(x, size(target, 1), 1);
    end
    assert(numel(x) == size(target, 1));
end

function [Pts, isNew, removed] = replacePointsByHull(Pts, P_target)
    Tri = delaunayTriangulation(P_target);
    keep = isnan(Tri.pointLocation(Pts));
    
    Pts = [Pts(keep, :); P_target];
    isNew = true(size(Pts, 1), 1);
    isNew(1:sum(keep)) = false;
    
    removed = ~keep;
end
