function G = makeLayeredGridNWM(g, p, varargin)
% Extrude 2D grid to layered 3D grid according the topology of 2D grid and
% provided surface point set. The surface points are given on all surfaces, 
% and topologically aligned in layered direction.
%
% SYNOPSIS:
%   G = makeLayeredGridNWM(g, p)
%   G = makeLayeredGridNWM(g, p, 'connectivity', t)
%
% PARAMETERS:
%   g  - The 2D grid to be extruded
%   p  - Points of all surfaces, topologically aligned in layered direction
%
% KEYWORD ARGUMENTS:
%  'connectivity'  - Connectivity list (nodes of cells) for g, 
%                    ncell_g x 1 cell. Note if the 2D grid are not on the 
%                    XY plane, the connectivity list of the 2D grid should
%                    be provided
%
% RETURNS:
%   G  - Valid 2.5D layerd grid structure
%
% EXAMPLE:
%   N     = 10;
%   N1    = 2*N;
%   N2    = 3*ceil(N/2)-2;
%   [X, Y] = meshgrid(0:1:N1, 0:1:N2);
%   X     = sqrt(3) / 2 * X;
%   Y(:,2:2:end)=Y(:,2:2:end)+0.5;
%   p     = [X(:), Y(:)];
%   t     = delaunayn(p);
%   g     = triangleGrid(p, t);
%   g     = pebi(g);
%   p     = g.nodes.coords;
%   z     = (0:10:50)';
%   pSurfs = arrayfun(@(z)[p-z/10, z*ones(size(p,1),1)], z, 'UniformOutput', false);
%   G = makeLayeredGridNWM(g, pSurfs);
%   figure, plotGrid(G), view(3)
% 
% SEE ALSO:
%   `makeLayeredGrid` `tessellationGrid` `layeredGrids`

    opt = struct('connectivity', []);
    opt = merge_options(opt, varargin{:});
    
    p = convertToColumn(p);
    assert(size(p, 2) == 1, 'Point cell array must be one dimensional!')
    
    % All sets of points should have 3 columns
    assert( all( cellfun(@(x)size(x, 2) == 3, p) ))

    % All sets of points should have same rows
    np = unique( cellfun(@(x)size(x, 1), p) );
    assert( length(np) == 1 );

    % Create a empty structure
    G = struct;

    % Cell number in layered direction
    nz = length(p)-1;

    % Cell number and face number of the 2D grid
    nc_g = g.cells.num;
    nf_g = g.faces.num;

    % G.cells -------------------------------------------------------------
    % num --------------------------------
    G.cells.num = nc_g * nz;

    % facePos ----------------------------
    ncf_g = diff(g.cells.facePos);
    ncf = repmat(ncf_g + 2, 1, nz);
    ncf = ncf(:);
    G.cells.facePos = cumsum([1; ncf]);

    % faces -------------------------------
    cf_g = arrayfunUniOut(@(c)g.cells.faces(g.cells.facePos(c):...
        g.cells.facePos(c+1)-1, 1),  (1:nc_g)');
    [cf, dir] = deal( cell(nc_g, nz) );
    for k = 1 : nz
        fXY = cellfunUniOut(@(x)x + (k-1) * nf_g, cf_g);
        fZ  = arrayfunUniOut(@(x)nf_g*nz + [(k-1)*nc_g + x; k*nc_g + x], ...
            (1:nc_g)');
        cf(:, k)  = cellfun(@(x ,y)[x;y], fXY, fZ, 'UniformOutput', false);
        dir(:, k) = arrayfunUniOut(@(x)[ones(x,1);5;6], cellfun(@length, fXY));
    end
    cf = cf(:);
    cf = cell2mat(cf);
    dir = dir(:);
    dir = cell2mat(dir);
    G.cells.faces = [cf, dir];

    % layers -----------------------------
    layers  = repmat((1:nz), nc_g, 1);
    G.cells.layers = layers(:);

    % G.faces -------------------------------------------------------------
    % num --------------------------------
    G.faces.num = nf_g * nz + (nz+1) * nc_g;
    assert( G.faces.num == max(G.cells.faces(:,1)) )

    % layers and surfaces ----------------
    G.faces.surfaces = zeros(G.faces.num, 1);
    surfaces  = repmat((1 : nz+1), nc_g, 1);
    G.faces.surfaces(nf_g * nz + 1 : end) = surfaces(:);
    nlayerF = nnz(G.faces.surfaces == 0);
    layerF = repmat((1:nz), nlayerF/nz, 1);
    layerF = [layerF(:); zeros(nnz(G.faces.surfaces>0),1)];
    assert(numel(layerF) == G.faces.num);
    G.faces.layers = layerF;
    
    % nodePos ----------------------------
    nfn1 = 4 * ones(nf_g * nz, 1);
    if isempty(opt.connectivity)
        cn_g = getNodesOfcell(g);
    else
        cn_g = opt.connectivity;
    end
    if isnumeric(cn_g)
        cn_g = mat2cell(cn_g, ones(size(cn_g, 1), 1));
    end

    nfn2 = cellfun(@length, cn_g);
    nfn2 = repmat(nfn2, 1, nz+1);
    nfn2 = nfn2(:);
    nfn  = [nfn1; nfn2];
    G.faces.nodePos = cumsum([1; nfn]);

    % nodes ------------------------------
    fn_g = arrayfunUniOut(@(f)g.faces.nodes(g.faces.nodePos(f):...
        g.faces.nodePos(f+1)-1), (1:nf_g)');

    fn1 = cell(nf_g, nz);
    for k = 1 : nz
        fn1(:, k) = cellfunUniOut(...
            @(x) [x([1,2]) + (k-1)*np; x([2,1]) + k*np], fn_g);
    end
    fn1 = fn1(:);
    fn1 = cell2mat(fn1);

    fn2 = cell(nc_g, nz+1);
    for k = 1 : nz + 1
        fn2(:, k) = cellfunUniOut(@(x) x + (k-1)*np, cn_g);
    end
    fn2 = fn2(:);
    fn2 = cellfunUniOut(@convertToColumn, fn2);
    fn2 = cell2mat(fn2);

    G.faces.nodes = [fn1; fn2];

    % G.nodes -------------------------------------------------------------
    G.nodes.num = np * (nz+1);
    G.nodes.coords = cell2mat(p);

    % G.type --------------------------------------------------------------
    G.type = [g.type, { mfilename }];

    % G.griddim -----------------------------------------------------------
    G.griddim    = 3;
    G.layers.num = nz;
    
    % G.surfGrid ----------------------------------------------------------
    G.surfGrid = g;
    
    % G.faces.neighbors ---------------------------------------------------
    G = computeGeometry(G, 'findNeighbors', true);
end

function cn = getNodesOfcell(G)
    [cn, pos] = gridCellNodes(G, (1:G.cells.num)');
    cn = arrayfunUniOut(@(c)cn(pos(c):pos(c+1)-1), (1:G.cells.num)');
    cn = sortPtsCounterClockWise(G.nodes.coords(:, [1,2]), cn);
end