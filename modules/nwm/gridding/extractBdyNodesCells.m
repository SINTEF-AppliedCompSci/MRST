function [bdNodes, bdCells] = extractBdyNodesCells(G, cI, varargin)
% An 2D version of 'VolumeOfInterest.getBoundaryInfoSingleSurface', which
% extracts the sorted boundary nodes and cells (in counterclockwise) of a
% inner continuous region what we name 'region of interest (ROI)' specified
% by cells 'cI'.
%
% SYNOPSIS:
%   [bdNodes, bdCells] = extractBdyNodesCells(G, cI)
%   [bdNodes, bdCells] = extractBdyNodesCells(G, cI, 'plotResults', false)
%
% PARAMETERS:
%  G     -   The 2D Cartesian grid
%  cI    -   Cells that specifies the ROI
%
% RETURNS:
%  bdNodes - The sorted boundary nodes, in counterclockwise
%  bdCells - The sorted boundary cells, in counterclockwise
%
% EXAMPLE:
%  G = cartGrid([25, 25], [200, 200]);
%  G = computeGeometry(G);
%  pbdy = [136,  150;
%          145,   95;
%           90,   30;
%           50,   50;
%           45,  105;
%           90,  160];
%  cCtro = G.cells.centroids;
%  in = inpolygon(cCtro(:,1), cCtro(:,2), pbdy(:,1), pbdy(:,2));
%  cI = find( in );
%  [bnv, bcv] = extractBdyNodesCells(G, cI);
% 
% SEE ALSO:
%  `VolumeOfInterest`

    opt = struct('plotResults', true);
    opt = merge_options(opt, varargin{:});

    % Build a local grid 'g' from nodes (connectivity list) of cI
    n = arrayfunUniOut(@(c)gridCellNodes(G, c), cI);
    n = sortPtsCounterClockWise(G.nodes.coords(:,1:2), n);
    assert(all(cellfun(@numel, n)==4))
    n = cell2mat(n);
    [nu, ~, ic] = unique(n);
    p = G.nodes.coords(nu, [1, 2]);
    t = reshape(ic, 4, [])';
    g = tessellationGrid(p, t);
    g = computeGeometry(g);
    
    % Get boundary faces of g, sorted, counterclockwise
    N = g.faces.neighbors;
    bf = find( ~all(N,2) );
    bf = sortPtsCounterClockWise(g.faces.centroids, {bf});
    bf = bf{1};

    % Get boundary nodes of g
    % Boundary nodes of bf
    [bfn, pos] = gridFaceNodes(g, bf);
    assert(all(diff(pos)==2))
    bfn = reshape(bfn, 2, [])';
    % Boundary nodes of g, sorted, counterclockwise
    bn = arrayfun(@(r)bfn(r, ~ismember(bfn(r,:), bfn(r-1,:))), ...
        (2:size(bfn,1)-1)');
    idx = ismember(bfn(1,:), bfn(2,:));
    bn  = [bfn(1,~idx); bfn(1,idx); bn];
    
    % Get boundary nodes of G in ROI
    bdNodes  = nu(bn); % sorted, counterclockwise
    
    % Get boundary cells of g
    bc = sum(N(bf, :),2);         % Sorted, counterclockwise
    bc = unique(bc, 'stable');    % Some cells appear twice
    % Insert the 'Z' cells
    N0 = [bc, [bc(2:end); bc(1)]];
    zc = zeros(size(bc));
    for ii = 1 : size(N0,1)
        c1 = N(any(N == N0(ii,1),2), :);
        c1 = unique(c1); 
        c1 = c1(c1~=0 & c1~=N0(ii,1));
        c2 = N(any(N == N0(ii,2),2), :);
        c2 = unique(c2); 
        c2 = c2(c2~=0 & c2~=N0(ii,2));
        if ~isempty(intersect(c1, c2))
            zc(ii) = intersect(c1, c2);
        end
    end
    bc = [bc, zc]';
    bc = bc(:);
    bc = bc(bc~=0);
    
    % Get boundary cells of G in ROI
    bdCells = cI(bc);
    if numel(bc) ~= unique(numel(bc))
        error(['Isolate boundary cells are detected, please redefine the', ...
            ' boundary polygon'])
    end
    
    if opt.plotResults
        % Plot the results
        figure, hold on, axis equal off
        plotGrid(G, 'facecolor', 'none')
        demoPlotPoly(G.nodes.coords(bdNodes,:), 'b^-', 'b', 4)
        demoPlotPoly(G.cells.centroids(bdCells,:), 'rs-', 'r', 4)
        demoPlotPoly(g.faces.centroids(bf,:), 'mp-', 'k', 4)
        legend('GC', 'Boundary-node polygon', 'Boundary-cell polygon', ...
            'Boundary-face polygon')
    end
end