function [bdNodes, bdCells] = demoGetBdyNodesCells(G, CI)
% An excerpt from 'VolumeOfInterest.getBoundaryInfoSingleSurface' to show
% how to get the sorted boundary nodes and cells (counterclockwise) of a
% inner continuous region what we name 'volume of interest (VOI)' specified
% by cells 'CI'.

    % 1. Build a local grid 'g' from nodes (connectivity list) of CI ------
    n = arrayfunUniOut(@(c)gridCellNodes(G, c), CI);
    n = sortPtsCounterClockWise(G.nodes.coords(:,1:2), n);
    assert(all(cellfun(@numel, n)==4))
    n = cell2mat(n);
    [nu, ~, ic] = unique(n);
    p = G.nodes.coords(nu, [1, 2]);
    t = reshape(ic, 4, [])';
    g = tessellationGrid(p, t);
    g = computeGeometry(g);
    
    % 2. Get boundary faces of g, sorted, counterclockwise ----------------
    N = g.faces.neighbors;
    bf = find( ~all(N,2) );
    bf = sortPtsCounterClockWise(g.faces.centroids, {bf});
    bf = bf{1};

    % 3. Get boundary nodes of g ------------------------------------------
    % Boundary nodes of bf
    [bfn, pos] = gridFaceNodes(g, bf);
    assert(all(diff(pos)==2))
    bfn = reshape(bfn, 2, [])';
    % Boundary nodes of g, sorted, counterclockwise
    bn = arrayfun(@(r)bfn(r, ~ismember(bfn(r,:), bfn(r-1,:))), ...
        (2:size(bfn,1)-1)');
    idx = ismember(bfn(1,:), bfn(2,:));
    bn  = [bfn(1,~idx); bfn(1,idx); bn];
    
    % 4. Get boundary nodes of G in VOI -----------------------------------
    bdNodes  = nu(bn); % sorted, counterclockwise
    
    % 5. Get boundary cells of g ------------------------------------------
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
    
    % 4. Get boundary cells of G in VOI -----------------------------------
    bdCells = CI(bc);
    if numel(bc) ~= unique(numel(bc))
        error(['Isolate boundary cells are detected, please redefine the', ...
            ' boundary polygon'])
    end
    
    % Plot the results
    figure, hold on, axis equal off
    plotGrid(G, 'facecolor', 'none')
    demoPlotPoly(G.nodes.coords(bdNodes,:), 'b^-', 'b', 4)
    demoPlotPoly(G.cells.centroids(bdCells,:), 'rs-', 'r', 4)
    demoPlotPoly(g.faces.centroids(bf,:), 'mp-', 'k', 4)
    legend('G', 'Boundary node polygon', 'Boundary cell polygon', ...
        'Boundary faces polygon')
end