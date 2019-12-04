function [G, t] = radCartHybridGrid(GC, CI, rW, rM, nR, pW)
% Build the hybrid grid by gluing the radial grid in the near-well region 
% to the Cartesian grid elsewhere in the reservoir
%
% SYNOPSIS:
%   G = radCartHybridGrid(GC, CI, rW, rM, nR, pW)
%
% PARAMETERS:
%  GC    -   The Cartesian grid structure
%  CI    -   Cells inside the well region
%  nR    -   Number of cells in radial direction
%  rW    -   The minimum radius (wellbore radius)
%  rM    -   The maximum radius
%  pW    -   The well point coordinates
%
% RETURNS:
%  G  -  Valid hybrid grid definition
%  t  -  Connectivity list of the hybrid grid
%
% EXAMPLE:
%   GC = cartGrid([20, 20], [200, 200]);
%   GC = computeGeometry(GC);
%   ij = gridLogicalIndices(GC);
%   idxI = ij{1} >= 10 & ij{1} <= 14 & ij{2} >= 10 & ij{2} <= 14;
%   CI = find(idxI); 
%   % Place the well at the region center
%   pCI = GC.cells.centroids(CI, :);
%   pW  = 0.5*[min(pCI(:,1)) + max(pCI(:,1)), min(pCI(:,2)) + max(pCI(:,2))];
%   [nR, rW, rM] = deal(10, 0.2, 16);
%   [G, t] = radCartHybridGrid(GC, CI, rW, rM, nR, pW);
%   figure, hold on; axis equal off, plotGrid(G)
%  
% SEE ALSO:
%  `tessellationGrid` `pebi` `buildRadialGrid` `glueRadCartGrids`

    % Get the sorted boundary nodes of the region (in counterclockwise)
    bn = extractBdyNodesCells(GC, CI, 'plotResults', false);
    
    % The angular dimension and grid angles are determined by the boundary
    % nodes to conform with the Cartesian grid
    nA = numel(bn);
    
    % Compute the angles
    pbn = GC.nodes.coords(bn, :);
    pbn0 = bsxfun(@minus, pbn, pW);
    th = cart2pol(pbn0(:,1), pbn0(:,2));
    
    % Get the grid radii
    r = logspace(log10(rW), log10(rM), nR+1);
    
    % Get the radial grid points
    [R, TH] = meshgrid(r, th);
    [px, py] = pol2cart(TH(:), R(:));
    pR = bsxfun(@plus, [px(:), py(:)], pW);
    
    % The boundary points are the outermost angular points. The total 
    % radial dimension is then nR+1
    pR = [pR; pbn];
    
    % Build the radial grid
    [GR, tR] = buildRadialGrid(pR, nA, nR+1);

    % Points and connectivity list of the radial grid
    pR = GR.nodes.coords;

    % Points and connectivity list of the Cartesian grid
    % Remove the cells inside the well region first
    [GC_Rem, ~, ~, mapn] = removeCells(GC, CI);
    pC = GC_Rem.nodes.coords;

    % Get the indices of boundary nodes in GC_Rem
    bn = arrayfun(@(n)find(mapn == n), bn);

    % Merge the common nodes (boundary nodes). The boundary node indices in 
    % pC are replaced by the ones in pR
    nNo  = (1:size(pC,1))';
    % The non-boundary nodes indices
    idx  = ~ismember(nNo, bn);
    nNo(idx) = (1:nnz(idx))' + size(pR,1);
    pC = pC(idx,:);
    % We already know the boundary nodes in pR: the last nA nodes
    nNo(bn)  = size(pR,1)+1 - (nA:-1:1)';

    % Map the connectivity list of GC_Rem
    [cnC, pos] = gridCellNodes(GC_Rem, (1:GC_Rem.cells.num));
    cnC = nNo(cnC);
    tC  = arrayfunUniOut(@(c)cnC(pos(c):pos(c+1)-1), (1:GC_Rem.cells.num)');

    % Assemble the points and connectivity lists
    p = [pR; pC];
    
    % Sort the tC to same direction with tR (in clockwise)
    tC = sortPtsClockWise(p, tC);
    t = [tR; tC];
    
    % Build the hybrid grid
    G = tessellationGrid(p, t);
    G = computeGeometry(G);
    G.subGrids = {GR; GC_Rem};
end