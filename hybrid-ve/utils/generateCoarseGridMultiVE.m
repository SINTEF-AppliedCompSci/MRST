function CG = generateCoarseGridMultiVE(G, VEGroups, varargin)
%Create a coarse hybrid VE grid from specified VE regions
%
% SYNOPSIS:
%   CG = generateCoarseGridMultiVE(G, VEGroups
%
% REQUIRED PARAMETERS:
%   G        - Grid structure
%
%   VEGroups - Indicator array with one entry per fine cell, containing the
%              index of the VE region where that cell belongs.
%
% RETURNS:
%   CG       - A modified MRST coarse grid with additional fields required
%              by the hybrid VE solvers.
%
% SEE ALSO:
%   convertToMultiVEModel

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

    opt = struct('useDepth', true);
    opt = merge_options(opt, varargin{:});
    [G.cells.topDepth, G.cells.bottomDepth, G.cells.height, G.cells.topCells] = getCellHeights(G, opt);
    
    require coarsegrid
    r = zeros(G.cells.num, 1);
    p = zeros(G.cells.num, 1);
    
    isFine = VEGroups == 0;
    r(isFine) = 1;
    
    [ii, jj, ~] = gridLogicalIndices(G);
    for i = 1:max(VEGroups)
        g = VEGroups == i;
        r(g) = i + 1;
        p(g) = max(p) + ii(g) + (jj(g) - 1).*G.cartDims(1) + 1;
    end
    p(isFine) = max(p) + (1:nnz(isFine))';
    p = processPartition(G, p);
    p = compressPartition(p);
    
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
    
    coarse = accumarray(p, r, [], @max);
    
    % Coarse cells, discretization type
    CG.cells.discretization = coarse;
    % Fine cells, discretization type
    CG.parent.cells.discretization = r;
    
    % Lateral connections between VE types should be set to zero
    
    % Vertical connections between types should be flagged
    
    % Define cells.height
    % Define cells.topDepth
    % Define cells.bottomDepth
    [CG.cells.topDepth, CG.cells.bottomDepth, CG.cells.height] = deal(zeros(CG.cells.num, 1));
    if opt.useDepth
        CG.cells.height = accumarray(CG.partition, CG.parent.cells.height, [CG.cells.num, 1], @sum);
        CG.cells.topDepth = accumarray(CG.partition, CG.parent.cells.topDepth, [CG.cells.num, 1], @min);
        CG.cells.bottomDepth = accumarray(CG.partition, CG.parent.cells.bottomDepth, [CG.cells.num, 1], @max);
    else
        for i = 1:CG.cells.num
            cells = CG.partition == i;
            CG.cells.height(i) = sum(CG.parent.cells.height(cells));
            % Alternate mode: We define depths from the top as metric
            % distance along the piecewise linear curve through the faces
            cells = find(cells);
            % Top and bottom is parametric distance from top
            td = min(CG.parent.cells.topDepth(cells));
            CG.cells.topDepth(i) = td;
            CG.cells.bottomDepth(i) = td + CG.cells.height(i);
            % Then fix top and bottom of cells as the cumulative sum along
            % the parametrization of the column
            d = CG.parent.cells.topDepth(cells);
            [~, sortIx] = sort(d);
            cells = cells(sortIx);
            b = cumsum(CG.parent.cells.height(cells));
            t = [0; b(1:end-1)];

            CG.parent.cells.topDepth(cells) = t + td;
            CG.parent.cells.bottomDepth(cells) = b + td;
            % Sort by depth
            % Assign new top and bottoms from this...
        end
    end
end

function [topz, bottomz, height, topcells] = getCellHeights(G, opt)
    faceNo = rldecode(1 : G.cells.num, ...
                    diff(G.cells.facePos), 2) .';
    
    hf_pt = G.faces.centroids(G.cells.faces(:, 1), :);
    
    hf_z = hf_pt(:, 3);
    hf_top = hf_z;
    hf_bottom = hf_z;
    if size(G.cells.faces, 2) > 1 && (any(strcmpi(G.type, 'cartGrid')) || any(strcmpi(G.type, 'processGRDECL')))
        ftype = G.cells.faces(:, 2);
    else
        fno = G.cells.faces(:, 1);
        ftype = ones(size(fno));
        cno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
        normals = abs(G.faces.normals(fno, :));
        nz = normals(:, 3);
        nxy = max(normals(:, 1), normals(:, 2));
        zc = G.cells.centroids(cno, 3);
        z = G.faces.centroids(fno, 3);
        
        ftype(nz > nxy & z < zc) = 5;
        ftype(nz > nxy & z > zc) = 6;
        warning('Grid is not corner point or cartGrid derived. Unable to guess column structure. Results may be inaccurate!');
    end
    % Avoid picking "wrong" type of face for columns
    maskTop    = ftype == 5;
    maskBottom = ftype == 6;

    hf_top(~maskTop) = inf;
    hf_bottom(~maskBottom) = -inf;

    topIx = accumarray(faceNo, hf_top, [G.cells.num, 1],  @minIndex);
    bottomIx = accumarray(faceNo, hf_bottom, [G.cells.num, 1],  @maxIndex);
    
    % We have local indexes, increment with the number of faces
    offsets = [0; cumsum(diff(G.cells.facePos(1:end-1)))];
    top = hf_pt(topIx + offsets, :);
    topSubs = topIx + offsets;
    topcells = G.faces.neighbors(G.cells.faces(topSubs, 1), :);
    cells = (1:G.cells.num)';
    for i = 1:2
        topcells(topcells(:, i) == cells, i) = 0;
    end
    topcells = sum(topcells, 2);
    topcells(topcells == 0) = cells(topcells == 0);
    
    bottom = hf_pt(bottomIx + offsets, :);
    topz = top(:, 3);
    bottomz = bottom(:, 3);
    if opt.useDepth
        height = bottomz - topz;
    else
        height = sqrt(sum((top - bottom).^2, 2));
    end
    assert(all(height >= 0))
end

function ix = minIndex(x)
    [~, ix] = min(x);
end

function ix = maxIndex(x)
    [~, ix] = max(x);
end