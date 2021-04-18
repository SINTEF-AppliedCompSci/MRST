function plotSingleFaceInfo(G, f, varargin)
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

opt = struct('rRel',     .1, ...
    'fromCell', [], ...
    'expandFac', 1, ...
    'plotInfo', true);

%G = computeGeometry(G);
[opt, extra] = merge_options(opt, varargin{:});

cf = G.faces.centroids(f,:);
if ~isempty(opt.fromCell)
    cc       = G.cells.centroids(opt.fromCell,:);
    cfVec    = cf-cc;
    cellView = true;
    assert(ismember(opt.fromCell, G.faces.neighbors(f,:)), ...
        'Cell %d is not a neighbor of face %d ', opt.fromCell, f);
else
    cfVec = zeros(1,3);
    cellView = false;
end
shiftVec = (opt.expandFac-1)*cfVec;


npos = G.faces.nodePos([f,(f+1)]);
nix  = G.faces.nodes(npos(1):(npos(2)-1));
seg  = [nix, [nix(2:end);nix(1)]];

coord = G.nodes.coords(nix,:);

[~, ia, ib] = uniquetol(coord, 'ByRows', true);

if numel(ia) ~= numel(ib)
    np = accumarray(ib, ones(size(ib)));
    fprintf('Face %d contains %d zero length segments\n', f, sum(np-1));
end

np   = numel(nix);

if ~cellView
    dv = bsxfun(@minus, coord, cf);
else
    dv = bsxfun(@minus, coord, cc);
end
d   = sqrt(sum(dv.^2,2));
dvu = bsxfun(@rdivide, dv, d);
r = opt.rRel*max(d);

ax = gca;
if ~ishold(ax)
    hold(ax, 'on');
end

% shift coordinates according to expVec
cf    = cf + shiftVec;
coord = bsxfun(@plus, coord, shiftVec);


patch('Faces', (1:np), 'Vertices', coord, 'FaceColor', 'y', 'FaceAlpha', .5, extra{:});

if ~opt.plotInfo
    return
end
plot3(cf(:,1), cf(:,2), cf(:,3), 'xk', 'LineWidth', 2);
plot3(coord(:,1), coord(:,2), coord(:,3), 'ok', 'LineWidth', 2);

bccol = [1 1 1];
for k = 1:np
    loc  = cf + dv(k,:)+ r*dvu(k,:);
    info = sprintf('node: %d, loc no: %d', nix(k), k);
    text(loc(1), loc(2), loc(3), info, 'FontSize', 10, 'BackgroundColor', bccol);
    lc = [coord(k,:); .8*loc + .2*coord(k,:)];
    plot3(lc(:,1), lc(:,2), lc(:,3), '--b');
end

neig = G.faces.neighbors(f,:);
% scale normal for plotting
n  = G.faces.normals(f,:);
n  = n/norm(n);
% check inner product of normal and cell-cell centroid
if nnz(neig) == 2
    cc = diff(G.cells.centroids(neig,:));
else
    cc = G.faces.centroids(f,:) - G.cells.centroids(neig(neig>0),:);
    cc = cc*(1-2*(neig(1)==0));
end
cc = cc/norm(cc);
if ~n*cc' > 0
    fail = true
else
    fail = false;
end

    
n  = n*max(d)/2;
sh = null(n);
sh = r*sh(:,1)';
for k = 1:2
    c     = neig(k);
    other = neig(mod(k,2)+1);
    if c ~= 0
        fpos  = G.cells.facePos([c,(c+1)]);
        fix   = G.cells.faces(fpos(1):(fpos(2)-1));
        order = find(fix==f);
    else
        order = nan;
    end
    
    if ~cellView || neig(k) == opt.fromCell
        if k == 1
            col = 'g';
            sgn = 1;
        else
            col = 'r';
            sgn = -1;
        end
        ns = sgn*n;
        quiver3(cf(1), cf(2), cf(3), ns(1), ns(2), ns(3), col, 'LineWidth', 2)
        ns = ns+sh;
        loc = cf + ns;
        info = sprintf('face: %d, loc no: %d\n cell: %d (%d), sign: %d', ...
            f,        order,      c, other,     sgn);
        text(loc(1), loc(2), loc(3), info, 'FontSize', 10, 'BackgroundColor', bccol);
    end
end
end
