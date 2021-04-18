function plotSingleCellInfo(G, c, varargin)
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

opt = struct('plotFaceInfo', []);
[opt, rest] = merge_options(opt, varargin{:});

fpos = G.cells.facePos([c,(c+1)]);
fix  = G.cells.faces(fpos(1):(fpos(2)-1));

plotInfo = opt.plotFaceInfo;
if ~isempty(plotInfo)
    if ~islogical(plotInfo)
        tmp = false(size(fix));
        tmp(plotInfo) = true;
        plotInfo = tmp;
    end
else
    plotInfo = true(size(fix));
end
    
ax = gca;
if ~ishold(ax)
    hold(ax, 'on');
end

for k = 1:numel(fix)
    plotSingleFaceInfo(G, fix(k), 'fromCell', c, 'plotInfo', plotInfo(k), rest{:})
end
    
tp = get(ax.Children, 'Type');
ix = strcmp(tp, 'text');
uistack(ax.Children(ix), 'top');
end
