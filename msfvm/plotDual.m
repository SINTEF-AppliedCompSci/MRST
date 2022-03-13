% Plot an implicitly defined dual grid for the multiscale finite volume method
%
% SYNOPSIS:
%   plotDual(G, dual)
%
% PARAMETERS:
%   G    - Grid structure
%   dual - Dual grid as defined by for example partitionUIdual.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


function plotDual(G, dual)

RGB = @(r,g,b,i) [r/(i*255) g/(i*255) b/(i*255)];
plotGrid(G, dual.nn, 'FaceColor', RGB(145,255,124,1), 'EdgeColor', 'k');

if ~isempty(dual.lineedge)
    plotGrid(G, dual.lineedge, 'FaceColor', RGB(0,0,255,1), 'EdgeColor', 'k', 'facea', .3, 'edgealpha',0);
else
    plotGrid(G, dual.ee(~ismember(dual.ee, dual.lineedge)), 'FaceColor', RGB(255,0,0,1), 'EdgeColor', 'k', 'facea', .3, 'edgealpha',0);
end
plotGrid(G, 'facea',0,'edgea',.1);
axis tight
end
