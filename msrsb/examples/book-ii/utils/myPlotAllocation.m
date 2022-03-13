function amax=myPlotAllocation(WP, WPc, names, amax)
%Plot well-allocation factors for two models of different solutions
% SYNOPSIS
%   amax=plotWellAllocationComparision(WP1, WP2, names, amax)
%
% PARAMETERS:
%   WP1, WP2 - data structure containing information about well pairs,
%              computed by a call to 'computeWellPairs'
%   names    - names of the wells
%   amax     - maximum allocation contained in these data
%
% See plotWellAllocationComparison for more information

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

colormap(.6*jet+.4);
cwp = cumsum(WPc.inj(2).alloc,1); amax = max([sum(cwp,2); amax]);
barh(WPc.inj(2).z, cwp, 'stacked', 'BarWidth', .98, 'EdgeColor','none');
hold on
cwp = cumsum(WP.inj(2).alloc,1);  amax = max([sum(cwp,2); amax]);
if size(cwp,1)<20
    barh(WP.inj(2).z, cwp, 'stacked','BarWidth', .98, 'FaceColor','none');
else
    c = cumsum(cwp,2);
    z =  WP.inj(2).z; dz = min(diff(z));
    n = size(c,2);
    stairs([zeros(1,n); c], repmat(z([1:end end]),1,n), '-k','LineWidth',1);
end
hold off, axis tight
lh=legend(names,'Location','SouthEast');
set(lh,'units','pixels','FontSize',8);
end
