function showWellCommunicationPercentage(d, ax, val)
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

val = val([1:end end],:);
val = double(val(:,[1:end end]));
% scale according to injectors
% val = bsxfun(@rdivide, val, sum(val, 2));
val(val==0) = nan;

pcolor(ax, val);
inames= arrayfun(@(x) x.label.String, d.WellPlot.injectors,'UniformOutput',false);
pnames= arrayfun(@(x) x.label.String, d.WellPlot.producers,'UniformOutput',false);
set(ax, 'XTick', 1.5:numel(pnames)+.5, ...
    'XTickLabel', pnames, 'YTick',1.5:numel(inames)+.5, ...
    'YTickLabel', inames, 'XTickLabelRotation',45, ...
    'XAxisLocation','top','Fontsize',8,'YDir','reverse');
axis(ax,'on','tight');
maxVal = max(max(val));
if maxVal == 0 || ~isfinite(maxVal)
    maxVal = 1;
end
caxis(ax, [0 maxVal]);
c = colorbar(ax,'southoutside');
c.TickLabels = (0:10:100);
end
