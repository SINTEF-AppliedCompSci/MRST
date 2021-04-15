function [] = plotSchedules(schedule, varargin)
% utility-function for analyseModel2D.m
% plot schedule, possibliy with gradient info
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
opt = struct('boxConst'    , [], ...
             'grad'         , [], ...
             'scaleGrad'    , true, ...
             'singlePlot'   , true);

opt = merge_options(opt, varargin{:});

boxConst    = opt.boxConst;
grad        = cell2mat(opt.grad);
doScale     = opt.scaleGrad;
singlePlot  = opt.singlePlot;
expandY     = .05;
numWells = numel(schedule.control(1).W);
if isempty(boxConst)
    boxConst = zeros(numWells, 0);
end

% control-times and midpoints
dt_c = accumarray(schedule.step.control, schedule.step.val);
t_c  = [0; cumsum(dt_c)/day]; 
mid  = t_c(1:end-1) + diff(t_c)/2;

for k =1:numWells
    if ~singlePlot
        figure(k);
    else
        subplot(numWells, 1, k);
    end
    hold on
    vals = arrayfun(@(x)x.W(k).val, schedule.control(:));
    tp   = schedule.control(1).W(k).type;
    bx   = boxConst(k,:);
    if strcmp(tp, 'bhp')
        vals = vals/barsa;
        bx   = bx/barsa;
    else
        vals = vals*day;
        bx   = bx*day;
    end
    stairs(t_c, [vals; vals(end)], 'LineWidth', 2)
    if ~isempty(bx)
        line([0,t_c(end)], [bx(1), bx(1)], 'color', 'r');
        line([0,t_c(end)], [bx(2), bx(2)], 'color', 'r');
    end
    set(gca, 'XTick', t_c);
    set(gca, 'XGrid', 'on')
    if ~isempty(grad)
        gk = grad(k,:);
        if doScale
            ax = axis(gca);
            dy = ax(4)-ax(3);
            gk = .25*gk*dy/max(abs(gk));
        end
        plotVertVecs(mid, vals, gk);
    end
    if expandY > 0
        axis tight;
        ax = axis(gca);
        dy = ax(4)-ax(3);
        ax = [ax(1) ax(2), ax(3)-dy*expandY, ax(4)+dy*expandY];
        axis(ax);
    end
    box on
    title( schedule.control(1).W(k).name)
    ylabel(tp)
end
end

function plotVertVecs(x, y, dy)
lw = 2; c = 'g'; mz = 10;
x = x(:); y = y(:); dy = dy(:); 
y2 = y+dy;
for k = 1:numel(x)
    if dy(k)>0, mk = '^'; else mk = 'v'; end 
    line([x(k) x(k)], [y(k) y2(k)], 'Color', c, 'LineWidth', lw);
    line(x(k),y2(k), 'Marker', mk, 'MarkerFaceColor', c, ...
                     'MarkerEdgeColor', c, 'MarkerSize', mz);
end
end
    
        
        
        
    
    