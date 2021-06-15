function [] = fracFlowContours(G, W, states, fluid, frac, varargin)
% utility-function for analyseModel2D.m
% plot contour of each state where fracFlow = frac
% only for 2D rectangular grids and oil/water systems!

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
assert(or(G.griddim == 2, G.cartDims(3)==1));
assert(prod(G.cartDims)==G.cells.num);

mn = min(G.nodes.coords(:,1:2));
mx = max(G.nodes.coords(:,1:2));
dd = (mx-mn)./G.cartDims(1:2);
x = (mn(1):dd(1):mx(1))';
x = .5*(x(1:end-1)+x(2:end));
y = (mn(2):dd(2):mx(2))';
y = .5*(y(1:end-1)+y(2:end));

ns = numel(states);

for k = 1:ns
    set(gca,'FontSize', 14)
    hold on
    plotCellData(G, states{k}.s(:,2))

    sw  = states{k}.s(:,1);
    krw = fluid.krW(sw);
    kro = fluid.krO(1-sw);

    muw = fluid.muW(states{k}.pressure);
    muo = fluid.muO(states{k}.pressure);
    bw = fluid.bW(200*barsa);   % FVF at reference pressure
    bo = fluid.bO(200*barsa); 
    fw = (bw*krw./muw)./(bw*krw./muw+bo*kro./muo);
    m = reshape(fw/frac, G.cartDims(1:2))';
    [c,h] = contour(x,y,m, [1 1], varargin{:});
    caxis([0 1])
    axis off, axis equal
    view([-1 -2 1.5])
    camproj perspective
    axis tight
    plotWell(G, W, 'Fontsize', 14, 'Color', 'k')
end
end
