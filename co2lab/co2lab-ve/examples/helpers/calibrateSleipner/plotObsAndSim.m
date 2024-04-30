function plotObsAndSim(Gt, obs, sim, varargin)
% Plot simulated and observed plume thicknesses
% obs and sim are cell arrays containing plume thicknesses (heights) "h"

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    assert(numel(sim) == 12) % 12 for 1999-2010
    assert(numel(obs) == 12)

    opt.edgealpha = 0.1;
    opt.edgecolor = 'none';
    opt = merge_options(opt, varargin{:});
    
    myplotGrid = @(G) plotGrid(G, 'facecolor','none','edgealpha',opt.edgealpha);
    myplotCellData = @(G,data,tol) plotCellData(G, data, data>tol, ...
        'edgealpha',opt.edgealpha, 'edgecolor',opt.edgecolor);

    %ty=[3,6,8,12];
    Tol = 0.01;

    figure,
    subplot(2,4,1)
    title({'Year: 2001';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{3}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,5)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{3}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,2)
    title({'Year: 2004';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{6}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,6)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{6}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,3)
    title({'Year: 2006';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{8}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,7)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{8}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,4)
    title({'Year: 2010';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{12}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,8)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{12}.h, Tol); colormap(jet); axis equal tight off; colorbar
end
