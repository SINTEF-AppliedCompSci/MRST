%% Compare grid-orientation effects for skew grids
% In this example we consider a symmetric well pattern on a skew grid to
% study the grid-orientation effects for the combined single-point upwind
% and TPFA/mimetic schemes. The script contains two different examples,
% both describing a 400 x 100 m^2 reservoir section. The first (exmpl=1)
% describes a horizontal reservoir cross-section in which water is injected
% from one well at the midpoint of the north perimeter and fluids are
% produced from two wells along the south perimeter, located 50 m from the
% SE and SW corners, respectively. The second case (exmpl=2) describes a
% vertical cross-section, where water is injected from two horizontal
% injectors at the bottom of the reservoir and fluids produced from a
% horizontal producer at the top of the reservoir. In both cases, the well
% operate under bottom-hole control.

mrstModule add incomp mimetic ad-core
gravity reset off;
exmpl = 1;        % type of cross-section: 1 - horizontal,  2 - vertical

% Final time and number of time steps
T = 1200*day;
nstep =120;

% Rectangular reservoir with a skew grid.
G = cartGrid([41,20],[2,1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 200*makeSkew(G.nodes.coords);
G.nodes.coords(:,2) = 100*G.nodes.coords(:,2);
G = computeGeometry(G);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv   = sum(poreVolume(G,rock));
hT   = computeTrans(G, rock);
IP   = computeMimeticIP(G, rock);

% Symmetric well pattern
wcells = findEnclosingCell(G,[200 97.5; 50 2.5; 350 2.5]);
rate   = sum(poreVolume(G,rock));
if exmpl==1
    W = addWell([], G, rock, wcells(1), 'Type', 'bhp', ...
        'Val', 200*barsa, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
    W = addWell(W, G, rock, wcells(2), 'Type', 'bhp', ...
        'Val', 100*barsa, 'name', 'P1', 'radius', .1, 'Comp_i', [0 1]);
    W = addWell(W, G, rock, wcells(3), 'Type', 'bhp', ...
        'Val', 100*barsa, 'name', 'P2', 'radius', .1, 'Comp_i', [0 1]);
else
    gravity(norm(gravity())*[0 -1 0]);
    W = addWell([], G, rock, wcells(1), 'Type', 'bhp', ...
        'Val', 100*barsa, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
    W = addWell(W, G, rock, wcells(2), 'Type', 'bhp', ...
        'Val', 200*barsa, 'name', 'I1', 'radius', .1, 'Comp_i', [1 0]);
    W = addWell(W, G, rock, wcells(3), 'Type', 'bhp', ...
        'Val', 200*barsa, 'name', 'I2', 'radius', .1, 'Comp_i', [1 0]);
end

% Two-phase fluid
fluid = initSimpleFluid('mu', [1 10].*centi*poise, ...
        'rho', [1000,  850].* kilogram/meter^3, 'n', [2, 2]);

%% Prepare plotting and allocate various objects
% Figure and figure settings
figure('Position',[400 460 900 350]);
pargs = {'EdgeColor','k','EdgeAlpha',.05};

[xtp,xmi] = deal(initState(G,W,100*barsa, [0 1]));

dt = repmat(T/nstep,1,nstep);
dt = [dt(1).*sort(repmat(2.^-[1:5 5],1,1)) dt(2:end)];
N  = numel(dt);
wellSols = cell(N,2);
t  = 0;
for n=1:N
    t = t + dt(n);
    
    xtp  = incompTPFA(xtp, G, hT, fluid, 'wells', W);
    xtp  = implicitTransport(xtp, G, dt(n), rock, fluid, 'wells', W);
    wellSols{n,1} = getWellSol(W, xtp, fluid);
    subplot(2,2,1), cla, 
       plotCellData(G, xtp.pressure/barsa, pargs{:});
       caxis([100 200]); title(sprintf('TPFA: %.0f days', t/day));
    subplot(2,2,3), cla,
       % plotGrid(G,'FaceColor','none');
       plotCellData(G, xtp.s(:,1),xtp.s(:,1)>.01,pargs{:});
       axis([0 400 0 100]); caxis([0 1]);

    xmi  = incompMimetic(xmi, G, IP, fluid, 'wells', W);
    xmi  = implicitTransport(xmi, G, dt(n), rock, fluid, 'wells', W);
    wellSols{n,2} = getWellSol(W, xmi, fluid);
    subplot(2,2,2), cla,
       plotCellData(G, xmi.pressure/barsa, pargs{:});
       caxis([100 200]); title(sprintf('Mimetic: %.0f days', t/day));
       colorbar('peer', gca, 'Position', [0.94 0.58 0.0175 0.345]);
    subplot(2,2,4), cla,
       % plotGrid(G,'FaceColor','none');
       plotCellData(G, xmi.s(:,1),xmi.s(:,1)>.01, pargs{:});
       axis([0 400 0 100]); caxis([0 1]);
       colorbar('peer', gca, 'Position', [0.94 0.105 0.0175 0.345]);
    
    for i=1:4,
        subplot(2,2,i), hold on
        plot(G.cells.centroids(wcells,1),G.cells.centroids(wcells,2),...
            'ok','MarkerSize',6,'MarkerFaceColor',[.6 .6 .6]);
        hold off
    end
    drawnow
end

plotWellSols({wellSols(:,1), wellSols(:,2)}, ...
    cumsum(dt), 'datasetnames', {'TPFA','MFD'});
%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
