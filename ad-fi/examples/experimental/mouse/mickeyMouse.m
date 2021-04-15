%% Read data
a = imread('mickey3.png');  a = sum(a(end:-1:1,:,:),3)';
k = imread('mickey2.jpeg'); K = sum(k(end:-1:1,:,:),3)'./255;
assert(all(size(a)==size(K)))

%% Construct grid
nlayers = 1;
G = cartGrid(size(a),size(a)*10*meter);
g = extractSubgrid(G, find(a<3*125));
G = makeLayeredGrid(g, nlayers);
G.nodes.coords(:,3) = 10*G.nodes.coords(:,3);
G = computeGeometry(G);
%% Construct petrophysical parameters (Carman-Kozeny relationship)
KG = K(g.cells.indexMap);
p = 0.35-.1*KG+.05;
k = p.^3.*(1.75e-5)^2./(0.81*72*(1-p).^2);

rock.poro = p;
rock.perm = k;
%% Visualize
clf, plotCellData(G, log10(k),'EdgeColor','k','EdgeAlpha',.1);
h=colorbar; caxis(log10([5 1000]*milli*darcy))
set(h,'XTick',1,'XTickLabel','mD');
set(h,'YTick',log10([10 100 1000].*milli*darcy),'YTickLabel',[10 100 1000]);
colormap(flipud(jet));

%%
findpt = @(pt) min(sum(abs(G.cells.centroids - repmat(pt, G.cells.num, 1)), 2));
[v, prod] = findpt([950, 500, 0]);
[v, inj1] = findpt([500, 1200, 0]);
[v, inj2] = findpt([1350, 1200, 0]);

plotGrid(G, [prod, inj1, inj2]);

%%
mrstVerbose true
mrstModule add deckformat ad-core ad-fi

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'mouse.data');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);

schedule = deck.SCHEDULE;

W = addWell([], G, rock, inj1, 'Type', 'rate', 'Val', 1, ...
            'InnerProduct', 'ip_tpf');
W = addWell(W, G, rock, inj2, 'Type', 'rate', 'Val', 1, ...
            'InnerProduct', 'ip_tpf');
W = addWell(W, G, rock, prod, 'Type', 'bhp', 'Val', 1, 'Sign', -1, ...
            'InnerProduct', 'ip_tpf');

W(1).poly = .01;
W(2).poly = .01;
%%


systemPolymer = initADISystem(deck, G, rock, fluid);
% systemPolymer = initADISystem({'Oil', 'Water'}, G, rock, fluid);
state0 = initResSol(G, 1, [.5, .5]);
state0.cmax = zeros(G.cells.num,1);
state0.c = state0.cmax;
state0.wellSol = initWellSolLocal(W, state0);


nt = 250;
dt = ones(nt,1)*day/10;
states = cell(nt,1);
state = state0;
for tstep = 1:nt
     state = solvefiADI(state, dt(tstep), W, G, systemPolymer);
     states{tstep} = state;
end

%% Plot oil
figure(1)
for i = 1:nt
    clf;
    plotCellData(G, states{i}.s(:,2));
    colorbar
    axis tight off
    pause(.5)
end
%% Plot polymer
h = figure(1);

for i = 1:nt
    clf;
    plotCellData(G, states{i}.c);
    colorbar
    axis tight off
    caxis([0 0.01])
    title(['Polymer concentration at timestep ' num2str(i)])
    print('-dpng', [num2str(i) '.png'])
%     pause(.1)
end

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
