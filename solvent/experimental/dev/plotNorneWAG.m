mrstModule add ad-core ad-eor ad-blackoil ad-props sequential matlab_bgl

gravity reset on

%%

pth = mrstPath('mrst-solvent');

load([pth, '/code/FourPhaseSolvent/examples/data/norne.mat']);

G = computeGeometry(G);

mrstModule add deckformat


%%


inj = [9, 15; ...
        26, 15; ...
       36, 80; ...
       10, 85; ...
       24, 30; ...
       14, 52; ...
       18, 80;
       23, 66];

nInj = size(inj,1);
W = [];
T = 5*year;
pv = 0.5*sum(poreVolume(G, rock));
rate = (pv/T)/nInj;
for i = 1:nInj
    W = verticalWell(W, G, rock, inj(i, 1), inj(i, 2), [], ...
        'comp_i', [0,0,0,1],...
        'type', 'rate', 'val', rate);
end
   
prod = [10, 66; ...
        12, 32; ...
        22, 49; ...
        13, 91; ...
        37, 95;
        35, 64];
nProd = size(prod,1);
for i = 1:nProd
    W = verticalWell(W, G, rock, prod(i, 1), prod(i, 2), [], ...
        'comp_i', [1,0,0,0], ...
        'type', 'bhp', 'val', 80*barsa);
end

%%

step = 500;

hs = ResultHandler('dataDirectory', '/media/strene/806AB4786AB46C92/mrst-solvent/norne/', ...
    'dataFolder', 'nStep2000_nCycl5', ...
    'dataprefix', 'state_step', 'cleardir', false);

n = 1625;
states = cell(n,1);
for i = 1:n
    states{i} = hs{i};
    
end

%%

mrstModule add mrst-gui
figure(1);
plotToolbar(G, states);

%%

pth = [mrstPath('mrst-solvent'), '/presentation/figures/norne/'];
savepng = @(name) print([pth, name], '-dpng', '-r300');
% savepng = @(name) [];
saveeps = @(name) print([pth, name], '-depsc');

close all

for i = 9:numel(W)
    W(i).sign = -1;
end

%%

close all

df = get(0, 'defaultfigureposition');
figure('position', [df(1:2), [1000, 800]]);

plotCellData(G, log10(rock.perm(:,1)));
ax = gca;
pos = ax.Position;
logColorbar('location', 'southoutside', 'position', [pos(1), pos(2)+0.15, 0.95*pos(3), 0.05*pos(4)]);
colormap(jet)
view([90,50]);
camlight(-20,0)
axis off

hold on
plotWell(G, W, 'fontsize', 0, 'color2', 'b')

savepng('perm');


df = get(0, 'defaultfigureposition');
figure('position', [df(1:2), [1000, 800]]);

plotCellData(G, rock.poro);
ax = gca;
pos = ax.Position;
colorbar('location', 'southoutside', 'position', [pos(1), pos(2)+0.15, 0.95*pos(3), 0.05*pos(4)]);
colormap(jet)
view([90,50]);
camlight(-20,0)
axis off

hold on
plotWell(G, W, 'fontsize', 0, 'color2', 'b')

savepng('poro');


%%

close all
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
azel = [90,75];
tol = 1e-2;

subplot(1,2,1)

wc = vertcat(W.cells);
M = struct('cdata',[],'colormap',[]);

ns = 1625;
% ns = 4;

for i = 1:ns

    clf;
    
    subplot(1,2,1)

    
    c = states{i}.s(:,4) > tol;
    c(wc) = false;
    plotFaces(G, boundaryFaces(G), 'edgealpha', 0.1, 'facecolor', 'none');
    plotCellData(G, states{i}.s(c,4), c)
    plotWell(G, W, 'fontsize', 0, 'color2', 'b')
    caxis([0,1]);
    colormap(jet)
    view(azel);
    axis tight off
    ax = gca;
    pos = ax.Position;
    colorbar('location', 'southoutside', 'position', [pos(1), pos(2)+0.07, pos(3), 0.05*pos(4)]);
    title('S_s');
    ax.FontSize = 20;
    
    
    subplot(1,2,2)
    plotCellData(G, states{i}.s(:,2))
    plotWell(G, W, 'fontsize', 0, 'color2', 'b')
    caxis([0,1]);
    colormap(jet)
    view(azel);
    axis tight off
    ax = gca;
    pos = ax.Position;
    colorbar('location', 'southoutside', 'position', [pos(1), pos(2)+0.07, pos(3), 0.05*pos(4)]);
    title('S_o');
    ax.FontSize = 20;
    
    drawnow
    
    ax = gca;
    ax.Units = 'pixels';
    
    m = 100; n = 90;
    rect = [m, n, 2000-m, 1000-n+20];
    
    M(i) = getframe(fig, rect);
    
    ax.Units = 'normalized';
    
    
end

name = [pth, 'norne'];
vo = VideoWriter(name);
vo.FrameRate = ns/30;
open(vo);

writeVideo(vo, M);

close(vo)

%%

step = 1500;
hs = ResultHandler('dataFolder', ['norneWAG_nStep', num2str(step), '_nCycles5_sOres_i0.38_sOres_m0.08_W'], ...
    'dataprefix', 'state_step', 'cleardir', false);

n = 1500;
sw = cell(n,1);
for i = 1:n
    sw{i} = hs{i};
end

%%

ns = 1625;
[qOWater, qOSolvent] = deal(zeros(ns,14));

for i = 1:ns
    qOr = horzcat(states{i}.wellSol.qOr);
    qOSolvent(i,:) = qOr;
    qOWater(i,:)   = qOr;
    if i > 500
        qOr = horzcat(sw{i-500+1}.wellSol.qOr);
        qOWater(i,:) = qOr;
    end
end

dT = 4*year/2000;
OS_cum = -cumsum(sum(qOSolvent,2)*dT,1);
OW_cum = -cumsum(sum(qOWater,2)*dT,1);

qOS = -sum(qOSolvent,2);
qOW = -sum(qOWater,2);

%%

close all;

pos = [500, 500, 2000,1000];
fig = figure('position', pos);

t = linspace(0,4*1625/2000,1625);
ll = lines(5);
l = ll([1,3],:);

hold on
tEq = round(1.225*2000/4);
tS = 501;
tE = 1625;
fill(t([tS:tEq,tEq:-1:tS])', [qOW(tS:tEq); qOS(tEq:-1:tS)], brighten(ll(2,:),0.2), 'edgecolor', 'none')
fill(t([tEq:tE,tE:-1:tEq])', [qOW(tEq:tE); qOS(tE:-1:tEq)], brighten(ll(5,:),0.2), 'edgecolor', 'none')

ind = [1,500, 600:100:1500,1625];
col = [1;repmat([2;1],5,1); 1];

hndl(1) = plot(t, qOW, 'color', l(1,:), 'lineWidth', 5);
for i = 1:numel(ind)-1
    ii = ind(i):ind(i+1)-1;
    hp = plot(t(ii), qOS(ii), 'color', l(col(i),:), 'lineWidth', 5);
    axis([t(1) t(end) 0 11.5])
    if i == 2
        hndl(2) = hp;
    end
end

for i = 2:2:numel(ind)
    line([t(ind(i)), t(ind(i))], [0, 11.5], 'color', 'k', 'linewidth', 1.5);
end

box on

xlabel('Time [years]');
ylabel('Oil production rate [m^3/s]');
leg = legend(hndl, 'Water', 'Solvent');

leg.FontSize = 45;
ax = gca;
ax.FontSize = 20;
ax.XLabel.FontSize = 45;
ax.YLabel.FontSize = 45;


savepng('norneProduction');
%%

close all
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
view([90,80]);

nf = diff(G.cells.facePos);

    
%     alpha = (s(nc(:,1)) + s(nc(:,2)))./2;

alpha = @(c,s) plotSaturationOpacity(G,c,s);

%     plotGrid(G, 'facecolor', 'b', 'facealpha', 'flat', 'faceVertexAlpha', alpha(G.cells.faces));

wc = vertcat(W.cells);
f = boundaryFaces(G);
plotFaces(G, f, 'facecolor', 'none', 'edgealpha', 0.1);
% ind = [2];
for i = 1:1625
    
    c = states{i}.s(:,2) < 0.9;
    c(wc) = false;
    
    ss = 1 - states{i}.s(:,2);
%     c = 1:G.cells.num;
%     plotSaturationOpacity(G, c, s{i}.s(:,4));
    
    plotGrid(G, c, 'facecolor', 'b', 'edgecolor', 'none', ...
         'facealpha', 'flat', 'faceVertexalphadata', alpha(c, ss));
% %     plotSaturationOpacity(G, c, s{i}.s(c,4));
    title(num2str(i))
    drawnow
    
    pause(0.001);
    
    
end

%%
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
