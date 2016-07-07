mrstModule add deckformat

if ~(makeNorneSubsetAvailable() && makeNorneGRDECL()),
   error('Unable to obtain simulation model subset');
end

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G(1));
rock = grdecl2Rock(grdecl, G.cells.indexMap);
%%
mrstModule add coarsegrid msrsb ad-core mrst-gui incomp
%%
totTime = 100*year;
N_step = 100;
dt = totTime/N_step;

pv = poreVolume(G, rock);

wells = [13, 88,  -1; ...
         18, 87,  -1; ...
         36, 90,  -1; ...
         10, 15,  -1; ...
         24, 32,  1; ...
         8,  45,  1; ...
         16, 55,  1];

W = [];
[inum, pnum] = deal(1);
for i = 1:size(wells, 1);
    W = verticalWell(W, G, rock, wells(i, 1), wells(i, 2), [], 'comp_i', [1, 0], 'type', 'bhp');
    if wells(i, 3) == 1
        if 0
            W(i).val = 250*barsa;
        else
            W(i).val = -sum(pv)/(totTime*sum(wells(:, 3) == 1));
            W(i).type = 'rate';
        end

        W(i).name = ['P', num2str(pnum)];
        W(i).sign = -1;
        pnum = pnum + 1;

    else
        W(i).val = 500*barsa;
        W(i).sign = 1;
        W(i).name = ['I', num2str(inum)];
        inum = inum + 1;
    end
end


close all
plotGrid(G, 'FaceColor', 'none', 'EdgeA', .2)
plotWell(G, W)
axis tight

plotGrid(G, vertcat(W.cells), 'FaceColor', 'none', 'EdgeColor', 'b')
% 
% figure;
% ijk = gridLogicalIndices(G);
% ijk = [ijk{:}];
% plotToolbar(G, ijk)
% axis tight off
%%

T = getFaceTransmissibility(G, rock);
% T(G.faces.tag > 0) = 0.1*T(G.faces.tag > 0);

% fluid = initSimpleFluid('mu', [1, 1]*centi*poise, 'n', [1, 1],...
%                         'rho', [750, 1000]*kilogram/meter^3);
fluid = initSimpleFluid('mu', [1, 5]*centi*poise, 'n', [2, 2],...
                        'rho', [750, 1000]*kilogram/meter^3);

%%

lim = inf
% lim = 5;
state0 = initResSol(G, 0, [0, 1]);
gravity reset off

psolve = @(state) incompTPFA(state, G, T, fluid, 'Wells', W, 'use_trans', true);
solver = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

state = psolve(state0);

states = state;
for i = 1:min(N_step, lim)
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state = psolve(states(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states = [states; solver(state)];
    fprintf('Ok!\n');
end

%%
global METISPATH
useMETIS = ~isempty(METISPATH);

G.cells.facePos = double(G.cells.facePos);

cdims = ceil(G.cartDims./[15, 10, 10]);
% p = partitionUniformPadded(G, cdims);
% p = partitionUI(G, cdims);

if useMETIS
    p = partitionMETIS(G, T, 250);
    p = processPartition(G, p);
else
    padded = partitionUniformPadded(G, [cdims(1:2), 1]);
    uni = partitionUI(G, [1, 1, cdims(3)]);
    p = padded.*uni;

    G_fault = makeInternalBoundary(G, find(G.faces.tag > 0));
    p = processPartition(G_fault, p);
end
p = compressPartition(p);

mrstModule add agglom
p0 = p;
fconn = ones(G.faces.num, 1);
fconn(G.faces.tag > 0) = 0;
p = mergeBlocksByConnections(G, p, fconn, 25);

p = processPartition(G, p);
p = compressPartition(p);

CG = generateCoarseGrid(G, p);
figure; plotCellData(G, mod(p, 13))
plotGrid(CG, 'facec', 'none', 'edgec', 'w', 'linewidth', 2)
%%
mrstModule add msrsb
CG = coarsenGeometry(CG);
CG = addCoarseCenterPoints(CG);
CG = setCentersByWells(CG, W);
           
CG = storeInteractionRegion(CG, 'ensureConnected', false);
%%
A = getIncomp1PhMatrix(G, T);

getBasis = @(A) getMultiscaleBasis(CG, A, 'type', 'rsb', 'useMex', true, 'iterations', 150);
basis0 = getBasis(A);

%%
basis = basis0;

updateBasis = true;
W_ms = W;
fn = getSmootherFunction('type', 'ilu');

psolve = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 0);

states_ms = psolve(state0, basis);
for i = 1:min(N_step, lim)
    state = states_ms(end);

    if updateBasis && mod(i, 10) == 0 && i > 1
        A = getIncomp1PhMatrix(G, T, state, fluid);
        basis = getBasis(A);
    end
    
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state, basis);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_ms = [states_ms; solver(state)];
    fprintf('Ok!\n');
end

%%
psolve = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 5, 'useGMRES', true);
basis = basis0;

states_it = psolve(state0, basis);
for i = 1:min(N_step, lim)
    state = states_it(end);
    A = getIncomp1PhMatrix(G, T, state, fluid);
    
    if updateBasis && mod(i, 10) == 0 && i > 1
        A = getIncomp1PhMatrix(G, T, state, fluid);
        basis = getBasis(A);
    end
    
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state, basis);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_it = [states_it; solver(state)];
    fprintf('Ok!\n');
end

%%
close all; plotToolbar(G, states);
axis equal tight off
daspect([1 1 0.1])
c = caxis();

figure; plotToolbar(G, states_ms);
axis equal tight off
daspect([1 1 0.1])
caxis(c);

figure; plotToolbar(G, states_it);
axis equal tight off
daspect([1 1 0.1])
caxis(c);

ws_ref = convertIncompWellSols(W, states, fluid);
ws_ms = convertIncompWellSols(W, states_ms, fluid);
ws_it = convertIncompWellSols(W, states_it, fluid);

ws = {ws_ref, ws_ms, ws_it};
wsn = {'Finescale', 'Multiscale', 'Multiscale (5 cycles)'};
shortname = {'ref', 'ms', 'it'};
plotWellSols(ws, 'datasetnames', wsn)


%%
close all

FZ = 15;
inj  = find(vertcat(W.sign) > 0);
prod = find(vertcat(W.sign) < 0);
markers = {'-', '--', '.-'};
figsize = [100, 100, 900, 450];
grafsize = [100, 100, 550, 350];
time = cumsum(repmat(dt, N_step, 1))/year;

%% Producer BHP
figure('position', grafsize, 'name', 'prod_bhp');
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    [welldata, wellnames, fldnames] = getWellOutput(ws{i}, 'bhp', prod);
    plot(time, welldata(2:end, :, 1)/barsa, markers{i})
    axis tight
    ylim([200, 400])
end

ylabel('Bottom hole pressure [bar]')
xlabel('Time [year]')

h = legend(wellnames, 'box', 'off', 'location', 'north', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
p = gcf;
tmpf = figure;
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    plot(1, 1, markers{i}, 'color', 'k', 'linewidth', 2, 'markersize', 20)
end
h2 = legend(wsn, 'location', 'north', 'box', 'off', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
copyobj(h2, p)
close(tmpf)
%% Producer oil rate
figure('position', grafsize, 'name', 'prod_oil');
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    [welldata, wellnames, fldnames] = getWellOutput(ws{i}, 'qOs', prod);
    d = abs(welldata(2:end, :, 1));
    d = convertTo(d, 1000*stb/day);
    plot(time, d, markers{i})
    axis tight
end

ylabel('Oil production [Mbbl/day]')
xlabel('Time [year]')

h = legend(wellnames, 'box', 'off', 'location', 'north', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
p = gcf;
tmpf = figure;
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    plot(1, 1, markers{i}, 'color', 'k', 'linewidth', 2, 'markersize', 20)
end
h2 = legend(wsn, 'location', 'north', 'box', 'off', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
copyobj(h2, p)
close(tmpf)

%% Producer water rate
figure('position', grafsize, 'name', 'prod_water');
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    [welldata, wellnames, fldnames] = getWellOutput(ws{i}, 'qWs', prod);
    d = abs(welldata(2:end, :, 1));
    d = convertTo(d, 1000*stb/day);
    plot(time, d, markers{i})
    axis tight
end

ylabel('Water production [Mbbl/day]')
xlabel('Time [year]')

h = legend(wellnames, 'box', 'off', 'location', 'north', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
p = gcf;
tmpf = figure;
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    plot(1, 1, markers{i}, 'color', 'k', 'linewidth', 2, 'markersize', 20)
end
h2 = legend(wsn, 'location', 'north', 'box', 'off', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
copyobj(h2, p)
close(tmpf)

%% Injector water
figure('position', grafsize, 'name', 'inj_water');
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    [welldata, wellnames, fldnames] = getWellOutput(ws{i}, 'qWs', inj);
    d = abs(welldata(2:end, :, 1));
    d = cumsum(d, 1);
    plot(time, d, markers{i})
    axis tight
end

ylabel('Cumulative injection [m^3]')
xlabel('Time [year]')

h = legend(wellnames, 'box', 'off', 'location', 'north', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
p = gcf;
tmpf = figure;
set(gca, 'fontsize', FZ)
hold on
for i = 1:numel(ws)
    plot(1, 1, markers{i}, 'color', 'k', 'linewidth', 2, 'markersize', 20)
end
h2 = legend(wsn, 'location', 'north', 'box', 'off', 'color', 'none', 'linewidth', 0, 'EdgeColor', [1 1 1]);
copyobj(h2, p)
close(tmpf)
%%
v = [-170, 15];
% v = [70, 20];

data = {states, states_ms, states_it};

for i = 1:numel(data)
    s = data{i};
    
    figure('position', figsize, 'name', [shortname{i}, '_sw']);
    plotCellData(G, s(end).s(:, 1), 'EdgeColor', 'w', 'EdgeAlpha', 0.25)
    view(v)
    axis tight off
    daspect([1 1 0.25])
    plotWell(G, W, 'color', 'k')
    if i == 1
        c = caxis();
    else
        caxis(c);
    end
end


for i = 1:numel(data)
    s = data{i};
    
    figure('position', figsize, 'name', [shortname{i}, '_pressure']);
    plotCellData(G, s(end).pressure, 'EdgeColor', 'w', 'EdgeAlpha', 0.25)
    view(v)
    axis tight off
    daspect([1 1 0.25])
    plotWell(G, W, 'color', 'k')
    if i == 1
        c = caxis();
    else
        caxis(c);
    end
end

%%
c = [min(log10(rock.perm(:))), max(log10(rock.perm(:)))];

snp = 'xyz';

for i = [1, 3]
    figure('position', figsize, 'name', ['perm', snp(i)]);
    plotCellData(G, log10(rock.perm(:, i)), 'EdgeColor', 'w', 'EdgeAlpha', 0.25)
    view(v)
    axis tight off
    daspect([1 1 0.25])
    plotWell(G, W, 'color', 'k')
    caxis(c)
end

figure('position', figsize, 'name', 'poro');
plotCellData(G, (rock.poro), 'EdgeColor', 'w', 'EdgeAlpha', 0.25)
view(v)
axis tight off
daspect([1 1 0.25])
plotWell(G, W, 'color', 'k')



%% Saving plots
mrstModule add export_fig

msdir = mrstPath('query', 'multiscale-devel');

fdir = fullfile(msdir, 'papers', 'mrsb', 'incomp', 'figures', 'norne');
if ~exist(fdir, 'dir')
    mkdir(fdir)
end
saveplot = @(h, ext) export_fig(fullfile(fdir, [get(h, 'Name'), '.', ext]), h, '-transparent');

ch = get(0, 'Children');
for i = 1:numel(ch)
    figure(ch(i));
    
    disp(['Saving plot ', num2str(i), ' of ', num2str(numel(ch))]);
    
    types = get(get(gca, 'Children'), 'Type');
    
    if any(strcmpi(types, 'patch'))
        saveplot(ch(i), 'png')
    else
        saveplot(ch(i), 'pdf')
    end
end
%% Debug stuff follows
return
%%
close all
plotGrid(G, 'facec', 'none', 'edgea', .05)
plotWellData(G, W, {states(1).wellSol.flux})
axis equal tight off
daspect([1 1 0.1])

figure
plotGrid(G, 'facec', 'none', 'edgea', .05)
plotWellData(G, W, {states_ms(1).wellSol.flux})
axis equal tight off
daspect([1 1 0.1])
%%
ms = states_ms(1).pressure;
rf = states(1).pressure;
it = states_it(1).pressure;
% ms = states_ms(end).s(:, 1);
% rf = states(end).s(:, 1);

N = 1
norm(ms - rf, N)/norm(rf, N)
norm(it - rf, N)/norm(rf, N)


%%
rf = sqrt(sum(faceFlux2cellVelocity(G, states(1).flux).^2, 2));
ms = sqrt(sum(faceFlux2cellVelocity(G, states_ms(1).flux).^2, 2));

fn = @(x) log10(x);

figure;
plotCellData(G, fn(rf), 'EdgeColor', 'w', 'EdgeAlpha', .2);
axis  tight off
view(90, 90)
c = caxis();
set(gcf, 'Name', 'v_ref');

figure;
plotCellData(G, fn(ms), 'EdgeColor', 'w', 'EdgeAlpha', .2);
axis  tight off
view(90, 90)
caxis(c);
set(gcf, 'Name', 'v_ms');


figure;
plotCellData(G, abs(rf-ms)./max(rf), 'EdgeColor', 'w', 'EdgeAlpha', .2);
axis  tight off
% caxis([0, 1])
view(90, 90)
set(gcf, 'Name', 'v_err');
%%
for i = 1:CG.cells.num
    figure(1); clf
    plotGrid(G, CG.cells.interaction{i}, 'facec', 'none')
    plotGrid(G, CG.partition == i)
    pause
end
%%
figure;
for i = 1:numel(W)
    plotCellData(G, CG.partition, ismember(CG.partition, CG.partition(W(i).cells)), 'facea', .2)
    
end
plotWell(G, W);
%%
bad = states_ms(1).pressure > max(states(1).pressure);
badc = unique(CG.partition(bad));

% for i = 1:CG.cells.num
for i = badc'
    figure(1); clf
    plotGrid(G, CG.cells.interaction{i}, 'facec', 'none')
    plotGrid(G, CG.partition == i)
    pause
end

%%
figure('position', figsize);
% Gplot = generateCoarseGrid(G, ones(G.cells.num, 1));
% plotGrid(Gplot, 'faceColor', 'w', 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FAcea', .1)
% plotFaces(G, 'faceColor', 'y', 'EdgeColor', 'none', 'outline', true, 'facea', .2)
% plotGrid(G, 'edgecolor', 'none', 'facea', .2)

fa = G.cells.faces(G.cells.faces(:, 2) == 6, 1);
d = log10(rock.perm(:, 3));

layers = 1:G.cells.num;
[ii, jj, kk] = gridLogicalIndices(G);
% layers = mod(kk, 5) == 0 | kk == max(kk) | kk == min(kk);
layers = kk == max(kk) | kk == min(kk);

plotFaces(G, fa(layers), d(layers), 'facea', .2, 'edgec', 'none', 'outline', true)

% set(gca, 'ZDir', 'normal')
% camlight(30, 60)
% set(gca, 'ZDir', 'reverse')

view(v)
axis tight off
daspect([1 1 0.25])
plotWell(G, W, 'color', 'k')
