mrstModule add deckformat ad-fi diagnostics spe10 internal/mrst-gui

clear rock
fn    = fullfile(ROOTDIR, 'modules', 'diagnostics', 'experimental', 'ad', 'trivial.data');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
% fluid_ad = initDeckADIFluid(deck);
% mu = [0.4 2 1]*centi*poise;
mu = [1 1 1]*centi*poise;
% n = [2 2 1];
n = [1 1 1];

fluid_ad = initSimpleADIFluid('mu', mu, 'n', n);


fluid = initSimpleFluid('mu' , mu(1:2), ...
                         'rho', [1, 1], ...
                         'n'  , n(1:2));

% offset [0 80 10], 60x60 and halfhalf works well
Nx = 60; %  60
Ny = 220; % 220
Nz = 1;
offset = [0 0 0];
dims = [Nx, Ny, Nz];
pdims = dims.*[20, 20, 2]*ft;
alpha = 1;
halfHalf = false;
spe = true;
totTime = 100*day;
qfs = true;
fullperf = true;
G = cartGrid(dims, pdims);
G = computeGeometry(G);

if spe
    rock = SPE10_rock((1:Nx) + offset(1), (1:Ny) + offset(2), (1:Nz) + offset(3));
    rock.perm = convertFrom(rock.perm, milli*darcy);

    if halfHalf
        rock2 = SPE10_rock((1:Nx) + offset(1), (1:Ny) + offset(2), (1:Nz) + offset(3) + 35);
        rock2.perm = convertFrom(rock2.perm, milli*darcy);

        ijk = gridLogicalIndices(G);
        subs = ijk{1} < round(max(ijk{1})/2);
        rock.perm(subs, :) = rock2.perm(subs, :);
        rock.poro(subs) = rock2.poro(subs);
    end

    minval = 0.01;
    zp = rock.poro < minval;
    rock.poro(zp) = minval;
else
    verbose = true;
    if 0
        grdecl  = simpleGrdecl(dims, 0.15);%, 'physDims', pdims);
        G       = processGRDECL(grdecl);
        for i = 1:3
            G.nodes.coords(:,i) = G.nodes.coords(:,i).*pdims(i);
        end
        G       = computeGeometry(G);
    end


    [rock.perm, L] = logNormLayers(dims, [10, 500, 25, 1000]);
    rock.perm = convertFrom(rock.perm, milli*darcy);
    rock.poro = rock.perm./max(rock.perm);
end

% rock.perm = 100*milli*darcy + 0*rock.perm;
% rock.poro = 0*rock.poro + 0.3;

pv = poreVolume(G, rock);
prorate = -.75*sum(pv)/totTime;


if fullperf
    v = @(i, n) [];
else
    v = @(i, n) ceil(Nz*(i+1)/n);
end
W = [];
if qfs
    %qfs
    Ninj = 4;
    injrate = -prorate/Ninj;
    W = verticalWell(W, G, rock, 1, 1,   v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
    W = verticalWell(W, G, rock, Nx, Ny, v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
    W = verticalWell(W, G, rock, Nx, 1,  v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
    W = verticalWell(W, G, rock, 1, Ny,  v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
else
    % half spot
    Ninj = 2;
    injrate = -prorate/Ninj;
    W = verticalWell(W, G, rock, 1, 1,   v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
    W = verticalWell(W, G, rock, Nx, Ny, v(numel(W), Ninj), 'Val', injrate, 'Type', 'rate');
end
targets = 1:numel(W);
% wlimit = mean(injrate)/2*ones(numel(W), 1);
wlimit = (injrate/2)*ones(numel(W), 1);
if 1
    W = verticalWell(W, G, rock, ceil(Nx/2), ceil(Ny/2), [], 'Val', prorate, 'Type', 'rate', 'Name', 'Injector');
%     W = verticalWell(W, G, rock, ceil(Nx/2), ceil(Ny/2), [], 'Val', 100*barsa, 'Type', 'bhp', 'Name', 'Injector');

else
    W = verticalWell(W, G, rock, ceil(Nx/2) - 5, ceil(Ny/2), Nz, 'Val', prorate/2, 'Type', 'rate', 'Name', 'Injector');
    W = verticalWell(W, G, rock, ceil(Nx/2) + 5, ceil(Ny/2), 1, 'Val', prorate/2, 'Type', 'rate', 'Name', 'Injector');
end
clf;
plotToolbar(G, rock, 'facea', .4)
plotWell(G, W)
axis tight off
%%
a = G.faces.neighbors;
f = a(:,1);
f(f==0) = a(f==0,2);
data = rock.poro(f, :);
% data = log10(rock.perm(f, 1));
% data = data./max(data);

clf
h = plotFaces(G, 1:G.faces.num, data, 'FaceVertexAlphaData', data.*((0.5*data>mean(data)) + 0.5), 'Edgecolor', 'none', 'FaceAlpha', 'flat');
alim([min(data), (Nz/2)*max(data)])
axis tight off
fastRotateButton

%%
T = computeTrans(G, rock);
objective = getObjectiveDiagnostics(G, rock, 'minlorenz', []);
% objective = getObjectiveDiagnostics(G, rock, 'mintargettof', []);

s = initADISystem(deck, G, rock, fluid_ad);


state0 = initResSol(G, 0*barsa, [0 1 0]);
state0.wellSol = initWellSol(W, 0);
state0.rs = zeros(G.cells.num, 1);

%%
figure(1); clf
% objective = getObjectiveDiagnostics(G, rock, 'minpvdiff')
[D_best W_best history] = optimizeTOF(G, W, fluid_ad, pv, T, s,...
                                     state0, wlimit, objective, ...
                                     'targets',targets, ...
                                     'alpha', alpha);
%%
clf;
subplot(1,3,1)
plotCellData(G, history.D(1).ipart);
axis tight off

plotWell(G, W)

subplot(1,3,2)
plotCellData(G, D_best.ipart);
plotWell(G, W)
axis tight off

subplot(1,3,3)
plotCellData(G, D_best.ipart, D_best.ipart ~= history.D(1).ipart);
plotWell(G, W)
axis tight off
%%




W_initial = W;
W_improved = W_best;
for i = 1:numel(W_initial)
    W_initial(i).compi = [1 0];
    W_improved(i).compi = [1 0];
end

psolve = @(state, w) incompTPFA(state, G, T, fluid, 'wells', w);

tsolve   = @(state, w, dT) implicitTransport(state, G, dT, rock, ...
                                                fluid, 'wells', w);

solve = @(state, w, dT) tsolve(psolve(state, w), w, dT);

%%
state_a = state0;
state_a.s = state_a.s(:,1:2);
state_b = state_a;
initial = [];
improved = [];

dt = 1*day;
% T = 50*day;
Nt = ceil(totTime/dt);
for i = 1:(Nt)
    state_a = solve(state_a, W_initial, dt);
    state_b = solve(state_b, W_improved, dt);

    initial = [initial; state_a];
    improved = [improved; state_b];
    i
%     if i > Nt/3
%         for j = 1:Ninj
%             W_initial(j).compi = [0 1];
%             W_improved(j).compi = [0 1];
%         end
%     end
end
%%
plotc = @(x) plotCellData(G, (x.s(:,1)));
plotd = @(x,y)plotCellData(G, x.s(:,1) - y.s(:,1));
for i = 1:1:numel(initial)
    clf;
    subplot(1,3,1);
    plotc(initial(i));
    title('Uniform')
    ca = caxis();
    colorbar
    axis tight off
    subplot(1,3,2);
    plotc(improved(i));
    title('Improved')
    caxis(ca);
    colorbar
    axis tight off
    subplot(1,3,3);
    plotd(improved(i), initial(i));
    title(sprintf('Difference (%d of %d)', i, numel(initial)))
    colorbar
    axis tight off
    drawnow
%     pause(0.01)
end
%%
clf
init = zeros(1, numel(initial));
impr = init;
c = W(end).cells;
% findMob = @(x) x.wellSol
for j = 1:numel(initial)
    init(:, j) =  sum( initial(j).s(c, 2));
    impr(:, j) =  sum(improved(j).s(c, 2));
end
pd = [init .' , impr .'];
% pd = cumsum(pd);
plot(pd)
legend('Initial', 'Optimized')
title('Oil cut')
grid on
%%
% sum(initial(end).s(:,2))
% sum(initial(end).s(:,2))

extracted = @(x) 1 - sum(x(end).s(:,2).*pv)./sum(pv);
extracted(initial)
extracted(improved)
%%
threshold = 0.1;
[sweep_initial sweep_improved] = deal(zeros(G.cells.num, 1));
for i = 1:numel(initial)
    sweep_initial = sweep_initial + double(initial(i).s(:,1) > threshold);
    sweep_improved = sweep_improved + double(improved(i).s(:,1) > threshold);
end

clf
subplot(1,2,1)
plotCellData(G, double(sweep_initial))
sum(sweep_initial)
axis tight off
subplot(1,2,2)
plotCellData(G, double(sweep_improved))
sum(sweep_improved)
axis tight off

%%
figure(1); clf;
plotToolbar(G, history.gradient)
plotWell(G, W)
axis tight off
figure(2); clf;
plotToolbar(G, rock)
plotWell(G, W)
axis tight off
%%
W_placement = optimizeWellPlacementDiagnostics(G, W, rock, objective, 1:4, history.D(1), wlimit(1), state0, fluid_ad, pv, T, s);

%%
W_placement = optimizeWellPlacementDiagnostics(G, W_placement, rock, objective, 5, history.D(1), W(5).val, state0, fluid_ad, pv, T, s);


%%


[state_placement, D_placement] = solveStationaryPressure(G, state0, s, W_placement, fluid_ad, pv, T);
[state_best, D_best] = solveStationaryPressure(G, state0, s, W_best, fluid_ad, pv, T);
objective(state_best, D_best)
objective(state_placement, D_placement)
%%

W_new = addGhostWellsEverywhere(G, W, rock);

%%
W_everywhere = [W; W_new];
state2 = state0;
state2.wellSol = initWellSol(W_everywhere, 0);

[state, D, gradnew] = solveStationaryPressure(G, state2, s, W_everywhere, fluid_ad, pv, T, 'objective', objective);
%%
tmp = zeros(G.cells.num, 1);
tmp([W_everywhere.cells]) = gradnew.well;

figure(1); clf;
plotToolbar(G, tmp)

figure(2); clf;
plotToolbar(G, history(end).gradient)
