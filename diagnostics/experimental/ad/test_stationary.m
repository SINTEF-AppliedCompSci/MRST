mrstModule add deckformat ad-fi diagnostics spe10 internal/mrst-gui
gravity off

bhpWells = false;
uniformRock = false;
grid = 'spe10';

clear rock


% fn    = fullfile(ROOTDIR, 'modules', 'diagnostics', 'experimental', 'ad', 'simple10x1x10.data');
% fn    = fullfile(ROOTDIR, 'modules', 'diagnostics', 'experimental', 'ad', 'trivial.data');



deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
fluid_ad = initDeckADIFluid(deck);
% fluid_ad = initSimpleADIFluid();

W = [];
alpha = 1;

min_well = 0.001/day;

switch grid
    case '1d_test'
        N = 150;
        mid = ceil(N/2);

        G = cartGrid([N, 1, 1]);
        G = computeGeometry(G);

        rock.poro = ones(G.cells.num, 1);

        rock.poro(1:mid) = .1;
        rock.poro(mid:end) = .5;
        rock.poro(mid) = 1;
        rock.perm = 1*milli*darcy*ones(G.cells.num, 1);

        W = verticalWell(W, G, rock, 1, 1, [], 'Val', 1/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, N, 1, [], 'Val', 1/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, mid, 1, [], 'Val', 0*barsa, 'Type', 'bhp');
        alpha = 1e-4;

    case 'norne_synthetic'
        prefix = '/data/norne/res/BC0407';
        mrstModule add deckformat ad-fi internal/mrst-gui diagnostics
        cdims = [46   112    22];
        tmp = readEclipseOutputFileUnFmt([prefix '.EGRID']);

        grdecl.COORD = tmp.COORD.values;
        grdecl.ZCORN = tmp.ZCORN.values;
        grdecl.ACTNUM = tmp.ACTNUM.values;
        grdecl.cartDims = cdims;

        G = processGRDECL(grdecl, 'SplitDisconnected', false);
        G = computeGeometry(G);

        badcells = [8897,11160,13423,15686,26019,28282,30545,32808,34944,37206,44347];
        ind = true(G.cells.num, 1);
        ind(badcells) = false;

        ijk = gridLogicalIndices(G);
        ind(ijk{3} < 6) = false;


        nornedeck = readEclipseDeck([prefix, '_3.DATA']);
        nornedeck = convertDeckUnits(nornedeck);

        rock  = initEclipseRock(nornedeck);
        rock  = compressRock(rock, G.cells.indexMap);
        fluid = initEclipseFluid(nornedeck);
        rock.perm = rock.perm(ind, :);
        rock.poro = rock.poro(ind);
        rock.ntg = rock.ntg(ind);
        G = extractSubgrid(G, find(ind));
        G = computeGeometry(G);
        G.cartDims = cdims;

        % inj
        W = verticalWell(W, G, rock, 6, 11, [], 'Val', 1.2/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, 29, 11, [], 'Val', 1/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, 41, 102, [], 'Val', 1.5/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, 16, 100, [], 'Val', 0.5/day, 'Type', 'rate');
        W = verticalWell(W, G, rock, 19, 90, [], 'Val', 0.5/day, 'Type', 'rate');
        % prod
        W = verticalWell(W, G, rock, 20, 33, [], 'Val', 100*barsa, 'Name', 'P1');
    %     W = verticalWell(W, G, rock, 12, 45, [], 'Val', 100*barsa, 'Name', 'P2');
    %     W = verticalWell(W, G, rock, 16, 46, [], 'Val', 100*barsa);
    %     W = verticalWell(W, G, rock, 10, 41, [], 'Val', 100*barsa);
        min_well = 0.1/day;
    case 'spe10'
        Nx = 60;
        Ny = 60;
        Nz = 1;
        offset = [0 0 36];
        dims = [Nx, Ny, Nz];
        pdims = dims.*[20, 20, 2]*ft;
        G = cartGrid(dims, pdims);
        G = computeGeometry(G);


        if ~uniformRock
            rock = SPE10_rock(1:Nx + offset(1), 1:Ny + offset(2), (1:Nz) + offset(3));
            rock.perm = convertFrom(rock.perm, milli*darcy);
            rock.poro = rock.poro + 0.01;
            zp = rock.poro < 0.01;
            rock.poro(zp) = 0.01;
        else
            rock.perm = repmat(0.1*darcy, G.cells.num, 1);
            rock.poro = 0.1 + 0*rock.perm(:,1);
        end


    %     rock.perm = repmat(100*milli*darcy, G.cells.num, 1);
    %     rock.poro = repmat(0.3, G.cells.num, 1);
        if bhpWells
            W = verticalWell(W, G, rock, 1, 1, 1, 'Val', 600*barsa);
            W = verticalWell(W, G, rock, Nx, Ny, 1, 'Val', 730*barsa);
            W = verticalWell(W, G, rock, Nx, 1, 1, 'Val', 550*barsa);
            W = verticalWell(W, G, rock, 1, Ny, 1, 'Val', 490*barsa);
            W = verticalWell(W, G, rock, ceil(Nx/2), ceil(Ny/2), 5, 'Val', 120*barsa);
            W = verticalWell(W, G, rock, ceil(Nx/2), ceil(Ny/2)+1, 5, 'Val', 100*barsa);
        else
            W = verticalWell(W, G, rock, 1, 1, 1, 'Val', 11/day, 'Type', 'rate');
            W = verticalWell(W, G, rock, Nx, Ny, ceil(2*Nz/4), 'Val', 13/day, 'Type', 'rate');
            W = verticalWell(W, G, rock, Nx, 1, ceil(3*Nz/4), 'Val', 25.5/day, 'Type', 'rate');
            W = verticalWell(W, G, rock, 1, Ny, ceil(4*Nz/4), 'Val', 1/day, 'Type', 'rate');
            W = verticalWell(W, G, rock, ceil(Nx/2), ceil(Ny/2), Nz, 'Val', 0*barsa, 'Type', 'bhp', 'Name', 'prod');
        end
    case 'norne'
        [G, rock, W] = Norne_setup([], true);
        min_well = 0.1/day;
end
W_incomp = W;
for i = 1:numel(W_incomp)
    W_incomp(i).compi = [1 0];
end

T = computeTrans(G, rock);
s = initADISystem(deck, G, rock, fluid_ad);
pv = poreVolume(G, rock);

state0 = initResSol(G, 0*barsa, [1 0 0]);
state0.wellSol = initWellSol(W, 0);
state0.rs = zeros(G.cells.num, 1);
%%
[state, D1, grad] = solveStationaryPressure(G, state0, s, W, fluid_ad, pv, T, 'objective', objective);

%%
figure(1); clf
plotToolbar(G, rock, 'FaceAlpha', 0.5)
plotWell(G, W)
axis tight off

%% Do optimization
objective = getObjectiveDiagnostics(G, rock, 'minlorenz', []);
% objective = getObjectiveDiagnostics(G, rock, 'minpvdiff')
[D_best W_best history] = optimizeTOF(G, W, fluid_ad, pv, T, s,...
                                     state0, 0.0001/day*ones(numel(D0.inj), 1), objective, ...
                                     'alpha', alpha);

%% Plot all the iterations where reduction in objective happened
clf
plotToolbar(G, [history.D])
axis tight off
plotWell(G, W);

%% Plot changes in injector partition
figure(2); clf;
subplot(1,3,1)
plotCellData(G, D_initial.ipart)
title('Initial partition')
view(0,90);
axis tight off

subplot(1,3,2)
plotCellData(G, D_best.ipart)
title('Optimized partition')
view(0,90);
axis tight off

subplot(1,3,3)
plotCellData(G, log10(sum(D_best.tof, 2) - sum(D_initial.tof, 2)), D_best.ipart ~= D_initial.ipart, 'edgea', .05)
plotGrid(G, 'facec', 'none', 'edgea', .1)
title('Changes in tof in difference in partition')
view(0,90);
axis tight off

fastRotateButton
%% Plot sweep diagrams
close all
pv = poreVolume(G, rock);

Fs = [];
Phis = [];
ltext = {};
hist = history.D;
for i = 1:numel(hist)
    [F, Phi] = computeFandPhi(pv, hist(i).tof);
    Fs = [Fs, F];
    Phis = [Phis, Phi];
    Lc = computeLorenz(F, Phi);
    ltext = {ltext{:}, sprintf('Iteration %d (L_c: %f)', i-1, Lc)};
end
plot(Phis, Fs)
legend(ltext, 'location', 'SouthEast')
title('Sweep diagrams')


%% Plot history over time
clf;
hold on
ci = cumsum(history.iterations);
ci = ci(1:end-1);

plot(history.values)
plot(ci, history.values(ci, :), '*')
