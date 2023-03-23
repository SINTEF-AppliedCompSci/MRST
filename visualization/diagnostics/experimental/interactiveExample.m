mrstModule add spe10 diagnostics

cartDims = [  60,  220,  15];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500,   500  ] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,  20, 40;
              1,    1,   220,  220, 130, 90];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1', 'I2'};


   rock = SPE10_rock(1:cartDims(end));
   rock.perm = convertFrom(rock.perm, milli*darcy);
   rock.poro = max(rock.poro, 1e-4);
   G  = cartGrid(cartDims, physDims);
   G  = computeGeometry(G);

W = [];
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w});
end
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);



%%
z = G.cells.centroids(:,3);
sat = (z - min(z))./(max(z) - min(z));
% rs = initState(G, W, 0);
rs = initResSol(G, 0, [sat 1-sat]);
T  = computeTrans(G, rock);
rs = incompTPFA(rs, G, T, fluid, 'wells', W);
%%
topP = 300*barsa;
sides = {'North', 'South', 'East', 'West'};
bc = [];
for i = 1:4
    bc = pside(bc, G, sides{i}, 1);
end
zc = G.cells.centroids(:, 3);

gravity reset
gravity on
grav = gravity();

p = @(z, top, rho) top + rho*grav(3)*(z-min(zc))./(max(zc) - min(zc));

bc.value = p(G.faces.centroids(bc.face, 3), topP, 1000*kilogram/meter^3);

src = addSource([], 1:G.cells.num, rock.poro);

rs2 = incompTPFA(rs, G, T, fluid, 'src', src, 'bc', bc);
%%
close all
interactiveDiagnostics(G, rock, W, rs2, 'state', rs)
view(-100, 70)
%%
close all
rs2.poro = rock.poro;
plotToolbar2(G, rs2)
%%
close all
h = figure;
plotToolbar2(G, [rock rock]);
editWells(G, W, rock, 'figure', h);
%%
close all
plotToolbar2(G, rock);
