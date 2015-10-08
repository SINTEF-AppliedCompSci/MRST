%% Lorenz coefficient for layers of SPE 10, Model 2
mrstModule add diagnostics spe10 incomp

%% Base model
% We set up a grid. Later, we will assign petrophysical properties from one
% single layer at the time to compute Lorenz coefficients inside the main
% loop
cartDims = [  60,  220,  1];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();
G  = cartGrid(cartDims, physDims);
G  = computeGeometry(G);

% Set parameters describing the 
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,  30;
              1,    1,   220,  220, 111];
wname    = {'P1', 'P2', 'P3', 'P4', 'I'};

% Set fluid model: Here, we use a simple fluid model with properties that
% are typical for water. Replace this by a multiphase fluid object if you
% also want to include fluid effects in the calculation
fluid = initSingleFluid('mu', 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);

%% Compute Lorenz coefficient for each layer
Lc = zeros(85,1);
h = waitbar(0,'Computing Lorenz coefficients ...');
for n=1:85
    
    % --- Set petrophysical data for this particular layer
    % To avoid problems with very small porosity values, we explicitly
    % impose a lower threshold of 1e-4
    rock = SPE10_rock(1:cartDims(1),1:cartDims(2),n);
    rock.perm = convertFrom(rock.perm, milli*darcy);
    rock.poro = max(rock.poro, 1e-4);
    
    % --- Set up well model
    % To ensure that we get the correct well index when updating the
    % petrophysical data, we simply regenerate the well objects.
    W = [];
    for w = 1 : numel(wtype),
        W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1, ...
            'Type', wtype{w}, 'Val', wtarget(w), ...
            'Radius', wrad(w), 'Name', wname{w}, ...
            'InnerProduct', 'ip_tpf');
    end
    
    % --- Initiate and solve flow problem 
    rS = initState(G, W, 0);
    T  = computeTrans(G, rock);
    rS = incompTPFA(rS, G, T, fluid, 'wells', W);
    
    % --- Compute flow diagnostics
    D       = computeTOFandTracer(rS, G, rock, 'wells', W, 'maxTOF', inf);
    [F,Phi] = computeFandPhi(poreVolume(G,rock), D.tof);
    Lc(n)   = computeLorenz(F,Phi);
    waitbar(n/85);
end
close(h);
clf; set(gcf,'Position',[470 420 960 400]);
h=bar(Lc,'Style','hist'); 
axis tight; set(h,'FaceColor',[.95 .95 1],'EdgeColor',[0 0 .7]);

% [nmin, nmax] = deal(22,61);

%% Look at case with lowest Lorenz coeffficient
minCase = false;
if minCase
    [~,n]=min(Lc);
else
    [~,n]=max(Lc);
end
rock = SPE10_rock(1:cartDims(1),1:cartDims(2),n);
rock.perm = convertFrom(rock.perm, milli*darcy);
rock.poro = max(rock.poro, 1e-4);
pv = poreVolume(G, rock);

nwloc = wloc;
for nstep=1:2
    W = [];
    for w = 1 : numel(wtype),
        W = verticalWell(W, G, rock, nwloc(1,w), nwloc(2,w), 1, ...
            'Type', wtype{w}, 'Val', wtarget(w), ...
            'Radius', wrad(w), 'Name', wname{w}, ...
            'InnerProduct', 'ip_tpf');
    end

    rS = initState(G, W, 0);
    T  = computeTrans(G, rock);
    rS = incompTPFA(rS, G, T, fluid, 'wells', W);
    D  = computeTOFandTracer(rS, G, rock, 'wells', W, 'maxTOF', inf);
    WP = computeWellPairs(rS, G, rock, W, D);
    [F,Phi] = computeFandPhi(pv, D.tof);
    computeLorenz(F,Phi)

    figure('Position',[710   420   720   400]);
    plotCellData(G,rock.poro,'EdgeColor','none');
    plotGrid(Gbb,'FaceColor','none');
    plotWell(G,W,'height',10,'LineWidth',4);
    plotWellPairConnections(G,WP,D,W,pv,1e-4);
    cmap=jet(128); colormap(.6*cmap + .4*ones(size(cmap))); clear cmap;
    view(0,70); set(gca,'DataAspectRatio',[1 1 .5]); axis off
    cax = caxis;

    pos = [.05 .025 .25 .35;  .725 .025 .25 .35; ...
        .05 .625 .25 .35; .725 .625 .25 .35];
    for i=1:4
        axes('position', pos(i,:));
        rad = sum(bsxfun(@minus,G.cells.centroids,....
            G.cells.centroids(W(i).cells,:)).^2,2);
        plotCellData(G,rock.poro,rad<3500,'EdgeColor','k','EdgeAlpha',.1);
        plotWell(G,W(i),'height',10,'LineWidth',10);
        caxis(cax);
        view(0,70); set(gca,'DataAspectRatio',[1 1 .5]); axis off tight;
    end
    
    if minCase,
        nwloc     = [  1,   53,     1,   59,  30;
                       1,    2,   218,  220, 111];
    else
        nwloc     = [  1,   60,     7,   60,  30;
                       6,   11,   220,  219, 111];
    end
end