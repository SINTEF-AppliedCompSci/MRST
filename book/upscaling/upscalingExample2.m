%% Example 3: Upscaling
% In this example we illustrate how allocation factors for well pairs can
% be used to assess the quality of upscaling.  As our example, we consider
% the Tarbert formation from the 10th SPE Comparative Solution Project
% with a different well pattern consisting of two central injectors and
% producers at each of the four corners. We will compare four different
% upscaling methods: 
% 
% # harmonic upscaling of K,
% # harmonic-arithmetic upscaling of K,
% # flow-based upscaling of K with axial pressure drop,
% # global upscaling of T and WI
%
% To assess the quality of the upscaling, we contrast the allocation
% factors for the injection wells computed on the fine and the coarse
% model. Ideally, bars on the negative axis that represent the allocation
% factors for the coarse model should be the mirror of the bars on the
% positive axis that represent the allocation factors for the fine model.
% To simplify the comparison, the fine-scale allocation factors are
% indicated by lines on top of those of the coarse scale.

%% Set up fine-scale problem
mrstModule add diagnostics spe10 coarsegrid upscaling AGMG incomp libgeometry

if ~exist('agmg', 'file') || ...
      norm(agmg(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) > 1.0e-8,
   error('This example requires the AGMG linear solver package');
end

fprintf(1,'Setting up fine-scale problem ...');
cartDims = [  60,  220, 36];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
rock = SPE10_rock([1 1:cartDims(end)-1]);
rock.perm = convertFrom(rock.perm, milli*darcy);
rock.poro = max(rock.poro, 1e-4);
G  = cartGrid(cartDims, physDims);
G  = mcomputeGeometry(G);
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500,   500  ] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,  20, 40;
              1,    1,   220,  220, 130, 90];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1', 'I2'};
W = [];
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'InnerProduct', 'ip_tpf');
end
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
fprintf(1,'done\n');

%% Solve flow problem and compute flow diagnostics
fprintf(1,'Solving fine-scale problem ...');
rS = initState(G, W, 0);
hT = computeTrans(G, rock);
rS = incompTPFA(rS, G, hT, fluid, 'wells', W, 'LinSolve', @(A,b) agmg(A,b,1));
D  = computeTOFandTracer(rS, G, rock, 'wells', W);
WP = computeWellPairs(rS, G, rock, W, D);
fprintf(1,'done\n');

for method=1:4
    %% Upscale petrophysical data
    fprintf(1,'Upscaling ...');
    cfac = [5 5 3];
    p  = partitionUI(G, cartDims./cfac);
    switch method
        case 1  % harmonic averaging
            for i=1:3;
                crock.perm(:,i) = accumarray(p,1)./accumarray(p,1./rock.perm(:,i));
            end
      
        case 2  % harmonic-arithmetic
            for i=1:3;
                coarse = cartDims./cfac;
                dims = G.cartDims; dims(i)=coarse(i);
                qq = partitionUI(G,dims);
                K = accumarray(qq,1)./accumarray(qq,1./rock.perm(:,i));
                crock.perm(:,i) = accumarray(p,K(qq))./accumarray(p,1);
            end
        
        case 3  % local, flow-based with axial pressure drop
            CG = generateCoarseGrid(G, p);
            crock.perm = upscalePerm(G, CG, rock, 'Verbose', true);
       
        case 4  % global, flow-based with upscaling of well index
            CG = coarsenGeometry(generateCoarseGrid(G, p));
            Tf = hT2T(G,hT);
            Wc = coarsenWells(CG, W);
            Wc = addDefaultWellFields(Wc);
            [Tc, Wc] = upscaleTransGlobal(CG, Wc, Tf, ...
                'GlobalFieldCases', 'revolving', ...
            'handleNegative', 'setToZero', ...
            'fluxThreshold', sqrt(eps), 'LinSolve', @(A,b) agmg(A,b,1));
        case 5
            CG = coarsenGeometry(generateCoarseGrid(G, p));
            Tf = 1./accumarray(G.cells.faces(:,1), 1./hT, [G.faces.num, 1]);
            [~, Tc] = upscaleTransNew(CG, Tf, 'match_method', 'max_flux', ...
                'bc_method', 'bc_simple'); % 'LinSolve', @(A,b) agmg(A,b,1));
            [~,~,Wc]   = upscaleTransNew(CG, Tf, 'match_method', 'lsq_flux', ...
                'bc_method', 'wells_simple', 'wells', {W}); %'LinSolve', @(A,b) agmg(A,b,1));
    end
    crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
    fprintf(1,'done\n');

    %% Setup the coarse-scale problem
    fprintf(1,'Setting up coarse-scale problem ...');
    Gc  = cartGrid(cartDims./cfac, physDims);
    Gc  = computeGeometry(Gc);
    if method<4
        Tc = computeTrans(Gc, crock);
        Tc = 1./accumarray(Gc.cells.faces(:,1), 1./Tc);
        cwloc(1,:) = ceil(wloc(1,:)/cfac(1));
        cwloc(2,:) = ceil(wloc(2,:)/cfac(2));
        Wc = [];
        for w = 1 : numel(wtype),
            Wc = verticalWell(Wc, Gc, crock, cwloc(1,w), cwloc(2,w), ...
                1 : (cartDims(end)/cfac(end)), ...
                'Type', wtype{w}, 'Val', wtarget(w), ...
                'Radius', wrad(w), 'Name', wname{w}, ...
                'InnerProduct', 'ip_tpf');
        end
    end
    fprintf(1,'done\n');

    %% Solve coarse-scale flow problem and compute flow diagnostics
    fprintf(1,'Solving coarse-scale problem ...');
    rSc = initState(Gc, Wc, 0);
    rSc = incompTPFA(rSc, Gc, Tc, fluid, 'wells', Wc, 'use_trans', true);
    Dc  = computeTOFandTracer(rSc, Gc, crock, 'wells', Wc);
    WPc = computeWellPairs(rSc, Gc, crock, Wc, Dc);
    fprintf(1,'done\n');

    %% Compare allocation factors
    figure; colormap(.6*jet+.4*ones(size(jet)));
    for i=1:numel(D.inj)
        subplot(1,2,i)
        barh(WP.inj(i).z, cumsum(WP.inj(i).alloc,1), ...
            'stacked','BarWidth', .98, 'EdgeColor','none');
        lh=legend(W(D.prod).name,4);
        hold on
        barh(WPc.inj(i).z, -cumsum(WPc.inj(i).alloc,1),...
            'stacked', 'BarWidth', .98);
        barh(WP.inj(i).z, -cumsum(WP.inj(i).alloc,1),'stacked',...
            'BarWidth', .98, 'FaceColor','none');
        hold off, axis tight
        set(lh,'units','pixels','FontSize',8);
        title(W(D.inj(i)).name);
    end
end
