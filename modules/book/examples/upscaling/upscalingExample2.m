%% Example 2: Upscaling of the SPE10 Data set
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
mrstModule add diagnostics spe10 coarsegrid upscaling incomp libgeometry linearsolvers

try
    callAMGCL(speye(5), ones(5,1));
    linsolver = @(A,b) callAMGCL(A,b);
catch
   warning(['Cannot use the algebraic multigrid solver from AMGCL.' ...
       ' Trying to proceed with the default backslash solver']);
   linsolver = @(A,b) A\b;
end

fprintf(1,'Setting up fine-scale problem ...');
cartDims  = [  60,  220, 36];
physDims  = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
rock      = getSPE10rock([1 1:cartDims(end)-1]);
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
for w = 1 : numel(wtype)
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'InnerProduct', 'ip_tpf', 'compi', 1);
end
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
fprintf(1,'done\n');

%% Solve flow problem and compute flow diagnostics
fprintf(1,'Solving fine-scale problem ...');
rS = initState(G, W, 0);
hT = computeTrans(G, rock);
rS = incompTPFA(rS, G, hT, fluid, 'wells', W, 'LinSolve', linsolver);
D  = computeTOFandTracer(rS, G, rock, 'wells', W);
WP = computeWellPairs(rS, G, rock, W, D);
fprintf(1,'done\n');

amax = -inf;
cfac = [10 10 3];
Gc   = cartGrid(cartDims./cfac, physDims);
Gc   = computeGeometry(Gc);

for method=1:4
   %% Upscale petrophysical data
   fprintf(1,'Upscaling ...');
   p  = partitionUI(G, cartDims./cfac);
   switch method
      case 1  % harmonic averaging
         tittel = 'Harmonic';
         for i=1:3
            crock.perm(:,i) = accumarray(p,1)./accumarray(p,1./rock.perm(:,i));
         end
      
      case 2  % harmonic-arithmetic
         tittel = 'Harmonic-arithmetic';
         for i=1:3
            coarse = cartDims./cfac;
            dims = G.cartDims; dims(i)=coarse(i);
            qq = partitionUI(G,dims);
            K = accumarray(qq,1)./accumarray(qq,1./rock.perm(:,i));
            crock.perm(:,i) = accumarray(p,K(qq))./accumarray(p,1);
         end
        
      case 3  % local, flow-based with axial pressure drop
         tittel = 'Flow-based, pressure drop';
         CG = generateCoarseGrid(G, p);
         crock.perm = upscalePerm(G, CG, rock, 'Verbose', true);
       
      case 4  % global, flow-based with upscaling of well index
         tittel = 'Global, flow-based';
         CG = coarsenGeometry(generateCoarseGrid(G, p));
         Tf = hT2T(G,hT);
         Wc = coarsenWells(CG, W);
         Wc = addDefaultWellFields(Wc);
         [Tc, Wc] = upscaleTransGlobal(CG, Wc, Tf, ...
            'GlobalFieldCases', 'revolving', ...
            'handleNegative', 'ignore', ...
            'fluxThreshold', sqrt(eps), 'LinSolve', linsolver);
         Gc = CG;
   end
   crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
   fprintf(1,'done\n');

   %% Setup the coarse-scale problem
   fprintf(1,'Setting up coarse-scale problem ...');
   if method<4
      Tc = computeTrans(Gc, crock);
      Tc = 1./accumarray(Gc.cells.faces(:,1), 1./Tc);
      cwloc(1,:) = ceil(wloc(1,:)/cfac(1));
      cwloc(2,:) = ceil(wloc(2,:)/cfac(2));
      Wc = [];
      for w = 1 : numel(wtype)
         Wc = verticalWell(Wc, Gc, crock, cwloc(1,w), cwloc(2,w), ...
            1 : (cartDims(end)/cfac(end)), ...
            'Type', wtype{w}, 'Val', wtarget(w), ...
            'Radius', wrad(w), 'Name', wname{w}, ...
            'InnerProduct', 'ip_tpf', 'compi', 1);
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
   subplot(2,2,method); colormap(.6*jet+.4*ones(size(jet)));
   cwp = cumsum(WPc.inj(2).alloc,1); amax = max([sum(cwp,2); amax]);
   h=barh(WPc.inj(2).z, cwp, 'stacked', 'BarWidth', .98, 'EdgeColor','none');
   hold on
   cwp = cumsum(WP.inj(2).alloc,1);  amax = max([sum(cwp,2); amax]);
   if size(cwp,1)<20
      barh(WP.inj(2).z, cwp, 'stacked','BarWidth', .98, 'FaceColor','none');
   else
      c = cumsum(cwp,2);
      z =  WP.inj(2).z; dz = min(diff(z));
      n = size(c,2);
      stairs([zeros(1,n); c], repmat(z([1:end end]),1,n), '-k','LineWidth',1);
   end
   hold off, axis tight
   lh=legend(W(D.prod).name,'Location','SouthEast');
   set(lh,'units','pixels','FontSize',8);
   title(tittel);
end
for i=1:4, subplot(2,2,i), set(gca,'XLim',[0 amax]); end
