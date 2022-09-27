%% Upscaling versus MsRSB: the Tarbert model
% This script is a continuation of upscalingExample2.m from the book
% module. We compare the MsRSB method with three upscaling methods for
% simulation of a two injectors and four producers in a subset of the SPE
% 10 benchmark. To assess the quality of the different methods, we contrast
% the allocation factors for the injection wells computed on the fine and
% the coarse model. Ideally, bars on the negative axis that represent the
% allocation factors for the coarse model should be the mirror of the bars
% on the positive axis that represent the allocation factors for the fine
% model. To simplify the comparison, the fine-scale allocation factors are
% indicated by lines on top of those of the coarse scale.

mrstModule add agmg linearsolvers book coarsegrid diagnostics 
mrstModule add incomp libgeometry msrsb spe10 upscaling 

%% Setup iterative linear solver
% When working with large models, one cannot use the standard MLDIVIDE
% ('\') solver in MATLAB. Here, we use the AGMG or AMGCL algebraic
% multigrid solver.
if ~norm(callAMGCL(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) < 1.0e-8
    linsolver = @(A,b) callAMGCL(A,b);
    disp('Using AMGCL as linear solver');
elseif exist('agmg', 'file') && ...
      norm(agmg(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) < 1.0e-8
    linsolver = @(A,b) agmg(A,b,1);
    disp('Using AGMG as linear solver');
else
    warning(['Cannot use an algebraic multigrid solver. ' ...
        'Trying to proceed with the default direct solver instead']);
    linsolver = @(A,b) A\b;
end

%% Set up fine-scale problem
fprintf(1,'Setting up fine-scale problem ...');
cartDims  = [  60,  220, 24];  % Number of z-layers must be multiple of 3.
physDims  = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
rock      = getSPE10rock([1 1:cartDims(end)-1]);
rock.poro = max(rock.poro, 1e-4);
G  = cartGrid(cartDims, physDims);
G  = mcomputeGeometry(G);
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500,   500  ] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
offset = 4;
wloc     = [  1+offset,   60-offset,     1+offset,   60-offset,  20+offset, 40+offset;
              1+offset,    1+offset,   220-offset,  220-offset, 130-offset, 90-offset];          
%wloc     = [  1,   60,     1,   60,  20, 40;
%              1,    1,   220,  220, 130, 90];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1', 'I2'};
W = [];
for w = 1 : numel(wtype)
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'InnerProduct', 'ip_tpf');
end
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
fprintf(1,'done\n');

%% Plot the setup
figure
Kx = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G,log10(Kx),'EdgeColor','none');
mrstColorbar(Kx,'south',true);
plotWell(G,W);
view(-75,20); set(gca,'DataAsp',[15 15 1]); axis tight off
colormap(colormap.^1.5)

%% Solve flow problem and compute flow diagnostics
fprintf(1,'Solving fine-scale problem ...');
rS = initState(G, W, 0);
hT = computeTrans(G, rock);
rS = incompTPFA(rS, G, hT, fluid, 'wells', W, 'LinSolve', @(A,b) linsolver(A,b));
D  = computeTOFandTracer(rS, G, rock, 'wells', W);
WP = computeWellPairs(rS, G, rock, W, D);
fprintf(1,'done\n');

amax = -inf;
cfac = [10 10 3];
Gc   = cartGrid(cartDims./cfac, physDims);
Gc   = computeGeometry(Gc);

pos = [1 3 5 2 4 6];
for method=1:4
   %% Upscale petrophysical data
   fprintf(1,'Upscaling ...');
   p  = partitionUI(G, cartDims./cfac);
   switch method
      case 1  % harmonic-arithmetic
         tittel = 'Harmonic-arithmetic';
         crock.perm = zeros(prod(cartDims./cfac),3);
         for i=1:3
            coarse = cartDims./cfac;
            dims = G.cartDims; dims(i)=coarse(i);
            qq = partitionUI(G,dims);
            K = accumarray(qq,1)./accumarray(qq,1./rock.perm(:,i));
            crock.perm(:,i) = accumarray(p,K(qq))./accumarray(p,1);
         end
        
      case 2  % local, flow-based with axial pressure drop
         tittel = 'Flow-based, pressure drop';
         CG = generateCoarseGrid(G, p);
         crock.perm = upscalePerm(G, CG, rock, 'Verbose', true);
       
      case 3  % global, flow-based with upscaling of well index
         tittel = 'Global, flow-based';
         CG = coarsenGeometry(generateCoarseGrid(G, p));
         Tf = hT2T(G,hT);
         Wc = coarsenWells(CG, W);
         Wc = addDefaultWellFields(Wc);
         [Tc, Wc] = upscaleTransGlobal(CG, Wc, Tf, ...
            'GlobalFieldCases', 'revolving', ...
            'handleNegative', 'ignore', ...
            'fluxThreshold', sqrt(eps), 'LinSolve', @(A,b) agmg(A,b,1));
         Gc = CG;

       case 4 % multiscale method with C-accelerated backend
           tittel = 'Multiscale: MsRSB only';
           A     = getIncomp1PhMatrix(G, hT, rS, fluid);
           CG    = generateCoarseGrid(G, p);
           CG    = coarsenGeometry(CG);
           CG    = storeInteractionRegionCart(CG);
           CG    = setupMexInteractionMapping(CG);
           basis = getMultiscaleBasis(CG, A, 'type', 'msrsb', ...
               'useMex',true, 'tolerance',1e-6,'iterations', 1000);
           Gc    = CG;
   end
   crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
   fprintf(1,'done\n');

   %% Setup the coarse-scale problem
   fprintf(1,'Setting up coarse-scale problem ...');
   if method<3
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
            'InnerProduct', 'ip_tpf');
      end
   end
   fprintf(1,'done\n');

   %% Solve coarse-scale flow problem and compute flow diagnostics
   if method<4
       fprintf(1,'Solving coarse-scale problem and computing flow diagnostics...');
       rSc = initState(Gc, Wc, 0);
       rSc = incompTPFA(rSc, Gc, Tc, fluid, 'wells', Wc, 'use_trans', true);
       Dc  = computeTOFandTracer(rSc, Gc, crock, 'wells', Wc);
       WPc = computeWellPairs(rSc, Gc, crock, Wc, Dc);
       fprintf(1,'done\n');
   else
       fprintf(1,'Solving multiscale problem and computing flow diagnostics...');
       rSc = initState(G, W, 0);
       rSc = incompMultiscale(rSc, Gc, hT, fluid, basis, 'wells', W);
       Dc  = computeTOFandTracer(rSc, G, rock, 'wells', W);
       WPc = computeWellPairs(rSc, G, rock, W, Dc);
       fprintf(1,'done\n');
   end
   
   %% Compare allocation factors
   subplot(3,2,pos(method),'FontSize',12);
   amax=myPlotAllocation(WP, WPc, {W(D.prod).name}, amax);
   title(tittel,'FontWeight','normal'); drawnow
end


%% Check iterative version
fn   = getSmootherFunction('type', 'ilu0');
rSc0 = initState(G, W, 0);

subplot(3,2,pos(5));
fprintf(1,'Solving iterative multiscale problem and computing flow diagnostics ..');
rSc = incompMultiscale(rSc, Gc, hT, fluid, basis, 'wells', W, ...
    'getSmoother',fn,'iterations', 3, 'useGMRES', true);
Di   = computeTOFandTracer(rSc, G, rock, 'wells', W);
WPi  = computeWellPairs(rSc, G, rock, W, Di);
amax = myPlotAllocation(WP, WPi, {W(D.prod).name}, amax);
title('Multiscale: 3 iterations','FontWeight','normal'); drawnow
fprintf(1,'done\n');

subplot(3,2,pos(6));
fprintf(1,'Solving iterative multiscale problem and computing flow diagnostics ..');
rSc = incompMultiscale(rSc, Gc, hT, fluid, basis, 'wells', W, ...
    'getSmoother',fn,'iterations', 7, 'useGMRES', true);
Di   = computeTOFandTracer(rSc, G, rock, 'wells', W);
WPi  = computeWellPairs(rSc, G, rock, W, Di);
amax = myPlotAllocation(WP, WPi, {W(D.prod).name}, amax);
title('Multiscale: 7 iterations','FontWeight','normal');
fprintf(1,'done\n');

for i=1:6, subplot(3,2,i), set(gca,'XLim',[0 amax]); end


%% Plot model
K = convertTo(rock.perm(:,1),milli*darcy);
figure, plotCellData(G,log10(K),'EdgeColor','none'); view(3);
set(gca,'DataAspectRatio',[15 15 1])
plotWell(G,W,'FontSize',12);
outlineCoarseGrid(G, p, 'EdgeColor','k');
mrstColorbar(K,'West',true);
axis tight off
colormap(.75*parula(32).^2+.25);

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
