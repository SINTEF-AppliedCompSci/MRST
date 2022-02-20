%% Upscaling a Subset of SPE 10
% In this example we illustrate how allocation factors for well pairs can
% be used to assess the quality of upscaling.  As our example, we consider
% a subsample of Model 2 from the 10th SPE Comparative Solution Project
% with a different well pattern consisting of two central injectors and
% producers at each of the four corners.

mrstModule add diagnostics spe10 coarsegrid agmg incomp
% Check for agmg
if ~exist('agmg', 'file') || ...
      norm(agmg(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) > 1.0e-8,
   error('This example requires the AGMG linear solver package');
end

%% Set up fine-scale problem
% We pick the fifteen topmost layers of the 3D model and use this as our
% test problem
fprintf(1,'Setting up fine-scale problem ...');

% Grid
cartDims = [  60,  220, 15];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m
G  = cartGrid(cartDims, physDims);
G  = computeGeometry(G);

% Rock
rock = getSPE10rock(1:cartDims(end));
rock.poro = max(rock.poro, 1e-4);

% Wells
W = [];
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500,   500  ] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,  20, 40;
              1,    1,   220,  220, 130, 90];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1', 'I2'};
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), 1 : cartDims(end), ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'InnerProduct', 'ip_tpf');
end

% Single-phase fluid
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
fprintf(1,'done\n');

%% Solve flow problem and compute flow diagnostics
fprintf(1,'Solving fine-scale problem with diagnostics ...');
rS = initState(G, W, 0);
T  = computeTrans(G, rock);
rS = incompTPFA(rS, G, T, fluid, 'wells', W, 'LinSolve', @(A,b) agmg(A,b,1));
D  = computeTOFandTracer(rS, G, rock, 'wells', W);
WP = computeWellPairs(rS, G, rock, W, D);
fprintf(1,'done\n');

%% Upscale petrophysical data
% Upscale the permeability using either simple harmonic averaging, which we
% expect will give quite poor results, or a standard flow-based method from
% the 'upscaling' module, which is expected to give reasonable results.
% Notice that the computational cost of the flow-based method may be quite
% high for large subsets of the SPE10 model. For the porosity, we use a
% simple average.
cfac      = [5 5 3];
flowbased = true;

fprintf(1,'Upscaling ...');
p  = partitionUI(G, cartDims./cfac);
if flowbased
   mrstModule add upscaling agglom coarsegrid
   CG = generateCoarseGrid(G, p);
   crock.perm = upscalePerm(G, CG, rock, 'Verbose',true);
else
   for i=1:3; %#ok<UNRCH>
      K = accumarray(p,1./rock.perm(:,i))./accumarray(p,1);
      crock.perm(:,i) = 1./K;
   end
end
crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
fprintf(1,'done\n');

%% Setup the coarse-scale problem
fprintf(1,'Setting up coarse-scale problem ...');
Gc  = cartGrid(cartDims./cfac, physDims);
Gc  = computeGeometry(Gc);
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
fprintf(1,'done\n');

%% Solve coarse-scale flow problem and compute flow diagnostics
fprintf(1,'Solving coarse-scale problem ...');
rSc = initState(Gc, Wc, 0);
Tc  = computeTrans(Gc, crock);
rSc = incompTPFA(rSc, Gc, Tc, fluid, 'wells', Wc);
Dc  = computeTOFandTracer(rSc, Gc, crock, 'wells', Wc);
WPc = computeWellPairs(rSc, Gc, crock, Wc, Dc);
fprintf(1,'done\n');

%% Compare allocation factors
% We contrast the allocation factors for the injection wells computed on
% the fine and the coarse model. In the plot, the colored bars represent
% allocation computed by the upscaled model, whereas the line bars come
% from the fine-scale model. Since the bars represent the cumulative
% influx from toe to heel of the well, similar to the plot from a
% production logging tool (PLT), the bars should match at the top of each
% coarse perforation if the upscaled model is able to reproduce the
% connections and flow patterns predicted by the fine-scale model
figure; set(gcf,'Position',[700 50 660 760]);
plotWellAllocationComparison(Dc,WPc,D,WP);

%% Copyright notice

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
% 
