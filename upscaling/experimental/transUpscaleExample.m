mrstModule add coarsegrid

%% Set up grid and create partition
g = computeGeometry(cartGrid([4, 1, 1]));

rock.perm = ones([g.cells.num, 1]);

% Calculate transmisibility on original grid

T = computeTrans(g, rock);

% Make partition for coarse grid
%p = partitionUI(g, floor(g.cartDims/2));

p = partitionUI(g, [g.cartDims(1), 1, 1]);
p = processPartition(g, p, 'Verbose', true);

%% Generate coarse grid and define wells
cg = generateCoarseGrid(g, p);
W  = [];

% Add four wells which will be used for upscaling in linearly independent
% combinations.
if (false)
   W = verticalWell(W, g, rock,  1,   1, [],     ...
                    'Type', 'bhp', 'Val', 300*barsa, ...
                    'Radius', 0.125, 'Name', 'P1');
   W = verticalWell(W, g, rock,  1,   g.cartDims(2), [], ...
                    'Type', 'bhp', 'Val', 300*barsa, ...
                    'Radius', 0.125, 'Name', 'P2');

   W = verticalWell(W, g, rock, 1, g.cartDims(2), [], ...
                    'Type', 'bhp', 'Val', 300*barsa, ...
                    'Radius', 0.125, 'Name', 'P3');

   W = verticalWell(W, g, rock, g.cartDims(1), g.cartDims(2), [], ...
                    'Type', 'bhp', 'Val', 300*barsa, ...
                    'Radius', 0.125, 'Name', 'P4');

   W = verticalWell(W, g, rock, ...
                    floor(g.cartDims(1) / 2) + 1, ...
                    floor(g.cartDims(2) / 2) + 1, ...
                    (1:floor(g.cartDims(3) / 2) + 1), ...
                    'Type', 'bhp', 'Val', 500*barsa, ...
                    'Radius', 0.125, 'Name', 'I4');

else

   W = verticalWell(W, g, rock,  1,   1, [],         ...
                    'Type', 'bhp', 'Val', 300*barsa, ...
                    'Radius', 0.125, 'Name', 'P1');

   W = verticalWell(W, g, rock, g.cartDims(1), g.cartDims(2), [], ...
                    'Type', 'bhp', 'Val', 500*barsa, ...
                    'Radius', 0.125, 'Name', 'P4');
end

clf
plotGrid(g, 'FaceAlpha', 0)
plotWell(g, W)

%% Upscale transmissibility and plot the result
[HT_cg, T_cg, Wcg, upscaled] = ...
   upscaleTrans(cg, T, 'match_method', 'lsq_flux', ...
                'bc_method', 'wells_simple', 'wells', W);
figure(1), clf
plot_press = ...
   @(x) plotCellData(cg.parent, ...
                     convertTo(x.pressure(cg.partition), barsa), ...
                     'EdgeColor', 'w', 'EdgeAlpha', 0.05, ...
                     'FaceAlpha', 0.85);

for i = 1 : numel(upscaled),
   subplot(numel(upscaled), 1, i)

   plot_press(upscaled(i).state); colorbar, view(3)
   outlineCoarseGrid(cg.parent, cg.partition, 'LineWidth', 3);
end

%%
% test options
% get an other well configuration

g  = cg.parent;
W2 = verticalWell([], g, rock, g.cartDims(1), g.cartDims(2), [], ...
                 'Type', 'bhp', 'Val', 300*barsa, ...
                  'Radius', 0.125, 'Name', 'P4');

W2 = verticalWell(W2, g, rock, ...
                  floor(g.cartDims(1)/2) + 1, ...
                  floor(g.cartDims(2)/2) + 1, ...
                  (1:floor(g.cartDims(3)/2) + 1), ...
                  'Type', 'bhp', 'Val', 500*barsa, ...
                  'Radius', 0.125, 'Name', 'I4');

% bc for setting outside
[pmin, pmax] = deal(100*barsa, 200*barsa);

bc    = cell([3, 1]);
bc{1} = pside([]   , cg.parent, 'West'  , pmin);
bc{1} = pside(bc{1}, cg.parent, 'East'  , pmax);
bc{2} = pside([]   , cg.parent, 'North' , pmin);
bc{2} = pside(bc{2}, cg.parent, 'South' , pmax);
bc{3} = pside([]   , cg.parent, 'Top'   , pmin);
bc{3} = pside(bc{3}, cg.parent, 'Bottom', pmax);

match_methods = {'lsq_flux', 'max_flux'};
bc_methods    = {'bc_simple', 'bc', 'wells_simple', 'wells', 'bc'   };
bcs           = {[]         , bc  , []            , []     , bc(1:2)};
wells         = {[]         , []  , {W}           , {W, W2}, []     };

%%
for m = reshape(match_methods, 1, []),
   for j = 1 : numel(bc_methods),
      fprintf('%s\n', repmat('-', [1, 50]));
      disp(['Testing ', m{1},' with ', bc_methods{j}]);

      [HT_cg, T_cg, Wcg, upscaled] = ...
         upscaleTransNew(cg, T, 'match_method', m{1}, 'bc', bcs{j}, ...
                         'bc_method', bc_methods{j}, 'wells', wells{j}, ...
                         'opt_trans_alg', 'none', 'use_trans', false);

      fprintf('Number of trans: %4d\n', numel(T_cg));
      fprintf('Number of negative trans: %d\n', sum(T_cg < 0));
   end
end
