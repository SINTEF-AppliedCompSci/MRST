%% Show upscaled relative permeabilities
% In this example, we show how to generate analytically upscaled fluid
% objects from an 1D aquifer model having small-scale caprock undulations.

mrstModule add co2lab;

%% Make aquifer models
% We generate two realizations of the aquifer model, one with small-scale
% caprock undulations and one without.
gravity on
avgAquif = makeAquiferModel('A', 0);
finAquif = makeAquiferModel('A', 2);

%% Determine the parameter that characterizes small-scale undulations
% First, we find the local trap height for the model with small-scale
% undulations. Because the large-scale shape of the model is the left half
% of an antiform, the trap height can be computed as a simple maximum
% operation.
z = finAquif.Gt.cells.z;
zt = max(z) * ones(size(z));
for i = 2:numel(zt) - 1
   zt(i) = max(z(i:end));
end
zt(end) = max(zt(end - 1),z(end));

% Convolve the local trap height with a Gaussian kernel to get a smooth
% variation
ht = zt - z;
ff = exp( - linspace(-25, 25, 501).^2); ff = ff' / sum(ff);
hts = filter2(ff, ht);

% Plot the local height and its filtered values
sz = get(0, 'ScreenSize');
figure('Position', [sz(3) - 900, sz(4) - 400 - 100, 900, 400]);
xc = avgAquif.Gt.cells.centroids(:, 1) / 1e3;
plot(xc, ht, 'k-')
hold on,
plot(xc, hts, 'k-', 'LineWidth', 2)
plot([20 21.5 21.5 20 20], [0 0 2 2 0], 'r', 'LineWidth', 2);
hold off
hf = gca;

% Inlet: zoom of the surfaces defining the local trap height
i = (xc >= 20) & (xc <= 21.5);
za = avgAquif.Gt.cells.z;
hi1 = axes('Position', get(hf, 'Position') * .35 + [.1 .6 0 0]);
plot(xc(i), [za(i) z(i) zt(i)], 'LineWidth', 1);
set(gca, 'YDir', 'reverse', 'YAxisLocation', 'right', 'XTick', []);
axis tight

% Inlet: zoom of the trap height and its filtered value
hi2 = axes('Position', get(hf, 'Position') * .35 + [.1 .3 0 0]);
plot(xc(i), ht(i), 'k-')
hold on, plot(xc(i), hts(i), 'k-', 'LineWidth', 2), hold off
set(gca, 'YAxisLocation', 'right');

%% Define upscaled fluids
% We will define four different fluid models: a reference model
% corresponding to the smooth aquifer model ('smooth'), an accretion layer
% model ('inf_rough'), and two models with square and sinusoidal subscale
% caprock variations ('square' and 'sinus'). As a starting point, we make a
% simple fluid model which we later will modify with the upscaled
% properties.
Gt = avgAquif.Gt;
surf_topos = {'smooth', 'inf_rough', 'square', 'sinus'};
fluid = cell(numel(surf_topos), 1);
for i = 1:numel(surf_topos)
   fluid{i} = makeVEFluid(Gt, avgAquif.rock, 'sharp interface', ...
                          'residual', [0 0], ...
                          'dissolution', false, ...
                          'top_trap', hts, ...
                          'surf_topo', surf_topos{i});
end

%% Plot the analytically and numerically upscaled models
% We extract and plot the relative permeabilities evaluated at x=20 km. The
% fluid objects are created so that they evaluate fluid parameters for the
% whole aquifer at a time. Therefore we insert a constant CO2 saturation in
% the whole model when evaluating the relative permeability for each height
% value.
f1 = figure; hold on
H = max(Gt.cells.H);
sGv = linspace(0, 4.0 / H, 100);
krG = zeros(size(sGv));
sGa = ones(Gt.cells.num, 1);
cno = floor(20 / 30 * Gt.cells.num);
p0 = Gt.cells.z(:) * norm(gravity) * fluid{1}.rhoWS;
markers = {'k', 'g', 'm', 'r', 'b'};
for k = 1:numel(fluid)
   for j = 1:numel(sGv)
      sGa(1:end) = sGv(j);
      buf = fluid{k}.krG(sGa, p0);
      krG(j) = buf(cno);
   end
   plot(sGv * 50, krG, markers{k}, 'LineWidth', 2)
end
%%
% The numerically upscaled models are loaded from a precomputed file. If it
% does not yet exist, compute it.
savefilename = fullfile(mrstPath('co2lab'), 'examples', 'papers', 'COMG-1', 'data', 'upscaled_relperm_theta');
if (exist([savefilename, '.mat'])==2 && ...
    ask_user('Saved result found.  Re-use? [y/n] '))%#ok
   fprintf('Re-using old result.\n');
else
   fprintf('Recomputing result.\n');
   upscaleRelPerms;
end
us = load([savefilename, '.mat']);

figure(f1);
plot(us.sat_mat*50, squeeze(us.krCO2(:,:,1)),'*-');

% Set legend and finish
lnames = {'fine scale', 'accretion layer', 'square', 'sinusoidal', ...
   'numerical: \theta=0','numerical: \theta=0.0162','numerical: \theta=0.03'};
legend(lnames, 'Location', 'NorthWest', 'FontSize', 12);
axis([0 4 0 0.08]);
