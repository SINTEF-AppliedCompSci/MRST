%% Example x: Compare the effects small-scale undulations
% In this example we consider a 1D antiform aquifer with a caprock given by
% the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% Injection of CO2 is first simulated using models with/without residual
% saturation for the case of a flat surface (A=0). Then, we include
% small-scale caprock undulations (A=2) and demonstrate that these will
% give a retardation effect on the plume. Finally, we show one of the
% solutions in more detail in physical space.

clear all;
moduleCheck('co2lab', 'ad-fi', 'ad-core')
gravity reset on

%% Check whether to recompute, or re-use old result
savefilename = 'data/upscalingExample1Data';
if (exist([savefilename, '.mat'])==2 && ...
    ask_user('Saved result found.  Re-use? [y/n] '))
   fprintf('Re-using old result.\n');
   recompute = false;
   load([savefilename, '.mat']);
else
   fprintf('Recomputing result.\n');
   results = cell(4,1);     
   recompute = true;
end

%% Time steps for injection and migration
Ti = 50 * year; 
dTi = 2 * year; 
istep = linspace(0.1 * year, dTi, 10)'; 
istep = [istep; ones(floor((Ti - sum(istep)) / dTi), 1) * dTi]; 
istep = [istep; Ti - sum(istep)]; 

Tm = 2000 * year; 
dTm = 20 * year; 
mstep = linspace(0.5 * year, dTm, 5)'; 
mstep = [mstep; ones(floor((Tm - sum(mstep)) / dTm), 1) * dTm]; 
mstep = [mstep; Tm - sum(mstep)]; 

legendtext = {'Fine scale', 'Accretion layer',...
              'Analytic: square', 'Analytic: sinus'};

%% Create model
upAquifer = makeAquiferModel_new('A', 0); 
fAquifer = makeAquiferModel_new('A', 2); 

z = fAquifer.Gt.cells.z; 
zt = max(z) * ones(size(z)); 
for i = 2:numel(zt) - 1
   zt(i) = max(z(i:end)); 
end
zt(end) = max(zt(end -1), z(end)); 
ht = zt - z; 
ff = exp(-linspace( -25, 25, 501).^2); ff = ff' / sum(ff); 
hts = filter2(ff, ht);

%% Main loop
figure; hold on
k = 1;
linetype = {'k-', 'b-','r-','g-'};
surf_topos = {'smooth', 'inf_rough', 'square', 'sinus'}; 
res_vals   = [.11, .21];

for i = 1:numel(surf_topos)
      
   %% Make fluid model
   surf_temp   = 12; % in Celsius
   temp_grad   = 30 / (kilo*meter); % thermal gradient
   p_range     = [1,  150] * mega * Pascal;
   t_range     = [12, 150] + 274;
   cw          = 4.3e-5 / barsa; % linear water compressibility

   if strcmp(surf_topos{i},'smooth')
      aquifer = fAquifer;
      temperature = aquifer.Gt.cells.z * temp_grad + (274 + surf_temp);

      fluid = makeVEFluid(aquifer.Gt, aquifer.rock2D, 'sharp interface' , ...
                          'fixedT'      , temperature           , ...
                          'co2_rho_pvt' , [p_range, t_range]    , ...
                          'wat_rho_pvt' , [cw, 100 * barsa]     , ...
                          'residual'    , res_vals);
   else
      aquifer = upAquifer;
      temperature = aquifer.Gt.cells.z * temp_grad + (274 + surf_temp);

      fluid = makeVEFluid(aquifer.Gt, aquifer.rock2D, 'sharp interface' , ...
                          'fixedT'      , temperature           , ...
                          'co2_rho_pvt' , [p_range, t_range]    , ...
                          'wat_rho_pvt' , [cw, 100 * barsa]     , ...
                          'residual'    , res_vals              , ...
                          'top_trap'    , hts                   , ...
                          'surf_topo'   , surf_topos{i});
   end
   Gt = aquifer.Gt;
   G = aquifer.G;

   %% Create well schedule and initial state object
   z  = G.cells.centroids(:,3);
   % Setting up wells (separate versions for injection and migration phase)
   Winj = aquifer.W;
   Winj(2).val = fluid.rhoWS * Gt.cells.z(Winj(2).cells)*norm(gravity);
   Wmig = Winj;
   Wmig(1).val = 0;

   % Initializing initial state object
   clear state;
   state.pressure = Winj(2).val +(z(:)-z(Winj(2).cells))*norm(gravity)*fluid.rhoWS;
   state.s = [ones(G.cells.num, 1), zeros(G.cells.num,1)];
   state.sGmax = state.s(:,2);

   % Defining schedule
   schedule.control = [struct('W', Winj), struct('W', Wmig)];
   schedule.step    = struct('control', [ones(size(istep)); 2 * ones(size(mstep))], ...
                             'val', [istep; mstep]);
   
   %% Run the schedule setup, if requested
   if recompute
      t2=tic;
      model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid);
      [wellSols, states] = simulateScheduleAD(state, model, schedule); 
      t2 = toc(t2); 
      xc = Gt.cells.centroids(:, 1)/1e3;
      results{k} = struct('states', {states}, 'ff', ff); 
   end
   
   %% Plot results
   state = results{k}.states{end - 70}; 
   sG = free_sg(state.s(:, 2), state.sGmax,...
                struct('res_gas', fluid.res_gas, 'res_water', fluid.res_water)); 
   hold on
   plot(xc, filter2(ff, sG .* Gt.columns.dz), linetype{k}, 'LineWidth', 2); 
   hold off
   drawnow; 
   
   results{k} = struct('states', {states}, 'ff', ff); 
   k = k+1;
end
axis tight
set(gca, 'YDir', 'reverse', 'FontSize', 12); 
h = legend(legendtext{:}, 4); set(h, 'FontSize', 14); 
set(gcf, 'Position', [680 580 800 420]); 
if recompute
   save(savefilename, 'xc', 'results')
end

% this data can be plotted later with showUpscalingExample1.m