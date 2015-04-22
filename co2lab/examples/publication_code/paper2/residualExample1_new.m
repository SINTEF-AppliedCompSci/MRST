%% Example 1: Compare the effects of residual and small-scale undulations
% In this example we consider a 1D antiform aquifer with a caprock given by
% the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% Injection of CO2 is first simulated using models with/without residual
% saturation for the case of a flat surface (A=0). Then, we include
% small-scale caprock undulations (A=2) and demonstrate that these will
% give a retardation effect on the plume. Finally, we show one of the
% solutions in more detail in physical space.

clear all;
moduleCheck('co2lab', 'ad-fi', 'ad-core');
gravity reset on

%% Check whether to recompute, or re-use old result
savefilename = 'data/residualExample1Data';
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

legendtext = {'No residual (A=0)', 'Residual (A=0)', ...
              'No residual (A=2)', 'Residual (A=2)'};

k = 1;
for n=1:2,
   %% Create model
   if n==1
      aquifer = makeAquiferModel_new('A', 0);
      linetype = {'b-', 'r-'}; 
      xc = aquifer.Gt.cells.centroids(:, 1) / 1e3; 
      ff = 1;
   else
      aquifer = makeAquiferModel_new('A', 2); 
      xx = xc(150:650); 
      ff = exp( -((xx - xc(400)) / (0.3)).^2); 
      ff = ff / sum(ff);
      linetype = {'b--', 'r--', 'b-', 'r-'};
      
      % print -depsc2 figs/ex1-fig1a.eps;
      set(get(gca,'Children'),'LineStyle', '--');
   end
      
   G  = aquifer.G;
   Gt = aquifer.Gt;
   
   for residual = [false, true]
            
      %% Make fluid model
      surf_temp   = 12; % in Celsius
      temp_grad   = 30 / (kilo*meter); % thermal gradient
      temperature = aquifer.Gt.cells.z * temp_grad + (274 + surf_temp);
      p_range     = [1,  150] * mega * Pascal;
      t_range     = [12, 150] + 274;
      res_vals    = [.11, .21] * residual;
      cw          = 4.3e-5 / barsa; % linear water compressibility

      fluid = makeVEFluid(Gt, aquifer.rock2D, 'sharp interface' , ...
                          'fixedT'      , temperature           , ...
                          'co2_rho_pvt' , [p_range, t_range]    , ...
                          'wat_rho_pvt' , [cw, 100 * barsa]     , ...
                          'residual'    , res_vals);

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
         t2 = tic; 
         model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid); 
         [wellSols, states] = simulateScheduleAD(state, model, schedule); 
         t2 = toc(t2); 
         xc = Gt.cells.centroids(:, 1)/1e3;
         results{k} = struct('states', {states}, 'ff', ff); 
      end
      
      %% Plot results
      state = results{k}.states{end - 70}; 
      sG = free_sg(state.s(:, 2), state.sGmax,... % @@ state.smax(:, 2) ??
           struct('res_gas', fluid.res_gas, 'res_water', fluid.res_water)); 
      hold on
      plot(xc, filter2(ff, sG .* Gt.columns.dz), linetype{k}, 'LineWidth', 2); 
      hold off
      drawnow;

      k = k + 1; 
   end
   axis tight
   set(gca, 'YDir', 'reverse', 'FontSize', 16); 
   legend(legendtext{:}, 4);
end

%% Plot solution in physical space
% To this end, we will use the second last case with no residual
% saturation, but with small-scale caprock undulations
figure
set(gcf, 'PaperPositionMode', 'auto')
p = get(gcf, 'Position'); set(gcf, 'Position', [p(1:2) 900 300]);

% Main plot
state = results{3}.states{end - 70}; 
ff = results{3}.ff; 
sG = free_sg(state.s(:, 2), state.sGmax, struct('res_gas', 0, 'res_water', 0)); 
z1 = sG .* Gt.columns.dz; 
z2 = filter2(ff, z1); 
hold on
plot(xc, z1, 'k-'); 
plot(xc, z2, 'k-', 'LineWidth', 1); 
plot([15 16 16 15 15]', [0.5 0.5 5 5 0.5], 'r', 'LineWidth', 1); 
hold off
h1 = gca; 
set(gca, 'YDir', 'reverse');

% Inlet: zoom of subscale solution and average
i = (xc >= 15) & (xc <= 16); 
axes('Position', get(h1, 'Position') * 0.3 + [0.6, 0.6, 0, 0]); 
plot(xc(i), z1(i), 'k - '); 
hold on; plot(xc(i), z2(i), 'k - ', 'LineWidth', 2); hold off
set(gca, 'YDir', 'reverse'); axis tight


% Inlet: zoom of the fluid distribution
x = xc(i); 
zt = Gt.cells.z(i); 
zc = zt + z1(i); 
h4 = axes('Position', get(h1, 'Position') * 0.3 + [0.6, 0.2, 0, 0]); 
patch([x; x(end:-1:1)], [zt; zt(end:-1:1) + 50], myCOColor(5)); 
patch([x; x(end:-1:1)], [zt; zc(end:-1:1)], myCOColor(2)); 
set(gca, 'YDir', 'reverse'); axis tight
set(gca, 'YLim', min(zt) + [0 25]); 
if recompute
   save(savefilename, 'xc', 'results')
end
% print -depsc2 figs/ex1-fig2.eps;