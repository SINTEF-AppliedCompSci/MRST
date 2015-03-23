%% Example x: Compare the effects small-scale undulations
% In this example we consider a 1D antiform aquifer with a caprock given by
% the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% Injection of CO2 is first simulated using models with/without residual
% saturation for the case of a flat surface (A=0). Then, we include
% small-scale caprock undulations (A=2) and demonstrate that these will
% give a retardation effect on the plume. Finally, we show one of the
% solutions in more detail in physical space.
try
   require co2lab ad-fi
catch %#ok<CTCH>
   mrstModule add co2lab ad-fi
end
gravity reset on

%% Time steps for injection and migration
Ti  =   50*year; 
dTi =  2*year;
istep = linspace(0.1*year, dTi, 10)';
istep = [istep; ones(floor((Ti-sum(istep))/dTi), 1)*dTi];
istep = [istep; Ti-sum(istep)];

Tm  = 2000*year; 
dTm = 20*year;
mstep = linspace(0.5*year, dTm, 5)';
mstep = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep = [mstep; Tm-sum(mstep)];

legendtext = {'Fine scale', 'Accretion layer', ...
              'Analytic: square', 'Analytic: sinus'};
results = cell(4,1);

%% Create model
upAquifer = makeAquiferModel('A',0);
fAquifer  = makeAquiferModel('A',2);

z  = fAquifer.Gt.cells.z;
zt = max(z)*ones(size(z));
for i=2:numel(zt)-1
   zt(i)=max(z(i:end));
end
zt(end) = max(zt(end-1),z(end));
ht  = zt - z;
ff = exp(-linspace(-25,25,501).^2); ff=ff'/sum(ff);
hts = filter2(ff, ht);

%% Main loop
figure; hold on
k = 1;
linetype = {'k-', 'b-','r-','g-'};
surf_topos = {'smooth','inf_rough','square','sinus'};
for i=1:numel(surf_topos)
      
   %% Make fluid model
   if strcmp(surf_topos{i},'smooth')
      aquifer = fAquifer;
      fluid = makeFluidModel(aquifer, 'residual', true, ...
         'dissolution', false, 'fluidType', 'sharp interface');
   else
      aquifer = upAquifer;
      fluid = makeFluidModel(upAquifer, 'residual', true, ...
         'dissolution', false, 'fluidType', 'sharp interface', ...
         'top_trap', hts, 'surf_topo',  surf_topos{i});
   end
   
   %% Setup system
   s = setupSimCompVe(aquifer.Gt, aquifer.rock2D);
   systemOG = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
      fluid, 'simComponents', s, 'VE', true);
   systemOG.nonlinear.linesearch    = false;
   systemOG.nonlinear.maxIterations = 10;
   systemOG.nonlinear.tol           = 1e-6;
      
   %% Create well schedule
   Gt = aquifer.Gt;
   z  = aquifer.G.cells.centroids(:,3);
   clear x0;
   W  = aquifer.W;
   W(2).val = fluid.rhoOS * Gt.cells.z(W(2).cells)*norm(gravity);
   x0.pressure = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
   x0.s(:,1)   = ones(Gt.cells.num,1);
   x0.s(:,2)   = zeros(Gt.cells.num,1);
   x0.rs       = ones(Gt.cells.num,1)*0.0;
   x0.smax     = x0.s;
   x0.smin     = x0.s;
   x0.sGmax    = x0.s(:,2);
   control = struct('W',[],'step',struct('val',[],'control',[]));
   control.W = {W, W(2)};
   control.step.val = [istep; mstep];
   control.step.control = [ones(size(istep));ones(size(mstep))*2];
      
   %% Run the schedule setup
   t2=tic;
   [wellSols, states] = runMrstADI(x0, Gt, systemOG, control, ...
      'force_step', false, 'dt_min', 0.5*year, 'report_all', false);
   t2=toc(t2);
   xc=Gt.cells.centroids(:,1)/1e3;
      
   %% Plot results
   state = states{end-70};
   sG = free_sg(state.s(:,2),state.smax(:,2), ...
      struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
   hold on
   plot(xc, filter2(ff,sG.*Gt.columns.dz), linetype{k}, 'LineWidth', 2);
   hold off
   drawnow;
      
   results{k}=struct('states',{states}, 'ff', ff);
   k = k+1;
end
axis tight
set(gca,'YDir','reverse','FontSize',12);
h=legend(legendtext{:}, 4); set(h,'FontSize',14);
set(gcf,'Position', [680 580 800 420]);
save('data/upscalingExample1Data','xc','results','control')
% this data can be plotted later with showUpscalingExample1.m