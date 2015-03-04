%% Example 1: Compare the effects of residual and small-scale undulations
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

legendtext = {'No residual (A=0)', 'Residual (A=0)', ...
              'No residual (A=2)', 'Residual (A=2)'};
results = cell(4,1);
for n=1:2,
   %% Create model
   if n==1
      aquifer = makeAquiferModel('A',0);
      linetype = {'b-', 'r-'}; k =1;
      xc = aquifer.Gt.cells.centroids(:,1)/1e3;
      ff = 1;
   else
      aquifer = makeAquiferModel('A',2);
      xx = xc(150:650);
      ff = exp(-((xx-xc(400))/(0.3)).^2);
      ff = ff/sum(ff);
      linetype = {'b--', 'r--', 'b-', 'r-'};
      
      % print -depsc2 figs/ex1-fig1a.eps;
      set(get(gca,'Children'),'LineStyle', '--');
   end
      
   G  = aquifer.G;
   Gt = aquifer.Gt;
   
   for residual= [false,true]
      
      %% Make fluid model
      fluid = makeFluidModel(aquifer, 'residual', residual, ...
         'dissolution', false, 'fluidType', 'sharp interface');
      
      %% Setup system
      s = setupSimCompVe(aquifer.Gt, aquifer.rock2D);
      systemOG = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
         fluid, 'simComponents', s, 'VE', true);
      systemOG.nonlinear.linesearch    = false;
      systemOG.nonlinear.maxIterations = 10;
      systemOG.nonlinear.tol           = 1e-6;
      
      %% Create well schedule
      z  = G.cells.centroids(:,3);
      clear x0;
      W  = aquifer.W;
      W(2).val = fluid.rhoOS * Gt.cells.z(W(2).cells)*norm(gravity);
      x0.pressure = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
      x0.s(:,1)   = ones(G.cells.num,1);
      x0.s(:,2)   = zeros(G.cells.num,1);
      x0.rs       = ones(G.cells.num,1)*0.0;
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
         struct('res_gas',fluid.res_gas, 'res_oil', fluid.res_oil));
      hold on
      plot(xc, filter2(ff,sG.*Gt.columns.dz), linetype{k}, 'LineWidth', 2);
      hold off
      drawnow;
      
      results{k}=struct('states',{states}, 'ff', ff);
      k = k+1;
   end
   axis tight
   set(gca,'YDir','reverse','FontSize',16);
   legend(legendtext{:}, 4);
end
%print -depsc2 figs/ex1-fig1b.eps;
%save('ex1D-1','results','xc');

%% Plot solution in physical space
% To this end, we will use the second last case with no residual
% saturation, but with small-scale caprock undulations
figure
set(gcf,'PaperPositionMode','auto')
p = get(gcf,'Position'); set(gcf,'Position',[p(1:2) 900 300]);

% Main plot
state = results{3}.states{end-70};
ff = results{3}.ff;
sG = free_sg(state.s(:,2),state.smax(:,2), struct('res_gas',0,'res_oil',0));
z1 = sG.*Gt.columns.dz;
z2 = filter2(ff,z1);
hold on
plot(xc,z1,'k-'); 
plot(xc,z2,'k-','LineWidth',1);
plot([15 16 16 15 15]',[0.5 0.5 5 5 0.5],'r','LineWidth',1);
hold off
h1 = gca;
set(gca,'YDir','reverse');

% Inlet: zoom of subscale solution and average
i = (xc>=15) & (xc<=16);
axes('Position',get(h1,'Position')*0.3+[0.6,0.6,0, 0]);
plot(xc(i), z1(i),'k-');
hold on; plot(xc(i), z2(i),'k-','LineWidth',2); hold off
set(gca,'YDir','reverse'); axis tight


% Inlet: zoom of the fluid distribution
x  = xc(i);
zt = Gt.cells.z(i);
zc = zt + z1(i);
h4 = axes('Position',get(h1,'Position')*0.3+[0.6,0.2,0, 0]);
patch([x; x(end:-1:1)], [zt; zt(end:-1:1)+50], myCOColor(5));
patch([x; x(end:-1:1)], [zt; zc(end:-1:1)], myCOColor(2));
set(gca,'YDir','reverse'); axis tight
set(gca,'YLim',min(zt)+[0 25]);
save('data/residualExample1Data','xc','results','control')
% print -depsc2 figs/ex1-fig2.eps;