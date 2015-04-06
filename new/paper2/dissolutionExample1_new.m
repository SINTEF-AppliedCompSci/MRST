%% Example x: Compare the effects small-scale undulations
% In this example we consider a 1D antiform aquifer with a caprock given by
% the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% Injection of CO2 is simulate using models with/without residual
% saturation. This is done for the case of a flat surface (A=0) and with
% the case with small-scale caprock undulations (A=2). For all this cases
% we simulate with and without disolution effects.  
% We show that the disolution of the residual saturation prevent the disolution 
% of the free CO2 left in the small traps in our model which first disolve the
% residual sone. Finaly we show the distribution of the CO2 in real space.
% The data from the simulations are also stored and can be inspected the
% figures can be generate by showDissolutionExample1.m

clear all;
moduleCheck('co2lab', 'ad-fi', 'ad-core');

gravity reset on
do_print = true; mkdir('figs')
depth = 1300;

%% Check whether to recompute or re-use old result
savefilename = ['data/disolutionExample1Data_',num2str(depth)];
if (exist([savefilename, '.mat'])==2 && ...
    ask_user('Saved result found.  Re-use? [y/n] '))
   fprintf('Re-using old result.\n');
   recompute = false;
   load([savefilename, '.mat']);
else
   fprintf('Recomputing result.\n');
   results = {};
   recompute = true;
end

%% Time steps for injection and migration
Ti    = 50 * year; 
dTi   = 2 * year; 
istep = linspace(0.1 * year, dTi, 10)'; 
istep = [istep; ones(floor((Ti - sum(istep)) / dTi), 1) * dTi]; 
istep = [istep; Ti - sum(istep)]; 

Tm    = 2000 * year; 
dTm   = 20 * year; 
mstep = linspace(0.5 * year, dTm, 5)'; 
mstep = [mstep; ones(floor((Tm - sum(mstep)) / dTm), 1) * dTm]; 
mstep = [mstep; Tm - sum(mstep)];

legendtext = {'No dissolution (A=0)' , ...
              'No dissolution (A=2)' , ...
              'Dissolution (A=0)'    , ...
              'Dissolution (A=2)'};

linetype = {'b--', 'r--', 'b-', 'r-'};

for residual= true; %@@[false,true] %residual saturation or not
    figure(),clf;
    k = 1;
    for n=1:2, % flat or non flat topsurface
        %% Create model
        if n==1
           aquifer = makeAquiferModel_new('A', 0, 'D', depth); 
           xc = aquifer.Gt.cells.centroids(:, 1) / 1e3; 
           ff = 1;
        else
           aquifer = makeAquiferModel_new('A', 2, 'D', depth); 
           xx = xc(150:650); 
           ff = exp( - ((xx - xc(400)) / (0.3)).^2); 
           ff = ff / sum(ff); 
           set(get(gca,'Children'),'LineStyle', '--');
        end
        
        G  = aquifer.G;
        Gt = aquifer.Gt;
        
        %% Make fluid model

        surf_temp   = 12; % in Celsius
        temp_grad   = 30 / (kilo*meter); % thermal gradient
        p_range     = [1,  150] * mega * Pascal;
        t_range     = [12, 150] + 274;
        res_vals    = [.11, .21] * residual;
        cw          = 4.3e-5 / barsa; % linear water compressibility
        temperature = Gt.cells.z * temp_grad + (274 + surf_temp);
        
        for dissolution = true%@@[false,true]

           fluid = makeVEFluid(aquifer.Gt, aquifer.rock2D, 'sharp interface' , ...
                               'fixedT'      , temperature        , ...
                               'co2_rho_pvt' , [p_range, t_range] , ...
                               'wat_rho_pvt' , [cw, 100 * barsa]  , ...
                               'dissolution' , dissolution        , ...
                               'residual'    , res_vals);
           %'dis_rate'    , 0, ... % @@@@
            
           % p_range(2) = 100000000 * 0.99;
           % fluid = makeFluidModel(aquifer, 'residual', residual, ...
           %                        'dissolution', dissolution, 'fluidType', 'sharp interface');
           % fluid.dis_rate = 0;
           
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
           state.s = [ones(Gt.cells.num, 1), zeros(Gt.cells.num,1)];
           state.sGmax = state.s(:,2);
           state.rs = zeros(Gt.cells.num, 1);
           
           % Defining schedule
           schedule.control = [struct('W', Winj), struct('W', Wmig)];
           schedule.step    = struct('control', [ones(size(istep)); 2 * ones(size(mstep))], ...
                                     'val', [istep; mstep]);

           %% Run the schedule setup, if requested
           if recompute
              t2 = tic; 
              model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid, ...
                                             'minimumPressure', p_range(1), ...
                                             'maximumPressure', p_range(2)); 
              
              % load('debugstate');
              % state = debugstate;
              % schedule.step.control = schedule.step.control(41:end);
              % schedule.step.val     = schedule.step.val(41:end);
              
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
           
           k = k +1;
        end
    end
    axis tight
    set(gca,'YDir','reverse','FontSize',16);
    legend(legendtext{:}, 4);
    if(do_print)
        if(~residual)
            print -depsc2 figs/ex1-fig3a.eps;
        else
            print -depsc2 figs/ex1-fig3b.eps;
        end
    end
end

if recompute
   save(savefilename,'xc','results');
end
