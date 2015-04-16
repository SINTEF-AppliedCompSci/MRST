function runStandardModel(save_filename, plot_routine, varargin)

   moduleCheck('co2lab', 'ad-fi', 'ad-core');
   gravity reset on;
   
   %% Check for presence of already-computed result; compute or re-use
   if (exist([save_filename, '.mat'], 'file') && ...
       ask_user('Saved result found.  Re-use? [y/n] '))
      fprintf('Re-using old result.\n');
      sim_outcome = load([savefilename, '.mat']);
   else
      fprintf('Recomputing result. \n');
      sim_outcome = run_standard_simulation(varargin{:});
      save(save_filename, sim_outcome);
   end
   
   %% Plot result to produce the figures
   plot_routine(sim_outcome);
   
end

% ----------------------------------------------------------------------------

function outcomes = run_standard_simulation(varargin)

   % Loop parameters
   opt.A             = 0;                   % magnitudes of subscale undulations
   opt.residual      = false;               % whether to enable residual saturation
   opt.relperm_types = {'sharp interface'}; % relperm model(s) to use
   opt.dis_types     = {'none'};            % dissol. types ('none'/'rate'/'instant')

   % Timestepping parameters
   opt.Ti   = 50   * year; % duration of injection phase
   opt.dTi  = 2    * year; % timestep size during injection
   opt.Tm   = 2000 * year; % duration of migration phase
   opt.dTm  = 20   * year; % timestep size during migration
   
   % fluid model parameters
   opt.surf_temp = 12;                        % in Celsius
   opt.temp_grad = 30 / (kilo*meter);         % thermal gradient
   opt.p_range   = [1,  150] * mega * Pascal; % fluid object supported pressure range
   opt.t_range   = [12, 150] + 274;           % fluid object supported temp. range
   opt.res_vals  = [.11, .21];                % residual saturation values (if enabled)
   opt.cw        = 4.3e-5 / barsa;            % linear water compressibility
   
   opt = merge_options(opt, varargin{:});
   
   simulation_count = 1; % global count of simulation runs
   total_count = numel(opt.A) * numel(opt.residual) * numel(opt.relperm_types) ...
       * numel(opt.dis_types);
   outcomes = cell(total_count, 1);
   
   %% Loop over grid types
   for A = opt.A

      % Constructing aquifer
      aquifer = makeAquiferModel_new('A', A);
      
      %% Loop over whether or not to use residual saturation
      for residual = opt.residual

         %% Loop over relperm models
         for relperm_type = opt.relperm_types
            
            %% Loop over dissolution types
            for dis_type = opt.dis_types

               % Make fluid model
               fluid = setup_fluid_model(opt, aquifer, residual, ...
                                              relperm_type{:}, dis_type{:});

               % Defining injection schedule
               [schedule, Winj] = setup_schedule(opt, fluid, aquifer.W, aquifer.Gt);

               % Defining initial state
               initState = setup_init_state(fluid, aquifer.G, Winj);
                  
               % Set up and run complete model
               model = CO2VEBlackOilTypeModel(aquifer.Gt, aquifer.rock2D, fluid);
               [wellSols, states] = simulateScheduleAD(initState, model, schedule);
               
               % Storing outcome
               outcomes{simulation_count} = struct('states', {states}, ...
                                                   'wellSols', wellSols);
               simulation_count = simulation_count + 1;
            end
         end
      end
   end
end

% ----------------------------------------------------------------------------

function [schedule, Winj, Wmig] = setup_schedule(opt, fluid, W, Gt)
   
   % Specify well structures
   Winj = W;
   Winj(2).val = fluid.rhoWS * Gt.cells.z(Winj(2).cells) * norm(gravity);
   Wmig = Winj;
   Wmig(1).val = 0;
   
   % Specify duration of individual injection timesteps
   istep = linspace(0.1 * year, opt.dTi, 10)'; 
   istep = [istep; ones(floor((opt.Ti - sum(istep)) / opt.dTi), 1) * opt.dTi]; 
   istep = [istep; opt.Ti - sum(istep)];
   
   % Specify duration of individual migration timesteps
   mstep = linspace(0.5 * year, opt.dTm, 5)';
   mstep = [mstep; ones(floor((opt.Tm - sum(mstep)) / opt.dTm), 1) * opt.dTm];
   mstep = [mstep; opt.Tm - sum(mstep)];
   
   % Put everything together in a schedule structure
   schedule = struct('control', [struct('W', Winj), struct('W', Wmig)]    , ...
                     'step'   , struct('control', [ones(size(istep));       ...      
                                                   2 * ones(size(mstep))] , ...
                                       'val'    , [istep; mstep]));
end

% ----------------------------------------------------------------------------

function initState = setup_init_state(fluid, G, Winj)
   
   W_val   = Winj(2).val;
   W_depth = G.cells.centroids(Winj(2).cells, 3);

   nc = G.cells.num;
   pfun = @(z) W_val + (z - W_depth) * norm(gravity) * fluid.rhoWS;
   
   initState = struct('pressure', pfun(G.cells.centroids(:, 3)), ...
                      's'       , [ones(nc, 1), zeros(nc, 1)], ...
                      'sGmax'   , zeros(nc, 1));
end

% ----------------------------------------------------------------------------

function fluid = setup_fluid_model(opt, aquifer, residual, relperm_type, dis_type)

   T = aquifer.Gt.cells.z * opt.temp_grad + (274 + opt.surf_temp); % add 274 to get Kelvin
   res_vals = opt.res_vals * residual; % becomes zero if 'residual' is false
   
   
   fluid = makeVEFluid(aquifer.Gt, aquifer.rock2D, relperm_type   , ...
                       'fixedT'      , T                          , ...
                       'co2_rho_pvt' , [opt.p_range, opt.t_range] , ...
                       'wat_rho_pvt' , [opt.cw, 100 * barsa]      , ...
                       'dissolution' , ~strcmpi(dis_type, 'none') , ...
                       'residual'    , res_vals);

   if strcmpi(dis_type, 'instant')
      fluid.dis_rate = 0; % value of zero indicates instant dissolution
   end
end
