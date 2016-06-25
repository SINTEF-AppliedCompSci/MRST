function sim_outcome = runStandardModel(save_filename, plot_routine, varargin)
% Script used to generate several of the published examples where different
% models are compared againt each other.  Arguments are:
%
% 'save_filename' - name of file where simulation results will be
%                   stored/retrieved
%
% 'plot_routine'  - handle to a function that takes a cell array of
%                   simulation results as argument, and does the
%                   corresponding analysis / generates plots.
%
% 'varargin'      - Various arguments that specify the specifics of the
%                   simulation(s) to be carried out, provided as 'key'/value
%                   pairs.  See documentation of options (fields of the 'opt'
%                   structure) in the function 'run_standard_simulation' below.

   moduleCheck('co2lab', 'ad-core');
   gravity reset on;

   % Check for presence of already-computed result; compute or re-use
   if (exist([save_filename, '.mat'], 'file') && ...
       ask_user('Saved result found.  Re-use? [y/n] '))
      fprintf('Re-using old result.\n');
      loaded_data = load([save_filename, '.mat']);
      sim_outcome = loaded_data.sim_outcome;
   else
      fprintf('Recomputing result. \n');
      sim_outcome = run_standard_simulation(varargin{:});
      ensure_path_exists(save_filename);
      save(save_filename, 'sim_outcome');
   end

   % Plot result to produce the figures
   plot_routine(sim_outcome);

end

% ----------------------------------------------------------------------------

function outcomes = run_standard_simulation(varargin)

   % Loop parameters
   opt.A              = 0;                   % magnitudes of subscale undulations
   opt.depth          = 2300;                % depth of aquifer
   opt.residual       = false;               % whether to enable residual saturation
   opt.subscale_types = {'smooth'};          % subscale geometry model(s) to use
   opt.dis_types      = {'none'};            % dissol. types ('none'/'rate'/'instant')
   opt.fluid_types    = {'sharp interface'}; % Fluid model to use (for
                                             % upscaled relperm/capillarity)
   % Timestepping parameters
   opt.Ti   = 50   * year; % duration of injection phase
   opt.dTi  = 2    * year; % timestep size during injection
   opt.Tm   = 2000 * year; % duration of migration phase
   opt.dTm  = 20   * year; % timestep size during migration

   % fluid model parameters
   opt.surf_temp = 12;                        % in Celsius
   opt.temp_grad = 30 / (kilo*meter);         % thermal gradient
   opt.p_range   = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   opt.t_range   = [  4, 250] + 274;           % CO2 default temperature range
   opt.res_vals  = [.11, .21];                % residual saturation values (if enabled)
   opt.cw        = 4.3e-5 / barsa;            % linear water compressibility

   opt.xres      = 1000;
   
   opt = merge_options(opt, varargin{:});

   simulation_count = 1; % global count of simulation runs
   total_count = numel(opt.A) * numel(opt.residual) * numel(opt.subscale_types) ...
       * numel(opt.dis_types) * numel(opt.fluid_types);
   outcomes = cell(total_count, 1);

   % Loop over grid types
   for A = opt.A

      % Loop over whether or not to use residual saturation
      for residual = opt.residual

         % Loop over fluid type
         for fluid_type = opt.fluid_types

            % Loop over topsurface type
            for subscale_type = opt.subscale_types

               % Loop over dissolution types
               for dis_type = opt.dis_types

                  % Constructing aquifer
                  if ~strcmpi(subscale_type{:}, 'smooth')
                     % We model upscaled caprock undulations implicitly, so the
                     % geometrical caprock model itself will be smooth
                     aquifer = makeAquiferModel('A', 0, 'D', opt.depth, ...
                                                'nx', opt.xres);
                  else
                     aquifer = makeAquiferModel('A', A, 'D', opt.depth, ...
                                                'nx', opt.xres);
                  end

                  % Make fluid model
                  [fluid, fluid_params] = ...
                      setup_fluid_model(opt, aquifer, residual, ...
                                             fluid_type{:}, ...
                                             toptrap(subscale_type{:}, A), ...
                                             subscale_type{:}, dis_type{:});

                  % Defining injection schedule
                  [schedule, Winj] = setup_schedule(opt, fluid, aquifer.W, aquifer.Gt);

                  % Defining initial state
                  initState = setup_init_state(fluid, aquifer.G, Winj, dis_type);

                  % Set up and run complete model
                  model = CO2VEBlackOilTypeModel(aquifer.Gt, aquifer.rock2D, fluid);
                  [wellSols, states] = ...
                      simulateScheduleAD(initState, model, schedule, ...
                                         'NonLinearSolver', ...
                                         NonLinearSolver('useRelaxation', true));

                  % Storing outcome
                  outcomes{simulation_count} = struct('states', {states}, ...
                                                      'wellSols', {wellSols}, ...
                                                      'Gt', {aquifer.Gt}, ...
                                                      'fluid_params', {fluid_params});
                  simulation_count = simulation_count + 1;
               end
            end
         end
      end
   end
end

% ----------------------------------------------------------------------------

function hts = toptrap(subscale_type, A)

   if strcmpi(subscale_type, 'smooth')
      hts = []; % no implicit subscale undulations - all geometry is
                % described explicitly by the model grid
   else
      % Create a model of the subscale undulations
      aquifer = makeAquiferModel('A', A);
      z = aquifer.Gt.cells.z;
      zt = max(z) * ones(size(z));
      for i = 2:numel(zt) - 1
         zt(i) = max(z(i:end));
      end
      zt(end) = max(zt(end -1), z(end));
      ht = zt - z;
      ff = exp(-linspace( -25, 25, 501).^2); ff = ff' / sum(ff);
      hts = filter2(ff, ht);
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
   last_dt = opt.Ti - sum(istep);
   if last_dt > 0
      istep = [istep; last_dt];
   end

   % Specify duration of individual migration timesteps
   mstep = linspace(0.5 * year, opt.dTm, 5)';
   mstep = [mstep; ones(floor((opt.Tm - sum(mstep)) / opt.dTm), 1) * opt.dTm];
   last_dt = opt.Tm - sum(mstep);
   if last_dt > 0
      mstep = [mstep; last_dt];
   end

   % Put everything together in a schedule structure
   schedule = struct('control', [struct('W', Winj), struct('W', Wmig)]    , ...
                     'step'   , struct('control', [ones(size(istep));       ...
                                                   2 * ones(size(mstep))] , ...
                                       'val'    , [istep; mstep]));
end

% ----------------------------------------------------------------------------

function initState = setup_init_state(fluid, G, Winj, dis_type)

   W_val   = Winj(2).val;
   W_depth = G.cells.centroids(Winj(2).cells, 3);

   nc = G.cells.num;
   pfun = @(z) W_val + (z - W_depth) * norm(gravity) * fluid.rhoWS;

   initState = struct('pressure', pfun(G.cells.centroids(:, 3)), ...
                      's'       , [ones(nc, 1), zeros(nc, 1)], ...
                      'sGmax'   , zeros(nc, 1));

   if ~strcmpi(dis_type, 'none')
      % We will model dissolution too
      initState.rs = zeros(nc, 1);
   end
end

% ----------------------------------------------------------------------------

function [fluid, params] = setup_fluid_model(opt, aquifer, residual, fluid_type, top_trap, ...
                                        subscale_type, dis_type)

   T = aquifer.Gt.cells.z * opt.temp_grad + (274 + opt.surf_temp); % add 274 to get Kelvin
   res_vals = opt.res_vals * residual; % becomes zero if 'residual' is false

   params.args = {aquifer.Gt, aquifer.rock2D, fluid_type, ...
                  'top_trap'    , top_trap, ...
                  'surf_topo'   , subscale_type, ...
                  'fixedT'      , T                          , ...
                  'co2_rho_pvt' , [opt.p_range, opt.t_range] , ...
                  'wat_rho_pvt' , [opt.cw, 100 * barsa]      , ...
                  'dissolution' , ~strcmpi(dis_type, 'none') , ...
                  'residual'    , res_vals};
   params.is_instant = strcmpi(dis_type, 'instant');
   
   fluid = makeVEFluid(params.args{:});
   
   if params.is_instant;
      fluid.dis_rate = 0; % value of zero indicates instant dissolution
   end
      
end

% ----------------------------------------------------------------------------
