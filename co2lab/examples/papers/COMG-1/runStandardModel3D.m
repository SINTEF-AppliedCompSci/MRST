function runStandardModel3D(save_filename, plot_routine, varargin)
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

   moduleCheck('co2lab', 'ad-fi', 'ad-core', 'ad-props');
   gravity reset on;

   %% Check for presence of already-computed result; compute or re-use
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

   %% Plot result to produce the figures
   plot_routine(sim_outcome);

end

% ----------------------------------------------------------------------------

function outcomes = run_standard_simulation(varargin)

   % Loop parameters
   opt.zres           = 20;                  % vertical resolution
   opt.A              = 0;                   % magnitudes of subscale undulations
   opt.depth          = 2300;                % depth of aquifer
   opt.residual       = false;               % whether to enable residual saturation
   opt.dis_types      = {'none'};            % dissol. types
                                             % ('none'/'rate'/'instant')
   
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
   opt.cap_press = 5 * kilo * Pascal;         % linear capillary pressure
   opt.n         = 1.5;                       % Corey exponent for rel. perms
   opt = merge_options(opt, varargin{:});

   simulation_count = 1; % global count of simulation runs
   total_count = numel(opt.A) * numel(opt.residual) * numel(opt.dis_types);
   outcomes = cell(total_count, 1);

   %% Loop over grid types
   for A = opt.A

      %% Loop over whether or not to use residual saturation
      for residual = opt.residual
   
         %% Loop over dissolution types
         for dis_type = opt.dis_types

            aquifer = makeAquiferModel('A', A, 'D', opt.depth, 'nz', opt.zres);

            % Make fluid model
            [fluid, fluid_params] = ...
                setup_fluid_model(opt, aquifer, residual, ...
                                       'sharp interface', ...
                                       [], ...
                                       'smooth', dis_type{:});
            
            % Defining injection schedule
            [schedule, Winj] = setup_schedule(opt, fluid, aquifer.W3D, aquifer.G);
            
            % Defining initial state
            initState = setup_init_state(fluid, aquifer.G, Winj, dis_type);
            
            % Set up and run complete model
            model = twoPhaseGasWaterModel(aquifer.G, aquifer.rock, fluid, ...
                                          opt.surf_temp + 274, ...
                                          opt.temp_grad * 1000);
            
            [wellSols, states] = simulateScheduleAD(initState, model, schedule);
            
            % Storing outcome
            outcomes{simulation_count} = struct('states', {states}, ...
                                                'wellSols', {wellSols}, ...
                                                'G', {aquifer.G}, ...
                                                'fluid_params', {fluid_params});
            simulation_count = simulation_count + 1;
         end
      end
   end
end
   
% ----------------------------------------------------------------------------

function [schedule, Winj, Wmig] = setup_schedule(opt, fluid, W, G)

   % Specify well structures
   Winj = W;
   Winj(2).val = fluid.rhoWS * max(G.cells.centroids((Winj(2).cells),3)) * norm(gravity);
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

function initState = setup_init_state(fluid, G, Winj, dis_type)

   W_val   = Winj(2).val;
   W_depth = max(G.cells.centroids(Winj(2).cells, 3));

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

   T = aquifer.G.cells.centroids(:,3) * opt.temp_grad + (274 + opt.surf_temp); % add 274 to get Kelvin
   res_vals = opt.res_vals * residual; % becomes zero if 'residual' is false

   params.args = {aquifer.Gt, aquifer.rock2D, fluid_type, ...
                  'top_trap'    , top_trap, ...
                  'surf_topo'   , subscale_type, ...
                  'fixedT'      , T                          , ...
                  'co2_rho_pvt' , [opt.p_range, opt.t_range] , ...
                  'wat_rho_pvt' , [opt.cw, 100 * barsa]      , ...
                  'dissolution' , ~strcmpi(dis_type, 'none') , ...
                  'residual'    , res_vals};
                  % 'co2_mu_ref'  , 1e-5 * Pascal * second, ...%6e-5 * Pascal * second, ...
                  % 'wat_mu_ref'  , 1e-5 * Pascal * second, ... %8e-4 * Pascal * second, ...
   params.is_instant = strcmpi(dis_type, 'instant');
   
   fluid = makeVEFluid(params.args{:});
   
   %% Converting to 3D fluid by changing relperm and cap. press. functions
   
   fluid = rmfield(fluid, 'krG');
   fluid = rmfield(fluid, 'krW');
   fluid = rmfield(fluid, 'pcWG');
   
   %fluid.pcGW = @(sg, p, varagin) 0;
   %fluid.pcGW = @(sg, p, varagin) 10 * kilo * Pascal * sg;
   fluid.pcGW = @(sg, p, varagin) opt.cap_press * sg;
   
   krW = coreyPhaseRelpermAD(opt.n, res_vals(2));
   krG = coreyPhaseRelpermAD(opt.n, res_vals(1));
   
   fluid.relPerm = @(sg) deal(krW(1-sg), krG(sg));
   
   if params.is_instant;
      fluid.dis_rate = 0; % value of zero indicates instant dissolution
   end
      
end

% ----------------------------------------------------------------------------
