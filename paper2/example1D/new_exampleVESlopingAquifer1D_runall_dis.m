function new_exampleVESlopingAquifer1D_runall_dis(varargin)

   opt = set_options(varargin{:});
   gravity on;
   mrstVerbose opt.verbose;
   
   % ensuring directory for saving data exists
   mkdir(opt.savedir_name);
   
   % Different capillarity models for use in simulations
   cap_type = {'sharp interface'  , ...
               'linear cap.'      , ...
               'P-scaled table'   , ...
               'P-K-scaled table' , ...
               'S table'};

   %% Generating helper structures common to all simulation runs
   rock      = setup_rock     (opt);
   schedule  = setup_schedule (opt);
   
   %% Running a series of simulations, with different parameters
   t1 = tic;
   for use_dis = [true, false]              % whether to include CO2 dissolution 
      for depth = [2300, 1300]              % depth of reservoir 
         for smooth = [false, true]         % use rough or smooth grid
            for res_sat = [true, false]     % use residual saturation
               for c_ix = 1:numel(cap_type) % capillarity model to use 
                  
                  Gt        = setup_grid  (opt, smooth, depth);
                  fluid     = makeADIFluid(cap_type{c_ix}, res_sat);
                  rock2D    = averageRock (rock, Gt);
                  initState = setup_state (opt);
                  model     = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, ...
                                                     'disgas', use_dis, ...
                                                     'disrate', opt.disrate);
                  [wellSols, states] = ...
                      simulateScheduleAD(initState, model, schedule);
                  
                  filename = make_filename(use_dis, depth, smooth, res_sat, ...
                                           cap_type(c_ix)); 
                  
                  save_result(opt.savedir_name, filename, Gt, result, @@)
               end
            end
         end
      end
   end

   toc(t1);
   
end

% ============================================================================

function initState = setup_state(opt, Gt)
   
end

% ----------------------------------------------------------------------------

function save_result(dirname, filename, Gt, @@)
   @@
end

% ----------------------------------------------------------------------------

function filename = make_filename(dis, depth, smooth, res_sat, ctype)
   
end

% ----------------------------------------------------------------------------

function schedule = setup_schedule(opt)
@@   
end

% ----------------------------------------------------------------------------

function rock = setup_rock(opt)
   
   rock = struct('perm', opt.perm * ones(opt.num_cells, 1), ...
                 'poro', opt.poro * ones(opt.num_cells, 1));
end

% ----------------------------------------------------------------------------
function Gt = setup_grid(opt, smooth, depth)

   % Constructing basic grid
   G = cartGrid([opt.num_cells, 1, 1], opt.physdim);
   
   % Setting grid to correct depth
   LL = opt.physdim(1) * 2 / 3;
   G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + ...
                          depth - LL * sin(x / LL) * tan(opt.slope);
   
   % Unless grid is supposed to be smooth, add fine-scale features
   if ~smooth
      G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 2 * sin(2 * pi * x / 0.3e3); 
   end
   
   Gt = topSurfaceGrid(computeGeometry(G));
end

% ----------------------------------------------------------------------------

function opt = set_options(varargin)

% Parameters for the simulation
   opt.verbose = true;
   
   % Grid-related parameters
   opt.num_cells = 1000;
   opt.physdim   = [30e3 10e3 50]; % physical dimension of reservoir
   opt.slope     = 0.03;

   % rock-related parameters
   opt.perm = 1000 * milli * darcy; % (isotropic) permeability
   opt.poro = 0.02;                 % porosity

   % fluid-related parameters
   opt.mu =  [6e-5 8e-4] * Pascal * second;   % viscosity [CO2, brine]
   opt.rho = [760 1100] * kilogram / meter^3; % density   [CO2, brine]
   opt.sr = 0.21; % residual saturation, CO2   (if option is active)
   opt.sw = 0.11; % residual saturation, water (if option is active)
   opt.disrate = 5e-11;
   
   
   % Time-related parameters
   opt.T_inj = 50 * year; 
   opt.dt_inj = 2 * year;
   opt.T_mig = 2000 * year;
   opt.dt_mig = 20 * year;
   
   % directory for saving results
   opt.savedir_name = 'data_all_results_after';
   
   opt = merge_options(opt, varargin{:});
   
end
