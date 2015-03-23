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
                  
                  % Setting up model
                  aquifer   = makeAquiferModel('D', depth, 'A', 2*(~smooth));
                  fluid     = makeFluidModel(struct('Gt'     , Gt           , ...
                                                    'rock'   , rock         , ...
                                                    'rock2D' , rock2D)      , ...
                                             'fluidType'   , cap_type{c_ix} , ...
                                             'residual'    , res_sat        , ...
                                             'dissolution' , use_dis);
                  initState = setup_state (opt);
                  model     = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, ...
                                                     'disgas', use_dis, ...
                                                     'disrate', opt.disrate);
                  
                  % Running simulation
                  [wellSols, states] = ...
                      simulateScheduleAD(initState, model, schedule);
                  
                  % Saving result
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

function opt = set_options(varargin)

   % Rate of dissolution
   opt.disrate = 5e-11;
   
   % Parameters for the simulation
   opt.verbose = true;
   
   % Time-related parameters
   opt.T_inj  = 50 * year; 
   opt.dt_inj = 2 * year;
   opt.T_mig  = 2000 * year;
   opt.dt_mig = 20 * year;
   
   % directory for saving results
   opt.savedir_name = 'data_all_results_after';
   
   opt = merge_options(opt, varargin{:});
   
ppend
