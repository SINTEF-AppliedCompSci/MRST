function runStandardModel(save_filename, plot_routine, varargin)

   %% Check for presence of already-computed result; compute or re-use
   if (exist([save_filename, '.mat']) == 2 &&
      ask_user('Saved result found.  Re-use? [y/n]' '))
      fprintf('Re-using old result.\n');
      sim_outcome = load([savefilename, '.mat']);
   else
      fprintf('Recomputing result. \n');
      sim_outcome = run_standard_simulation(varargin{:});
      save(save_filename, sim_outcome);
   end
   
   %% Plot result to produce the figures
   plot_routine(results);
   
end

% ----------------------------------------------------------------------------

function outcomes = run_standard_simulation(varargin)

   opt.A = [0]; % magnitudes of subscale undulations
   opt.residual = [false]; % whether to consider residual saturation
   opt.relperm_types = {'sharp interface'}; % relperm model(s) to use
   opt.dis_types = {'none'}; % dissolution types ('none', 'rate' or 'instant')
   
   
   opt = merge_options(opt, varargin{:});
   
   simulation_count = 1; % global count of simulation runs
   
   %% Loop over grid types
   for A = opt.A

      % Constructing aquifer
      aquifer = makeAquiferModel_new('A', A);
      G = aquifer.G;
      Gt = aquifer.Gt;
      
      % Defining injection schedule
      schedule = ;
      
      % Defining initial state
      initState = ;
   
      %% Loop over whether or not to use residual saturation
      for residual = opt.residual

         %% Loop over relperm models
         for relperm_type = opt.relperm_types
            
            %% Loop over dissolution types
            for dis_type = opt.dis_types

               % Make fluid model
               fluid = ;
               
               % Set up and run complete model
               model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid);
               [wellSols, states] = simulateSchedule(state, model, schedule);
               
               % Storing outcome
               outcomes{simulation_count} = struct('states', {states});
               simulation_count = simulation_count + 1;
            end
         end
      end
   end
end
