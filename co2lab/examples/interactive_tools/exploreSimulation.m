function exploreSimulation(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
   
   gravity on;
   mrstModule add ad-core;

   rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density
   
   opt.grid_coarsening   = 4;
   opt.default_formation = 'Utsirafm';
   opt.window_size       = [1024 768];
   opt.seafloor_depth    = 100 * meter;
   opt.seafloor_temp     =  7; % in Celsius
   opt.temp_gradient     = 35.6; % degrees per kilometer
   opt.water_density     = 1000; % kg per m3
   opt.dis_max           = (53 * kilogram / meter^3) / rhoCref; % value from CO2store
   opt.max_num_wells     = 10;
   opt.default_rate      = 1 * mega * 1e3 / year / rhoCref; % default injection rate
   opt.max_rate          = 10 * mega * 1e3 / year / rhoCref; % maximum allowed injection rate
   opt.seafloor_depth    = 100 * meter;
   opt.seafloor_temp     =  7; % in Celsius
   opt.temp_gradient     = 35.6; % degrees per kilometer
   opt.water_compr_val   = 4.3e-5/barsa;
   opt.pvMult            = 1e-5/barsa; % pore volume multiplier
   opt.water_residual    = 0.11;
   opt.co2_residual      = 0.21;
   opt.inj_time          = 50 * year;
   opt.inj_steps         = 10;
   opt.mig_time          = 3000 * year;
   opt.mig_steps         = 30;
   opt.well_radius       = 0.3;
   opt.C                 = 0.1;
   opt.subtrap_file      = 'utsira_subtrap_function_3.mat';
   opt.outside_distance  = 100 * kilo * meter; % used to adjust
                                               % transmissibilities of
                                               % semiopen faces
   opt.savefile = []; % default is that we do not save results
   
   opt = merge_options(opt, varargin{:});

   var.Gt                = []; % to be set by the call to 'set_formation'
   var.rock2D            = []; % rock properties of current formation
   var.ta                = []; % trapping analysis of current formation
   var.current_formation = ''; % name of formation currently loaded
   var.data              = []; % current face data to plot
   var.loops             = []; % one or more loops of boundary faces
   var.loops_bc          = []; % bc for all faces (0: closed, 1: semi-open, 2: open)
   var.co2               = CO2props('sharp_phase_boundary', false);% @@ line can probably be removed
   var.wells             = reset_wells(opt.max_num_wells);
   % Temporary variables used by different functions
   temps = [];
   
   
   
   set_formation(opt.default_formation, false);

   % Setting up interactive interface
   
   var.h = figure('Color',[.95 .95 .97]);%'KeyPressFcn', @(obj, e) parse_keypress(e.Key));
   
   % Graphical window
   var.ax = axes('parent', var.h, 'position', [0.05, 0.12, 0.5, 0.87]);

   % Formation selection
   names = formation_names();
   fsel = uicontrol('parent', var.h, ...
                    'style', 'popup', ...
                    'units', 'normalized', ...
                    'position', [0.05, 0.0 0.3 0.05], ...
                    'string', listify(names), ...
                    'value', formation_number(var.current_formation));
   set (fsel, 'Callback', @(es, ed) set_formation(names{get(es, 'Value')}, ...
                                                  true));
   
   % Group for setting boundaries
   bc_group = setup_bc_group([.56 .84 .19 .15]);

   % Group for setting interaction mode
   imode_group = setup_imode_group([.77 .84 .19 .15]);
   
   % Group for setting timesteps and durations
   time_entries = setup_time_group([.56, .72, .4, .1]);
   
   % Group for displaying and modifying wells
   [well_group, well_entries] = setup_well_group([.56 .20, .4, .50]);%#ok

   % group for toggling highlighting of traps
   trap_choice = setup_trap_group([.76, .01, .20, .05]);
   
   % launch button
   launch_button = uicontrol('parent', var.h, ...
                             'style', 'pushbutton', ...
                             'units', 'normalized', ...
                             'position', [.58 .08 .16 .09], ...
                             'string', 'Launch new simulation!', ...
                             'backgroundcolor', [0.9 0.95 1], ...
                             'callback', @(varargin) launch_simulation()); %#ok
                       
   
   % Other options
   [opt_group, opt_choices] = setup_opt_group([.76, .07, .20, .11]);%#ok
   
   % Launching by calling redraw function
   redraw();
   set_size(var.h, opt.window_size); % doing this twice       
   set_size(var.h, opt.window_size); % occasionally prevents  
                                     % some strange placements
   
   % ============================= LOCAL HELPER FUNCTIONS =============================

   function launch_simulation()

      if isempty(var.wells(1).pos)
         % no wells present
         msgbox('Add at least one well before running the simulation', ...
                      'No wells present', 'error', 'modal');
         return;
      end
      
      
      use_dissolution = logical(get(opt_choices.dissolution, 'value'));
      use_trapping    = logical(get(opt_choices.subscale,    'value'));
      use_cap_fringe  = logical(get(opt_choices.cap_fringe,  'value'));
      
      dh = [];
      topo = 'smooth';
      if use_trapping
         dh = computeToptraps(load(opt.subtrap_file), var.Gt, true);
         topo = 'inf_rough';
      end
      
      % Set up input parameters
      initState = setup_initstate();
      ref_p     = mean(initState.pressure); % use mean pressure as reference
                                            % pressure for linear compressibilities
      fluid     = makeVEFluid(var.Gt, var.rock2D, ...
                              if_else(use_cap_fringe, 'P-scaled table', 'sharp interface') , ...
                              'C'             , opt.C                                      , ...
                              'fixedT'        , caprock_temperature()                      , ...
                              'wat_rho_pvt'   , [opt.water_compr_val, ref_p]               , ...
                              'wat_rho_ref'   , opt.water_density                          , ...
                              'pvMult_p_ref'  , ref_p                                      , ...
                              'pvMult_fac'    , opt.pvMult                                 , ...
                              'residual'      , [opt.water_residual,  opt.co2_residual]    , ...
                              'dissolution'   , use_dissolution                            , ...
                              'dis_max'       , opt.dis_max                                , ...
                              'surf_topo'     , topo                                       , ...
                              'top_trap'      , dh);
                              
      model     = CO2VEBlackOilTypeModel(var.Gt, var.rock2D, fluid);
      schedule  = setup_schedule();
      if isempty(schedule)
         return; % something went wrong
      end

      semiopen_faces = get_bfaces_of_type(1);
      if ~isempty(semiopen_faces)
         % modifying transmissibilities for semi-open boundary faces

         semiopen_cells = sum(var.Gt.faces.neighbors(semiopen_faces,:), 2); 
         d = var.Gt.cells.centroids(semiopen_cells) - var.Gt.faces.centroids(semiopen_faces);
         d = sqrt(sum(d.^2, 2)); % norm of distance
         
         model.operators.T_all(semiopen_faces) = ...
             model.operators.T_all(semiopen_faces) .* d ./ (d + opt.outside_distance);
      end
      
      % spawn simulation window 
      
      visualSimulation(initState, model, schedule, 'rhoCref', rhoCref, ...
                       'trapstruct', var.ta, 'dh', dh, 'savefile', opt.savefile, ...
                       'trap_outlines', if_else(do_outline_traps(), trapfaces(), []));
   end

   % ----------------------------------------------------------------------------

   function res = get_bfaces_of_type(type)
      
      res = [];
      for i = 1:numel(var.loops)
         faces = var.loops{i};
         bcs   = var.loops_bc{i};
         res = [res; faces(bcs==type)];%#ok
      end
   end
   
   
   function schedule = setup_schedule()

      schedule = [];
      % Read specified time information
      inj_time = read_number_from_edit(time_entries.inj_time_edit, 'injection time');
      if isempty(inj_time) 
         return
      end
      mig_time = read_number_from_edit(time_entries.mig_time_edit, 'migration time');
      if isempty(mig_time) 
         return
      end
      inj_steps = read_number_from_edit(time_entries.inj_steps_edit, 'injection steps'); 
      if isempty(inj_steps)
         return
      else
         inj_steps = floor(inj_steps);
      end
      mig_steps = read_number_from_edit(time_entries.mig_steps_edit, 'migration steps'); 
      if isempty(mig_steps)
         return
      elseif mig_time == 0
         mig_steps = 0;
      else
         mig_steps = floor(mig_steps);
      end
   
      % Create wells 
      W = [];
      for i = 1:opt.max_num_wells

         if ~isempty(var.wells(i).pos)
            wcell_ix = closest_cell(var.Gt, [var.wells(i).pos,0], 1:var.Gt.cells.num);
            W = addWellVE(W, var.Gt, var.rock2D, wcell_ix , ...
                          'type'   , 'rate'               , ...
                          'val'    , var.wells(i).rate    , ...
                          'radius' , opt.well_radius      , ...
                          'comp_i' , [0 1]                , ...
                          'name'   , ['I', num2str(i)]);
         end
      end
      W_shut = W;
      for i = 1:numel(W_shut)
         W_shut(i).val = 0;
      end
      
      schedule.control(1).W = W;
      schedule.control(2).W = W_shut;
      
      % Define boundary conditions
      open_faces = [get_bfaces_of_type(1); get_bfaces_of_type(2)];

      schedule.control(1).bc = addBC([], open_faces, ...
                                     'pressure', ...
                                     var.Gt.faces.z(open_faces) * ...
                                     opt.water_density * norm(gravity), ...
                                     'sat', [1 0]);
      schedule.control(2).bc = schedule.control(1).bc;
      
      
      dTi = (inj_time * year) / inj_steps; %opt.inj_time / opt.inj_steps;
      dTm = (mig_time * year) / mig_steps; %opt.mig_time / opt.mig_steps;
      istepvec = ones(inj_steps, 1) * dTi;
      mstepvec = ones(mig_steps, 1) * dTm;
      
      schedule.step.val = [istepvec; mstepvec];
      schedule.step.control = [ones(inj_steps, 1); ones(mig_steps, 1) * 2];
      
   end

   % ----------------------------------------------------------------------------
   
   function num = read_number_from_edit(edit, name)
      num = str2num(get(edit, 'string'));%#ok
      if isempty(num) || ~isreal(num)
         msgbox(['Non-numeric value entered for ' name '. Please fix and try again'], ...
                'Could not parse value', 'error', 'modal');
         num = [];
      end
   end
   
   
   % ----------------------------------------------------------------------------
   
   function T = caprock_temperature()
      % Return temperature in Kelvin
      T = 273.15 + ...
          opt.seafloor_temp + ...
          (var.Gt.cells.z - opt.seafloor_depth) / 1e3 * opt.temp_gradient;
   end
   
   % ----------------------------------------------------------------------------
   
   function state = setup_initstate()
   
      state.pressure = var.Gt.cells.z * norm(gravity) * opt.water_density;
      state.s = repmat([1 0], var.Gt.cells.num, 1);
      state.sGmax = state.s(:,2);
      
      % If dissolution is activated, we need to add a field for that too
      if logical(get(opt_choices.dissolution, 'value'))
         state.rs = 0 * state.sGmax;
      end
      
   end
      
   % ----------------------------------------------------------------------------
   
   function wells = reset_wells(num)
      wells = repmat(struct('pos', [], 'rate', 0), num, 1);
   end
      
   % ----------------------------------------------------------------------------
   
   function res = get_interaction_type()
      res = get(get(imode_group, 'selectedobject'), 'string');
   end
   
   % ----------------------------------------------------------------------------
   
   function res = get_active_bc_type()
      switch(get(get(bc_group, 'selectedobject'), 'string'))
        case 'Closed'
          res = 0;
        case 'Semi-open'
          res = 1;
        case 'Open'
          res = 2;
        otherwise
          error('missing case');
      end
   end
   
   % ----------------------------------------------------------------------------
   
   function set_uniform_bc_callback(varargin)
      bc_type = get_active_bc_type();
      for l = 1:numel(var.loops_bc)
         var.loops_bc{l} = var.loops_bc{l} * 0 + bc_type;
      end
      redraw();
   end
   
   % ----------------------------------------------------------------------------
   
   function set_rotate_state_callback()
      sel = get(get(imode_group, 'selectedobject'), 'string');
      if strcmpi(sel, 'Rotate model')
         rotate3d(var.ax, 'on'); 
      else
         rotate3d(var.ax, 'off');
      end
   end
   
   % ----------------------------------------------------------------------------
   
   function [group, choices] =  setup_opt_group(pos)
      
      % Create group
      group = uipanel('Visible', 'off',...
                      'units', 'normalized', ...
                      'BackgroundColor', [1  .9  .9], ...
                      'position', pos);
      
      % create widgets (dissolution)
      choices.dissolution = uicontrol('parent', group, ...
                                      'style', 'checkbox', ...
                                      'units', 'normalized', ...
                                      'BackgroundColor', [1  .9  .9], ...
                                      'position' , [.1, .1, .2, .2]);

      label = uicontrol('parent', group, ...
                        'style', 'text', ...
                        'units', 'normalized', ...
                        'horizontalalignment', 'left', ...
                        'position', [.2, .08, .55, .2], ...
                        'BackgroundColor', [1  .9  .9], ...
                        'string', 'Use dissolution');%#ok

      % create widgets (subscale trapping)      
      choices.subscale = uicontrol('parent', group, ...
                                      'style', 'checkbox', ...
                                      'units', 'normalized', ...
                                      'BackgroundColor', [1  .9  .9], ...
                                      'position' , [.1, .4, .2, .2]);

      label2 = uicontrol('parent', group, ...
                         'style', 'text', ...
                         'units', 'normalized', ...
                         'horizontalalignment', 'left', ...
                         'position', [.2, .38, .80, .2], ...
                         'BackgroundColor', [1  .9  .9], ...
                         'string', 'Use subscale trapping');%#ok

      % create widgets (capillary fringe)      
      choices.cap_fringe = uicontrol('parent', group, ...
                                     'style', 'checkbox', ...
                                     'units', 'normalized', ...
                                     'BackgroundColor', [1  .9  .9], ...
                                     'position', [.1, .70, .2, .2]);

      label3 = uicontrol('parent', group, ...
                         'style', 'text', ...
                         'units', 'normalized', ...
                         'horizontalalignment', 'left', ...
                         'position', [.2, .68, .80, .2], ...
                         'BackgroundColor', [1  .9  .9], ...
                         'string', 'Use capillary fringe');%#ok
      
      
      set(group, 'visible', 'on');      
   end
   
   % ----------------------------------------------------------------------------
   
   function choice = setup_trap_group(pos)
      group = uipanel('visible', 'off', 'units', 'normalized', ...
                      'position',pos, 'BackgroundColor', [1 1 .9]);
      

      choice.display_traps = uicontrol('parent', group, ...
                                       'style', 'checkbox', ...
                                       'units', 'normalized', ...
                                       'backgroundcolor', [1 1 .9], ...
                                       'callback', @(varargin) redraw(), ...
                                       'position' , [.1, .05, .2, .9]);

      label = uicontrol('parent', group, ...
                        'style', 'text', ...
                        'units', 'normalized', ...
                        'horizontalalignment', 'left', ...
                        'backgroundcolor', [1 1 .9], ...
                        'position', [.2, .00, .55, .7], ...
                        'string', 'Outline traps');%#ok

      set(group, 'visible', 'on');
   end
   
   % ----------------------------------------------------------------------------
   
   function entries = setup_time_group(pos)
      group = uipanel('visible', 'off', 'units', 'normalized', ...
          'position', pos, 'BackgroundColor', [.9 .95 1]);
      set(group, 'visible', 'on');
      
      inj_time_label = uicontrol('parent', group, ...
                                 'style', 'text', ...
                                 'units', 'normalized', ...
                                 'horizontalalignment', 'left', ...
                                 'position', [.02 .45, .37, .4], ...
                                 'BackgroundColor', [.9 .95 1], ...
                                 'string', 'Injection time (years):');%#ok
      inj_step_label = uicontrol('parent', group, ...                
                                 'style', 'text', ...                
                                 'units', 'normalized', ...          
                                 'horizontalalignment', 'left', ...  
                                 'position', [.52 .45, .37, .4], ...  
                                 'BackgroundColor', [.9 .95 1], ...
                                 'string', '# injection timesteps:');%#ok
      mig_time_label = uicontrol('parent', group, ...
                                 'style', 'text', ...
                                 'units', 'normalized', ...
                                 'horizontalalignment', 'left', ...
                                 'position', [.02 .05, .37, .4], ...
                                 'BackgroundColor', [.9 .95 1], ...
                                 'string', 'Migration time (years):');%#ok
      mig_step_label = uicontrol('parent', group, ...                
                                 'style', 'text', ...                
                                 'units', 'normalized', ...          
                                 'horizontalalignment', 'left', ...  
                                 'position', [.52 .05, .37, .4], ...  
                                 'BackgroundColor', [.9 .95 1], ...
                                 'string', '# migration timesteps:');%#ok
      
      entries.inj_time_edit = uicontrol('parent', group, ...
                                'style', 'edit', ...
                                'units', 'normalized', ...
                                'horizontalalignment', 'left', ...
                                'position', [.4, .55, .09, .35], ...
                                'fontsize', 8, ...
                                'string', num2str(opt.inj_time/year));
      entries.inj_steps_edit = uicontrol('parent', group, ...
                                 'style', 'edit', ...
                                 'units', 'normalized', ...
                                 'horizontalalignment', 'left', ...
                                 'position', [.85, .55, .11, .35], ...
                                 'fontsize', 8, ...
                                 'string', num2str(opt.inj_steps));
      entries.mig_time_edit = uicontrol('parent', group, ...
                                'style', 'edit', ...
                                'units', 'normalized', ...
                                'horizontalalignment', 'left', ...
                                'position', [.4, .05, .09, .35], ...
                                'fontsize', 8, ...
                                'string', num2str(opt.mig_time/year));
      entries.mig_steps_edit = uicontrol('parent', group, ...
                                 'style', 'edit', ...
                                 'units', 'normalized', ...
                                 'horizontalalignment', 'left', ...
                                 'position', [.85, .05, .11, .35], ...
                                 'fontsize', 8, ...
                                 'string', num2str(opt.mig_steps));
      
   end
   % ----------------------------------------------------------------------------
   
   function group = setup_imode_group(pos)
      
      % Create group
      group = uibuttongroup('Visible', 'off',...
                            'units', 'normalized', ...
                            'position', pos, ...
                            'backgroundcolor', [1 0.95 0.9], ...
                            'selectionchangefcn', @(varargin) set_rotate_state_callback());
      % Create radiobuttons
      b1 = uicontrol(group, 'style', 'radiobutton', ...
                            'string', 'Edit boundaries', ...
                            'units', 'normalized', ...
                            'backgroundcolor', [1 0.95 0.9], ...
                            'position', [.1 .1 .9 .3]);%#ok
      b2 = uicontrol(group, 'style', 'radiobutton', ...
                            'string', 'Select wellsites', ...
                            'units', 'normalized', ...
                            'backgroundcolor', [1 0.95 0.9], ...
                            'position', [.1 .38 .9 .3]);%#ok
      b3 = uicontrol(group, 'style', 'radiobutton', ...
                            'string', 'Rotate model', ...
                            'units', 'normalized', ...
                            'backgroundcolor', [1 0.95 0.9], ...
                            'position', [.1 .66 .9 .3]);%#ok
                            
      set(group, 'visible', 'on');      
   end
   
   % ----------------------------------------------------------------------------

   function [group, entries] = setup_well_group(pos)
      % Create group
      group = uipanel('Visible', 'off',...
                      'units', 'normalized', ...
                      'BackgroundColor', [.925 .925 .925], ...
                      'position', pos);
      
      entries = [];
      
      yheight = (1 / opt.max_num_wells) * 0.95;
      for i = 1:opt.max_num_wells
         ypos    = ((opt.max_num_wells - (i)) / opt.max_num_wells) * 0.95 + .025;
         entries = [entries; add_well_entry(group, [0.05, ypos, 0.9, yheight], i)];%#ok
      end
      set(group, 'visible', 'on');
   end
   
   % ----------------------------------------------------------------------------
   
   function we = add_well_entry(group, pos, index)
      we.name = uicontrol('parent', group, ...
                          'style', 'edit',...
                          'string', sprintf('Well %i:', index), ...
                          'horizontalalignment', 'left', ...
                          'units', 'normalized', ...
                          'enable', 'inactive', ...
                          'fontsize', 8, ...
                          'handlevisibility', 'off', ...
                          'position', [pos(1), pos(2), pos(3)*0.1, pos(4)]);

      we.status = uicontrol('parent', group, ...
                            'style', 'edit',... % 'edit' rather than 'text' for vertical alignment
                            'string', '<none>', ...
                            'enable', 'inactive', ...
                            'horizontalalignment', 'center', ...
                            'units', 'normalized', ...
                            'fontsize', 8, ...
                            'handlevisibility', 'off', ...
                            'position', [pos(1) + pos(3)*0.1, pos(2), pos(3)*0.3, pos(4)]);
      
      we.delete = uicontrol('parent', group, ...
                            'style', 'pushbutton', ...
                            'string', 'X', ...
                            'units', 'normalized', ...
                            'position', [pos(1) + pos(3)*0.4, pos(2), pos(3)*0.1, pos(4)], ...
                            'handlevisibility', 'off', ...
                            'callback', @(varargin) clear_well_callback(index));
      we.rate = uicontrol('parent', group, ...
                          'style', 'slider', ...
                          'units','normalized', ...
                          'position', [pos(1) + pos(3)*0.51, pos(2), pos(3) * 0.35, pos(4)], ...
                          'value', opt.default_rate, ...
                          'min', 0, ...
                          'max', opt.max_rate, ...
                          'callback', @(varargin) set_new_rate_callback(index));
      we.rate_view = uicontrol('parent', group, ...
                            'style', 'edit',... % 'edit' rather than 'text' for vertical alignment
                            'string', sprintf('%3.1f Mt', opt.default_rate * year * rhoCref/1e9), ...
                            'enable', 'inactive', ...
                            'horizontalalignment', 'left', ...
                            'units', 'normalized', ...
                            'fontsize', 8, ...
                            'handlevisibility', 'off', ...
                            'position', [pos(1) + pos(3)*0.87, pos(2), pos(3)*0.11, pos(4)]);
   end

   % ----------------------------------------------------------------------------
   
   function set_new_rate_callback(ix)
      if isempty(var.wells(ix).pos)
         % well is inactive, keep rate to default value
         set(well_entries(ix).rate, 'value', opt.default_rate);
      else
         var.wells(ix).rate = get(well_entries(ix).rate, 'value');
      end      
      redraw();
   end
   
   % ----------------------------------------------------------------------------
   
   function clear_well_callback(ix)
      
      % shifting remaining wells up
      for i = ix:(numel(var.wells)-1)
         var.wells(i) = var.wells(i+1);
      end
      var.wells(end) = struct('pos', [], 'rate', 0);

      % redraw with new well information
      redraw();
   end
   
   % ----------------------------------------------------------------------------
   
   function group = setup_bc_group(pos)

      % create radiobutton group
      group = uibuttongroup('Visible', 'off',...
                            'units', 'normalized', ...
                            'backgroundcolor', [.95 .975 .95], ...
                            'position', pos);
      % create radiobuttons
      b1 = uicontrol(group, 'style', 'radiobutton', ...
                            'string', 'Closed', ...
                            'units', 'normalized', ...
                            'backgroundcolor', [1 .7 .7], ...
                            'position', [.1 .075 .4 .27], ...
                            'HandleVisibility', 'off');%#ok
      b2 = uicontrol(group, 'style', 'radiobutton', ... 
                            'string', 'Semi-open', ...
                            'units', 'normalized', ...
                            'position', [.1 .375 .4 .27], ...
                            'backgroundcolor', [1 1 .7], ...
                            'HandleVisibility', 'off');%#ok
      b3 = uicontrol(group, 'style', 'radiobutton', ...
                            'string', 'Open', ...
                            'units', 'normalized', ...
                            'position', [.1 .675 .4 .27], ...
                            'backgroundcolor', [.7 1 .7], ...
                            'HandleVisibility', 'off');%#ok
      pb1 = uicontrol(group, 'style', 'pushbutton', ...
                             'string', 'Set all', ...
                             'units', 'normalized', ...
                             'position', [.6 .7 .35 .2], ...
                             'handlevisibility', 'off', ...
                             'callback', @set_uniform_bc_callback);%#ok
      
      set(group, 'visible', 'on');
   end
   
   % ----------------------------------------------------------------------------
   
   function click_handler(varargin)

      pt = get(gca,'CurrentPoint'); pt = pt(end,:); 
      fn = [];
      
      switch get_interaction_type()
        case 'Edit boundaries'
          if ~isfield(temps, 'bc_segment_start')
             [temps.bc_segment_start, temps.bc_segment_loop_ix] = ...
                 closest_bface(var.Gt, pt, var.loops);
             fn = @() plot(pt(1), pt(2), '*r');
          else
             bc_segment_end = closest_bface(var.Gt, pt, var.loops, temps.bc_segment_loop_ix);
             
             cur_loop = var.loops{temps.bc_segment_loop_ix};
             ix1 = find(cur_loop == temps.bc_segment_start, 1);
             ix2 = find(cur_loop ==  bc_segment_end, 1);
             if ix1 > ix2 % ensure ix2 > ix1
                tmp = ix1;
                ix1 = ix2;
                ix2 = tmp;
             end
             if ix2-ix1 < numel(cur_loop)/2
                seq = ix1:ix2;
             else
                seq = [ix2:numel(cur_loop), 1:ix1];
             end

             var.loops_bc{temps.bc_segment_loop_ix}(seq) = get_active_bc_type();
             
             temps = rmfield(temps, 'bc_segment_start');
             temps = rmfield(temps, 'bc_segment_loop_ix');
          end
        case 'Select wellsites'

          % find first empty slot
          for i = 1:numel(var.wells)
             if isempty(var.wells(i).pos)
                break;
             end
          end
          full = (i==numel(var.wells) && ~isempty(var.wells(i).pos));
          if full
             % Discard oldest well, shift others downwards
             var.wells(1:end-1) = var.wells(2:end);
          end
          var.wells(i).pos = pt(1:2);
          var.wells(i).rate = opt.default_rate;
        case 'Rotate model'
          % do nothing - the radio button itself has toggled on rotate when activated
          % rotate3d(var.ax, 'on');
        otherwise
          disp('unimplemented');
          return;
      end
      redraw(fn);
      
      
      % %disp(pts);
      % ix = closest_cell(var.Gt, pts(end,:), 1:var.Gt.cells.num);
      % var.data(ix) = 1;
      % redraw();
   end

   % ----------------------------------------------------------------------------
   function res = do_outline_traps()

     res = logical(get(trap_choice.display_traps, 'value'));
      
   end
   
   % ----------------------------------------------------------------------------

   function res = trapfaces()
      trapcells = find(var.ta.traps);
      half_faces = var.Gt.cells.faces(mcolon(var.Gt.cells.facePos(trapcells), ...
                                             var.Gt.cells.facePos(trapcells+1)-1))';
      sorted_half_faces = sort(half_faces);
      repeats = diff(sorted_half_faces)==0;
      interior_faces = sorted_half_faces([repeats;false]);
      unique_faces = unique(half_faces);
      exterior_faces = setdiff(unique_faces, interior_faces);
      
      res = exterior_faces;
   end
   
   % ----------------------------------------------------------------------------
   
   function redraw(post_fnx)
      axes(var.ax); cla;
      axis auto;
      cla;
      
      % Draw current field (depth)
      plotCellData(var.Gt, var.data, 'EdgeAlpha',.1, 'buttondownfcn', @click_handler);
      
      if do_outline_traps()
         plotFaces(var.Gt, trapfaces(), 'edgecolor', 'y', 'linewidth', 2, 'edgealpha', 1);
      end
         
      % Draw boundary conditions
      for i = 1:numel(var.loops)
         loop = var.loops{i};
         loop_bc = var.loops_bc{i};
         plotFaces(var.Gt, loop(loop_bc==0), 'edgecolor', 'r', 'linewidth', 4); % closed
         plotFaces(var.Gt, loop(loop_bc==1), 'edgecolor', 'y', 'linewidth', 4); % closed
         plotFaces(var.Gt, loop(loop_bc==2), 'edgecolor', 'g', 'linewidth', 4); % closed
      end
      
      % Draw wells
      for i = 1:opt.max_num_wells
         w = var.wells(i);
         if ~isempty(w.pos)
            hold on;
            lon = w.pos(1);
            lat = w.pos(2);
            wellcell = closest_cell(var.Gt, [w.pos,0], 1:var.Gt.cells.num);
            plotWell(var.Gt.parent, ...
                     addWell([], var.Gt.parent, var.rock2D, wellcell, 'name', sprintf('W%i', i)), ...
                     'color', 'k', 'fontsize', 14);
            plot3(lon, lat, var.Gt.cells.z(wellcell)*0.98, 'ro', 'markersize', 8, ...
                  'MarkerFaceColor',[0 0 0]);
            set(well_entries(i).status, 'string', sprintf('(%4.2e, %4.2e)', lon, lat));
         else
            set(well_entries(i).status, 'string', '<none>');
         end
         annual_rate = var.wells(i).rate * year * rhoCref/1e9;
         set(well_entries(i).rate_view, 'string', sprintf('%3.1f Mt', annual_rate));
         set(well_entries(i).rate, 'value', var.wells(i).rate);
      end
      
      % Call optional function to complete redraw process
      if nargin>0 && ~isempty(post_fnx)
         hold on;
         post_fnx();
      end
      
      view(0, 90);
   end

   % ----------------------------------------------------------------------------
   
   function set_formation(name, do_redraw)
   
      % Default values, in case values are lacking in model file.
      default_perm = 200 * milli * darcy;
      default_poro = 0.2;
      
      var.current_formation = name;
      
      % Load grid and rock, and assure rock values are valid
      [var.Gt, var.rock2D] = getFormationTopGrid(name, opt.grid_coarsening);
      
      if any(isnan(var.rock2D.poro))
         warning('Replacing missing porosity value with default value.');
         var.rock2D.poro = default_poro * ones(size(var.rock2D.poro));
      end
      if any(isnan(var.rock2D.perm))
         warning('Replacing missing permeability value with default value.');
         var.rock2D.perm = default_perm * ones(size(var.rock2D.perm));
      end
      
      % Run trapping analysis (we need this information to compute
      % inventories)
      var.ta = trapAnalysis(var.Gt,false);
      
      var.data = var.Gt.cells.z; 
      var.loops = find_boundary_loops(var.Gt);
      % Setting all boundary conditions to open (2)
      var.loops_bc = cellfun(@(x) 0 * x + 2, var.loops, 'uniformoutput', false);

      var.wells = reset_wells(opt.max_num_wells);
      temps = []; % reset all temporary variables
      
      % Call 'redraw' if requested
      if do_redraw
         redraw();
      end
      
   end
   
end

% ======================= INDEPENDENT HELPER FUNCTIONS =======================

function num = formation_number(name)
   num = find(cellfun(@(x) strcmpi(x, name), formation_names()), 1);
end

% ----------------------------------------------------------------------------

function str = listify(names)

   str = cellfun(@(x) [x,'|'], names, 'uniformoutput', false);
   str = [str{:}];
   str = str(1:end-1);
end

% ----------------------------------------------------------------------------

function names = formation_names()

   names = {'Brentgrp', 'Brynefm', 'Fensfjordfm', 'Gassumfm', 'Huginfmeast', ...
            'Huginfmwest', 'Johansenfm', 'Krossfjordfm', 'Pliocenesand', ...
            'Sandnesfm', 'Skadefm', 'Sleipnerfm', 'Sognefjordfm', 'Statfjordfm', ...
            'Ulafm', 'Utsirafm', ...
            'Stofm', 'Nordmelafm', 'Tubaenfm', ...
            'Bjarmelandfm', ...
            'Arefm', 'Garnfm', 'Ilefm', 'Tiljefm'};
end

% ----------------------------------------------------------------------------

function h = set_size(h, res)
% Utility function to resize a graphical window
   
   pos  = get(h, 'position');
   screensize = get(0,'screensize');
   res = min([res; screensize(3:4)-[0 85]]);
   set(h, 'position', [min([pos(1:2); screensize(3:4)-res-[0 85]]) res]);
   
end

% ----------------------------------------------------------------------------

function ix = closest_cell(Gt, pt, candidates)
   
   d = bsxfun(@minus, [Gt.cells.centroids(candidates,:), Gt.cells.z(candidates)], pt);
   d = sum(d.^2, 2);
   [~, ix] = min(d);
   ix = candidates(ix);
end

% ----------------------------------------------------------------------------

function [cix, fix] = boundary_cells_and_faces(Gt)
   
   fix = find(prod(Gt.faces.neighbors, 2) ==0);
   cix = unique(sum(Gt.faces.neighbors(fix, :), 2));
end

% ----------------------------------------------------------------------------

function loops = find_boundary_loops(Gt)

   [~, fix] = boundary_cells_and_faces(Gt); % boundary face indices

   tmp = [fix, Gt.faces.nodes(Gt.faces.nodePos(fix));
          fix, Gt.faces.nodes(Gt.faces.nodePos(fix)+1)];
   tmp = sortrows(tmp, 2);
   tmp = reshape(tmp(:,1), 2, []); % columns now express face neighborships
   
   % defining connectivity matrix
   M = sparse(tmp(1,:), tmp(2,:), 1, Gt.faces.num, Gt.faces.num);
   M = spones(M+M');
   num_loops = 0;
   while nnz(M) > 0
      num_loops = num_loops + 1;
      [loop,~] = find(M, 1);
      next = find(M(loop, :), 1); % find one of the two neighbor faces
      while ~isempty(next)
         M(loop(end), next) = 0; %#ok
         M(next, loop(end)) = 0; %#ok
         loop = [loop; next]; %#ok
         next = find(M(next, :)); 
         assert(numel(next) <= 1);
      end
      assert(loop(1) == loop(end));
      loops{num_loops} = loop(1:end-1); %#ok
   end
end

% ----------------------------------------------------------------------------

function [face_ix, loop_ix] = closest_bface(Gt, pt, loops, imposed_loop)

   if (nargin > 3)
      loops_ixs = imposed_loop;
   else
      loops_ixs = 1:numel(loops);
   end
   
   closest_face = [];
   
   for l = loops_ixs
      % computing closest face for this loop
      loop = loops{l};
      loop_coords = [Gt.faces.centroids(loop,:), Gt.faces.z(loop,:)];
      dist = bsxfun(@minus, loop_coords, pt);
      dist = sum(dist.^2, 2);
      [dmin, num] = min(dist);
      closest_face = [closest_face; [dmin, num]]; %#ok
   end
   
   [~, i] = min(closest_face(:,1));
   loop_ix = loops_ixs(i);
   face_ix = loops{loop_ix}(closest_face(i, 2));
   
end
% ----------------------------------------------------------------------------

function val = if_else(cond, yes, no)
   
   if cond
      val = yes;
   else
      val = no;
   end
   
end
