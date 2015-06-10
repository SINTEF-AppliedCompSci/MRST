function visualSimulation(initState, model, schedule, varargin)

   gravity on;
   moduleCheck('ad-core', 'ad-props');

   opt.h_min_threshold = 1e-2; % 1 cm
   opt.rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density
   opt.window_size = [700 900];
   
   opt = merge_options(opt, varargin{:});
   
      
   %% Setting up interactive interface
   iface.h = figure();
   set_size(iface.h, opt.window_size(1), opt.window_size(2));

   % Graphical window
   iface.ax = axes('parent', iface.h, 'position', [.05 .2 .90 .78]);
   
   % Slider
   iface.slider = uicontrol('parent', iface.h, ...
                            'style', 'slider', ...
                            'units','normalized', ...
                            'position', [.05 .1 .80 .03], ...
                            'value', 0, ...
                            'min', 0, ...
                            'max', sum(schedule.step.val)/year, ...
                            'callback', @(varargin) ...
                            slider_change_callback());
   iface.slider_display = uicontrol('parent', iface.h, ...
                                    'style', 'edit', ...
                                    'enable', 'inactive', ...
                                    'horizontalalignment', 'left', ...
                                    'units', 'normalized', ...
                                    'fontsize', 8, ...
                                    'handlevisibility', 'off', ...
                                    'position', [.86 .1 .1 .03]);
   slider_change_callback(); % call directly here to set display value
   
   % Variable picker
   iface.picker = uicontrol('parent', iface.h, ...
                            'style', 'popup', ...
                            'units', 'normalized', ...
                            'position', [.05 .05 .4 .03], ...
                            'string', variable_string(initState), ...
                            'value', 2, ... % 2 typically means saturation
                            'callback', @(varargin) picker_change_callback());
      
   %% Run simulation

   [wellSols, states] = simulateScheduleAD(initState, model, schedule, ...
                                           'afterStepFn', @after_step_callback);
   
   
   %% Post-simulation viewing opportunity
   
   
   % for i = 1:20
   %    if ~ishandle(iface.h)
   %       disp('Window closed.  Aborting simulation');
   %       return;
   %    end
   %    pause(1)
   % end
   % keyboard;
   
   %% ================================= Callbacks =================================
   
   function after_step_callback(model, states, reports, solver, schedule, simtime)
      
      % Identify last computed state
      state_ix = find(cellfun(@(x) ~isempty(x), states), 1, 'last');
         
      % Determine what to display
      field = fieldnames(initState);
      field = field(get(iface.picker, 'value'));
      data = states{state_ix}.(field);
      if strcmpi(field, 's')
         % keep only co2 saturation, and convert to height
         data = data(:,2); 
         data = data .* model.Gt.cells.H;
         data(h < opt.h_min_threshold) = NaN;
      end
      
      axes(var.ax); axis auto; cla;

      plotCellData(model.Gt, data);
      
      % Update slider
      set(iface.slider, 'value', simtime/year);
      
   end
   
   % ----------------------------------------------------------------------------
   
   function slider_change_callback()
      set(iface.slider_display, 'string', sprintf('%3.1f year', get(iface.slider, 'value')));
   end
   
   % ----------------------------------------------------------------------------
   function picker_change_callback()
      % str = fieldnames(initState);
      % str = str{get(iface.picker, 'value')};
      % disp(str);
   end
   
   
end
% ============================================================================

function str = variable_string(state)

   str = cellfun(@(x) [x,'|'], fieldnames(state),'uniformoutput', false);
   str = [str{:}];
   str = str(1:end-1);
   
end

% ----------------------------------------------------------------------------

function h = set_size(h, res_x, res_y)
% Utility function to resize a graphical window
   
   pos = get(h, 'position');
   set(h, 'position', [pos(1:2), res_x, res_y]);
   
end
