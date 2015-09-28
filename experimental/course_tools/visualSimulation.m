function visualSimulation(initState, model, schedule, varargin)

   gravity on;
   moduleCheck('ad-core', 'ad-props');

   opt.h_min_threshold = 1e-2; % 1 cm
   opt.s_min_threshold = 1e-4;
   opt.rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density
   opt.window_size = [700 900];
   opt.trapstruct = []; % trapping structure as computed by trapAnalysis
   opt.dh = []; % subscale trapping structure, if available
   opt.savefile = []; % if nonempty string, save result to the specified filename
   
   opt = merge_options(opt, varargin{:});
   
   if isempty(opt.trapstruct)
      opt.trapstruct = trapAnalysis(model.G, false);
   end
      
   %% Setting up interactive interface
   iface.h = figure();
   set_size(iface.h, opt.window_size(1), opt.window_size(2));

   % Graphical window
   iface.ax = axes('parent', iface.h, 'position', [.05 .15 .90 .78]);
   
   % Slider
   step = 1/(numel(schedule.step.val)-1);
   iface.slider = uicontrol('parent', iface.h, ...
                            'style', 'slider', ...
                            'units','normalized', ...
                            'position', [.05 .1 .80 .03], ...
                            'enable', 'off', ...
                            'value', 1, ...
                            'min', 1, ...
                            'max', numel(schedule.step.val), ...
                            'sliderstep', [step 5*step], ...
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
   update_slider_display(); % call directly here to set display value
   
   % Variable picker
   iface.picker = uicontrol('parent', iface.h, ...
                            'style', 'popup', ...
                            'enable', 'inactive', ...
                            'units', 'normalized', ...
                            'position', [.05 .05 .4 .03], ...
                            'string', variable_string(initState), ...
                            'value', 2, ... % 2 typically means saturation
                            'callback', @(varargin) picker_change_callback());

   drawnow;
   %% Run simulation

   [wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule, ...
                                                     'afterStepFn', @after_step_callback);
   
   % check that context still exists
   if ~ishandle(iface.h)
      disp('Window closed.  Aborting simulation');
      return;
   end

   % save outcome if requested
   if ~isempty(opt.savefile)
      save(opt.savefile, 'model', 'schedule', 'wellSols', 'states', 'sim_report');
   end
   
   % Plot the inventory in separate window
   h2 = figure; plot(1); ax = get(h2, 'currentaxes');

   reports = makeReports(model.G, {initState, states{:}}, model.rock, model.fluid, schedule, ...
                         [model.fluid.res_water, model.fluid.res_gas], ...
                         opt.trapstruct, opt.dh);
   directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
   fsize = 24;
   set(get(gca, 'xlabel'), 'fontsize', fsize)
   set(get(gca, 'ylabel'), 'fontsize', fsize)
   set(gca,'fontsize', fsize);
   set(gcf, 'position', [1, 1, 850, 850]);
   
   %% Post-simulation viewing opportunity
   rotate3d on;
   set(iface.h, 'toolbar','figure');
   set(iface.slider, 'enable', 'on');
   set(iface.picker, 'enable', 'on');
   slider_change_callback(); % call directly here to set display value
   
   
   %% ================================= Callbacks =================================
   
   function [model, states, reports, solver, ok] = after_step_callback(model, states, reports, solver, schedule, simtime)
      
      % If htis fails, the window was likely closed and we should abort the simulation
      try 
         
         % Identify last computed state
         state_ix = find(cellfun(@(x) ~isempty(x), states), 1, 'last');
         
         % Determine what to display
         field = fieldnames(initState);
         field{end+1} = 'h'; % adding height
         field{end+1} = 'plume_s_avg';
         field = field(get(iface.picker, 'value'));
         %data = states{state_ix}.(field{:});
         switch lower(field{:})
           case 's'
             % keep only co2 saturation
             data = states{state_ix}.s(:,2);
             data(data < opt.s_min_threshold) = NaN;
             % data = data(:,2); 
             % data = data .* model.G.cells.H;
             % data(data < opt.h_min_threshold) = NaN;
           case 'pressure'
             data = states{state_ix}.pressure - initState.pressure; % visualize overpressure
           case {'h', 'plume_s_avg'}
             p = states{state_ix}.pressure;
             sG = states{state_ix}.s(:,2);
             sG_max = states{state_ix}.sGmax;
             drho = (model.fluid.rhoWS .* model.fluid.bW(p) - model.fluid.rhoGS * model.fluid.bG(p));
             pc = model.fluid.pcWG(sG, p, 'sGmax', sG_max);
             h = pc ./ (drho*norm(gravity));
             h(h<opt.h_min_threshold) = NaN;
             if strcmpi(field{:}, 'h')
                data = h;
             else % free plume s-average
                sfree = free_sg(sG, sG_max, struct('res_water', model.fluid.res_water, ...
                                                   'res_gas', model.fluid.res_gas));
                data = sfree.*model.G.cells.H./h./(1-model.fluid.res_water);
             end
           otherwise
             data = states{state_ix}.(field{:});
         end
         
         axes(iface.ax); axis auto; cla;
         % plot wells
         hold on;
         W = schedule.control(schedule.step.control(state_ix)).W;
         for i = 1:numel(W)
            wellcell = W(i).cells;
            pt = model.G.cells.centroids(wellcell,:);
            z  = model.G.cells.z(wellcell);
            plot3(pt(1), pt(2), z *0.98, 'ro', 'markersize', 8, ...
                  'MarkerFaceColor',[0 0 0]);
         end
         
         plotCellData(model.G, data); colorbar;
         %rotate3d on;
            
         % Update slider
         set(iface.slider, 'value', state_ix);
         update_slider_display();
         ok = true;
      catch
         ok=false;
      end
      
   end
   
   % ----------------------------------------------------------------------------
   
   function update_slider_display()
      tvals = cumsum(schedule.step.val);
      t = tvals(floor(get(iface.slider, 'value'))) /year;
      set(iface.slider_display, 'string', sprintf('%3.1f year', t));
                                                        
   end
   

   % ----------------------------------------------------------------------------
   
   function slider_change_callback()
      update_slider_display();
      ix = floor(get(iface.slider, 'value'));
      
      after_step_callback(model, states(1:ix), [], [], schedule, []);

   end
   
   % ----------------------------------------------------------------------------
   function picker_change_callback()
      
      after_step_callback(model, states(1:floor(get(iface.slider, 'value'))), [], [], schedule, []);

   end
   
   
end
% ============================================================================

function str = variable_string(state)

   % Basic state fields
   str = cellfun(@(x) [x,'|'], fieldnames(state),'uniformoutput', false);
   str = [str{:}];
   str = str(1:end-1);
   
   % Additional fields
   str = [str, '|h'];
   str = [str, '|plume_s_avg'];
   
end

% ----------------------------------------------------------------------------

function h = set_size(h, res_x, res_y)
% Utility function to resize a graphical window
   
   pos = get(h, 'position');
   set(h, 'position', [pos(1:2), res_x, res_y]);
   
end
