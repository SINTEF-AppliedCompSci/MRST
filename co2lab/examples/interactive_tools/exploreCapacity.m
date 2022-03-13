function exploreCapacity(varargin)
   
% opt - structure containing variables that can be overridden by user at
%       command line
% var - structure containing variables that user cannot override at command
%       line 

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
   % Initializing 'opt'
   rhoCref = 760 * kilogram / meter ^3;
   
   opt.grid_coarsening = 2;
   opt.default_formation = 'Utsirafm';
   opt.window_size = [1024 768];
   opt.seafloor_depth = 100 * meter;
   opt.seafloor_temp  =  7; % in Celsius
   opt.temp_gradient  = 35.6; % degrees per kilometer
   opt.water_density  = 1000; % kg per m3
   opt.press_deviation = 0; % pressure devation from hydrostatic (percent)
   opt.res_sat_co2 = 0.21; 
   opt.res_sat_wat = 0.11;
   opt.dis_max = (53 * kilogram / meter^3) / rhoCref; % value from CO2store
   
   opt = merge_options(opt, varargin{:});
   
   % Initializing 'var'
   var.Gt_cached = struct();
   var.current_formation = '';
   var.co2 = CO2props('sharp_phase_boundary', true);
   set_formation(opt.default_formation, false);
   
   
   % Setting up interactive interface
   
   var.h = figure('Color',[.95 .95 .97],'KeyPressFcn', @(obj, e) parse_keypress(e.Key));
   set_size(var.h, opt.window_size);
   
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
   set (fsel, 'Callback', @(es, ed) set_formation(names{get(es, 'Value')}, true));
          
   % Radiobuttons for selection of data to map
   buttons = setup_button_group([.625, .60, .3, .38]);
   
   % Set sliders
   [s1, l1, v1] = add_slider('Temperature gradient (deg. C per km)', 'temp_gradient', 10, 50, [.59, .5, .35, .03]);%#ok
   [s2, l2, v2] = add_slider('Pressure deviation from hydrostatic (%)', 'press_deviation', -50, 50, [.59, .43, .35, .03]);%#ok
   [s3, l3, v3] = add_slider('Residual CO2 saturation', 'res_sat_co2', 0, 0.5, [.59, .36, .35, .03]);%#ok
   [s4, l4, v4] = add_slider('Residual brine saturation', 'res_sat_wat', 0, 0.5, [.59, .29, .35, .03]);%#ok
   [s5, l5, v5] = add_slider('Max dissolution', 'dis_max', 0, (100*kilogram/meter^3)/rhoCref, [.59, .22, .35, .03]);%#ok
   
   % Total figure reporting window
   rep_area = uicontrol('parent', var.h, ...
                        'style', 'text', ...
                        'units', 'normalized', ...
                        'position', [0.59, 0.02, 0.35, 0.16], ...
                        'horizontalalignment', 'left', ...
                        'fontsize', 12, ...
                        'backgroundcolor', [.95 .95 .97],...
                        'string', 'uninitialized');
   
   
   % Launching by calling redraw function
   redraw();
   
   % ============================= HELPER FUNCTIONS =============================

   function p = compute_pressure()
      Gt = var.Gt_cached.(var.current_formation).Gt;
      p = (Gt.cells.z * opt.water_density * norm(gravity)).* (1+opt.press_deviation/100);
   end
   % ----------------------------------------------------------------------------
   
   function t = compute_temperature()
      % NB: returns value in Celsius, not Kelvin
      Gt = var.Gt_cached.(var.current_formation).Gt;
      t = opt.seafloor_temp + (Gt.cells.z - opt.seafloor_depth) .* opt.temp_gradient ./ 1000;
   end

   % ----------------------------------------------------------------------------
   
   function cum_cap = compute_cumul_trapcap()
      
      [~, strap] = recompute_trapcap();
      
      Gt = var.Gt_cached.(var.current_formation).Gt;
      ta = var.Gt_cached.(var.current_formation).ta;
      trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, 'removeOverlap', false);
      tvols = [trees.value];%#ok
      int_tr = find(ta.trap_regions); %#ok ixs of cells spilling into interior trap
      [dummy, reindex] = sort([trees.root], 'ascend');%#ok

      cum_cap = zeros(Gt.cells.num, 1);   
      for i = 1:numel(ta.trap_z) % loop over each trap
       
         % ix of cells spilling directly into this trap
         cix = find(ta.trap_regions == i); 
         
         % cell indices of all cells of this trap, and its upstream traps
         aix = find(ismember(ta.traps, [trees(reindex(i)).traps]));
       
         % compute total structural trap capacity (in mass terms) of these
         % cells, and store result 
         cum_cap(cix) = sum(strap(aix));%#ok
      end
   end      
      
   % ----------------------------------------------------------------------------
   
   function [H1, strap, btrap_res, btrap_dis] = recompute_trapcap()
      p = compute_pressure();
      t = compute_temperature() + 273.15;
      
      Gt   = var.Gt_cached.(var.current_formation).Gt;
      ta   = var.Gt_cached.(var.current_formation).ta;
      poro = var.Gt_cached.(var.current_formation).poro;
      
      % computing structural trap heights (H1) for each cell
      trapcells     = find(ta.traps);
      H1            = zeros(Gt.cells.num, 1);
      if ~isempty(trapcells)
         H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
      end
      H1=min(H1,Gt.cells.H);
      H2 = Gt.cells.H - H1;
      assert(all(H1<=Gt.cells.H));
      
      % Computing total trapping volume in structural traps (dissolved and
      % structurally trapped
      strap_pvol_tot       = Gt.cells.volumes .* H1 .* poro;
      strap_pvol_co2_plume = strap_pvol_tot .* (1 - opt.res_sat_wat);
      strap_pvol_co2_diss  = strap_pvol_tot .* opt.res_sat_wat .* opt.dis_max;
      
      strap = strap_pvol_co2_plume .* var.co2.rho(p,t) + ...
              strap_pvol_co2_diss  .* rhoCref;

      % Computing total trapping volume below structural traps (dissolved and
      % residually trapped
      btrap_pvol_tot          = Gt.cells.volumes .* H2 .* poro;
      btrap_pvol_co2_residual = btrap_pvol_tot .* opt.res_sat_co2;
      btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-opt.res_sat_co2) .* opt.dis_max;
   
      btrap_res = btrap_pvol_co2_residual .* var.co2.rho(p,t);
      btrap_dis = btrap_pvol_co2_dissol   .* rhoCref;
   end
   
   % ----------------------------------------------------------------------------
   
   function [slider, label, vdisplay] = add_slider(label, field, minval, maxval, pos)
   
      slider = uicontrol('parent', var.h, ...
                         'style', 'slider', ...
                         'units', 'normalized', ...
                         'position', pos, ...
                         'value', opt.(field), ...
                         'backgroundcolor', [.97 .97 .95], ...
                         'min', minval, ...
                         'max', maxval);
      label = uicontrol('parent', var.h, ...
                        'style', 'text', ...
                        'units', 'normalized', ...
                        'position', [pos(1), pos(2) + pos(4), pos(3), pos(4)], ...
                        'horizontalalignment', 'left', ...
                        'string', label);
      set(label, 'BackgroundColor', get(gcf, 'color'));
      
      vdisplay = uicontrol('parent', var.h, ...
                           'style', 'text', ...
                           'units', 'normalized', ...
                           'position', [pos(1) + pos(3) + 0.01, pos(2) - 0.005, 0.05, pos(4)], ...
                           'horizontalalignment', 'left', ...
                           'string', num2str(opt.(field)));
      set(vdisplay, 'BackgroundColor', get(gcf, 'color'));
      
      set(slider, 'Callback', @(es, ed) update_after_slider_change(get(es, 'Value'), vdisplay, field));
      
   end
   
   % ----------------------------------------------------------------------------
   
   function update_after_slider_change(val, vdisp, field)
   
      opt.(field) = val;
      set(vdisp, 'string', sprintf('%6.2f',val));
      redraw();
   end
   
   % ----------------------------------------------------------------------------
   
   function bgroup = setup_button_group(pos)
      
      % Create radiobutton group 
      bgroup = uibuttongroup('Visible', 'off', ...
                             'units', 'normalized', ...
                             'position', pos, ...
                             'SelectionChangeFcn', @redraw);
      
      % Create radiobuttons
      r1 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'caprock depth (m)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.07, 0.8 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r2 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'formation thickness (m)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.16, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r3 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'caprock temperature (C)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.25, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r4 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'caprock pressure (MPa)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.34, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r5 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'caprock topography', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.43, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r6 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'caprock co2 density (kg/m3)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.52, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r7 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'total capacity (tonnes/m2)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.61, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r8 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'structural trap capacity (Mt)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.70, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
      r9 = uicontrol(bgroup, 'style',  'radiobutton', ...
                             'string', 'reachable structural capacity (Mt)', ...
                             'units', 'normalized', ...
                             'position', [0.1, 0.79, 0.8, 0.1], ...
                             'HandleVisibility', 'off'); %#ok
       
      set(bgroup, 'Visible', 'on');
   end

   
   
   % ----------------------------------------------------------------------------
   
   function parse_keypress(key)
   % Handling of keypresses should go in here.  Dummy for now.
      switch key
        case 'uparrow'
        case 'downarrow'
        case 'leftarrow'
        case 'rightarrow'
      end         
   end

   % ----------------------------------------------------------------------------

   function set_formation(name, do_redraw)
      
      default_poro = 0.2;
      var.current_formation = name;
      
      if ~isfield(var.Gt_cached, var.current_formation)
         wb = waitbar(0,strcat("Loading ",name));
         [Gt, rock] = getFormationTopGrid(var.current_formation, opt.grid_coarsening);
         waitbar(0.20,wb);
         var.Gt_cached.(var.current_formation).Gt = Gt;
         waitbar(0.40,wb);
         var.Gt_cached.(var.current_formation).poro = rock.poro;
         waitbar(0.60,wb);
         var.Gt_cached.(var.current_formation).ta = trapAnalysis(Gt, false);
         waitbar(0.80,wb);
      end
      
      if any(isnan(var.Gt_cached.(var.current_formation).poro))
         assert(all(isnan(var.Gt_cached.(var.current_formation).poro)));
         warning(['NaN porosity value encountered.  Replacing with default ' ...
                  'value']);
         var.Gt_cached.(var.current_formation).poro = default_poro;
      end
      
      if do_redraw
         redraw();
      end
      if exist('wb','var')
        waitbar(1,wb);
        close(wb);
      end
   end
   
   % ----------------------------------------------------------------------------
   
   function str = capacity_report()
      [~, strap, btrap_res, btrap_dis] = recompute_trapcap();
      
      tot_trap_capa = strap + btrap_res + btrap_dis;
      sum_tot = sum(tot_trap_capa);
      
      str = sprintf('Total trapping capacity: %.2f Gtons\n\n', sum_tot / giga / 1e3);
      str = [str, sprintf('Breakdown:\n')];
      str = [str, sprintf('Structural: %.2f Gtons (%5.2f%%)\n', sum(strap) / giga / 1e3, (sum(strap) / sum_tot) * 100)];
      str = [str, sprintf('Residual: %.2f Gtons (%5.2f%%)\n', sum(btrap_res)  / giga / 1e3, sum(btrap_res / sum_tot) * 100)];
      str = [str, sprintf('Dissolved: %.2f Gtons (%5.2f%%)\n', sum(btrap_dis) / giga / 1e3, sum(btrap_dis / sum_tot) * 100)];
      
   end
      
   % ----------------------------------------------------------------------------

   function redraw(varargin) 
      axes(var.ax); cla;
      axis auto;
      Gt = var.Gt_cached.(var.current_formation).Gt;
      ta = var.Gt_cached.(var.current_formation).ta;
      
      % Report total trapping figures
      set(rep_area, 'string', capacity_report());
      
      ealpha = 0.05; % edge alpha value
      
      % Determine what to plot
      switch (get(get(buttons, 'selectedobject'), 'String'))
      
        case 'caprock depth (m)'
          plotCellData(Gt.parent, Gt.cells.z, 'edgealpha', ealpha);
          colorbar; rotate3d on;
          
        case 'formation thickness (m)'
          plotCellData(Gt.parent, Gt.cells.H, 'edgealpha', ealpha);
          colorbar; rotate3d on;          
          
        case 'caprock temperature (C)'
          plotCellData(Gt.parent, compute_temperature(), 'edgealpha', ealpha);
          colorbar; rotate3d on;          
          
        case 'caprock pressure (MPa)'
          plotCellData(Gt.parent, compute_pressure() / 1e6, 'edgealpha', ealpha);
          colorbar; rotate3d on;          
          
        case 'caprock topography'
          colorbar('off'); rotate3d off; view(0, 90);
          mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
          
        case 'caprock co2 density (kg/m3)'
          t = opt.seafloor_temp + (Gt.cells.z - opt.seafloor_depth) .* opt.temp_gradient ./ 1000;
          p = (Gt.cells.z * opt.water_density * norm(gravity)) .* (1+opt.press_deviation/100);
          plotCellData(Gt.parent, var.co2.rho(p, t + 273.15), 'edgealpha', ealpha);
          colorbar; rotate3d on;
        case 'total capacity (tonnes/m2)'
          [~, strap, btrap_res, btrap_dis] = recompute_trapcap();
          tot_cap = (strap + btrap_res + btrap_dis) ./ Gt.cells.volumes ./ 1e3;
          plotCellData(Gt.parent, tot_cap, 'edgealpha', ealpha);
          colorbar; rotate3d on;
        case 'structural trap capacity (Mt)'
          trapcells = ta.traps~=0;
          [~, strap] = recompute_trapcap();
          trapcaps = accumarray(ta.traps(trapcells), strap(trapcells));
          trapcap_tot = ones(size(ta.traps)) * NaN;
          trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));
          
          plotGrid(Gt.parent, 'facecolor','none', 'edgealpha', ealpha);
          plotCellData(Gt.parent, trapcap_tot/1e3/1e6, 'edgecolor','none'); 
          colorbar; rotate3d on;

        case 'reachable structural capacity (Mt)'
          cumul_trap = compute_cumul_trapcap();
          plotCellData(Gt.parent, cumul_trap/1e9, 'edgealpha', ealpha);
          colorbar; rotate3d on;
      end
      %      plotGrid(var.Gt_cached.(var.current_formation).Gt);
   end

end

% ============================================================================

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
            'Ulafm', 'Utsirafm', 'Stofm', 'Nordmelafm', 'Tubaenfm', ...
            'Bjarmelandfm', 'Arefm', 'Garnfm', 'Ilefm', 'Tiljefm'};
end

% ----------------------------------------------------------------------------

function h = set_size(h, res)
% Utility function to resize a graphical window
   
   pos  = get(h, 'position');
   screensize = get(0,'screensize');
   res = min([res; screensize(3:4)-[0 85]]);
   set(h, 'position', [min([pos(1:2); screensize(3:4)-res-[0 85]]) res]);
   
end
