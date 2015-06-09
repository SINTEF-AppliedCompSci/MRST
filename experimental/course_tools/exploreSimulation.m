function exploreSimulation(varargin)
   
   gravity on;
   moduleCheck('ad-core');

   rhoCref = 760 * kilogram / meter ^3; % an (arbitrary) reference density
   
   opt.grid_coarsening = 4;
   opt.default_formation = 'Utsirafm';
   opt.window_size = [1200 900];
   opt.seafloor_depth = 100 * meter;
   opt.seafloor_temp  =  7; % in Celsius
   opt.temp_gradient  = 35.6; % degrees per kilometer
   opt.water_density  = 1000; % kg per m3
   opt.press_deviation = 0; % pressure devation from hydrostatic (percent)
   opt.res_sat_co2 = 0.21; 
   opt.res_sat_wat = 0.11;
   opt.dis_max = (53 * kilogram / meter^3) / rhoCref; % value from CO2store
   
   opt = merge_options(opt, varargin{:});

   var.Gt                = []; % to be set by the call to 'set_formation'
   var.rock2D            = [];
   var.ta                = [];
   var.current_formation = '';
   var.data              = [];
   var.co2               = CO2props('sharp_phase_boundary', false, 'rhofile', 'rho_demo');
   var.interaction_mode  = 'set bc'; % valids are 'set bc', 'set wells' and 'none'
   
   set_formation(opt.default_formation, false);

   %% Setting up interactive interface
   
   var.h = figure();%'KeyPressFcn', @(obj, e) parse_keypress(e.Key));
   set_size(var.h, opt.window_size(1), opt.window_size(2));
   
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
   
   %% Launching by calling redraw function
   redraw();
   
   
   % ============================= LOCAL HELPER FUNCTIONS =============================
   
   function click_handler(varargin)

      persistent segment_start = [];
      
      pt = get(gca,'CurrentPoint'); pt = pt(end,:); 

      switch var.interaction_mode
        case 'set bc'
        ix = closest_cell(var.Gt, pt, boundary_cells_and_faces(var.Gt));
        if isempty(segment_start)
           segment_start = ix;
        else
           face_ixs = shortest_path(segment_start, ix);
           segment_start = 0;
           
        
        otherwise
          disp('unimplemented');
          return;
      end
      redraw();
      
      
      % %disp(pts);
      % ix = closest_cell(var.Gt, pts(end,:), 1:var.Gt.cells.num);
      % var.data(ix) = 1;
      % redraw();
   end
   
   function redraw()
      axes(var.ax); cla;
      axis auto;
      cla;
      plotCellData(var.Gt, var.data, 'buttondownfcn', @click_handler);
      view(0, 90);
   end
      
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
      var.ta = trapAnalysis(var.Gt, false);
      
      % Call 'redraw' if requested
      if do_redraw
         redraw();
      end
      
      var.data = zeros(var.Gt.cells.num, 1);
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
            'Ulafm', 'Utsirafm'};     
end

% ----------------------------------------------------------------------------

function h = set_size(h, res_x, res_y)
% Utility function to resize a graphical window
   
   pos = get(h, 'position');
   set(h, 'position', [pos(1:2), res_x, res_y]);
   
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

function neighbor_mat = find_boundary_face_neighbors(Gt)

   [~, fix] = boundary_cells_and_faces(Gt); % boundary face indices
   nix = unique(Gt.faces.neighbors(fix,:)); % boundary node indices

   tmp = [(1:numel(fix))', Gt.faces.neighbors(fix, 1); ...
          (1:numel(fix))', Gt.faces.neighbors(fix, 2)];
   tmp = sortrows(tmp, 2);
   tmp = reshape(tmp(:,1), 2, []); % columns now express face neighborships
   
   neighbor_mat = sparse(tmp(1,:), tmp(2,:), 1, Gt.faces.num, Gt.faces.num);
   
end
