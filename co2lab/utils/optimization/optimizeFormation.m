function [Gt, optim, init, history, other] = optimizeFormation(varargin)
%
% Optimize a injection scenario for a formation from the CO2-atlas.  Well
% positions are selected based on structural trapping potential, and rates
% are chosen to maximise storage while minimizing leakage.
% 
% SYNOPSIS:
%   function [Gt, optim, init, history, other] = optimizeFormation(varargin)
%
% PARAMETERS:
%   varargin - a number of paired arguments on the form: 'option', value.
%              Options include:
%              - modelname: name of the formation (default: 'utsirafm')
%              - schedule:     proposed initial schedule (leave empty for
%                              automatic specification of schedule with well
%                              placement and initial rates based on top
%                              surface geometry).
%              - coarse_level: downsampling level of formation grid (to
%                              lower computational cost)
%              - num_wells:    desired number of wells (default 10)
%              - subtrap_file: file containing information on estimated
%                              subscale trapping for the given formation
%                              (empty for no subscale trapping).
%              - sr:           residual CO2 saturation
%              - sw:           residual brine saturation
%              - isteps:       number of injection timesteps
%              - msteps:       number of migration timesteps
%              - itime:        total duration of injection phase
%              - mtime:        total duration of migration phase
%              - maxmise_boundary_distance (true/false): 
%                              if true, wells will be placed as far as
%                              possible from boundary within the chosen
%                              catchment region
%              - well_buffer_dist: 
%                              if 'maximize_boundary_distance' is fasle,
%                              wells should (if possible) be within this
%                              distnace of the relevant spill region's
%                              boundary. 
%              - report_basedir: directory for saving results
% RETURNS:
%   Gt      - Top surface grid of the chosen formation
%   optim   - Data structure containing all information about the final,
%             optimized scenario (including wells and the optimized schedule).
%   init    - Data structure containing all information about the starting
%             point scenario (before optimization of rates).
%   history - Information about the rate optimization process
%   other   - Data structure containing additional information needed to
%             rerun the optimized scenario (fluid and rock objects, residual
%             saturation, trapping structure, subscale trapping, and initial
%             state) 
% EXAMPLE:
%
%   For an example, refer to the sample script 'optimizeUtsira'.
%
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


   moduleCheck('ad-core');
   
   gravity on;
   opt = merge_options(opt_defaults(), varargin{:});
   
   %% Physical grid and rock
   [Gt, rock2D, ~] = getFormationTopGrid(opt.modelname, opt.coarse_level);
   
   %% Spill-point analysis objec
   ta = trapAnalysis(Gt, false);
   
   %% CO2 property object
   co2 = CO2props();
   
   %% Load subscale trapping function, if present
   dh = [];
   if ~isempty(opt.trapfile_name)
      if ~exist(opt.trapfile_name, 'file')
         warning('Subtrap datafile not found.  Proceeding without subtraps.');
      else
         dh = computeToptraps(load(opt.trapfile_name), Gt, true);
      end
   end

   %% Defining fluid
   
   % computing fixed temperature field
   T_ref = opt.ref_temp + opt.temp_grad * (Gt.cells.z - opt.ref_depth) / 1000;   
   
   fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                                , ...
                       'fixedT'       , T_ref                                       , ...
                       'wat_rho_pvt'  , [4.3e-5/barsa  , opt.refPress]              , ...
                       'residual'     , [opt.sw        , opt.sr]                    , ...
                       'wat_rho_ref'  , opt.rhoW                                    , ...
                       'co2_rho_ref'  , opt.refRhoCO2                               , ... 
                       'wat_rho_pvt'  , [opt.c_water   , opt.refPress]              , ...
                       'co2_rho_pvt'  , [opt.p_range   , opt.t_range]               , ...
                       'co2_mu_pvt'   , [opt.p_range   , opt.t_range]               , ...
                       'wat_mu_ref'   , opt.muBrine                                 , ...
                       'pvMult_fac'   , opt.pvMult                                  , ...
                       'dissolution'  , false                                       , ...
                       'pvMult_p_ref' , opt.refPress                                , ...
                       'surf_topo'    , if_else(isempty(dh), 'smooth', 'inf_rough') , ...
                       'top_trap'     , dh);
   
   %% Construct schedule
   use_default_schedule = isempty(opt.schedule);
   if use_default_schedule
      [wc, qt] = pick_wellsites(Gt, rock2D, co2, ta, opt.num_wells, opt.rhoW, ...
                                opt.sw, opt.ref_temp, opt.ref_depth, opt.temp_grad, ...
                                opt.well_buffer_dist, opt.maximise_boundary_distance); 
      opt.schedule = setSchedule(Gt, rock2D, wc, qt/opt.refRhoCO2, opt.isteps, ...
                                     opt.itime, opt.msteps , opt.mtime, true, ...
                                     'minval', sqrt(eps));
   else
      tmp = load(opt.schedule);
      name = fields(tmp);
      opt.schedule = tmp.(name{1});
   end
   
   % Add constant pressure boundary conditions
   bfaces = identifyBoundaryFaces(Gt);
   for i = 1:numel(opt.schedule.control)
      opt.schedule.control(i).bc = addBC([], bfaces, 'pressure', ...
                                         Gt.faces.z(bfaces) * fluid.rhoWS * ...
                                         norm(gravity) + opt.surface_pressure, ...
                                         'sat', [1 0]);
   end
   
   %% Defining initial state
   initState.pressure = compute_hydrostatic_pressure(Gt, fluid.rhoWS, opt.surface_pressure);
   initState.s = repmat([1 0], Gt.cells.num, 1);
   initState.sGmax = initState.s(:,2);
   
   
   other.fluid = fluid;
   other.rock = rock2D;
   other.residual = [opt.sw, opt.sr];
   other.traps = ta;
   other.dh = dh;
   other.initState = initState;

   %% Set up model and run optimization
   model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
   min_rates = sqrt(eps) * ones(opt.num_wells, 1);
   max_rates = 3 * max([opt.schedule.control(1).W.val]) * ones(opt.num_wells, 1);
   [optim, init, history] = ...
       optimizeControls(initState, model, opt.schedule, min_rates, max_rates, ...
                     'last_control_is_migration', true);  
end

% ----------------------------------------------------------------------------
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end

% ----------------------------------------------------------------------------
function wc = select_wellcell(Gt, candidates, buffer)

   % Identify cells at the boundary of the selected region
   num    = Gt.cells.num;
   n_rels = Gt.faces.neighbors;
   n_rels = double(n_rels(~any(n_rels == 0,2),:)); % remove exterior face relations
   adj    = sparse(n_rels(:,1), n_rels(:,2), 1, num, num, 2 * size(n_rels,1));
   adj    = adj + adj';
   boundary_candidates = candidates(sum(adj(candidates, candidates)) < 4);
   interior_candidates = setdiff(candidates, boundary_candidates);
   
   % compute smallest distance for each interior candidate to boundary
   distances = inf(numel(interior_candidates), 1);
   for bix = boundary_candidates(:)'
      dvec = bsxfun(@minus, ...
                    Gt.cells.centroids(bix, 1:2), ...
                    Gt.cells.centroids(interior_candidates, 1:2));
      dist = sqrt(sum(dvec.^2, 2));
      distances = min(distances, dist);
   end
   
   inside_buffer_candidates = interior_candidates(distances > buffer);

   if isempty(inside_buffer_candidates)
      % none of the candidates were far enough from boundary.  We will just
      % choose the highest point
      candidate_z = Gt.cells.z(candidates);
      [~, min_ix] = min(candidate_z);
      wc = candidates(min_ix);
   else
   
      inside_buffer_z = Gt.cells.z(inside_buffer_candidates);
      [~, max_z_ix] = max(inside_buffer_z);
      
      wc = inside_buffer_candidates(max_z_ix);
   end
end

% ----------------------------------------------------------------------------
function [wc, qt] = pick_wellsites(Gt, rock2D, co2, ta, num_wells, ...
                                   rhoW, sw, seafloor_temp, seafloor_depth, ...
                                   tgrad, buf_dist, maximize_boundary_distance)

    % Computing local pressure and temperature conditions (in order to
    % compute CO2 densities) 
    P = rhoW .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure
    T = seafloor_temp + (Gt.cells.z - seafloor_depth) .* tgrad ./ 1000;
    
    % Computing local CO2 densities
    rhoCO2 = co2.rho(P, T);

    % Computing trap volumes
    tcells               = find(ta.traps);
    trap_heights         = zeros(Gt.cells.num, 1);
    trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
    trap_heights         = min(trap_heights,Gt.cells.H);
    strap_vol            = Gt.cells.volumes .* trap_heights .* rock2D.poro .* (1-sw);

    assert(all(trap_heights <= Gt.cells.H));
    
    % Computing trap capacities in mass terms and suggesting wellsites
    trapcap = accumarray(ta.traps(tcells), strap_vol(tcells) .* rhoCO2(tcells));
    
    % Suggest wellsites based on maximum cumulative trap capacity reached

    bdist = sqr_distance_from_boundary(Gt); % used below to discourage placement
    bdist = bdist/max(bdist(:));            % of wells close to boundary

    wc = zeros(num_wells, 1); 
    qt = wc;

    for i = 1:num_wells
        cumcap = [0; cumulative_capacity(trapcap, ta.trap_adj)];
        field  = cumcap(ta.trap_regions + 1);

        if maximize_boundary_distance
           % Pick well site as far as possible away from boundary.
        
           % Adding an insignificant quantity so that if everything else is
           % equal, the cell farthest from the boundary is chosen
           grid_eps = 10 * eps(max(cumcap(:)));
           field    = field + grid_eps * bdist;
           [qt(i), wc(i)] = max(field);
        else
           % Pick well site as far downslope as possible (while keeping a
           % buffer distance from boundary
           
           qt(i) = max(field); 
           candidate_cells = find(field == qt(i));
           wc(i) = select_wellcell(Gt, candidate_cells, buf_dist);
        end
        
        % Setting available capacity of used traps to zero
        trapcap = set_used_to_zero(trapcap, ta.trap_adj, ta.trap_regions(wc(i)));
    end

end

% ----------------------------------------------------------------------------

function cumcap = cumulative_capacity(trapcap, adj)
    % Compute cumulative capacity of each trap and those it is connected to upstream
    cumcap = trapcap;
    adj_orig = adj;
    while nnz(adj)
        cumcap = cumcap + adj * trapcap;
        adj = spones(adj * adj_orig);
    end
end


% ----------------------------------------------------------------------------

function trapcap = set_used_to_zero(trapcap, adj, trap_ix)
    % Set capacity of trap with index 'trap_ix' and all its ustream traps to zero
    trapcap(trap_ix) = 0;
    ustr_traps = find(adj(trap_ix,:));
    while ~isempty(ustr_traps)
        trapcap(ustr_traps) = 0;
        ustr_traps = find(sum(adj(ustr_traps, :), 1));
    end
end

% ----------------------------------------------------------------------------

function opt = opt_defaults()

    opt.modelname = 'utsirafm';
    opt.schedule = [];
    opt.coarse_level = 3;
    opt.num_wells = 10;
    opt.subtrap_file = 'utsira_subtrap_function_3.mat';

    opt.surface_pressure = 1 * atm;
    
    % Default fluid properties
    [rho_default, mu_default, sr, sw] = getValuesSPE134891();
    opt.refRhoCO2 = rho_default(2);
    opt.rhoW      = rho_default(1);
    opt.muBrine   = mu_default(1);
    opt.muCO2     = 0; % zero value will trigger use of variable viscosity 
    opt.pvMult    = 1e-5/barsa;
    opt.refPress  = 100 * barsa;
    opt.c_water   = 4.3e-5/barsa; % water compressibility
    opt.p_range   = [0.1, 400] * mega * Pascal; % pressure span of sampled property table
    opt.t_range   = [4 250] + 274; % temperature span of sampled property table 

    % residual saturations
    opt.sr = sr; % gas
    opt.sw = sw; % brine
    
    % number of timesteps (injection and migration)
    opt.isteps = 10;
    opt.msteps = 31;

    % durations of injection an migration phases
    opt.itime  = 50   * year;
    opt.mtime  = 3000 * year;
    
    % Temperature regime control parameters
    opt.ref_temp = 7 + 273.15;
    opt.ref_depth = 80; 
    opt.temp_grad = 35.6;

    % dissolution parameters (default zero - suggested nonzero values in comments)
    opt.dis_rate = 0;  % %0.44 * kilogram / rho / poro / (meter^2) / year = 8.6924e-11;
    opt.dis_max = 0;   % 53/760 = 0.07

    % factor for setting upper limit of injection rates
    opt.lim_fac = 10;

    
    % Well selection strategy
    opt.maximise_boundary_distance = false; % If true, wells will be placed
                                            % as far as possible from
                                            % boundary within the chosen
                                            % catchment region
    opt.well_buffer_dist = 5 * kilo * meter; % If 'maximise_boundary_distance' is
                                             % false, wells should if possible
                                             % be within this distance of spill
                                             % region boundary

    % directory for saving reslts
    opt.report_basedir = './simulateUtsira_results/';
    opt.trapfile_name = ''; % 'utsira_subtrap_function_3'
end

% ----------------------------------------------------------------------------

function p = compute_hydrostatic_pressure(Gt, rho_water, surface_pressure)
    
    p = rho_water * norm(gravity()) * Gt.cells.z + surface_pressure;
end

% ----------------------------------------------------------------------------

function [field, closest_bnd_edge] = sqr_distance_from_boundary(Gt)

    field = inf * ones(Gt.cells.num, 1);
    closest_bnd_edge = zeros(Gt.cells.num, 1);
    
    % Looping over boundary edges
    bnd_edges = find(prod(Gt.faces.neighbors, 2) == 0);
    for b = bnd_edges'
        dx = Gt.cells.centroids(:,1) - Gt.faces.centroids(b,1);
        dy = Gt.cells.centroids(:,2) - Gt.faces.centroids(b,2);
        
        d2 = dx.*dx + dy.*dy;
        replace = (d2 < field);
        field(replace) = d2(replace);
        closest_bnd_edge(replace) = b;
    end
    
 end
 
 % ----------------------------------------------------------------------------
 function res = if_else(cond, yes, no)
    if cond
       res = yes;
    else
       res = no;
    end
 end
 