function [ will_stay, will_leak ] = massAtInfinity( Gt, rock, p, sG, sGmax, sF, rs, fluid, ta, dh, varargin )
% Forecast amount of co2 (in kg) to remain in formation by time infinity.
%
% SYNOPSIS:
%   future_mass = massAtInfinity(Gt, rock, p, sG, sGmax, sF, rs, fluid, ta, dh)
%   future_mass = massAtInfinity(Gt, rock, p, sG, sGmax, sF, rs, fluid, ta, dh, ...
%                  'p_future', pf, 'surface_pressure', ps, 'plotsOn', true)
%
% DESCRIPTION:
%   This calculation of co2 mass at time infinity is based on the
%   formations' co2 mass at a point in time in which it is safe to assume
%   the flow dynamics are gravity-dominated. Thus, the co2 mass in a
%   particular spill tree at time infinity is based on that spill tree's
%   capacity, and any extra co2 mass is assumed to have been leaked by time
%   infinity.
%
%   NB: this routine was developed to handle ADI variables (p, sG, sGmax,
%   rs, will_stay, etc.), thus the syntax used here respects ADI variables,
%   and in some instances a cell array of ADI variables.
%
%
% INPUTS    - Gt, top surface grid
%           - rock, includes porosity (and possibly ntg) data
%           - p, current pressure (possible ADI var)
%           - sG, sGmax, sF (possible ADI vars)
%           - fluid structure containing:
%               - res_water, res_gas, rhoGS, rhoWS, bG(p), pvMultR(p)
%               - (NB: fluid is assumed to be same structure used to get
%               state results by simulateSchedule)
%           - ta, trapping structure obtained using cell-based method
%           - dh, sub-trapping (if applicable, otherwise can be empty, [])
%
% (optional) - plotsOn, true or false for plotting
%            - surface_pressure (used in calculation of p_future)
%            - p_future (default is hydrostatic pressure)
%
% RETURNS   - co2 mass in terms of:
%               1. amount forecast to remain at time infinity (will_stay)
%               2. amount forecast to leak at time infinity (will_leak)

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

   moduleCheck('ad-core', 'ad-props', 'ad-blackoil')

   opt.plotsOn             = false;    % controls whether plots will be made
   opt.h                   = 48;       % figure handle for plotting
   opt.surface_pressure    = 1 * atm;
   opt.p_future            = [];       % hydrostatic unless otherwise specified
   opt.report_trap_regions = [];       % specify specific catchments numbers
                                       % for which to print reports, e.g.,
                                       % [37, 45, 48, 53, 52]; 
   opt.tot_inj             = [];       % total injected (in kg)
   opt = merge_options(opt, varargin{:});

   if isempty(opt.p_future)
      opt.p_future = fluid.rhoWS * norm(gravity()) * Gt.cells.z + opt.surface_pressure;
   end
   %opt.p_future = state.pressure; % @@ for de-bugging

   % Renaming variables:
   sw      = fluid.res_water;
   sr      = fluid.res_gas;
   poro    = rock.poro;
   ntg     = ones(Gt.cells.num,1);
   if isfield(rock,'ntg')
      ntg = rock.ntg;
   end
   pvMult_future = fluid.pvMultR(opt.p_future);
   bG_future     = fluid.bG(opt.p_future);
   

   %% 1. Get masses (totals and per cell)
   % NB: use current pressure, p

   [h, h_max]         = upscaledSat2height(sG, sGmax, Gt, 'resSat', [sw sr]);  % can be ADI
   [masses, masses_0] = massTrappingDistributionVEADI(Gt, ...
                                                     p, sG, sF, h, h_max, rock, ...
                                                     fluid, ta, dh, 'rs', rs);
   % Renaming:
   resDis_0    = masses_0{1};
   resStruc_0  = masses_0{2};
   resTrap_0   = masses_0{3};
   freeRes_0   = masses_0{4};
   freeStruc_0 = masses_0{5};
   subtrap_0   = masses_0{6};
   freeMov_0   = masses_0{7};
   
   % Report for testing purposes:
   tmp = resDis_0 + resStruc_0 + resTrap_0 + freeRes_0 + freeStruc_0 + subtrap_0 + freeMov_0;
   dispif(mrstVerbose, '\n------------------\n')
   dispif(mrstVerbose,   'Initial inventory:' )
   dispif(mrstVerbose, '\n------------------\n')
   dispif(mrstVerbose, '\nThere is a total of %5.3f Mt CO2 in domain.\n', sum(tmp)/1e9 )
   dispif(mrstVerbose, '\nThe boundary catchments contain %5.3f Mt CO2 (%3.3f %% of tot).\n', ...
          sum( tmp(ta.trap_regions==0) )/1e9, sum( tmp(ta.trap_regions==0) )./sum(tmp)*100 )
   dispif(mrstVerbose, '\nThe spill-path catchments contain %5.3f Mt (%3.3f %% of tot) CO2,\n', ...
          sum( tmp(ta.trap_regions~=0) )/1e9, sum( tmp(ta.trap_regions~=0) )./sum(tmp)*100 )
   dispif(mrstVerbose, 'of which %5.3f Mt (%3.3f %% of tot) is in traps.\n\n', ...
          sum( tmp(ta.traps~=0) )/1e9, sum( tmp(ta.traps~=0) )./sum(tmp)*100 )
   
   % Reports of specific catchments (and boundary catchment by default):
   if ~isempty(opt.report_trap_regions)
      treg = opt.report_trap_regions;
      if isempty(opt.tot_inj)
         for i=1:numel(treg)
            mass_in_region = sum( tmp(ta.trap_regions == treg(i)) )/1e9; % Mt
            fprintf('\nThere is %5.3f Mt in region %3.0f.\n', mass_in_region, treg(i) );
         end
         mass_in_bdry_region = sum( tmp(ta.trap_regions==0) )/1e9; % Mt
         fprintf('\nThere is %5.3f Mt in boundary catchments.\n', mass_in_bdry_region );
      else
         % give report as percentage of tot inj
         TINJ = opt.tot_inj; % kg
         fprintf('\nA total of %5.3f Mt CO2 was injected.\n', TINJ/1e9 )
         fprintf('\nA total of %5.3f Mt (%3.3f %% of tot inj.) CO2 exists in domain now.\n', ...
                 sum(tmp)/1e9, sum(tmp)/TINJ*100 )
         fprintf('\n     The boundary catchments contain %5.3f Mt (%3.3f %% of tot inj.).\n', ...
                 sum( tmp(ta.trap_regions==0) )/1e9, sum( tmp(ta.trap_regions==0) )/TINJ*100 )
         fprintf('\n     The spill-path catchments contain %5.3f Mt (%3.3f %% of tot inj.),\n', ...
                 sum( tmp(ta.trap_regions~=0) )/1e9, sum( tmp(ta.trap_regions~=0) )/TINJ*100 )
         fprintf('       of which %5.3f Mt (%3.3f %% of tot inj.) is in traps.\n', ...
                 sum( tmp(ta.traps~=0) )/1e9, sum( tmp(ta.traps~=0) )/TINJ*100 )
         for i=1:numel(treg)
            mass_in_region = sum( tmp(ta.trap_regions == treg(i)) ); % kg
            fprintf('\nThere is %5.3f Mt (%3.3f %% of tot inj.) in region %3.0f.\n', ...
                    mass_in_region/1e9, mass_in_region/TINJ*100, treg(i) );
         end
         mass_in_bdry_region = sum( tmp(ta.trap_regions==0) ); % kg
         fprintf('\nThere is %5.3f Mt (%3.3f %% of tot inj.) in boundary catchments.\n', ...
                 mass_in_bdry_region/1e9, mass_in_bdry_region/TINJ*100 );
      end
   end
   

   %% 2. Get capacities
   % use pvMult_future and bG_future, which are functions of future pressure

   % Computing structural trap heights (H1) for each cell
   tcells               = find(ta.traps);
   trap_heights         = zeros(Gt.cells.num, 1);
   trap_heights(tcells) = ta.trap_z(ta.traps(tcells)) - Gt.cells.z(tcells);
   trap_heights         = min(trap_heights, Gt.cells.H);
   H2                   = Gt.cells.H - trap_heights;
   assert(all(trap_heights <= Gt.cells.H));

   % Computing structural trap pore vol (m3)
   strap_pvol           = Gt.cells.volumes .* trap_heights .* poro .* ntg .* pvMult_future;

   % Computing co2 capacity in structural traps (m3 co2, at ref. depth)
   strap_co2_vol        = strap_pvol .* (1-sw) .* bG_future;               % per cell
   num_traps            = max(unique(ta.traps));
   trapcap_vol          = cell(num_traps,1);
   for i = 1:num_traps
       trapcap_vol{i}   = sum(strap_co2_vol(ta.traps == i));               % per trap
   end

   % Computing pore vol outside of straps (m3)
   btrap_pvol         = Gt.cells.volumes .* H2 .* poro .* ntg .* pvMult_future;

   % computing capacity for residual co2 (m3 at ref. depth) below straps
   btrap_res_co2_vol  = btrap_pvol .* sr .* bG_future;


   %% 2. Initial accounting of masses

   % The mass currently in structural traps is: freeStruc_0 + resStruc_0

   % Available capacity in straps to take on more co2
   avail_trapcap = (strap_co2_vol .* fluid.rhoGS) - freeStruc_0 - resStruc_0; % kg
   avail_trapcap_per_trap = cell(num_traps,1);
   for i = 1:num_traps
      avail_trapcap_per_trap{i} = sum(avail_trapcap(ta.traps == i));
   end

   % Movable masses below straps and within catchments/trap regions
   spillable = freeMov_0;
   spillable( ta.trap_regions == 0 ) = 0;
   spillable_per_trap   = cell(num_traps,1);
   for i = 1:num_traps
       spillable_per_trap{i} = sum(spillable(ta.trap_regions == i));       % per trap region
   end
   
   % Movable masses outside catchments/trap regions
   leaked_from_outside_catchments = freeMov_0;
   leaked_from_outside_catchments( ta.trap_regions ~=0 ) = 0;
   
   assert( all((value(freeMov_0) - value(leaked_from_outside_catchments) - value(spillable)) == 0) )


   %% 3. Spill masses along trees
   % We do this by spilling any masses a trap cannot retain down the tree,
   % and out of the domain.
   
   % Get the trap index of the root trap of each tree
   
   trees = maximizeTrapping(Gt, 'res', ta); % 'calculateAll', 'removeOverlap' @@
   treeRoots = [trees.root];
   clear trees
   
   % Compare "avail_trapcap" to "spillable", keeping track of:
   % a) mass_added_per_trap_region
   % b) mass_leaked_per_tree
   
   % Initialize
   mass_leaked_per_tree = cell(numel(treeRoots),1);
   mass_added_per_trap_region = spillable_per_trap;                 % can be ADI
   
   dispif(mrstVerbose, '\n----------------------------\n');
   dispif(mrstVerbose,   'Spilling masses along trees:' );
   dispif(mrstVerbose, '\n----------------------------\n');
   dispif(mrstVerbose, '\nThere are %d trees.\n', numel(treeRoots));
   for i = 1:numel(treeRoots)
      dispif(mrstVerbose, '\nLooking at tree %d... \n', i);

      % Get the trap indexes of all traps downstream from a tree
      % root, including the trap index of the root itself. (NB: While
      % trap_vols is not required, we use its value to compare
      % against the trapcap computed above. We pass in poro * ntg *
      % pvMult_future as the porosity, to account for the impact of
      % pressure and possible ntg data.)
      [trap_ixs, trap_vols] = ...
          downstreamTraps(Gt, poro.*ntg.*pvMult_future.*bG_future, ta, treeRoots(i));
      
      % Compare trap_vols (m3 pore space of traps) with trapcap_vol
      % (m3 reservoir co2). NB: we must account for residual water (sw).
      trap_co2_vols = trap_vols' .* (1-sw);
      trapcap_tmp = cellfun( @(x) value(x), {trapcap_vol{trap_ixs}})';
      assert( all( abs( trap_co2_vols - trapcap_tmp )./abs(trap_co2_vols) < 1e-8 ) ) % syntax handles cell array trapcap
      
      % Spill masses along a tree...
      dispif(mrstVerbose, 'The trap indexes in this tree are: %s \n\n', num2str(trap_ixs));
      for j = 1:numel(trap_ixs)
         dispif(mrstVerbose, ' Looking at trap %d... \n', trap_ixs(j));
         dispif(mrstVerbose, '  Trap %d has an available capacity of %d (kg co2). \n', ...
                trap_ixs(j), value(avail_trapcap_per_trap{trap_ixs(j)}) );
         dispif(mrstVerbose, '  Trap region %d currently has %d (kg co2) spillable. \n', ...
                trap_ixs(j), value(mass_added_per_trap_region{trap_ixs(j)}) );
         
         
         % Computing the mass that will spill out of trap_ixs(j),
         spill_mass = max(0, mass_added_per_trap_region{trap_ixs(j)} - ...
                          avail_trapcap_per_trap{trap_ixs(j)}); % may be an ADI var 
         assert( spill_mass >= 0 )
         
         % and removing this mass from trap_ixs(j),
         mass_added_per_trap_region{trap_ixs(j)} = mass_added_per_trap_region{trap_ixs(j)} - spill_mass; % may be an ADI var
         
         % and spilling this "spill mass" into the next trap
         % downstream or out of the tree:
         if j < numel(trap_ixs)
            % Haven't reach end of spill-tree yet, so add this
            % spill mass to the next trap down the spill-tree
            % (i.e., trap_ixs(j+1))
            dispif(mrstVerbose, '  Thus %d (kg co2) spills out of trap %d, and into trap %d.\n', ...
                   value(spill_mass), trap_ixs(j), trap_ixs(j+1));
            mass_added_per_trap_region{trap_ixs(j+1)} = mass_added_per_trap_region{trap_ixs(j+1)} + spill_mass; % may be an ADI var
            dispif(mrstVerbose, '  Trap %d now has %d (kg co2). \n', ...
                   trap_ixs(j+1), value(mass_added_per_trap_region{trap_ixs(j+1)}) );
         else
            % We've reached the end of the spill-tree, thus the
            % spill_mass goes out of the spill-tree (i.e., out of
            % the domain)
            dispif(mrstVerbose, '  Thus %d (kg co2) spills out of trap %d, and out of tree %d.\n', ...
                   value(spill_mass), trap_ixs(j), i);
            mass_leaked_per_tree{i} = spill_mass;                     % may be an ADI var
            
         end
      end
   end
   
   assert( abs( sum_cellArray_of_ADI(spillable_per_trap) - ...
                sum_cellArray_of_ADI(mass_added_per_trap_region) - ...
                sum_cellArray_of_ADI(mass_leaked_per_tree) ) < 1e-1 ) 
   
   % CO2 ultimately remaining in the structural traps
   % (summation taken in the following way to preserve possible ADI)
   mass_added = mass_added_per_trap_region{1};
   for i = 2:numel(mass_added_per_trap_region)
      mass_added = mass_added + mass_added_per_trap_region{i}; % can be ADI
   end
   
   assert( abs( sum_cellArray_of_ADI(spillable_per_trap) - ...
                (mass_added + sum_cellArray_of_ADI(mass_leaked_per_tree)) ) < 1e-1 )


   %% 4. Final accounting

   % Predicted mass of co2 in straps:
   stay_straps = mass_added + sum(freeStruc_0) + sum(resStruc_0);

   % Conservative estimate of co2 residually trapped:
   stay_resid = sum(resTrap_0 + freeRes_0);

   % Predicted amount of leaked co2:
   % (ADI not preserved)
   will_leak = sum(value(leaked_from_outside_catchments)) + ...
               sum_cellArray_of_ADI(mass_leaked_per_tree); 

   % Varargout:
   will_stay = stay_straps + stay_resid; % kg
   
   % we ensure balance of mass inventory within tolerance of 0.1% of total
   % mass
   assert( abs( sum(masses) - ( value(will_stay) + will_leak ) ) < 1e-3 * sum(masses), ...
       'There is a mismatch between mass calculations')
   
   %% 5. Reporting (optional)

   dispif(mrstVerbose, '\n----------------\n');
   dispif(mrstVerbose,   'Final inventory:'  );
   dispif(mrstVerbose, '\n----------------\n');
   dispif(mrstVerbose, ...
          '\nThe total structural trapping capacity is %5.3f Mt co2. \n', ...
          sum(strap_co2_vol)*fluid.rhoGS/1e9); 
   dispif(mrstVerbose, ...
          'The total residual trapping capacity is %5.3f Mt co2. \n\n', ...
          sum(btrap_res_co2_vol)*fluid.rhoGS/1e9);
   
   dispif(mrstVerbose, ...
          '%5.3f Mt co2 is currently within formation. \n', sum(masses)/1e9);
   dispif(mrstVerbose, ...
          '   --- %5.3f Mt will remain: \n', (value(stay_straps) + value(stay_resid))/1e9)
   dispif(mrstVerbose, ...
          '       --- %5.3f Mt structurally (%5.3f perc. of strap capacity) \n', ...
          value(stay_straps)/1e9, value(stay_straps)/(sum(strap_co2_vol)*fluid.rhoGS) * 100);
   dispif(mrstVerbose, ...
          '       --- %5.3f Mt residually (%5.3f perc. of resid capacity) \n', ...
          value(stay_resid)/1e9, value(stay_resid)/(sum(btrap_res_co2_vol)*fluid.rhoGS) * 100);
   dispif(mrstVerbose, ...
          '   --- %5.3f Mt will leak: \n', will_leak/1e9);
   dispif(mrstVerbose, ...
          '       --- %5.3f Mt from outside catchments \n', sum(value(leaked_from_outside_catchments))/1e9);
   dispif(mrstVerbose, ...
          '       --- %5.3f Mt from trees \n', sum_cellArray_of_ADI(mass_leaked_per_tree)/1e9);
   
   
   %% 6. Plotting (optional)

   if opt.plotsOn
      % Pie to show breakdown of leaked and remaining:
      will_leak_out_catch = sum(value(leaked_from_outside_catchments));
      will_leak_via_trees = sum_cellArray_of_ADI(mass_leaked_per_tree);
      if ~ishandle(opt.h)
         figure(opt.h);
      else
         % Avoid stealing focus if figure already exists
         set(0, 'CurrentFigure', opt.h); clf(opt.h)
      end
      % NB: max(X,eps) is used as hack to make pie chart plot "0" values
      ph = pie(max([  will_leak_out_catch, ...
                      will_leak_via_trees, ...
                      value(stay_straps), ...
                      value(stay_resid)       ], eps));
      % A work-around way to display with decimal point precision
      tots = will_leak_out_catch + will_leak_via_trees + ...
             value(stay_straps) + value(stay_resid);
      chld = get(gca,'Children');
      set(chld(7), 'String', [num2str(will_leak_out_catch/tots*100,'%2.1f'), '%']);
      if will_leak_via_trees/tots*100 < 0.1
         set(chld(5), 'String', ' ');
      else
         set(chld(5), 'String', [num2str(will_leak_via_trees/tots*100,'%2.1f'), '%']);
      end
      set(chld(3), 'String', [num2str(value(stay_straps)/tots*100,'%2.1f'), '%']);
      set(chld(1), 'String', [num2str(value(stay_resid)/tots*100,'%2.1f'), '%']);
      set(chld(8),'FaceColor','m');
      set(chld(6),'FaceColor','r');
      set(chld(4),'FaceColor','y');
      set(chld(2),'FaceColor','g');
      lh = legend('will leak outside catchments','will leak along trees',...
                  'will remain (in straps)','will remain (residually)');
      set(lh,'Location','SouthEast')
      title(['Future of ',num2str(sum(masses)/1e9),' Mt co2'])
      %drawnow;
   end
   
   
end


function total = sum_cellArray_of_ADI( ADI_cellArray )
   total = 0;
   for i = 1:numel(ADI_cellArray)
      total = total + value(ADI_cellArray{i});
   end
end
