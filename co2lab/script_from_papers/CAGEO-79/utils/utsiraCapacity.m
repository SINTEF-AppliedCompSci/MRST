function utsiraCapacity(varargin)
% Compute total theoretical trapping capacity of the formation, taking structural,
% residual and dissolution trapping into account.

   moduleCheck('mex', 'libgeometry');
    
   % Setting default values
   [default_rho, dummy, default_sr, default_sw] = getValuesSPE134891();%#ok
   
   opt.rhoCfun        = @(p,T) default_rho(2);
   opt.rhoCref        = default_rho(2); % reference rho used for dissolution computations
   opt.rhoW           = default_rho(1);
   opt.sr             = default_sr;
   opt.sw             = default_sw;
   opt.dis_max        = (53 * kilogram / meter^3) / opt.rhoCref; % value from CO2store
   opt.tgrad          = 35.6; % degrees per kilometer
   opt.seafloor_temp  = 7; % in Celsius
   opt.seafloor_depth  = 80; 
   opt.path_mul        = 0;
   opt.well_num        = 10;
   opt.coarse_level    = 1;
   opt.boundary_margin = 20e3;
   opt.do_plot         = true;
   opt.wellpos         = [438516, 6471210];

   % Overriding with user parameters
   opt = merge_options(opt, varargin{:});
   
   % Loading Utsira grid and performing structural trapping analysis
   fprintf('Loading grid and petrophysical parameters ... ');
   [Gt, rock] = getFormationTopGrid('utsirafm', opt.coarse_level); 
   ta         = trapAnalysis(Gt, false);
   fprintf('done\n');
   
   % computing pressure and temperature fields (used for estimating densities)
   fprintf('Computing pressure and temperature ... ');
   gravity on;
   p = opt.rhoW .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure 
   t = 273.15 + opt.seafloor_temp + ...
       (Gt.cells.z - opt.seafloor_depth) .* opt.tgrad ./ 1000;
   fprintf('done\n');
   
   % computing structural trap heights (H1) for each cell
   fprintf('Computing structural trap heights ... ');
   trapcells     = find(ta.traps);
   H1            = zeros(Gt.cells.num, 1);
   H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
   H1=min(H1,Gt.cells.H);
   assert(all(H1<=Gt.cells.H));
   fprintf('done\n');
   
   % computing the height of the part of the column not in a structural trap
   H2 = Gt.cells.H - H1;
   
   % Computing total trapping volume in structural traps (dissolved and
   % structurally trapped
   fprintf('Computing total volume in structural traps ... ');
   strap_pvol_tot       = Gt.cells.volumes .* H1 .* rock.poro;
   strap_pvol_co2_plume = strap_pvol_tot .* (1 - opt.sw);
   strap_pvol_co2_diss  = strap_pvol_tot .* opt.sw .* opt.dis_max;
   
   strap_mass_co2 = strap_pvol_co2_plume .* opt.rhoCfun(p,t) + ...
                    strap_pvol_co2_diss  .* opt.rhoCref;
   fprintf('done\n');
   
   % Computing total trapping volume below structural traps (dissolved and
   % residually trapped
   fprintf('Computing trapping volume below structural traps ... ');
   btrap_pvol_tot          = Gt.cells.volumes .* H2 .* rock.poro;
   btrap_pvol_co2_residual = btrap_pvol_tot .* opt.sr;
   btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-opt.sr) .* opt.dis_max;
   
   btrap_mass_co2_res = btrap_pvol_co2_residual .* opt.rhoCfun(p,t);
   btrap_mass_co2_dis = btrap_pvol_co2_dissol   .* opt.rhoCref;
   fprintf('done\n');
   
   % Computing total trapping capacity per cell
   tot_trap_capa = strap_mass_co2 + btrap_mass_co2_res + btrap_mass_co2_dis;
   
   % Reporting trapping capacities:
   fprintf('Total trapping capacity: %f Gtons\n', sum(tot_trap_capa) / giga / 1e3);
   fprintf('Breakdown:\n');
   fprintf('Structural: %f Gtons\n', sum(strap_mass_co2)    / giga / 1e3);
   fprintf('Residual: %f Gtons\n', sum(btrap_mass_co2_res)  / giga / 1e3);
   fprintf('Dissolved: %f Gtons\n', sum(btrap_mass_co2_dis) / giga / 1e3);
   
   % Computing cumulative structural trap volume reached when injecting at a
   % given cell
   fprintf('Computing cumulative trap volume for given cell ... ');
   trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, 'removeOverlap', false);
   tvols = [trees.value];%#ok
   int_tr = find(ta.trap_regions); %#ok ixs of cells spilling into interior trap
   [dummy, reindex] = sort([trees.root], 'ascend');%#ok

   structural_mass_reached = zeros(Gt.cells.num, 1);   
   for i = 1:numel(ta.trap_z) % loop over each trap
       
       % ix of cells spilling directly into this trap
       cix = find(ta.trap_regions == i); 
        
       % cell indices of all cells of this trap, and its upstream traps
       aix = find(ismember(ta.traps, [trees(reindex(i)).traps]));
       
       % compute total structural trap capacity (in mass terms) of these
       % cells, and store result 
       structural_mass_reached(cix) = sum(strap_mass_co2(aix));%#ok
       
   end
   fprintf('done');

   if opt.do_plot
       structural_trapping_plot(Gt, ta, strap_pvol_tot, strap_mass_co2, opt.wellpos);
       reachable_multiplot(Gt, tot_trap_capa./Gt.cells.volumes, ...
                           structural_mass_reached, opt.wellpos);
   end
end

function plot_well(h, well, size, facecolor)
    hold on;
    if ~isempty(well)
       plot(well(1), well(2), 'ok', ...
            'markersize'      , 8         , ...
            'markerfacecolor' , facecolor , ...
            'markeredgecolor' , 'y');
    end
end

% ----------------------------------------------------------------------------
function structural_trapping_plot(Gt, ta, pvol, mass, wellpos)
    
    % plotting trapping structure
    figure;
    mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
    plot_well(gcf, wellpos, 12, 'k');
    set(gcf, 'position', [1 1 560 1048]);
    axis off					
    set(gcf, 'color', [1 1 1]);
    title('Traps and spill paths');		
    set(get(gca, 'title'), 'FontSize', 20);	
    set(get(gca, 'title'), 'FontWeight', 'bold');

    % plotting trapping capacity in mass terms
    figure;
    trapcells = ta.traps~=0;
    trapcaps = accumarray(ta.traps(trapcells), mass(trapcells));
    trapcap_tot = ones(size(ta.traps)) * NaN;
    trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));
    
    plotCellData(Gt, trapcap_tot/1e3/1e6, 'EdgeColor','none'); 
    cbh = colorbar;
    plotGrid(Gt, 'facecolor','none','edgealpha', 0.05);
    plot_well(gcf, wellpos, 8, 'k');
    title('Trap capacity (megatonnes)');
    axis off;
    set(cbh, 'fontSize', 20);
    set(gcf,'color', [1 1 1]); % white background
    set(gcf, 'position', [1 1 560 1048]);
    set(get(gca, 'title'), 'FontSize', 20);
    set(get(gca, 'title'), 'FontWeight', 'bold');
    set(gca, 'fontSize', 18);
    
end

function reachable_multiplot(Gt, tot, reachable_struct, wellpos)

    % plotting total capacity
    figure;
    plotCellData(Gt, tot/1e3, 'EdgeColor','none');
    plot_well(gcf, wellpos, 8, 'k');
    title({'Pillar capacity','(tonnes/m^2)'});
    axis off;
    set(gcf,'color', [1 1 1]); % white background
    set(gcf, 'position', [1 1 560 1048]);
    set(get(gca, 'title'), 'FontSize', 20);
    set(get(gca, 'title'), 'FontWeight', 'bold');
    set(gca, 'fontSize', 18);
    cbh = colorbar;
    set(cbh, 'fontSize', 20);
    
    % plotting total reachable structural capacity
    figure;
    plotCellData(Gt, reachable_struct/1e3/1e6,'EdgeColor','none');
    plot_well(gcf, wellpos, 8, 'w');
    title({'Reachable structural', 'capacity (megatonnes)'});
    axis off;
    set(gcf,'color', [1 1 1]); % white background
    set(gcf, 'position', [1 1 560 1048]);
    set(get(gca, 'title'), 'FontSize', 20);
    set(get(gca, 'title'), 'FontWeight', 'bold');
    set(gca, 'fontSize', 18);
    cbh = colorbar;
    set(cbh, 'fontSize', 20);
end

function c_ix = closest_gridcell(Gt, coords)

    dvec = bsxfun(@minus, Gt.cells.centroids, coords);
    dist = sum((dvec.^2), 2);
    [~, c_ix] = min(dist);
    
end
