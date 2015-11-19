function [ res ] = getTrappingInfo(name, coarsening, varargin)
% Computes structural trapping info and makes plots (optional).
%
% SYNOPSIS
%   res = getTrappingInfo('Stofm',3)
%   res = getTrappingInfo('Stofm',3, 'plotsOn',false)
%
% DESCRIPTION
%   Parameter values (seafloor_depth, etc.) are obtained using
%   getSeaInfo(), which returns values that correspond to the specified
%   formation. Then upper bound (theoretical) trapping capacities are
%   computed and returned. Plotting is optional.
%
% See also
%   exploreCapacity.m
    

    opt.plotsOn = true;
    opt = merge_options(opt, varargin{:});
    
    
    %% Arbitrary reference co2 density
    rhoCref = 760 * kilogram / meter ^3;
    
    
    %% Get formation grid and info
    [var.Gt, var.rock2D] = getFormationTopGrid(name, coarsening);
    var.ta               = trapAnalysis(var.Gt,'false');
    var.co2              = CO2props('sharp_phase_boundary', true, ...
                                    'rhofile', 'rho_demo');
    info                 = getSeaInfo(name, rhoCref);
    
    
    %% make plots (optional) and get trapping capacity results (res)
    draw();
    
    res.Gt      = var.Gt;
    res.rock2D  = var.rock2D;
    res.ta      = var.ta;

    res.info    = info;
    res.co2     = var.co2; % NB: co2.rho(p,t + 273.15) where t is Celsius
    


    % ---------------------------------------------------------------------
    % ------------------------- Helper Functions --------------------------
    % ---------------------------------------------------------------------

    function p = compute_pressure()
        Gt = var.Gt;
        p = (Gt.cells.z * info.water_density * norm(gravity)) ...
            .* (1+info.press_deviation/100);
    end
    % ---------------------------------------------------------------------

    function t = compute_temperature()
        % NB: returns value in Celsius, not Kelvin
        Gt = var.Gt;
        t = info.seafloor_temp + (Gt.cells.z - info.seafloor_depth) ...
            .* info.temp_gradient ./ 1000;
    end

    % ---------------------------------------------------------------------

    function cum_cap = compute_cumul_trapcap()
        % returns kg

        [~, strap] = recompute_trapcap(); % kg

        Gt = var.Gt;
        ta = var.ta;
        trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, ...
                                'removeOverlap', false);
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

    % ---------------------------------------------------------------------

    function [H1, strap, btrap_res, btrap_dis] = recompute_trapcap()
        p = compute_pressure();
        t = compute_temperature() + 273.15;

        Gt   = var.Gt;
        ta   = var.ta;
        poro = var.rock2D.poro;

        % computing structural trap heights (H1) for each cell
        trapcells     = find(ta.traps);
        H1            = zeros(Gt.cells.num, 1);
        if ~isempty(trapcells)
         H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
        end
        H1=min(H1,Gt.cells.H);
        H2 = Gt.cells.H - H1;
        assert(all(H1<=Gt.cells.H));

        % Computing total trapping volume in structural traps (dissolved
        % and structurally trapped
        strap_pvol_tot       = Gt.cells.volumes .* H1 .* poro;
        strap_pvol_co2_plume = strap_pvol_tot .* (1 - info.res_sat_wat);
        strap_pvol_co2_diss  = strap_pvol_tot .* info.res_sat_wat .* info.dis_max;

        strap = strap_pvol_co2_plume .* var.co2.rho(p,t) + ...
              strap_pvol_co2_diss  .* rhoCref;

        % Computing total trapping volume below structural traps (dissolved
        % and residually trapped
        btrap_pvol_tot          = Gt.cells.volumes .* H2 .* poro;
        btrap_pvol_co2_residual = btrap_pvol_tot .* info.res_sat_co2;
        btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-info.res_sat_co2) .* info.dis_max;

        btrap_res = btrap_pvol_co2_residual .* var.co2.rho(p,t);
        btrap_dis = btrap_pvol_co2_dissol   .* rhoCref;
    end

    % ---------------------------------------------------------------------

    function str = capacity_report()
        [~, strap, btrap_res, btrap_dis] = recompute_trapcap();

        tot_trap_capa = strap + btrap_res + btrap_dis;
        sum_tot = sum(tot_trap_capa);

        str = sprintf('Total trapping capacity: %.2f Gtons\n\n', ...
            sum_tot / giga / 1e3);
        str = [str, sprintf('Breakdown:\n')];
        str = [str, sprintf('Structural: %.2f Gtons (%5.2f%%)\n', ...
            sum(strap) / giga / 1e3, (sum(strap) / sum_tot) * 100)];
        str = [str, sprintf('Residual: %.2f Gtons (%5.2f%%)\n', ...
            sum(btrap_res)  / giga / 1e3, sum(btrap_res / sum_tot) * 100)];
        str = [str, sprintf('Dissolved: %.2f Gtons (%5.2f%%)\n', ...
            sum(btrap_dis) / giga / 1e3, sum(btrap_dis / sum_tot) * 100)];

    end

    % ---------------------------------------------------------------------

    function draw(varargin) 
        if opt.plotsOn
          figure; set(gcf,'Position',[1 1 1402 1144])
          [nr,nc] = deal(3);
        end
        Gt = var.Gt;
        ta = var.ta;


        % Calculations:

        % reservoir conditions
        t = compute_temperature();% Celsius
        p = compute_pressure();   % Pascals

        % total capacity (kg/m2): (volume is top-surface cell area)
        [~, strap, btrap_res, btrap_dis] = recompute_trapcap();
        tot_cap = (strap + btrap_res + btrap_dis) ./ Gt.cells.volumes; % kg/m2
        
        % structural trap capacity (kg)
        trapcells = ta.traps~=0;
        [~, strap] = recompute_trapcap(); % redundant?
        trapcaps = accumarray(ta.traps(trapcells), strap(trapcells));
        trapcap_tot = ones(size(ta.traps)) * NaN;
        trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells)); % kg

        % reachable structural capacity (kg)
        cumul_trap = compute_cumul_trapcap();


        % Store output
        res.tot_cap       = tot_cap;       % kg/m2
        res.trapcap_tot   = trapcap_tot;   % kg
        res.cumul_trap    = cumul_trap;    % kg
        % breakdown totals
        [~, strap, btrap_res, btrap_dis] = recompute_trapcap(); % redundant?
        tot_trap_capa = strap + btrap_res + btrap_dis;
        sum_tot = sum(tot_trap_capa);
        res.breakdown.total_trapping_capacity      = sum_tot / giga / 1e3;        % Gt
        res.breakdown.structural_trapping_capacity = sum(strap) / giga / 1e3;     % Gt
        res.breakdown.residual_trapping_capacity   = sum(btrap_res) / giga / 1e3; % Gt
        res.breakdown.dissolved_trapping_capacity  = sum(btrap_dis) / giga / 1e3; % Gt


        % Report total trapping values
        str = capacity_report();
        disp(str)


        % Plots (Optional)
        if opt.plotsOn

            subplot(nr,nc,1)
            title('caprock depth (m)')
            plotCellData(Gt.parent, Gt.cells.z, 'edgealpha', 0.2);
            colorbar; rotate3d on;

            subplot(nr,nc,2)
            title('formation thickness (m)')
            plotCellData(Gt.parent, Gt.cells.H, 'edgealpha', 0.2);
            colorbar; rotate3d on;            

            subplot(nr,nc,3)
            title('caprock temperature (C)')
            plotCellData(Gt.parent, compute_temperature(), 'edgealpha', 0.2);
            colorbar; rotate3d on;  

            subplot(nr,nc,4)
            title('caprock pressure (MPa)')
            plotCellData(Gt.parent, compute_pressure() / 1e6, 'edgealpha', 0.2);
            colorbar; rotate3d on;          

            subplot(nr,nc,5)
            title('caprock topography')
            colorbar('off'); rotate3d off; view(0, 90);
            mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);

            subplot(nr,nc,6)
            title('caprock co2 density (kg/m3)')
            plotCellData(Gt.parent, var.co2.rho(p, t + 273.15), 'edgealpha', 0.2);
            colorbar; rotate3d on;

            subplot(nr,nc,7)
            title('total capacity (tonnes/m2)')
            plotCellData(Gt.parent, tot_cap./1e3, 'edgealpha', 0.2);
            colorbar; rotate3d on;

            subplot(nr,nc,8)
            title('structural trap capacity (Mt)')
            plotGrid(Gt.parent, 'facecolor','none', 'edgealpha', 0.1);
            plotCellData(Gt.parent, trapcap_tot/1e3/1e6, 'edgecolor','none'); 
            colorbar; rotate3d on;

            subplot(nr,nc,9)
            title('reachable structural capacity (Mt)')
            plotCellData(Gt.parent, cumul_trap./1e3./1e6, 'edgealpha', 0.2);
            colorbar; rotate3d on;

        end

        %plotGrid(Gt);

    end

    % -----------------------  End of Helper Functions  -------------------

end