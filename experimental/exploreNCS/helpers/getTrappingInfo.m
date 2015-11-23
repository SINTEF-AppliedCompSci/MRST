function [ varargout ] = getTrappingInfo(name, coarsening, varargin)
% Computes structural trapping info and makes plots (optional).
%
% SYNOPSIS
%   res = getTrappingInfo('Stofm',3)
%   res = getTrappingInfo('Stofm',3, 'plotsOn',false)
%   res = getTrappingInfo('Stofm',3, 'cells',trapcells)
%
% DESCRIPTION
%   Parameter values (seafloor_depth, etc.) are obtained using
%   getSeaInfo(), which returns values that correspond to the specified
%   formation. Then upper bound (theoretical) trapping capacities are
%   computed and returned. Plotting is optional.
%
%   Formations may contain net-to-gross rock property. If so, the pore
%   volume of the formation is updated using:
%       poro = poro.*ntg,
%   in order to reflect any potentially reduced pore volumes. NB:
%       Bulk rock volume = (Gt.cells.volume.*Gt.cells.H)
%       Bulk rock volume * NTG = Net rock volume
%       Net rock volume * poro = volume which fluid can occupy
%
%   If varargin cells is specified, trapping info is returned for those
%   specific cells. If not specified, default is cells = 1:Gt.cells.num
%
% See also
%   exploreCapacity.m
    

    opt.plotsOn     = true;
    opt.cells       = [];
    opt.trapName    = {'default'};
    opt = merge_options(opt, varargin{:});
    
    
    %% Arbitrary reference co2 density
    rhoCref = 760 * kilogram / meter ^3;
    
    
    %% Get formation grid and info
    % entire formation grid
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
    res.co2     = var.co2; % NB: co2.rho(p,t) where t is in Kelvin
    res.caprock_pressure    = compute_pressure();               % Pascals
    res.caprock_temperature = compute_temperature() + 273.15;   % Kelvin
    
    
    %% trapping info for specific cells
    % 'cells' is intended to belong to the same structural trap (whole or
    % portion of trap)
    if ~isempty(opt.cells)
        assert( unique(var.ta.traps(opt.cells)) ~= 0 , ...
            'Cells are not part of a structural trap.')
        assert( numel( unique(var.ta.traps(opt.cells)) ) == 1 , ...
            'Cells contain more than one unique trap ID.')
        
        [output] = compute_trapcap(opt.cells, opt.trapName);
    end

    %% function returns:
    varargout{1} = res;
    if exist('output','var')
       varargout{2} = output;
    end

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
        t = compute_temperature() + 273.15; % Kelvin

        Gt   = var.Gt;
        ta   = var.ta;
        poro = var.rock2D.poro;
        if isfield(var.rock2D,'ntg')
            ntg = var.rock2D.ntg;
            fprintf(['\n Note: net-to-gross rock property exists.',...
               ' Net rock volume is (on avg) %d percent of bulk rock volume. \n'], mean(ntg)*100);
        else
            ntg = ones(numel(poro),1);
        end

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
        strap_pvol_tot       = Gt.cells.volumes .* H1 .* poro .* ntg;
        strap_pvol_co2_plume = strap_pvol_tot .* (1 - info.res_sat_wat);
        strap_pvol_co2_diss  = strap_pvol_tot .* info.res_sat_wat .* info.dis_max;

        strap = strap_pvol_co2_plume .* var.co2.rho(p,t) + ...
              strap_pvol_co2_diss  .* rhoCref;
        %fprintf('\nAverage strap is %d kg.\n', mean(strap))

        % Computing total trapping volume below structural traps (dissolved
        % and residually trapped
        btrap_pvol_tot          = Gt.cells.volumes .* H2 .* poro .* ntg;
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
          figure; set(gcf,'Position',[1 1 1102 1261])
          [nr,nc] = deal(3);
        end
        Gt = var.Gt;
        ta = var.ta;


        % Calculations:

        % reservoir conditions
        t = compute_temperature() + 273.15; % Kelvin
        p = compute_pressure();             % Pascals

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
        fprintf('\n\n -------- %s ------- \n\n', name)
        disp(str)


        % Plots (Optional)
        if opt.plotsOn
            
            myplotCellData = @(g,d) plotCellData(g,d,'EdgeColor','none'); % 'EdgeAlpha',0.2
            
            subplot(nr,nc,1)
            title({'caprock';'depth (m)'})
            plotCellData(Gt.parent, Gt.cells.z, 'EdgeColor','none');
            colorbar; rotate3d on;

            subplot(nr,nc,2)
            title({'formation';'thickness (m)'})
            plotCellData(Gt.parent, Gt.cells.H, 'EdgeColor','none');
            colorbar; rotate3d on;            

            subplot(nr,nc,3)
            title({'caprock';'temperature (C)'})
            plotCellData(Gt.parent, compute_temperature(), 'EdgeColor','none');
            colorbar; rotate3d on;  

            subplot(nr,nc,4)
            title({'caprock';'pressure (MPa)'})
            plotCellData(Gt.parent, compute_pressure() / 1e6, 'EdgeColor','none');
            colorbar; rotate3d on;          

            subplot(nr,nc,5)
            title({'caprock';'topography'})
            colorbar('off'); rotate3d off; view(0, 90);
            mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);

            subplot(nr,nc,6)
            title({'caprock';'CO_2 density (kg/m^3)'})
            plotCellData(Gt.parent, var.co2.rho(p, t), 'EdgeColor','none');
            colorbar; rotate3d on;

            subplot(nr,nc,7)
            title({'total capacity';'(tonnes/m^2)'})
            plotCellData(Gt.parent, tot_cap./1e3, 'EdgeColor','none');
            colorbar; rotate3d on;

            subplot(nr,nc,8)
            title({'structural trap';'capacity (Mt)'})
            plotGrid(Gt.parent, 'facecolor','none', 'edgealpha', 0.1);
            plotCellData(Gt.parent, trapcap_tot/1e3/1e6, 'EdgeColor','none'); 
            colorbar; rotate3d on;

            subplot(nr,nc,9)
            title({'reachable structural';'capacity (Mt)'})
            plotCellData(Gt.parent, cumul_trap./1e3./1e6, 'EdgeColor','none');
            colorbar; rotate3d on;
            
            hfig = gcf;
            set(findobj(hfig.Children,'Type','axes'),'Fontsize',16,'box','on')
            axis(findobj(hfig.Children,'Type','axes'),'equal','tight','off')

        end

        %plotGrid(Gt);

    end

    % ---------------------------------------------------------------------

    function [output, H1, strap, btrap_res, btrap_dis] = compute_trapcap(cells, trapName)
    % computes trapping capacity for specific grid cells
    % 'cells' is intended to belong to the same structural trap

        %assert()

        p = compute_pressure();
        t = compute_temperature() + 273.15; % Kelvin

        Gt   = var.Gt;
        ta   = var.ta;
        poro = var.rock2D.poro;
        perm = var.rock2D.perm;
        if isfield(var.rock2D,'ntg')
            ntg = var.rock2D.ntg;
            fprintf(['\n Note: net-to-gross rock property exists.',...
               ' Net rock volume is (on avg) %d percent of bulk rock volume. \n'], mean(ntg)*100);
        else
            ntg = ones(numel(poro),1);
        end

        % computing structural trap heights (H1) for each cell
        trapcells     = cells;
        H1            = zeros(Gt.cells.num, 1);
        if ~isempty(trapcells)
         H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
        end
        H1=min(H1,Gt.cells.H);
        H2 = Gt.cells.H - H1;
        assert(all(H1<=Gt.cells.H));

        % Computing total trapping volume in structural traps (dissolved
        % and structurally trapped
        strap_pvol_tot       = Gt.cells.volumes .* H1 .* poro .* ntg;
        strap_pvol_co2_plume = strap_pvol_tot .* (1 - info.res_sat_wat);
        strap_pvol_co2_diss  = strap_pvol_tot .* info.res_sat_wat .* info.dis_max;

        strap = strap_pvol_co2_plume .* var.co2.rho(p,t) + ...
              strap_pvol_co2_diss  .* rhoCref;
        %fprintf('\nAverage strap is %d kg.\n', mean(strap))

        % Computing total trapping volume below structural traps (dissolved
        % and residually trapped
        btrap_pvol_tot          = Gt.cells.volumes .* H2 .* poro .* ntg;
        btrap_pvol_co2_residual = btrap_pvol_tot .* info.res_sat_co2;
        btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-info.res_sat_co2) .* info.dis_max;

        btrap_res = btrap_pvol_co2_residual .* var.co2.rho(p,t);
        btrap_dis = btrap_pvol_co2_dissol   .* rhoCref;


        % Other info about The Structural Trap that might be useful for
        % comparions later:
        output.numtrapcells     = numel(cells);
        output.rockVol_m3       = sum(Gt.cells.volumes .* H1);
        output.netVol_m3        = sum(Gt.cells.volumes .* H1 .* ntg);
        output.poreVol_m3       = sum(Gt.cells.volumes .* H1 .* ntg .* poro);
                                % cubic meters, total in structural trap

        % av ntg, poro, perm, and their ranges:
        output.ntg_min = min(ntg(trapcells));
        output.ntg_max = max(ntg(trapcells));
        output.ntg_av  = mean(ntg(trapcells));

        output.perm_min = min(perm(trapcells));
        output.perm_max = max(perm(trapcells));
        output.perm_av  = mean(perm(trapcells));

        output.poro_min = min(poro(trapcells));
        output.poro_max = max(poro(trapcells));
        output.poro_av  = mean(poro(trapcells));

        % min/max depth of structural trap (not formation), and average:
        output.depth_min = min(Gt.cells.z(trapcells));
        output.depth_max = max(Gt.cells.z(trapcells));
        output.depth_av  = mean(Gt.cells.z(trapcells));


        % storage capacity in structure, in volume (m3), and in mass (kg,
        % Mt) assuming 100% of pore volume is capable of storing CO2 (no
        % storage efficiency factor used here):
        output.storageCap_m3 = output.poreVol_m3;
        output.storageCap_kg = sum( Gt.cells.volumes(trapcells) ...
            .* H1(trapcells) .* ntg(trapcells) .* poro(trapcells) ...
            .* var.co2.rho(p(trapcells),t(trapcells)) ); % m3 * kg/m3 = kg
        output.storageCap_Mt = output.storageCap_kg / 1e9; % 10^9 kg = 1 Mt


        % corresponding storage efficiency in structural trap:
        output.Seff = output.storageCap_m3 / output.poreVol_m3;


        output.T = table(  [output.numtrapcells;        ...
                    output.rockVol_m3;          ...
                    output.netVol_m3;           ...
                    output.poreVol_m3;          ...
                    output.ntg_min;             ...
                    output.ntg_max;             ...
                    output.ntg_av;              ...
                    output.perm_min;            ...
                    output.perm_max;            ...
                    output.perm_av;             ...
                    output.poro_min;            ...
                    output.poro_max;            ...
                    output.poro_av;            ...
                    output.depth_min;           ...
                    output.depth_max;           ...
                    output.depth_av;            ...
                    output.storageCap_Mt],       ...
            'RowNames', {'NumTrapCells' ...
            'RockVolm3' 'NetVolm3' 'PoreVolm3' ...
            'ntgmin' 'ntgmax' 'ntgav' ...
            'permmin' 'permmax' 'permav' ...
            'poromin' 'poromax' 'poroav' ...
            'depthmin' 'depthmax' 'depthav' ...
            'StorageCapacityMt'}, ...
            'VariableNames', trapName  );

        output.trapName = trapName;

    end

    % -----------------------  End of Helper Functions  -------------------

end