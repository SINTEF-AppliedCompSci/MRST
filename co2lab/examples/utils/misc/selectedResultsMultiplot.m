function selectedResultsMultiplot(Gt, reports, plot_steps, varargin)

    % background variables can be: h|h_max|pressure|rs|sGmax|s|smax|totalCO2

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
    
    % Massaging input arguments
    opt.plot_distrib = true;
    opt.plot_wells = true;
    opt.plot_well_numbering = false;
    opt.plot_plume = true;
    opt.plot_traps = false;
    opt.plot_rivers = false;
    opt.plume_threshold = 0.3;
    opt.background = '';
    opt.background_threshold = [];
    opt.backgroundalpha = 0.7;
    opt.background_use_log = false;
    opt.quickclip = true;
    opt.ref_temp = 7 + 273.15; % used for 'totalCO2' (value from SPE paper)
    opt.ref_depth = 80;        % used for 'totalCO2' (value from SPE paper)
    opt.temp_grad = 35.6;      % used for 'totalCO2' (value from SPE paper)
    opt.poro = 0.36;           % used for 'totalCO2' (value from SPE paper)
    opt.ntg = 1;
    opt.ymax = inf;
    opt.maplines = 40;
    opt.init_state = [];
    opt.trapmethod = false;     % cell vs node based
    opt.ta = [];                % computed if required and empty
    opt.rivercolor = [1 0 0];
    opt.trapcolor = [1 0 0];
    
    rho = getValuesSPE134891();
    opt.rhoRef = rho(2);
    
    opt = merge_options(opt, varargin{:});

    % cut away irrelevant part of grid
    kept_cells = Gt.cells.centroids(:,2) < opt.ymax;
    [GtSub, ~, ~, nix] = extractSubgrid(Gt, kept_cells);
    GtSub.cells.z = Gt.cells.z(kept_cells);
    GtSub.nodes.z = Gt.nodes.z(nix);
    
    % plotting
    h = figure;
    num_tsteps = numel(plot_steps);

    % mymap = colormap;
    % mymap(1,:) = [1 1 1];
    
    for i = 1:numel(plot_steps);
       ix = plot_steps(i);
        subplot(1, num_tsteps, i);
        
        plume = []; 
        if opt.plot_plume % && isfield(states{i}, 'h')
           plume = reports(ix).sol.h(kept_cells);
           %plume = stepdata(i).report.sol.h(kept_cells);
        end
        field = [];
        if ~isempty(opt.background)
            if ischar(opt.background)
                % 'background' gives us a variable name, or 'totalCO2'
                if strcmpi(opt.background, 'totalCO2')
                    field = compute_total_co2(Gt, reports(ix).sol , ...
                                              opt.ref_temp,         ...
                                              opt.ref_depth,        ...
                                              opt.temp_grad,        ...
                                              opt.poro,             ...
                                              opt.ntg,              ...
                                              opt.rhoRef);
                elseif strcmpi(opt.background, 'overpressure')
                   if isempty(opt.init_state)
                      error(['Plotting of overpressure requires that init_state ' ...
                             'is provided.']);
                   end
                   field = (reports(ix).sol.pressure - opt.init_state.pressure) / 1e5;
                    fprintf('Max overpressure at step %i: %f bar\n', ix, max(field(:)));
                else
                    field = reports(ix).sol.(opt.background);
                end
                if size(field, 2) > 1
                    % this happens for saturations, in which case we only want to
                    % keep the second column (CO2 saturation)
                    field = field(:,2);
                end
            elseif iscell(opt.background)
                % 'background' provides us directly with data, one cell per timestep
                field = opt.background{ix};
            else
                % 'background' provides us directly with data, same for all timesteps
                field = opt.background;
            end
            field = field(kept_cells);
        end
        
        if opt.background_use_log
           field = log10(field);
        end
                
        wellcells = [];
        if opt.plot_wells
            wellcells = [reports(ix).W.cells];
            % keeping only wells within area that is not cut away
            tmp = zeros(Gt.cells.num,1);
            tmp(wellcells) = 1;
            tmp = find(tmp(kept_cells)); % indices of kept well cells
            order = arrayfun(@(x) find(wellcells == x, 1, 'first'), tmp);
            [~, b] = sort(order);
            wellcells = tmp(b);
        end
    
        traps = [];
        rivers = [];
        if opt.plot_traps
            if isempty(opt.ta)
                opt.ta = trapAnalysis(GtSub, opt.trapmethod);
            end
            traps = opt.ta.traps;
            if opt.plot_rivers
                rivers = opt.ta.cell_lines;
            end
        end
        
        mapPlot(h, GtSub, ...
                'traps', traps, ...
                'trapcolor', opt.trapcolor, ...
                'rivers', rivers, ...
                'rivercolor', opt.rivercolor, ...
                'title', sprintf('Year: %i', ceil(reports(ix).t/year)), ...
                'plumes', plume, ...
                'wellcells', wellcells, ...
                'well_numbering', opt.plot_well_numbering, ...
                'plume_h_threshold', opt.plume_threshold, ...
                'background', field, ...
                'backgroundalpha', opt.backgroundalpha,...
                'background_threshold', opt.background_threshold, ...
                'quickclip', opt.quickclip,...
                'maplines', opt.maplines, ...
                'colorbarposition', 'southoutside');
        
        % Tweaking result
        axis off;
        set(gcf, 'color', [1 1 1]);
        %colormap(h, mymap);
        %set(gcf, 'position', [1, 1, 1920, 850]);
        %set(gcf, 'position', [1, 1, 1600, 960]);
        %set(get(gca, 'title'), 'fontsize', 30);
        set(gcf, 'position', [30 30 660 930]);
        %colorbar off;
    end
    % adjusting height of map plots to be the same (not necessarily the case
    % initially, due to a bug(??) in Matlab)
    ypos = [inf, -inf];
    for i = 1:num_tsteps
        subplot(1, num_tsteps, i);
        gcapos = get(gca, 'position');
        ypos(1) = min(ypos(1), gcapos(2));
        ypos(2) = max(ypos(2), gcapos(4));
    end
    for i = 1:num_tsteps
        subplot(1, num_tsteps, i);
        gcapos = get(gca, 'position');
        gcapos(2) = ypos(1);
        gcapos(4) = ypos(2);
        set(gca, 'position', gcapos);
    end
    
    if opt.plot_distrib
        h = figure; plot(1);
        ax = get(h, 'currentaxes');
        
        % load all timesteps up to last plotted one (a bit of a hack)
        plotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
        %fsize = 24;
        %set(get(gca, 'xlabel'), 'fontsize', fsize)
        %set(get(gca, 'ylabel'), 'fontsize', fsize)
        %set(gca,'fontsize', fsize);
        set(gcf, 'position', [1, 1, 850, 850]);
    end
end
% ----------------------------------------------------------------------------

function field = compute_total_co2(Gt, state, reftemp, refdepth, temp_grad, ...
                                   poro, ntg, rhoRef)
    
    co2 = CO2props;
    P   = state.pressure;
    T   = reftemp + (Gt.cells.z - refdepth) * temp_grad / 1000;
    % 'field' represents tons of mass of CO2 per lateral square meter
    field = state.s(:,2) .* co2.rho(P, T) .* poro .* ntg .* Gt.cells.H / 1e3; 
                                               
    % Adding contribution from dissolution
    if isfield(state, 'rs')
        diss = state.s(:,1) .* rhoRef .* state.rs .* poro .* ntg .* Gt.cells.H / 1e3;
        field = field + diss;
    end
end
