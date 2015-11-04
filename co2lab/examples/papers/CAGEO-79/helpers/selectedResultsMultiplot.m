function selectedResultsMultiplot(basedir, casename, tsteps, varargin)

    % background variables can be: h|h_max|pressure|rs|sGmax|s|smax|totalCO2
    
    % Massaging input arguments
    opt.plot_distrib = true;
    opt.plot_wells = true;
    opt.plot_well_numbering = false;
    opt.plot_plume = true;
    opt.plot_traps = false;
    opt.plume_threshold = 0.3;
    opt.background = '';
    opt.backgroundalpha = 0.7;
    opt.background_use_log = false;
    opt.quickclip = true;
    opt.ref_temp = 7 + 273.15; % used for 'totalCO2' (value from SPE paper)
    opt.ref_depth = 80;        % used for 'totalCO2' (value from SPE paper)
    opt.temp_grad = 35.6;      % used for 'totalCO2' (value from SPE paper)
    opt.poro = 0.36;           % used for 'totalCO2' (value from SPE paper)
    opt.ymax = inf;

    rho = getValuesSPE134891();
    opt.rhoRef = rho(2);
    
    opt = merge_options(opt, varargin{:});

    % Loading data
    [Gt, stepdata] = load_data(fullfile(basedir, casename, 'report'), tsteps);
    
    % cut away irrelevant part of grid
    kept_cells = Gt.cells.centroids(:,2) < opt.ymax;
    [GtSub, ~, ~, nix] = extractSubgrid(Gt, kept_cells);
    GtSub.cells.z = Gt.cells.z(kept_cells);
    GtSub.nodes.z = Gt.nodes.z(nix);
    
    % plotting
    h = figure;
    num_tsteps = numel(tsteps);
    
    for i = 1:num_tsteps
        subplot(1, num_tsteps, i);
        
        plume = []; 
        if opt.plot_plume && isfield(stepdata(i).report.sol, 'h')
            plume = stepdata(i).report.sol.h(kept_cells);
        end
        field = [];
        if ~isempty(opt.background)
            if ischar(opt.background)
                % 'background' gives us a variable name, or 'totalCO2'
                if strcmpi(opt.background, 'totalCO2')
                    field = compute_total_co2(Gt, stepdata(i).report, ...
                                              opt.ref_temp,           ...
                                              opt.ref_depth,          ...
                                              opt.temp_grad,          ...
                                              opt.poro,               ...
                                              opt.rhoRef);
                elseif strcmpi(opt.background, 'overpressure')
                    [~, initstate] = load_data(fullfile(basedir, casename, 'report'), 0);
                    field = (stepdata(i).report.sol.pressure - ...
                            initstate.report.sol.pressure) / 1e5;
                    fprintf('Max overpressure at step %i: %f bar\n', i, max(field(:)));
                else
                    field = stepdata(i).report.sol.(opt.background);
                end
                if size(field, 2) > 1
                    % this happens for saturations, in which case we only want to
                    % keep the second column (CO2 saturation)
                    field = field(:,2);
                    %field(field <0.0005) = NaN;
                end
            elseif iscell(opt.background)
                % 'background' provides us directly with data, one cell per timestep
                field = opt.background{i};
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
            wellcells = [stepdata(i).report.W.cells];
            % keeping only wells within area that is not cut away
            tmp = zeros(Gt.cells.num,1);
            tmp(wellcells) = 1;
            tmp = find(tmp(kept_cells)); % indices of kept well cells
            order = arrayfun(@(x) find(wellcells == x, 1, 'first'), tmp);
            [a, b] = sort(order);
            wellcells = tmp(b);
            %wellcells = find(tmp);
        end
    
        traps = [];
        if opt.plot_traps
            ta = trapAnalysis(GtSub, false);
            traps = ta.traps;
        end
        
        mapPlot(h, GtSub, ...
                'traps', traps, ...
                'title', sprintf('Year: %i', ceil(stepdata(i).report.t/year)), ...
                'plumes', plume, ...
                'wellcells', wellcells, ...
                'well_numbering', opt.plot_well_numbering, ...
                'plume_h_threshold', opt.plume_threshold, ...
                'background', field, ...
                'backgroundalpha', opt.backgroundalpha,...
                'quickclip', opt.quickclip,...
                'colorbarposition', 'southoutside');
        
        % Tweaking result
        axis off;
        set(gcf, 'color', [1 1 1]);
        %set(gcf, 'position', [1, 1, 1920, 850]);
        set(gcf, 'position', [1, 1, 1600, 960]);
        set(get(gca, 'title'), 'fontsize', 30);
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
        [~, allsteps] = load_data(fullfile(basedir, casename, 'report'), 1:tsteps(end));
        plotTrappingDistribution(ax, allsteps, 'legend_location', 'northwest');
        fsize = 24;
        set(get(gca, 'xlabel'), 'fontsize', fsize)
        set(get(gca, 'ylabel'), 'fontsize', fsize)
        set(gca,'fontsize', fsize);
        set(gcf, 'position', [1, 1, 850, 850]);
    end
end
% ----------------------------------------------------------------------------

function [Gt, stepdata] = load_data(dir, tsteps)
    
    % loading grid
    tmp = load(fullfile(dir, 'simulation_grid'), 'Gt');
    Gt = tmp.Gt;
    
    % loading timesteps
    for i = 1:numel(tsteps)
        stepdata(i) = load(fullfile(dir, sprintf('report_%i', tsteps(i))));
    end
end
% ----------------------------------------------------------------------------

function field = compute_total_co2(Gt, data, reftemp, refdepth, temp_grad, ...
                                   poro, rhoRef)
    
    co2 = CO2props;
    P   = data.sol.pressure;
    T   = reftemp + (Gt.cells.z - refdepth) * temp_grad / 1000;
    % 'field' represents tons of mass of CO2 per lateral square meter
    field = data.sol.s(:,2) .* co2.rho(P, T) .* poro .* Gt.cells.H / 1e3; 
                                               
    % Adding contribution from dissolution
    if isfield(data.sol, 'rs')
        diss = data.sol.s(:,1) .* rhoRef .* data.sol.rs .* poro .* Gt.cells.H / 1e3;
        field = field + diss;
    end
end
