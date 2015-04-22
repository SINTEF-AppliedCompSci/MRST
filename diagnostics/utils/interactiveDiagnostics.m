function interactiveDiagnostics(G, rock, W, varargin)
%Launch an interactive diagnostics session
%
% SYNOPSIS:
%   interactiveDiagnostics(G, rock, W);
%   interactiveDiagnostics(G, rock, W, 'state', state)
%   interactiveDiagnostics(G, rock, W, dataset, 'state', state)
%
% DESCRIPTION:
%   This function launches an interactive session for doing flow
%   diagnostics. The functionality differs slightly based on the input
%   arguments given:
%     - If a dataset is given, this set of cellwise data will be available
%     for visualization.
%     - If a state is given, this state will allow for visualization of
%     the component ratios inside a drainage volume. The flux from this
%     state can also be used to calculate time of flight if "computeFlux"
%     is disabled, for instance if the user has some external means of
%     computing fluxes.
%
%   Once the initialization is complete, two windows will be produced:
%     - A plotting window, showing the reservoir along with the wells and
%     visualized quantitites.
%         In the plotting window, it is possible to click wells to get
%         additional information, such as the allocation factors per
%         perforation in the well, the phase distribution grouped by time
%         of flight and the corresonding pore volumes swept/drained.
%
%     - A controller window which is used to alter the state of the
%     plotting window:
%     Wells can be selected for display (if a well is selected, the
%     corresponding drainage (producer) or flooding (injector) volumes will
%     be visualized. A simple playback function can be used to show
%     propagation of time of flight. Different quantitites can be
%     visualized to get a better understanding of the system.
%
%     A set of buttons provide (experimental) well editing, access to
%     visualization of the well pair connections, Phi/F diagram with
%     Lorenz' coefficient etc.
%
%
% REQUIRED PARAMETERS:
%
%   G    - Valid grid structure.
%
%   rock - Rock with valid permeability and porosity fields.
%
%   W    - A set of wells which are compatible with the incompTPFA solver.
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%  'state' - Reservoir state containing fluid saturations and optionally
%            flux and pressure (if 'computeFlux' is false)
%
%  'tracerfluid' - Fluid used for tracer computation. This defaults to a
%            trivial fluid.
%
%  'fluid'  - Fluid used for mobility calculations if 'useMobilityArrival'
%             is enabled.
%
%  'LinSolve' - Linear solver for pressure systems. Defaults to mldivide
%  (backslash)
%
%  'computeFlux' - If set to false, fluxes are extracted from the provided
%                  state keyword argument. This requires a state to be
%                  provided and can be used if the fluxes are computed
%                  externally (for instance from a expensive full-physics
%                  simulation)
%
%  'useMobilityArrival' - If the well plot showing nearby saturations
%                         should plot mobility instead of saturations. This
%                         may be interesting in some cases because the
%                         mobile fluids are more likely to be extracted.
%                         However, this plot is often dominated by very
%                         mobile gas regions.
%
% RETURNS:
%
%   Nothing. Creates two figures.
%
% EXAMPLE:
%     G = computeGeometry(cartGrid([10, 10, 2]));
%     rock = struct('poro', ones(G.cells.num, 1), 'perm', ones(G.cells.num, 1)*darcy)
%
%     W = verticalWell([], G, rock, 1,  1, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  1,  10, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  10, 1, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  10, 10, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  5,  5, [], 'Val', 0)
%     interactiveDiagnostics(G, rock, W);
%
% SEE ALSO:
%   Diagnostics examples

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
    try
      require mrst-gui
    catch %#ok<CTCH>
       mrstModule add mrst-gui
    end

    if nargin > 3 && mod(numel(varargin), 2) == 1
        dsname = inputname(4);
        datasets = varargin{1};
        varargin = varargin(2 : end);
    else
        dsname = '';
        datasets = [];
    end
    water    = initSingleFluid('mu' , 1*centi*poise, 'rho', 1000*kilogram/meter^3);
    oilwater = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                               'rho', [1000, 859]*kilogram/meter^3, ...
                               'n'  , [   2,   2]);

    oilwater.names = {'Water', 'Oil', 'Gas'};
    
    window_title = 'Interactive Diagnostics';

    opt = struct('state',               [],...
                 'tracerfluid',         water, ...
                 'LinSolve',            @mldivide, ...
                 'computeFlux',         true, ...
                 'useMobilityArrival',  false,...
                 'fluid',               oilwater, ...
                 'window_title',        window_title ...
    );

    opt = merge_options(opt, varargin{:});

    assert(opt.computeFlux || ~isempty(opt.state),...
        'If computeFlux is off a state must be provided!')

    if isempty(datasets) && ~isempty(opt.state)
        dsname = 'State';
        datasets = opt.state;
    end

    % Currently playing back
    playback = false;

    % Main scope variables for tracers and tof
    [D, WP] = deal([]);


    fig_main = figure('name', opt.window_title);
    
    % Check for which version of handle graphics MATLAB is currently
    % running
    isHG1 = isnumeric(fig_main);
    
    axis tight off

    % Main scope variables for control panel
    [fig_ctrl,...
     wellPlot,...
     draintoggle,...
     floodtoggle,...
     speedsh,...
     mtofsh,...
     mtofeh,...
     Mtofsh,...
     Mtofeh,...
     alfash,...
     tofext,...
     hdataset,...
     hset_op,...
     ni, np] = deal(NaN);
    
    % Setup "persistent" variables
    pm_htop = [];
    pm_htext = [];
    pm_hs = [];
    pm_hline = [];
    pm_mainph = [];
    pm_outlineph = [];
    phiPlot = []; % Used in plotPhi
    pwc_wch = []; % Used in plotWellConnections

    % Precompute TOF etc.
    pv = poreVolume(G, rock);
    computeValues();

    % Create control panel
    createMainControl();

    % Plot initial setup
    plotMain()

    % Private helpers
    function createMainControl()
        % Set up figure handles
        if ~ishandle(fig_ctrl)
            pos = get(fig_main, 'OuterPosition');
            size_xy = [475 550];
            pos_xy = pos(1:2) + pos(3:4) - size_xy - [400 0];
            fig_ctrl = figure('Position',[pos_xy,  size_xy], 'Toolbar','none', 'MenuBar', 'none');
            set(fig_ctrl, 'Name', ['Controller ', opt.window_title]);
        else
            clf(fig_ctrl)
        end
        % Drainage / flooding controls
        gpw = .5;
        gph = .55;

        drainp = uibuttongroup('Parent', fig_ctrl, 'Title', 'Drainage volumes', 'Position', [0 .45 gpw gph]);
        floodp = uibuttongroup('Parent', fig_ctrl, 'Title', 'Flooding volumes', 'Position', [.5 .45 gpw gph]);

        pp = unique(D.ppart);
        ip = unique(D.ipart);

        np = numel(pp);
        ni = numel(ip);

        draintoggle = zeros(np, 1);
        floodtoggle = zeros(ni, 1);

        draintoggle = uicontrol(drainp, 'Style', 'listbox', ...
                                        'Units', 'normalized',...
                                        'Max', 2, ...
                                        'Min', 0, ...
                                        'String', {W(D.prod).name},...
                                        'Value', [], ...
                                        'Callback', @plotMain, ...
                                        'Position', [0 0 1 1]);

        floodtoggle = uicontrol(floodp, 'Style', 'listbox', ...
                                        'Units', 'normalized',...
                                        'Max', 2, ...
                                        'Min', 0, ...
                                        'String', {W(D.inj).name},...
                                        'Value', 1:numel(D.inj), ...
                                        'Callback', @plotMain, ...
                                        'Position', [0 0 1 1]);
        % TOF adjustment
        tofp = uipanel('Parent', fig_ctrl, 'Title', 'Time of flight', 'Position', [0 .2 1 .25]);
        tof = D.tof(:);
        
        isNeg = tof <= 0;
        tof(isNeg) = min(tof(~isNeg));
        
        tofext = convertTo(([min(tof), 5*10^(mean(log10(tof)))]), year);
        tofext(2) = max(tofext(2), 15);

        tof_N = 5;
        tof_h = 1/tof_N;
        [speedsh, speedeh] = linkedSlider(tofp, [0 1*tof_h 1 tof_h], .15, [1 1000], 50, 'Resolution');
        [mtofsh, mtofeh]   = linkedSlider(tofp, [0 2*tof_h 1 tof_h], .15, tofext, tofext(1), 'Min TOF');
        [Mtofsh, Mtofeh]   = linkedSlider(tofp, [0 3*tof_h 1 tof_h], .15, tofext, tofext(2), 'Max TOF');
        [alfash, alfaeh]   = linkedSlider(tofp, [0 4*tof_h 1 tof_h], .15, [0 1], 1, 'Alpha');

        uicontrol(tofp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0, .2 1/tof_N],...
                   'String', 'Play',...
                   'Callback', @(src, event) playBackTof(src, event)...
                   );
        uicontrol(tofp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [.2 0, .2 1/tof_N],...
                   'String', 'Stop',...
                   'Callback', @(src, event) stopPlayBackTof(src, event)...
                   );
        % General config
        confp = uipanel('Parent', fig_ctrl, 'Title', 'Configuration', 'Position', [0 0 1 .2]);

        cellfields = getStructFields(G, datasets, dsname);
        perm = strcat({'X', 'Y', 'Z'}, ' permeability');
        
        hdataset = uicontrol(confp, 'Style', 'popupmenu',...
                                    'Units', 'normalized',...
                                    'Position', [.0 0 .475 1],...
                                    'Callback', @plotMain,...
                                    'String', {'Forward TOF',...
                                               'Backward TOF',...
                                               'Sum of TOFs',...
                                               'Drainage region',...
                                               'Flooding region',...
                                               'Porosity', ...
                                               perm{1:size(rock.perm, 2)},...
                                               cellfields{:}}...
                                   );
        hset_op = uicontrol(confp, 'Style', 'popupmenu',...
                                   'Units', 'normalized',...
                                   'Position', [.525 0 .475 1],...
                                   'Callback', @plotMain,...
                                   'String', {'Union (Flood, Drain)',...
                                              'Intersection (Flood, Drain)',...
                                              'Flood',...
                                              'Drain'}...
                                   );
        bheight = .05;
        bwidth = .2;

        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 bheight bwidth .5],...
                   'String', 'Apply',...
                   'Callback', @(src, event) plotMain(src, event)...
                   );


        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 + bwidth*1 bheight bwidth .5],...
                   'String', 'Phi/F diagram',...
                   'Callback', @(src, event) plotPhi(src, event)...
                   );

        if opt.computeFlux
            uicontrol(confp, 'Style', 'pushbutton',...
                       'Units', 'normalized',...
                       'Position', [0 + bwidth*2 bheight bwidth .5],...
                       'String', 'Edit wells',...
                       'Callback', @(src, event) changeWells()...
                       );
        end

        uicontrol(confp, 'Style', 'togglebutton',...
                   'Units', 'normalized',...
                   'Position', [0 + bwidth*3 bheight bwidth .5],...
                   'String', 'Wellpairs',...
                   'Callback', @plotWellConnections...
                   );

    end

    function plotMain(src, event)

        [cdata clim cmap] = selectDataset();

        if ishandle(fig_main)
            set(0, 'CurrentFigure', fig_main);
        else
            figure(fig_main);
        end

        if any(ishandle(pm_mainph))
            delete(pm_mainph);
        end

        if ~any(ishandle(pm_outlineph))
            if isHG1
                pm_outlineph = plotGrid(G, 'facec', 'none', 'edgea', .05, 'edgec', 'black');
            else
                pm_outlineph = plotGrid(G, 'facec', 'none', 'edgec', 'black');
            end
            set(pm_outlineph, 'UserData', 'gridoutline');
        end

        % Limit dataset based on tof
        min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
        max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);


        alpha   = get(alfash, 'Value');

        % Drainage volumes
        tmp = get(draintoggle, 'Value');
        psubs = ismember(D.ppart, tmp);
        data = D.tof(:,2);
        psubset = data >= min_tof & data <= max_tof;

        % Flooding volumes
        tmp = get(floodtoggle, 'Value');
        isubs = ismember(D.ipart, tmp);
        data = D.tof(:,1);
        isubset = data >= min_tof & data <= max_tof;

        % Find the selection
        tmp = get(hset_op, 'Value');
        set_op_strs = get(hset_op, 'String');
        set_op_str = set_op_strs{tmp};
        selection = [];
        if (any(regexpi(set_op_str, 'union')))
            selection = (isubset & isubs) | (psubset & psubs);
        elseif (any(regexpi(set_op_str, 'intersection')))
            selection = (isubset & isubs) & (psubset & psubs);
        elseif (any(strcmpi(set_op_str, 'flood')))
            selection = (isubset & isubs);
        elseif (any(strcmpi(set_op_str, 'drain')))
            selection = (psubset & psubs);
        end
        
        % Plot selection
        if any(selection)
            pm_mainph = plotCellData(G, cdata, selection, 'EdgeColor', 'none', 'FaceAlpha', alpha);
        end


        fastRotateButton();

        % Plot wells
        if ~all(ishandle([pm_htop, pm_htext, pm_hs])) || isempty([pm_htop, pm_htext, pm_hs])

            [pm_htop, pm_htext, pm_hs, pm_hline] = plotWell(G, W,  'color', 'red', 'height',  0);
            for i = 1:numel(W)
                color = colorizeWell('global', i, D);
                set([pm_htop(i) pm_htext(i) pm_hs(i)],    'ButtonDownFcn', @(src, event) onClickWell(src, event, i));
                set([pm_htop(i) pm_hs(i)], 'FaceColor', color, 'EdgeColor', color)
                set([pm_htext(i) pm_hline(i)], 'Color', color)
                set(pm_htext(i), 'FontWeight', 'bold', 'Interpreter', 'none')
            end

        end
         axis tight
        if all(isfinite(clim)) && clim(2) > clim(1)
            caxis(clim);
        end
    end

    function playBackTof(src, event)
        if playback
            return
        end
        N = round(get(speedsh, 'Value'));
        t_0 = get(mtofsh, 'Value');
        t_end = get(Mtofsh, 'Value');
        if t_0 >= t_end
            disp 'Minimum TOF is larger than maximum TOF...'
            return
        end
        playback = true;
        % Pop figure to front explicitly before loop
        figure(fig_main);
        for i = 0:N
            if ~playback
                return
            end
            set(Mtofsh, 'Value', min(t_0 + i*(t_end - t_0)/N, tofext(2)));
            f = get(Mtofsh, 'Callback');
            % Fire event
            f(Mtofsh, []);
            timer = tic();
            plotMain();
            drawnow
            % The whole event should take a minimum of 10 seconds
            pause(min(0, 10/N - toc(timer)));
        end
        playback = false;
    end

    function stopPlayBackTof(src, event)
        playback = false;
    end

    function onClickWell(src, event, wk)

        winj = repmat(D.inj,[1,np]);
        wpro = rldecode(D.prod,ni,2);

        isInj = ismember(wk, D.inj);
        if isInj
            ik = find(D.inj == wk);
            sub = winj == wk;

            % well pair stuff
            wp = @(x) WP.inj(x);
            otherNames = {W(D.prod).name};

            % set plots to match piecharts
            v = find(strcmpi(get(hdataset, 'String'), 'drainage region'));
            set(hdataset, 'Value', v(1))
            set(floodtoggle, 'Value', ik)
            set(draintoggle, 'Value', []);
        else
            ik = find(D.prod == wk);
            sub = wpro == wk;

            % well pair stuff
            wp = @(x) WP.prod(x);
            otherNames = {W(D.inj).name};

            v = find(strcmpi(get(hdataset, 'String'), 'flooding region'));
            set(hdataset, 'Value', v(1))
            set(draintoggle, 'Value', ik)
            set(floodtoggle, 'Value', []);
        end

        if isempty(wellPlot) || ~ishandle(wellPlot)
            wellPlot = figure;
        else
            set(0, 'CurrentFigure', wellPlot); clf
        end

        plotArrival = ~isempty(opt.state) && ~isInj;

        if plotArrival;
            subplot(2, 2, 1);
        else
            subplot(2, 2,[1 3])
        end

        pie(max(WP.vols(sub), eps), ones(size(WP.vols(sub))))
        title('Pore volumes')
        if plotArrival
            subplot(2,2,3);  cla;
            plotTOFArrival(opt.state, W, pv, opt.fluid, find(D.prod == wk), D, opt.useMobilityArrival)
        end

        subplot(2,2,[2 4])
        tmp = wp(ik);
        if numel(tmp.z) > 1
            % Allocation factors by depth does not make sense for only one
            % perforation!
            [z, zind] = sort(tmp.z, 'descend');
            alloc = tmp.alloc(zind, :);
            % Always some allocation - avoid division by zero
            alloc(alloc == 0) = eps;
            walloc = bsxfun(@rdivide, cumsum(alloc,1), sum(alloc(:)));
            area(z, walloc, eps); axis tight;
            hold on
            plot(z, zeros(numel(z, 1)), '>', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
            % Flip it around, xlabel is really ylabel
            view(90, -90);
            set(gca, 'XDir', 'reverse')
            legend(otherNames, 'Location', 'EastOutside');
            xlabel('Depth')
            title(['Allocation factors by depth for ', W(wk).name]);
        end

        plotMain();
    end

    function plotPhi(src, event)
        if isempty(phiPlot) || ~ishandle(phiPlot)
            phiPlot = figure;
        else
            set(0, 'CurrentFigure', phiPlot);
        end
        [F,Phi] = computeFandPhi(pv,D.tof);
        plot(Phi,F,'.');
        title(sprintf('Lorenz coefficient: %f\n', computeLorenz(F,Phi)));
    end

    function plotWellConnections(src, event)

        if get(src, 'Value')
            if ishandle(fig_main)
                set(0, 'CurrentFigure', fig_main);
            else
                fig_main = figure;
                axis tight off
                plotMain();
            end
            if ~any(ishandle(pwc_wch))
                pwc_wch = plotWellPairConnections(G, WP, D, W, pv);
            end
        else

            delete(pwc_wch)
            pwc_wch = [];
        end

    end

    function [cdata clim cmap] = selectDataset()
        dataind = get(hdataset, 'Value');
        datanames = get(hdataset, 'String');
        clim = [];
        cmap = 'jet';
        switch lower(datanames{dataind})
            case 'forward tof'
                cdata = log10(D.tof(:,1));
                clim = log10(convertFrom(tofext, year));
            case 'backward tof'
                cdata = log10(D.tof(:,2));
                clim = log10(convertFrom(tofext, year));
            case 'sum of tofs'
                data = sum(D.tof(:,1:2), 2);
                cdata = log10(data);
                clim = log10(convertFrom([min(data), tofext(2)], year));
            case 'drainage region'
                cmap = 'gray';
                cdata = D.ppart;
                clim = [1, max(D.ppart)];
            case 'flooding region'
                cdata = D.ipart;
                clim = [1, max(D.ipart)];
            case 'porosity'
                cdata = rock.poro;
            case 'x permeability'
                cdata = log10(rock.perm(:,1));
            case 'y permeability'
                cdata = log10(rock.perm(:,2));
            case 'z permeability'
                cdata = log10(rock.perm(:,3));
            otherwise
                cdata = readStructField(datasets, datanames{dataind});
        end
        if isempty(clim)
            m = min(cdata(:));
            M = max(cdata(:));

            clim = [m - eps(m), M + eps(m)];
        end
    end

    function computeValues()
        if opt.computeFlux
            rS = initState(G, W, 0);
            T  = computeTrans(G, rock);
            rS = incompTPFA(rS, G, T, opt.tracerfluid, 'wells', W, 'LinSolve', opt.LinSolve);
        else
            rS = initState(G, W, 0);
            if isfield(opt.state, 'wellSol')
                rS.wellSol = opt.state.wellSol;
            end
            rS.flux     = opt.state.flux;
            rS.pressure = opt.state.pressure;
        end
        D = computeTOFandTracer(rS, G, rock, 'wells', W);
        D.itracer(isnan(D.itracer)) = 0;
        D.ptracer(isnan(D.ptracer)) = 0;
        % Cap tof to maximum tof for unreachable areas for the time being
        tf = D.tof(:);
        D.tof(isinf(D.tof)) = max(tf(isfinite(tf)));

        WP = computeWellPairs(rS, G, rock, W, D);
    end

    function changeWells()
        W = editWells(G, W, rock);
        recomputeValues();
        createMainControl();
    end

    function recomputeValues()
        computeValues()
        if ishandle(fig_main)
            close(fig_main)
        end
        plotMain();
    end
end



function [sliderhandle, edithandle] = linkedSlider(parent, pos, fieldsize, ext, defaultval, title)
    x = pos(1);
    y = pos(2);
    minval = abs(ext(1));
    maxval = abs(ext(2));
    dims = pos(3:4);
    defaultval = abs(defaultval);
    
    uicontrol(parent, 'Style', 'text', 'Units', 'normalized', 'Position', [x, y, fieldsize*dims(1) dims(2)], 'string', title)

    cap = @(x) max(minval, min(x, maxval));

    edithandle = uicontrol(parent, 'Style', 'edit', 'Units', 'normalized', 'Value', defaultval, 'Position', [x + (1-fieldsize)*dims(1), y, fieldsize*dims(1) dims(2)], 'String', sprintf('%.1f', defaultval));
    fun = @(src, event) set(edithandle, 'String', sprintf('%.1f', get(src, 'Value')));

    sliderhandle = uicontrol(parent, 'Style', 'slider', 'Units', 'normalized', 'Position', [x + fieldsize*dims(1), y, (1-2*fieldsize)*dims(1) dims(2)], 'Min', minval, 'Max', maxval, 'Value', defaultval, 'Callback', fun);
    fun2 = @(src, event) set(sliderhandle, 'Value', cap(sscanf(get(src, 'String'), '%f')));
    set(edithandle, 'Callback', fun2);
end

