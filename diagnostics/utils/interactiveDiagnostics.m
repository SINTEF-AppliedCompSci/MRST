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
%     A player controller (experimental), if a time-series of wells (and
%     states) are supplied.
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
%  'computeGrid' - Use G for plotting, but computeGrid for computing time
%                  of flight etc.
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
%  'daspect' - Data aspect ratio, in a format understood by daspect()
%
%  'name' - Name to use for windows
%
%  'leaveOpenOnClose' - Default false. Leaves all figures open when closing
%                       the controller.
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
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
       require mrst-gui incomp
    catch %#ok<CTCH>
       mrstModule add mrst-gui incomp
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
    
    opt = struct('state',               [],...
                 'computeGrid',         [],...
                 'tracerfluid',         water, ...
                 'LinSolve',            @mldivide, ...
                 'computeFlux',         true, ...
                 'useMobilityArrival',  false,...
                 'fluid',               oilwater, ...
                 'name',                [], ...
                 'daspect',             [], ...
                 'leaveOpenOnClose',    false, ...
                 'lineWells',           false, ...
                 'maxTOF',              [] ...
    );

    opt = merge_options(opt, varargin{:});
    
    if (isempty(opt.name))
        opt.name = {'Interactive Diagnostics'};
    elseif (isa(opt.name, 'char'))
        opt.name = {opt.name};
    end

    assert(opt.computeFlux || ~isempty(opt.state),...
        'If computeFlux is off a state must be provided!')
    
    if (isstruct(W))
        W = {W};
    end
    
    %Create anonymous function for remapping
    %values from the compute grid to fit the plotting grid G
    cdataToPlotGrid = @(f) f;
    wellsToPlotGrid = @(f) f;
    computeGrid = G;
    if (~isempty(opt.computeGrid))
        assert(isfield(G.cells, 'eMap'), ...
            'G must have an eMap field when using a separate computeGrid');
        cdataToPlotGrid = @(cdata) cdata(G.cells.eMap);
        wellsToPlotGrid = @(wells) remapWells(G, wells);
        computeGrid = opt.computeGrid;
    end
        
    state = [];
    state_idx = 1;
    if (~isempty(opt.state))        
        if (numel(opt.state) == 1)
            state = {opt.state};
        else
            state = opt.state;
        end
        
        assert(numel(W) == numel(state), ...
            'W and state must have equal number of elements');
    end
    

    if (~isempty(datasets))
        if (numel(datasets) == 1)
            datasets = {datasets};
        end
        
        assert(numel(W) == numel(datasets), ...
            'W and datasets must have equal number of elements');
    elseif(isempty(datasets) && ~isempty(state))
        dsname = 'State';
        datasets = state;
    else
        dsname = '';
        datasets = cell(1,1);
    end
    


    % Main scope variables for tracers and tof
    [D, D_int, WP] = deal([]);

    % Main figures
    [fig_main,... 
     fig_ctrl,...
     fig_well,...
     fig_phi,...
     fig_well_alloc,...
    ] = deal(NaN);

    %Gui elements
    [ctrl_drain_vols,...
     ctrl_flood_vols,...
     ctrl_well_conn,...
     ctrl_use_average,...
     ctrl_show_grid,...
     speedsh,...
     tofsh,...
     nwtofsh,...
     mtofsh,...
     mtofeh,...
     Mtofsh,...
     Mtofeh,...
     alfash,...
     tofext,...
     hdataset,...
     hset_op,...
     mrst_ds,...
     ds_panel] = deal(NaN);
    
    %Handles for well plot in main figure
    fig_main_wells = {};
    [fig_main_wells.htop,...
     fig_main_wells.htext,...
     fig_main_wells.hs,...
     fig_main_wells.hline] = deal([]); 
    fig_main_wells.dirty = true;
    
    %Handles for cell data in main figure
    fig_main_cell_data = []; % "Grid with cell values"
    fig_main_grid = []; %Grid outline
    
    %Well connection pairs in main figure
    fig_main_well_conns = []; 
    
    %Handles that should be enabled / disabled
    %according to if we're using averaged or single 
    %timestep values
    average_gui_handles = [];
    non_average_gui_handles = [];
    
    %"GUI State" variables
    window_name = opt.name{state_idx};
    compute_average = false;
    playback = false; %Play back TOF?
    time_varying = (~isempty(state) && numel(state) > 1) || (~isempty(W) && numel(W) > 1);
    selection = []; % Selection of cells
    
    
    
    %Create main figure
    fig_main = figure('Name', window_name);    
    
    %Set scaling etc.
    if (~isempty(opt.daspect))
        daspect(opt.daspect);
    else
        max_G = max(G.cells.centroids);
        min_G = min(G.cells.centroids);
        da = max_G - min_G;
        da(da == 0) = 1;
        daspect(da);
    end
    % Create control panel
    % Precompute TOF etc. for creating main control
	pv = poreVolume(computeGrid, rock);
    [D, WP] = computeTOFAndTracerAndWellPairs(W{state_idx}, state{state_idx});
    tofext = getTOFRange(D, opt);
    createMainControl();
    
    % Trigger plot initial setup
    plotMain();
    axis tight off
    view(3);

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Only function definitions from here on %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Clamp real: Clamps data to real valued numbers (non-inf)
    function val = clamp_real(in_val)
        val = max(-realmax, min(realmax, in_val));
    end

    %
    function tofext = getTOFRange(D, opt)
        if (D.isvalid)
            if ~isempty(opt.maxTOF)
                tofext = convertTo([0, opt.maxTOF], year);
            else
            tof = D.tof(:);
            isNeg = (tof <= 0);
            tof(isNeg) = min(tof(~isNeg));

            tofext = convertTo(([min(tof), 5*10^(mean(log10(tof)))]), year);
            tofext(2) = clamp_real(max(tofext(2), 15));
            end
        else
            tofext = [0 1];
        end
    end
        
    function setTOFSliderExtents()
        %Update min/max and current value
        min_val = tofext(1);
        max_val = tofext(2);
        cur_m_val = str2double(get(mtofeh, 'String'));
        cur_m_val = max(min(cur_m_val, max_val), min_val);

        cur_M_val = str2double(get(Mtofeh, 'String'));
        cur_M_val = max(min(cur_M_val, max_val), min_val);

        set(mtofsh, 'Min', tofext(1));
        set(mtofsh, 'Max', tofext(2));
        set(mtofsh, 'Value', cur_m_val);
        %set(mtofeh, 'String', num2str(cur_m_val));
        mtofs_callback([], []);

        set(Mtofsh, 'Min', tofext(1));
        set(Mtofsh, 'Max', tofext(2));
        set(Mtofsh, 'Value', cur_M_val);
        %set(Mtofeh, 'String', num2str(cur_M_val));
        Mtofs_callback([], []);
    end

    % Set special functions for the min/max time of flight slider
    % handle
    function mtofs_callback(src, event)
        min_val = get(mtofsh, 'Min');
        max_val = get(mtofsh, 'Max');
        cur_val = get(mtofsh, 'Value');
        if (cur_val <= min_val)
            set(mtofsh, 'Value', min_val);
            set(mtofeh, 'String', sprintf('%.1f', 0));
        elseif (cur_val >= max_val)
            set(mtofsh, 'Value', max_val);
            set(mtofeh, 'String', sprintf('%.1f', Inf));
        else
            set(mtofeh, 'String', sprintf('%.1f', cur_val));
        end

        %Only plot if src is nonempty (true callback)
        if (~isempty(src))
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
    end

    function Mtofs_callback(src, event)
        min_val = get(Mtofsh, 'Min');
        max_val = get(Mtofsh, 'Max');
        cur_val = get(Mtofsh, 'Value');
        if (cur_val <= min_val)
            set(Mtofsh, 'Value', min_val);
            set(Mtofeh, 'String', sprintf('%.1f', -Inf));
        elseif (cur_val >= max_val)
            set(Mtofsh, 'Value', max_val);
            set(Mtofeh, 'String', sprintf('%.1f', Inf));
        else
            set(Mtofeh, 'String', sprintf('%.1f', cur_val));
        end

        %Only plot if src is nonempty (true callback)
        if (~isempty(src))
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
    end


    function addVisualizationControls(parent, pos)
        vis_N = 3;
        vis_H = 1 / vis_N;
        vis_h = vis_H*0.9;
        bwidth = .2;
        
        vis_control = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'None');
        
        [alfash, ~, ~]   = linkedSlider(vis_control, [0 2*vis_H 1 vis_h], .15, [0 1], 1, 'Alpha', @plotMain);
        
        
               
        ctrl_show_grid = uicontrol(vis_control, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [0 1*vis_H bwidth vis_h],...
                   'String', 'Show grid',...
                   'Callback', @plotMain, ...
                   'Value', 0 ...
                   );
               
        ctrl_well_conn = uicontrol(vis_control, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [bwidth 1*vis_H bwidth vis_h],...
                   'String', 'Well pairs',...
                   'Callback', @plotWellConnections ...
                   );
               
               
        confp = uipanel(vis_control, 'Units', 'normalized', 'Position', [0 0*vis_H 1 vis_h], 'BorderType', 'None');
        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0 bwidth 1],...
                   'String', 'Phi/F diagram',...
                   'Callback', @(src, event) plotPhi(src, event)...
                   );
               
        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [1*bwidth 0 bwidth 1],...
                   'String', 'Well allocations',...
                   'Callback', @plotWellAllocations ...
                   ); 

        if opt.computeFlux
            uicontrol(confp, 'Style', 'pushbutton',...
                       'Units', 'normalized',...
                       'Position', [2*bwidth 0 bwidth 1],...
                       'String', 'Edit wells',...
                       'Callback', @(src, event) changeWells()...
                       );
        end
    end


    function addTOFControls(parent, pos)
        %Function that forces update of time integral TOF
        function dirtyTOFPlotMain(src, event)
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
        
        % Panel for Time of flight selection parameters
        tof_panel = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'none');
        tof_N = 7;
        tof_H = 1/tof_N;
        tof_h = tof_H*0.9;
        
        cellfields = getStructFields(G, datasets{state_idx}, dsname);
        perm = strcat({'X', 'Y', 'Z'}, ' permeability');
        uicontrol(tof_panel, 'Style', 'text', ...
                             'string', 'Display:', ...
                             'HorizontalAlignment', 'left', ...
                             'Units', 'normalized',...
                             'Position', [.0 6*tof_H .25 tof_h]);
                
        hdataset = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                    'Units', 'normalized',...
                                    'Position', [.2 6*tof_H .33 tof_h],...
                                    'Callback', @plotMain,...
                                    'String', {'Forward TOF',...
                                               'Backward TOF',...
                                               'Sum of TOFs',...
                                               'Drainage region',...
                                               'Flooding region',...
                                               'Drainage blend',...
                                               'Flooding blend',...
                                               'Porosity', ...
                                               perm{1:size(rock.perm, 2)},...
                                               cellfields{:}}...
                                   );
        hset_op = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                   'Units', 'normalized',...
                                   'Position', [.66 6*tof_H .33 tof_h],...
                                   'Callback', @plotMain,...
                                   'String', {'Union {Flood, Drain} volumes',...
                                              'Intersection {Flood, Drain} volumes',...
                                              'Flood volumes',...
                                              'Drain volumes'}...
                                   );
        
        [nwtofsh, ~, ~]   = linkedSlider(tof_panel, [0 5*tof_H 1 tof_h], .15, [0 1], 0.2, 'Near well TOF', @dirtyTOFPlotMain);
        [Mtofsh, Mtofeh, ~]   = linkedSlider(tof_panel, [0 4*tof_H 1 tof_h], .15, tofext, tofext(2), 'Max TOF', @dirtyTOFPlotMain);
        [mtofsh, mtofeh, ~]   = linkedSlider(tof_panel, [0 3*tof_H 1 tof_h], .15, tofext, tofext(1), 'Min TOF', @dirtyTOFPlotMain);
        set(mtofsh, 'Callback', @mtofs_callback);
        set(mtofeh, 'String', num2str(-Inf))
        set(Mtofsh, 'Callback', @Mtofs_callback);
        set(Mtofeh, 'String', num2str(+Inf))
        
        if (time_varying)
            [tofsh, tofeh, ~]   = linkedSlider(tof_panel, [0 2*tof_H 1 tof_h], .15, [0 1], 0.25, 'TOF Freq.', @plotMain);
            average_gui_handles = [average_gui_handles(:); tofsh; tofeh];
        end
        
        
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 1*tof_H .2 tof_h],...
                   'String', 'Play TOF',...
                   'Callback', @playBackTof...
                   );
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [.2 1*tof_H .2 tof_h],...
                   'String', 'Stop TOF',...
                   'Callback', @stopPlayBackTof...
                   );
        [speedsh, ~, ~] = linkedSlider(tof_panel, [.45 1*tof_H .55 tof_h], .15, [1 1000], 50, 'Steps', []);  
        
        
        
        % Function which saves variables to the workspace on demand
        function saveDataToWorkSpace(src, event)
            indices = find(selection);
            
            labels = {'Sector model', 'Sector model (index map)'};
            varnames = labels;
            vars = {indices, computeGrid.cells.indexMap(indices)};
            
            %Clean the varnames of spaces etc
            clean = @(x) lower(regexprep(x,'[^a-zA-Z0-9]','_'));
            varnames = cellfun(@(x) clean(x), varnames, 'UniformOutput', false);

            export2wsdlg(labels, varnames, vars)
        end
        
        function saveDataToFile(src, event)
            title = 'Export selection to FLUXNUM';
            num_lines = 1;
            prompt = {'Enter filename:'};
            default = {'FLUXNUM.PROP'};
            
            while (true)
                answer = inputdlg(prompt, title, num_lines, default);
                if (numel(answer) == 0)
                    return;
                end

                filename = answer{1};
                if (exist(filename, 'file'))
                    choice = questdlg('File already exists. Do you want to overwrite?', ...
                        'File exists', ...
                        'Yes','No', 'No');
                    if (strcmp(choice, 'Yes'))
                        break;
                    end
                else
                    break;
                end
            end
            
            [fid, msg] = fopen(filename, 'wt');
            if fid < 0
                errordlg(['Failed to open output file ', msg, '. Aborting.']);
                return
            end
            
            fprintf(fid, 'FLUXNUM\n');
            var = zeros(prod(computeGrid.cartDims), 1);
            var( computeGrid.cells.indexMap(selection) ) = 1;
            out_format = [repmat('%d\t', 1, 20), '\n'];
            fprintf(fid, out_format, var);
            fprintf(fid, '/\n');
            fclose(fid);
        end
        
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0*tof_H 0.3 tof_h],...
                   'String', 'Save to workspace', ...
                   'Callback', @saveDataToWorkSpace...
                   );
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0.3 0*tof_H 0.3 tof_h],...
                   'String', 'Save to file', ...
                   'Callback', @saveDataToFile...
                   );
    end


    function addTimestepSelector(parent, pos)
        ds_panel = uipanel('Parent', parent,...
                'Units', 'normalized',...
                'Position', pos,...
                'BorderType', 'none');
            
        uicontrol(ds_panel,'Style','text',...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'String','Simulation timestep',...
            'Position', [0 0.5 0.4 .5]);

        selector_datasets = [];
        if (numel(state) > 1)
            selector_datasets = state;
        elseif(numel(W) > 1)
            selector_datasets = W;
        else
            error('Programmer error');
        end
        
        function selectTimeStep(datas, index, fullname)
            %Get ID of currently selected wells
            drain_wells = get(ctrl_drain_vols, 'Value');
            flood_wells = get(ctrl_flood_vols, 'Value');
            drain_selection = D.prod(drain_wells);
            flood_selection = D.inj(flood_wells);
            well_selection = [drain_selection, flood_selection];
            fig_main_wells.dirty = true;
            
            state_idx = index;
            computeValues(state_idx);
            
            %Update main plot
            plotMain();
            
            %Update producers / injectors selection box text
            set(ctrl_drain_vols, 'String', {W{state_idx}(D.prod).name});
            set(ctrl_flood_vols, 'String', {W{state_idx}(D.inj).name});
            
            %Find out which items we should be selecting
            drain_wells = find(ismember(D.prod, well_selection));
            flood_wells = find(ismember(D.inj, well_selection));
            
            %Set selections
            set(ctrl_drain_vols, 'Value', drain_wells);
            set(ctrl_flood_vols, 'Value', flood_wells);
            
            %Update aux plots
            if (ishandle(fig_phi))
                plotPhi([], []);                
            end
            if (ishandle(fig_well_alloc))
                plotWellAllocations([], [])
            end
        end

        datasetSelector(G, selector_datasets, 'Parent', ds_panel, 'Location', ...
            [0.4, 0.5, 0.6, .5], 'Callback', @selectTimeStep, ...
            'active', 1, 'Nofields', true,...
            'Tag', 'datasetselector');
        mrst_ds = findobj('Tag', 'datasetselector', '-and', 'Parent', ds_panel);
        set(mrst_ds, 'BorderType', 'None');
        mrst_ds_children = allchild(mrst_ds);
        non_average_gui_handles = [non_average_gui_handles(:); mrst_ds_children];

        %Control for showing / computing average over time of e.g., TOF
        ctrl_use_average = uicontrol(ds_panel, 'Style', 'checkbox',...
                                   'Units', 'normalized',...
                                   'Position', [0 0 1 0.5],...
                                   'Callback', @handleIntegralChange,...
                                   'HorizontalAlignment', 'right',...
                                   'String', 'Average over timesteps'...
                                   );                               
                               
        function handleIntegralChange(src, event)
            %Options specific when using the average over timesteps of the TOF
            time_integral_options = {
                'Forward TOF frequency',...
                'Backward TOF frequency',...
                'Sum TOF frequency'};

            %Check if we are to compute the average over time
            %or for an instantaneous snapshot.
            if (ishandle(ctrl_use_average))
                compute_average = logical(get(ctrl_use_average, 'Value'));
                if (compute_average)
                    set(average_gui_handles, 'Enable', 'on');
                    set(non_average_gui_handles, 'Enable', 'off');
                    
                    %Add new integral options
                    old_options = get(hdataset, 'String');
                    new_options = {old_options{:}, time_integral_options{:}};
                    set(hdataset, 'String', new_options);
                    
                    computeAverage();
                else
                    set(average_gui_handles, 'Enable', 'off');
                    set(non_average_gui_handles, 'Enable', 'on');

                    %Remove integral options
                    options = get(hdataset, 'String');
                    for idx=1:numel(time_integral_options)
                        mask = ~ismember(options, time_integral_options{idx});
                        options = options(mask);
                    end
                    set(hdataset, 'String', options);
                    curr_val = get(hdataset, 'Value');
                    if (curr_val > numel(options))
                        set(hdataset, 'Value', 1);
                    else
                        set(hdataset, 'Value', curr_val);
                    end
                    
                    computeValues(state_idx);
                end
                
            end
            
            %Only plot if src is nonempty (true callback)
            if (~isempty(src))
                plotMain(src, event);
            end
        end
    end



    function createMainControl()
        % Set up figure handles
        if ~ishandle(fig_ctrl)
            pos = get(fig_main, 'OuterPosition');
            size_xy = [475 550];
            pos_xy = pos(1:2) + pos(3:4) - size_xy - [400 0];
            fig_ctrl = figure('Position',[pos_xy,  size_xy], ...
                'Toolbar','none',...
                'MenuBar', 'none',...
                'CloseRequestFcn', @close_all_func);
            set(fig_ctrl, 'Name', ['Controller ', window_name]);
        else
            clf(fig_ctrl)
        end

        % Drainage / flooding controls
        np = numel(D.prod);
        ni = numel(D.inj);

        ctrl_drain_vols = zeros(np, 1);
        ctrl_flood_vols = zeros(ni, 1);
        
        gpy = 0.45;
        gph = .55;
        if (ishandle(mrst_ds))
            gpy = gpy + 0.15;
            gph = gph - 0.15;
        end

        drainp = uipanel(fig_ctrl, 'Title', 'Drainage volumes',...
                                   'Units', 'normalized',...
                                   'Position', [0 gpy .5 gph]);
        floodp = uipanel(fig_ctrl, 'Title', 'Flooding volumes',...
                                   'Units', 'normalized',...
                                   'Position', [.5 gpy .5 gph]);

        function dirtyWellsPlotMain(src, event)
            fig_main_wells.dirty = true;
            plotMain(src, event);
        end
        
        ctrl_drain_vols = uicontrol(drainp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.prod).name},...
                                          'Value', [], ...
                                          'Callback', @dirtyWellsPlotMain, ...
                                          'Position', [0 0 1 1]);

        ctrl_flood_vols = uicontrol(floodp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.inj).name},...
                                          'Value', [], ...
                                          'Callback', @dirtyWellsPlotMain, ...
                                          'Position', [0 0 1 1]);
        
                                      
        tabgp = uitabgroup(fig_ctrl,'Position', [0 0 1 .4]);
        
        tof_controls_tab = uitab(tabgp, 'Title', 'Time of flight region');
        addTOFControls(tof_controls_tab, [0 .1 1 .9]); 
        
        vis_controls_tab = uitab(tabgp, 'Title', 'Visualization config');
        addVisualizationControls(vis_controls_tab, [0 .6 1 .4]); 
        
        %Add time-varying dataset slider
        if (time_varying)
            data_controls_tab = uitab(tabgp, 'Title', 'Dataset selection');
            addTimestepSelector(data_controls_tab, [0 0.7 1 0.3]);
        end
        
        %Enable/disable the correct controls by firing callback on
        %time_integral button
        if (ishandle(ctrl_use_average))
            f = get(ctrl_use_average, 'Callback');
            f([], []);
        end
    end

    function computeAverage()
        if (~isempty(D_int))
            D = D_int;
        else
            extra_args = {};
            if (isfield(state{1}, 'time'))
                times = cellfun(@(x) x.time, state);
                times = [0; times];
                dt = (times(2:end) - times(1:end-1)) / times(end);
                extra_args = {extra_args{:}, 'dt', dt};
            end
            
            min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
            max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);

            D_int = computeTOFandTracerAverage(state, computeGrid, rock, ...
                'wells', W, ...
                'max_tof', max_tof,...
                'min_tof', min_tof,...
                extra_args{:});
                
            D = D_int;
        end
    end

    function selection = getTOFSelection()
        %Get options for TOF region selection
        drain_wells = get(ctrl_drain_vols, 'Value');
        flood_wells = get(ctrl_flood_vols, 'Value');
        near_well_max_tof = convertFrom(get(nwtofsh, 'Value'), year);
        min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
        max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);
        switch(get(hset_op, 'Value'))
            case 1
                set_op = 'union';
            case 2
                set_op = 'intersection';
            case 3
                set_op = 'flood';
            case 4
            	set_op = 'drain';
        end
        
        %Get TOF region
        if (compute_average)
            % Get current TOF frequency threshold
            tof_freq = get(tofsh, 'Value');
            
            isubset = D.isubset;
            psubset = D.psubset;
            
            frequency = (psubset+isubset) / 2.0;
            psubset = psubset & (frequency >= tof_freq);
            isubset = isubset & (frequency >= tof_freq);
            
            selection = selectTOFRegion(D, max_tof, min_tof, ...
                                        'drain_wells', drain_wells, ...
                                        'flood_wells', flood_wells, ...
                                        'set_op', set_op, ...
                                        'near_well_max_tof', near_well_max_tof, ...
                                        'psubset', psubset, ...
                                        'isubset', isubset);
        else
            selection = selectTOFRegion(D, max_tof, min_tof, ...
                                        'drain_wells', drain_wells, ...
                                        'flood_wells', flood_wells, ...
                                        'set_op', set_op, ...
                                        'near_well_max_tof', near_well_max_tof);
        end
    end

    function plotMain(src, event)
        if (compute_average)
            window_name = 'Average over timesteps';
        else
            if (numel(opt.name) == 1)
                window_name = [opt.name{1}, ' timestep ', num2str(state_idx)];
            else
                window_name = opt.name{state_idx};
            end
        end
        
        if(ishandle(fig_ctrl))
            set(fig_ctrl, 'Name', ['Controller ', window_name]);
        end            
            
        if ishandle(fig_main)
            set(0, 'CurrentFigure', fig_main);
        else
            fig_main = figure();
        end
        set(fig_main, 'Name', window_name);

        plot_grid = get(ctrl_show_grid, 'Value');
        if (plot_grid)
            if ~any(ishandle(fig_main_grid))
                fig_main_grid = plotGrid(G, 'facec', 'none', 'edgea', .05, 'edgec', 'black');
                set(fig_main_grid, 'UserData', 'gridoutline');
            end
        else
            if (any(ishandle(fig_main_grid)))
                delete(fig_main_grid);
            end
        end
        
        tofext = getTOFRange(D, opt);
        setTOFSliderExtents();
        selection = getTOFSelection();
        
        dataind = get(hdataset, 'Value');
        [cdata, clim, ~] = selectDataset(dataind);

        alpha   = get(alfash, 'Value');

        % Plot current selection
        if any(ishandle(fig_main_cell_data))
            delete(fig_main_cell_data);
        end
        if any(selection)
            if (dataind == 6 || dataind == 7)
                fig_main_cell_data = plotTracerBlend(G, ...
                    cdataToPlotGrid(cdata{1}), ...
                    cdataToPlotGrid(cdata{2}), ...
                    'cells', cdataToPlotGrid(selection), ...
                    'FaceAlpha', alpha);
            else
                fig_main_cell_data = plotCellData(G, ...
                    cdataToPlotGrid(cdata), cdataToPlotGrid(selection), ...
                    'EdgeColor', 'none', 'FaceAlpha', alpha);
            end
        end

        fastRotateButton();

        plotWells();
        
        if (clim(2) > clim(1))
            caxis(clim);
        end
        
        %Update plot of wells (if enabled)
        plotWellConnections();
    end

    function playBackTof(src, event)
        if playback
            return
        end
        N = round(get(speedsh, 'Value'));
        
        if (compute_average)
            t_0 = 1;
            t_end = get(tofsh, 'Value');
            callback = get(tofsh, 'Callback');
        else
            t_0 = get(mtofsh, 'Value');
            t_end = get(Mtofsh, 'Value');
            callback = get(Mtofsh, 'Callback');
            if t_0 >= t_end
                disp 'Minimum TOF is larger than maximum TOF...'
                return
            end
        end
        playback = true;
        % Pop figure to front explicitly before loop
        figure(fig_main);
        for i = 0:N
            if ~playback
                return
            end
            
            if (compute_average)
                set(tofsh, 'Value', t_0 + i*(t_end - t_0)/N);
                callback(tofsh, []);
            else
                set(Mtofsh, 'Value', min(t_0 + i*(t_end - t_0)/N, tofext(2)));
                callback(Mtofsh, []);
            end
            
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
        np = numel(D.prod);
        ni = numel(D.inj);
            
        winj = repmat(D.inj,[1,np]);
        wpro = rldecode(D.prod,ni,2);
        
        wp = [];

        isInj = ismember(wk, D.inj);
        if isInj
            ik = find(D.inj == wk);
            sub = winj == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.inj(x);
            end
            otherNames = {W{state_idx}(D.prod).name};

            % set plots to match piecharts
            v = find(strcmpi(get(hdataset, 'String'), 'drainage region'));
            set(hdataset, 'Value', v(1))
            set(ctrl_flood_vols, 'Value', ik)
            set(ctrl_drain_vols, 'Value', []);
        else
            ik = find(D.prod == wk);
            sub = wpro == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.prod(x);
            end
            otherNames = {W{state_idx}(D.inj).name};

            v = find(strcmpi(get(hdataset, 'String'), 'flooding region'));
            set(hdataset, 'Value', v(1))
            set(ctrl_drain_vols, 'Value', ik)
            set(ctrl_flood_vols, 'Value', []);
        end

        if isempty(fig_well) || ~ishandle(fig_well)
            fig_well = figure();
        else
            set(0, 'CurrentFigure', fig_well); clf
        end
        set(fig_well, 'name', ['Well ', W{state_idx}(wk).name]);

        plotArrival = ~isempty(state) && ~isInj;

        if plotArrival;
            subplot(2, 2, 1);
        else
            subplot(2, 2,[1 3])
        end

        if (numel(WP) > 0)
            tmp = wp(ik);
            salloc = abs(sum(tmp.alloc,1));
            [sa, tmp_ix] = sort(salloc, 'descend');
            ssa = sum(sa);
            ix = min(5, find(sa > .01*ssa, 1, 'last'));
            [sa, tmp_ix] = deal(sa(ix:-1:1), tmp_ix(ix:-1:1));
            pie(sa/ssa, otherNames(tmp_ix))
            %pie(max(WP.vols(sub), eps), ones(size(WP.vols(sub))))
            title('Well allocation factors')
            if isInj
               set(ctrl_drain_vols, 'Value', tmp_ix);
            else
               set(ctrl_flood_vols, 'Value', tmp_ix); 
            end
            fig_main_wells.dirty = true;
        end
        if plotArrival
            subplot(2,2,3);  cla;
            plotTOFArrival(state{state_idx}, W{state_idx}, pv, opt.fluid, find(D.prod == wk), D, opt.useMobilityArrival, 'inj_ix', tmp_ix)
            
        end

        
        if (numel(WP) > 0)
            if numel(tmp.z) > 1
                subplot(2,2,[2 4])
                % Allocation factors by connection does not make sense for only one
                % perforation!
                %[z, zind] = sort(tmp.z, 'descend');
                alloc  = abs(tmp.alloc(:, tmp_ix));
                na = size(alloc,1);
                % reverse cumsum
                calloc = cumsum(alloc, 1, 'reverse');
                % Use standard units
                calloc = convertTo(calloc, 1/day);
                area((1:na)', calloc); axis tight;
                hold on
                plot((1:na)', zeros(na,1), '>', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
                % Flip it around, xlabel is really ylabel
                view(90, -90);
                set(gca, 'XDir', 'reverse')
                %legend(otherNames, 'Location', 'EastOutside');
                xlabel('Connection #')
                ylabel('Accumulated flux [m^3/day]')
                title('Allocation by connection');
            end
        end

        plotMain();
    end

    function plotPhi(src, event)
        if (~D.isvalid)
            warndlg('It appears that the time of flight is invalid. Aborting.');
            return
        end
        
        if isempty(fig_phi) || ~ishandle(fig_phi)
            fig_phi = figure();
        else
            set(0, 'CurrentFigure', fig_phi);
        end
        set(fig_phi, 'Name', ['Phi/F diagram (', window_name, ')']);
        [F,Phi] = computeFandPhi(pv,D.tof);
        plot(Phi,F,'.');
        title(sprintf('Lorenz coefficient: %f\n', computeLorenz(F,Phi)));
    end

    function plotWells()
        % Plot wells
        if (fig_main_wells.dirty)
            if (any(ishandle(fig_main_wells.htop)))
                delete(fig_main_wells.htop);
            end
            if (any(ishandle(fig_main_wells.htext)))
                delete(fig_main_wells.htext);
            end
            if (any(ishandle(fig_main_wells.hs)))
                delete(fig_main_wells.hs);
            end
            if (any(ishandle(fig_main_wells.hline)))
                delete(fig_main_wells.hline);
            end
    
            %Get the selected wells only
            drain_wells = D.prod(get(ctrl_drain_vols, 'Value'));
            flood_wells = D.inj(get(ctrl_flood_vols, 'Value'));
            sel_wells = [drain_wells, flood_wells];
            
            if (numel(sel_wells) > 0)
                W_sel = W{state_idx}(sel_wells);

                [fig_main_wells.htop, ...
                    fig_main_wells.htext, ...
                    fig_main_wells.hs, ...
                    fig_main_wells.hline] = plotWell(G, wellsToPlotGrid(W_sel), ...
                                                    'color', 'red', 'height',  0);

                for j = 1:numel(W_sel)
                    i = sel_wells(j);
                    color = colorizeWell('global', i, D);
                    set([fig_main_wells.htop(j) fig_main_wells.htext(j) fig_main_wells.hs(j)], 'ButtonDownFcn', @(src, event) onClickWell(src, event, i));
                    set([fig_main_wells.htop(j) fig_main_wells.hs(j)], 'FaceColor', color, 'EdgeColor', color)
                    set([fig_main_wells.htext(j) fig_main_wells.hline(j)], 'Color', color)
                    set(fig_main_wells.htext(j), 'FontWeight', 'bold', 'Interpreter', 'none')
                end
            end
            
            fig_main_wells.dirty = false;
        end
    end

    function plotWellAllocations(~, ~)
        if (~D.isvalid)
            warndlg('It appears that the time of flight is invalid. Aborting');
            return
        end
        
        if isempty(fig_well_alloc) || ~ishandle(fig_well_alloc)
            fig_well_alloc = figure();
        else
            set(0, 'CurrentFigure', fig_well_alloc);
        end
        set(fig_well_alloc, 'Name', ['Well allocations (', window_name, ')']);
        
        %Filter well pairs to only include selected wells
        drain_wells = D.prod(get(ctrl_drain_vols, 'Value'));
        flood_wells = D.inj(get(ctrl_flood_vols, 'Value'));
        sel_wells = [drain_wells, flood_wells];
        
        plotWellAllocationComparison(D, WP, [], [], 'plotOnly', sel_wells);
    end

    function plotWellConnections(~, ~)
        if get(ctrl_well_conn, 'Value')
            if ishandle(fig_main)
                set(0, 'CurrentFigure', fig_main);
            else
                fig_main = figure;
                axis tight off
                plotMain();
            end
            
            if any(ishandle(fig_main_well_conns))
                delete(fig_main_well_conns)
                fig_main_well_conns = [];
            end
            
            if (D.isvalid)
                fig_main_well_conns = plotWellPairConnections(computeGrid, WP, D, W{state_idx}, pv);
            else
                return
            end
        else
            delete(fig_main_well_conns)
            fig_main_well_conns = [];
        end
    end

    function [cdata, clim, cmap] = selectDataset(dataind)
        datanames = get(hdataset, 'String');
        clim = [];
        cmap = 'jet';
        
        switch lower(datanames{dataind})
            case 'forward tof'
                cdata = log10(D.tof(:,1));
                clim = log10(clamp_real(convertFrom(tofext, year)));
            case 'backward tof'
                cdata = log10(D.tof(:,2));
                clim = log10(clamp_real(convertFrom(tofext, year)));
            case 'sum of tofs'
                data = sum(D.tof(:,1:2), 2);
                cdata = log10(data);
                clim = log10(clamp_real(convertFrom([min(data), tofext(2)], year)));
            case 'drainage region'
                cmap = 'gray';
                cdata = D.ppart;
                clim = clamp_real([1, max(D.ppart)]);
            case 'flooding region'
                cdata = D.ipart;
                clim = clamp_real([1, max(D.ipart)]);
            case 'drainage blend'
                cdata = {D.ppart, max(D.ptracer, [], 2)};
                clim = [1, max(D.ipart)];
            case 'flooding blend'
                cdata = {D.ipart, max(D.itracer, [], 2)};
                clim = [1, max(D.ipart)];
            case 'porosity'
                cdata = rock.poro;
            case 'x permeability'
                cdata = log10(rock.perm(:,1));
            case 'y permeability'
                cdata = log10(rock.perm(:,2));
            case 'z permeability'
                cdata = log10(rock.perm(:,3));
            case 'forward tof frequency'
                cdata = D.psubset;
            case 'backward tof frequency'
                cdata = D.isubset;
            case 'sum tof frequency'
                cdata = D.isubset + D.psubset;
            otherwise
                assert(isfield(datasets{state_idx}, datanames{dataind}), 'Trying to access non-existent field');
                cdata = readStructField(datasets{state_idx}, datanames{dataind});
        end
        
        if isempty(clim)
            m = min(cdata(:));
            M = max(cdata(:));

            clim = [m - eps(m), M + eps(m)];
        end
    end

    function [D, WP] = computeTOFAndTracerAndWellPairs(W_step, state_step)
        if opt.computeFlux
            rS = initState(G, W_step, 0);
            T  = computeTrans(G, rock);
            rS = incompTPFA(rS, G, T, opt.tracerfluid, 'wells', W_step, 'LinSolve', opt.LinSolve);
        else
            rS = initState(G, W_step, 0);
            if isfield(state_step, 'wellSol')
                rS.wellSol = state_step.wellSol;
            end
            rS.flux     = state_step.flux;
            rS.pressure = state_step.pressure;
        end
        
        D = computeTOFandTracer(rS, computeGrid, rock, 'wells', W_step, 'processCycles', true);
        D.itracer(isnan(D.itracer)) = 0;
        D.ptracer(isnan(D.ptracer)) = 0;
        
        % Cap tof to maximum tof for unreachable areas for the time being
        tf = D.tof(:);
        if (all(isinf(tf)))
            D.isvalid = false;
            WP = [];
        else
            D.tof(isinf(D.tof)) = max(tf(isfinite(tf)));
            D.isvalid = true;
            WP = computeWellPairs(rS, computeGrid, rock, W_step, D);
        end
    end

    function computeValues(idx)
        [D, WP] = computeTOFAndTracerAndWellPairs(W{idx}, state{idx});

        if (~D.isvalid || isempty(WP))
            warning('Time of flight returned inf. Are there both active injectors and producers present?')
        end
    end

    function changeWells()
        W{state_idx} = editWells(G, W{state_idx}, rock);
        computeValues();
        createMainControl();
        
        fig_main_wells.dirty = true;
        
        plotMain();
    end



    function [sliderhandle, edithandle, grouphandle] = linkedSlider(parent, pos, fieldsize, ext, defaultval, title, callback)
        minval = abs(ext(1));
        maxval = abs(ext(2));
        defaultval = abs(defaultval);

        grouphandle = uipanel(parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'None');
        uicontrol(grouphandle, 'Style', 'text', 'Units', 'normalized', 'Position', [0, 0, fieldsize 1], 'string', title)

        cap = @(x) max(minval, min(x, maxval));

        edithandle = uicontrol(grouphandle, 'Style', 'edit', 'Units', 'normalized', 'Value', defaultval, 'Position', [0 + (1-fieldsize), 0, fieldsize 1], 'String', sprintf('%.1f', defaultval));
        sliderhandle = uicontrol(grouphandle, 'Style', 'slider', 'Units', 'normalized', 'Position', [fieldsize, 0, (1-2*fieldsize) 1], 'Min', minval, 'Max', maxval, 'Value', defaultval);

        slidercallback2 = @(src, event) set(edithandle, 'String', sprintf('%.1f', get(src, 'Value')));
        editcallback2 = @(src, event) set(sliderhandle, 'Value', cap(sscanf(get(src, 'String'), '%f')));

        function slidercallback3(src, event)
            callback(src, event);
            slidercallback2(src, event);
        end
        function editcallback3(src, event)
            callback(src, event);
            editcallback2(src, event);
        end

        if (~isempty(callback))
            set(sliderhandle, 'Callback', @slidercallback3);
        else 
            set(sliderhandle, 'Callback', slidercallback2);
        end

        if (~isempty(callback))
            set(edithandle, 'Callback', @editcallback3);
        else
            set(edithandle, 'Callback', editcallback2);
        end
    end


    % Private helpers
    function close_all_func(src, event)
        if (~opt.leaveOpenOnClose)
            if any(ishandle(fig_ctrl))
                delete(fig_ctrl);
            end
            if any(ishandle(fig_well_alloc))
                delete(fig_well_alloc);
            end
            if any(ishandle(fig_phi))
                delete(fig_phi);
            end
            if any(ishandle(fig_main))
                delete(fig_main);
            end
            if any(ishandle(fig_well))
                delete(fig_well);
            end
        end
    end

end