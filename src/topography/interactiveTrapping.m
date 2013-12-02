function varargout = interactiveTrapping(inp, varargin)
%Create an interactive figure showing trapping structure for a top surface grid
%
% For a detailed list of functionality and controls, please see the
% showTrapsInteractively example.
%
% SYNOPSIS:
%       interactiveTrapping('Johansenfm')
%       interactiveTrapping(G_top)
%
% PARAMETERS:
%   inp     - Either a valid top surface grid as defined by
%             topSurfaceGrid(G) or a string which is valid input for
%             getAtlasGrid.
%
%   'pn'/pv - List of optional property names/property values:
%                   
%               - coarsening: If inp is the name of a Atlas grid, this
%                 argument will be passed onto getAtlasGrid to coarsen the
%                 surface. Default: 1 for the full grid.
%
%               - light: toggle the lighting of the surface grid.
%                 Default: false
%
%               - spillregions: toggle the display of spill regions on the
%                 surface grid
%                 Default: false
%
%               - traps: toggle the display of traps on the surface grid
%                 Default: true
%
%               - colorpath: toggle red/gray color scheme to mark the traps
%                 encountered along a spill path; same color is applied on
%                 the histogram of trapping volumes.
%                 Default: true
%
%               - method: Choose cell or node based algorithm. Valid
%                 inputs: 'node', 'cell'. Default: 'cell'.
%
%               - injpt: Choose this cell number as injection point at
%                 startup. If zero, no injection point is selected.
%                 Default: zero
%
% RETURNS:
%   h  - Handle to resulting figure object.
%
%
% SEE ALSO:
%   trapAnalysis, showTrappingStructure, showTrapsInteractively

%{
#COPYRIGHT#
%}
   opt = struct('coarsening',    1,     ...
                'light',        false,  ...
                'spillregions', false,  ...
                'traps',        true,   ...
                'method',       'cell', ...
                'colorpath',    true,   ...
                'injpt',        0       ...
                );
   
   global selectedVECell;
   selectedVECell = 1;
   
   if nargout > 0
       varargout{1} = NaN;
   end
   
   if strcmpi(opt.method, 'cell')
       require coarsegrid gridtools
   end
   
   opt = merge_options(opt, varargin{:});

    if isfield(inp, 'cells')
        % Assume input is a top surface grid
        Gt = inp;
        assert(isfield(Gt, 'columns'), 'Please provide a top surface grid');
        top = [];

    elseif ischar(inp)
        % Assume input is the name of a CO2 Storage Atlas grid

        [gr data petroinfo] = getAtlasGrid(inp, 'coarsening', opt.coarsening , 'nz', 1);

        % Grab topgrid
        top = data(cellfun(@(x) strcmpi(x.variant, 'top'), data));
        if numel(top)>0
            top = top{1};
        else
            top = [];
        end
        % Process the grid
        G = processGRDECL(gr{1});

        % The coarsening may disconnect small components, take only the largest and
        % first grid produced by processGRDECL and add geometry data to it.
        try
            mlist = mrstModule();
            mrstModule add mex
            G = mcomputeGeometry(G(1));
            mrstModule('reset', mlist{:})
        catch %#ok
            G = computeGeometry(G(1));
        end

        Gt = topSurfaceGrid(G);
        % Pass it silently along with the grid
        Gt.petroinfo = petroinfo{1};
    end

    res = trapAnalysis(Gt, strcmpi(opt.method, 'cell'));

    if isempty(res.trap_adj) || all(res.trap_adj(:)==0)
        disp('No trap adjacency graph produced. Does the grid have any traps?')
        return
    end
    
    res.Gflat = flattenTraps(Gt, res);
    % Store a list of traps to avoid slicing sparse matrices
    % Commented out because it wasn't needed. May be used for very large
    % grids.
%     [i, j] = find(res.trap_adj);
%     T = sortrows([i,j]);
%     i = T(:,1);
%     j = T(:,2);
%     p = cumsum([1; accumarray(i, 1)]);
%     res.p = p;
% 
%     res.i = i;
%     res.j = j;


    defaultpos = get(0, 'DefaultFigurePosition');
    % Try to be sort of not stupid when creating a large figure by going via
    % the default setting.

    h = figure('Position', [defaultpos(1:2) - [0 100], 1000 500]);
    bf = boundaryFaces(Gt);

    % Create a toolbar for toggling different displays
    ht = uitoolbar(h);
    
    uitoggletool(ht, 'TooltipString', 'Plot traps',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('traps'),...
                     'Tag', 'toggleTraps', ...
                     'State', boolToToggle(opt.traps));

    uitoggletool(ht, 'TooltipString', 'Plot spill regions',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('spill'),...
                     'Tag', 'toggleSpillRegions',...
                     'State', boolToToggle(opt.spillregions));
    
                 
    uitoggletool(ht, 'TooltipString', 'Lighting',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top), ...
                     'cdata', geticon('light'),...
                     'Tag', 'toggleLight', 'State', boolToToggle(opt.light));

    uitoggletool(ht, 'TooltipString', 'Toggle contourlines',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('contour'),...
                     'Tag', 'toggleContour');
                  
    uitoggletool(ht, 'TooltipString', 'Red/gray color along spill paths', ...
                      'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top), ...
                      'cdata', geticon('colorpath'), ...
                      'Tag', 'toggleColorpath', 'State', boolToToggle(opt.colorpath));

    uipushtool(ht, 'TooltipString', 'Reset view',...
                   'cdata', geticon('reset'),...
                   'ClickedCallback', @(src, event) view(0, 90));

    uipushtool(ht, 'TooltipString', 'Setup VE simulation',...
                   'cdata', geticon('simulate'),...
                   'ClickedCallback', @(src, event) runSimulation(Gt, res, src, event));
               
               
               
    ax = plotMain(Gt, res, bf, top);
    if opt.injpt > 0,
       clickHandler(Gt, res, bf, [], [], [], opt.injpt);
    end
    axes(ax); %#ok<MAXES>

    if nargout > 0
       % Return a handle to the figure if it was asked for.
       varargout{1} = h;
    end
end

function clickHandler(Gt, res, bf, data, src, event, flag)                 %#ok
    % Get data
    if ~flag,
       pts = get(gca, 'CurrentPoint');
       [c f] = nearestCellLine(Gt, 1:Gt.faces.num, pts);                      %#ok
       c = c(1);
    else
       c = flag;
    end
    global selectedVECell;
    
    if strcmpi(get(gcf, 'SelectionType'), 'normal');
        selectedVECell  = c;
    end
    
    reptext = findobj('Tag', 'gridReport');
    if ~isempty(reptext)
        set(reptext(1), 'string', getReport(Gt, selectedVECell));
    end
    
    if isfield(res, 'trap_regions')
        % We use spill region to categorize cells
        trap = res.trap_regions(c);
    else
        % Fail if not directly clicking a trap
        trap = res.traps(c);
    end
    if trap == 0
        plotMain(Gt, res, bf, data);
        plotGrid(Gt, c, 'facec', 'yellow')
        disp(['Current positition is not downstream from any trap, '...
              'or spillpoint data was not produced.'])
        return
    end

    if strcmpi(get(gcf, 'SelectionType'), 'alt')
        figure(302); clf;
        
        currTrap = res.traps == trap;
        
        volume = volumesOfTraps(Gt, res, trap);
        title(sprintf('Trap %d, containing %1.3g m^3 volume', trap, volume));
        gc = Gt.cells.centroids(currTrap, :);
        
        m = min(gc);
        M = max(gc);

        % Find the general area around the trap
        currRegion = findRegion(Gt, m, M, .3);

        % Actual plot
        
        % Plot bottom of trap
        plotGrid(res.Gflat, currTrap, 'facec', 'red')
        % Plot ceiling of trap
        plotGrid(Gt, currTrap & currRegion, 'facec', 'blue', 'facea', .8, 'edgec', 'none')
        plotCellData(Gt, Gt.cells.z, ~currTrap & currRegion, 'facea', .3);
        
        axis off tight
        view(-20, 45)
        legend({'Trap floor', 'Trap ceiling'})
        cbar = colorbar('South');
        set(get(cbar, 'YLabel'), 'String', 'Depth')
        return
    end
    
    % Plot several things
    % - The grid with the traps highlighted
    % - Amount trapped in a pie chart
    % - Log-bar-plot of total volume of each trap reached
    if strcmpi(get(gcf, 'SelectionType'), 'extend')
        outline = 'green';
        A = res.trap_adj';
    else
        outline = 'blue';
        A = res.trap_adj;
    end
    subplot('position', [.025 .025 .65 .95]);
    % Store view
    [az, el] = view();
    
    cla;
    %%% Plot path of migration
    subplot('position', [.025 .025 .65 .95]); cla
    axmain = gca;
    plotMain(Gt, res, bf, data);
    fn = @(src, event) clickHandler(Gt, res, bf, data, src, event, 0);
    plotGrid(Gt, res.traps == trap, 'edgec', 'none', 'facec', 'red', 'ButtonDownFcn', fn);
    
    
    subt = getMigrationTree(Gt, A, trap, 0, fn);
    
    t = [trap subt];
    colorizePath = strcmpi(get(findobj(gcf, 'tag', 'toggleColorpath'), 'State'), 'on');
    for i = 1:numel(t)
        ci = res.traps == t(i);
        
        plotGrid(Gt, ci, 'facec', colorizeDepth(i, numel(t), colorizePath),...
           'edgec', 'none', 'ButtonDownFcn', fn);
    end
    partition = res.traps;
    partition(~ismember(partition, t)) = max(t)+1;
    
    Gtmp = Gt;
    Gtmp.nodes.z = Gt.nodes.z - (max(Gt.cells.z) - min(Gt.cells.z))/1000;
    plotPartitionOutlineTopsurface(Gtmp, partition, 'LineWidth', 1.5, 'EdgeColor', outline)
    % Highlight current selection
    plotGrid(Gt, c, 'facec', 'yellow', 'ButtonDownFcn', fn)
    
    view(az, el);
    
    if isfield(res, 'cell_lines')
        [x y z] = getPlotRivers(Gt, [res.cell_lines{[unique(subt) trap]}], res);
    end
    
    if strcmpi(get(findobj(gcf, 'tag', 'toggleLight'), 'State'), 'on')
        lighting phong
    end
    hold on
    plot3(x, y, z, 'k', 'LineWidth', 2);
    
    %%% Plot piechart
    subplot('position', [.75 .525 .225 .425]); cla reset; hold off;
    subtraps = unique(subt);
    subtraps(subtraps == trap) = [];
    
    % Cap above zero since pie refuses to plot zero values
 
    vprimary   = volumesOfTraps(Gt, res, trap);
    vsecondary = volumesOfTraps(Gt, res, subtraps);
    if isempty(subtraps) || sum(vsecondary) < 1;
        % To get nice plot, set zero values to 1 m^3 volume. This is small
        % at the scales we are operating on, and will give a better plot.
        vsecondary = 1;
    end
    
    trapped = vprimary + sum(vsecondary);
    
    total = sum(volumesOfTraps(Gt, res, 1:max(res.traps))) - trapped;
    
    pie([vprimary, sum(vsecondary), total], [1,0,0], {'Primary', 'Migration', 'Not filled'})

    %%%  Plot barchart
    subplot('position', [.775 .025 .2 .45]); cla
    % While setting the axis to logarithmic would be the best idea, this
    % unfortunately breaks opengl support for ALL subplots. We instead
    % manually create a logarithmic axes
    data = [vprimary vsecondary(end:-1:1)];
    ldata = log10(data);
    ldata(ldata < 0) = 0;
    % Transform and plot the data to the unit axis
    % trick: By using hist we get tighter spacing *and* the bars are for
    % some reason drawn as patches.
    hbar = bar(ldata./max(ldata), 'hist');
    N = numel(subtraps) + 1;
    set(hbar, 'FaceVertexCData', colorizeDepth(1:N, N, colorizePath));
    set(gca, 'YTickLabel', sprintf('%2.3g|', 10.^(get(gca, 'YTick')*max(ldata))))
    set(gca, 'YTickMode', 'manual')
    axis tight
    xlabel('Trap graph distance')
    ylabel('Trapped volume (m^3)')
    
    % Set axis to the main plot to ensure that immediate changes to view
    % etc. go to the actual 3D axis.
    axis(axmain);
end

function ax = plotMain(Gt, res, bf, atlasdata)
% Main plotting of grid
    
    isOn = @(tag) strcmpi(get(findobj(gcf, 'tag', tag), 'State'), 'on');
    
    showTraps = isOn('toggleTraps');
    showSpill = isOn('toggleSpillRegions');
    showLight = isOn('toggleLight');  
    showContour = isOn('toggleContour');  
    
    fn = @(src, event) clickHandler(Gt, res, bf, atlasdata, src, event,0);
    subplot('position', [.025 .025 .65 .95]);

    cla; ax = gca;

    zdata = Gt.cells.z;
    zrange = (zdata - min(zdata))./(max(zdata) - min(zdata));
    
    if isfield(res, 'trap_regions') && showSpill
        % Plot spill regions
        map = greyMask(jet(max(res.traps) + 1));
        map(1,:) = get(gcf, 'Color');
        tmp = res.trap_regions;
        tmp(tmp>max(res.traps)) = max(res.traps);
        
        h = plotCellData(Gt, ones(Gt.cells.num, 1), 'ButtonDownFcn', fn, 'EdgeColor', 'none');
        set(h, 'FaceVertexCData', map(tmp + 1, :))
    else
        % Plot heightdata
        
        h = plotCellData(Gt, zdata, 'ButtonDownFcn', fn);
        map = (jet());
        v = interp1(linspace(0,1,64), map, zrange );
        set(h, 'FaceVertexCData', v);
    end
    
    if showTraps
        plotCellData(Gt, res.traps + zrange./5, res.traps ~= 0, 'ButtonDownFcn', fn, 'edgecolor', [.3 .3 .3], 'edgea', .05)
        
        plotPartitionOutlineTopsurface(Gt, res.traps, 'LineWidth', 1.5)
        caxis([0 max(res.traps)]);
    end
    
    colormap jet
    axis tight off

    if ~isempty(atlasdata) && showContour
        contourAtlas(atlasdata, 10, 1)
    end
    
    if showLight
        light;
        lighting phong
        light('Position',[max(Gt.cells.centroids) 4*max(Gt.cells.z)],'Style','infinite');
    end
end

function subtraps = getMigrationTree(G, A, trap, depth, fn)
% Recursively traverse and find the full migration tree
    subtraps = find(A(trap, :));
    tmp = [];
    for i = 1:numel(subtraps);
        trp = getMigrationTree(G, A, subtraps(i), depth + 1, fn);
        tmp = [tmp trp]; %#ok
    end
    subtraps = [reshape(subtraps, 1, []) reshape(tmp, 1, [])];
end

function [c f] = nearestCellLine(G, bf, pts)
% Find a point on the surface based on a click.
    N = numel(bf);
    x1 = repmat(pts(1,:), N, 1);
    x2 = repmat(pts(2,:), N, 1);
    x0 = G.faces.centroids(bf,:);
    if size(x0, 2) == 2
        x0 = [x0 G.faces.z(bf)];
    end
    % Find distance projected on line with a fudge factor for strange grids
    linedist = sum(cross(x0-x1, x0-x2 ,2).^2, 2);
%     ptdist = sum((x0-x1).^2, 2);
%     dist = sqrt(ptdist - 10*linedist);
    
    dist = linedist;
    
    [val ind] = min(dist); %#ok backwards compatability
    f = bf(ind);
    c = G.faces.neighbors(f, :);
    c = c(c~=0);
end


function c = findRegion(G, m, M, radius)
% Find a square region in some radius given minimum and maximum values for
% x and y.
    nc = G.cells.num;
    r = radius*(M-m);
    
    c = true(nc, 1);
    for i = 1:2
        c = c & G.cells.centroids(:,i) < M(i) + r(i);
        c = c & G.cells.centroids(:,i) > m(i) - r(i);
    end
end

function [x y z] = getPlotRivers(Gt, cell_lines, res)
% Plot rivers
    cc=[Gt.cells.centroids, Gt.cells.z];
    [x,y,z] = deal([]);
    
    for i=1:numel(cell_lines);
       cl = unique(cell_lines{i});
       cl = cl(res.traps(cl) == 0);
       x = [x; cc(double(cl),1); NaN]; %#ok
       y = [y; cc(double(cl),2); NaN]; %#ok
       z = [z; cc(double(cl),3); NaN]; %#ok
    end

end

function c = colorizeDepth(depth, N, colorize)
% Colorize from bright red to nearly black
    depth = reshape(depth, [], 1);
    c= max(1-eps - (depth-1)/(N+1), 0);
    if colorize
       c = [.25+.75*c .15*ones(length(c),2)];
    else
       c = .25 + .5*repmat(c, 1, 3);
    end
end

function cdata = geticon(name)
    n = mfilename('fullpath');
    icon = fullfile(n(1:end-length(mfilename)), 'icons', [name '.gif']);
    [cdata, map] = imread(icon);
    map(find(map(:,1)+map(:,2)+map(:,3)==3)) = NaN;     %#ok trivial lookup
    cdata = ind2rgb(cdata,map);
end

function c = greyMask(c)
    c = (c + repmat(get(gcf, 'Color'), size(c, 1), 1))./2;
end

function plotPartitionOutlineTopsurface(Gt, p, varargin)
    % Hack together a sort of 3d-grid to plot z correctly
    Gtmp = Gt;
    Gtmp.nodes.coords = [Gtmp.nodes.coords, Gtmp.nodes.z];
    Gtmp.faces.centroids = [Gtmp.faces.centroids, Gtmp.faces.z];
    Gtmp.cells.centroids = [Gtmp.cells.centroids, Gtmp.cells.z];
    Gtmp.cartDims = [Gtmp.cartDims 1];
    Gtmp.nodes = rmfield(Gtmp.nodes, 'z');
    Gtmp.griddim = 3;
    outlineCoarseGrid(Gtmp, p, 'FaceColor', 'none', varargin{:})
end

function runSimulation(Gt, res, src, event) %#ok<INUSD>
    global selectedVECell veSimAborted
    persistent fi
    
    
    pos = get(gcf, 'OuterPosition');
    if isempty(fi) || ~ishandle(fi) 
        fi = figure('Position',[pos(1:2)-[300 0],  [510 400]], 'Toolbar','none', 'MenuBar', 'none');
    else
        figure(fi);
        clf;
    end
    
    set(fi, 'Name', 'Set up simulation');
    
    bottom = 50;

    ph = uipanel(fi, 'Position', [0 0 1 1], 'Title', 'Set up VE simulation');
    
    % Petrophysical stuff
    poro = .3;
    perm = 300;
    if isfield(Gt, 'petroinfo')
        if ~isnan(Gt.petroinfo.avgperm)
            perm = convertTo(Gt.petroinfo.avgperm, milli*darcy);
        end
        if ~isnan(Gt.petroinfo.avgporo)
            poro = Gt.petroinfo.avgporo;
        end
    end
    sporo = linkedSlider(ph, [10 bottom], 0, 1, poro, 'Porosity');
    
    sperm = linkedSlider(ph, [10 bottom + 1*25], min(.1, perm), max(1000, perm), perm, 'Permeability (mD)');
    
    % Timesteps etc
    sinj  = linkedSlider(ph, [10 bottom + 3*25 + 10], 0, 200, 50, 'Injection time (years)');
    stni  = linkedSlider(ph, [10 bottom + 2*25 + 10], 1, 400, 25, 'Time steps');
    smig  = linkedSlider(ph, [10 bottom + 5*25 + 10], 0, 10000, 2450, 'Migration time (years)');
    stnm  = linkedSlider(ph, [10 bottom + 4*25 + 10], 1, 1000, 245, 'Time steps');
    spres = linkedSlider(ph, [10 bottom + 6*25 + 20], 300, 1000, 300, 'Pressure (bar)');

    % Rough CO2 volume estimate for the reservoir to ensure that the
    % scrollbar shows something (numerically) sensibile
    estco2vol = 500*.3*sum(Gt.cells.volumes.*Gt.cells.H)/1e9;
    
    sco2 = linkedSlider(ph, [10 bottom + 7*25 + 20], 0, estco2vol/50, estco2vol/10000, 'CO2 rate (Mt/year)');
    
    uicontrol(ph, 'Style', 'pushbutton', 'Position', [100, 5, 100, 40], 'string', 'Simulate!', 'callback', @(src, event) uiresume(fi))
    uicontrol(ph, 'Style', 'pushbutton', 'Position', [240, 5, 100, 40], 'string', 'Abort', 'callback', @(src, event) abortSim()) 

    advancedph = uicontrol(ph, 'Style', 'checkbox', 'Position', [350, 5, 120, 40], 'string', 'Advanced plot', 'Value', true);
    
    uicontrol(ph, 'Style', 'text', 'Tag', 'gridReport', 'Position', [10 bottom + 7*25 + 50, 455, 100], 'string', getReport(Gt, selectedVECell))
    
    while ishandle(fi)
        uiwait(fi);

        if ~ishandle(fi)
            % User closed window
            return
        end
        
        petrodata.avgperm = get(sperm, 'Value')*milli*darcy;
        petrodata.avgporo = get(sporo, 'Value');
        
        veSimAborted = false;
        migrateInjection(Gt, res, petrodata, selectedVECell,...
                         'amount',      get(sco2, 'Value'),...
                         'T_injection', get(sinj, 'Value')*year,...
                         'T_migration', get(smig, 'Value')*year,...
                         'topPressure', get(spres, 'Value')*barsa,...
                         'Ni',          round(get(stni, 'Value')),...
                         'Nm',          round(get(stnm, 'Value')), ...-
                         'plotPanel',   get(advancedph, 'Value') == 1);
    end
end

function [s, e] = linkedSlider(fi, pos, md, Md, val, title)
    uicontrol(fi, 'Style', 'text', 'Position', [pos(1:2) 150, 25], 'string', title)
    size1 = [pos(1)+150, pos(2), 250, 25];
    sizee1 = [size1(1:2) + [size1(3) + 10, 0], [50 25]];
    cap = @(x) max(md, min(x, Md));
    
    e = uicontrol(fi, 'Style', 'edit', 'Value', val, 'Position', sizee1, 'String', sprintf('%.2f', val));
    fun = @(src, event) set(e, 'String', sprintf('%.2f', get(src, 'Value')));
    s = uicontrol(fi, 'Style', 'slider', 'Position', size1, 'Min', md, 'Max', Md, 'Value', val, 'Callback', fun);
    fun2 = @(src, event) set(s, 'Value', cap(sscanf(get(src, 'String'), '%f')));
    set(e, 'Callback', fun2);
end

function abortSim()
    global veSimAborted
    veSimAborted = true;
end

function v = boolToToggle(boolvar)
    if boolvar
        v = 'on';
    else
        v = 'off';
    end
end

function s = getReport(Gt, c)
    [i j] = ind2sub(Gt.cartDims, c);
    s = '';
    s = [s sprintf('Reservoir with geometric volume of %2.4g m^3 and %d cells\n', sum(Gt.cells.volumes.*Gt.cells.H), Gt.cells.num)];
    s = [s sprintf('Injection site located at cell %d at \n\tx: %1.2f   y: \t%1.2f      i: %d   j: %d', c, Gt.cells.centroids(c, 1:2), i, j)];
end
