function plotResult(results, timesteps, includes)

    if ~exist('includes') || isempty(includes)
        includes = [1 1 1 1];
    end
    
    % loading results (presumably from same simulation grid and with same
    % timesteps)
    num_models = numel(results);
    Gt = [];
    states = cell(num_models, 1);
    for i = 1:num_models
        states{i} = load(results{i});
        if isempty(Gt) 
            if isfield(states{i}, 'Gt')
                Gt = states{i}.Gt;
            elseif isfield(states{i}, 'G')
                Gt = topSurfaceGrid(states{i}.G);
            end
        end
    end
    assert(~isempty(Gt))
    yrow = ceil(Gt.cartDims(2)/2);
    ixs = sub2ind(Gt.cartDims, 1:Gt.cartDims(1), yrow*ones(1, Gt.cartDims(1)));
       
    z = Gt.cells.z(ixs);
    x = Gt.cells.centroids(ixs,1);
    
    %% plotting graphs
    figure('Color', [1 1 1]);
    style1 = {'r', 'g', 'b', 'k'}
    style2 = {'r--', 'g--', 'b--', 'k--'};
    tnum = numel(timesteps);
    rows = sum(includes);
    for t_ix = 1:tnum
        t = timesteps(t_ix);
        
        % Plotting height profiles
        if includes(1)
            subplot(rows, tnum, t_ix);
            heights = [];
            for m = 1:num_models
                heights = [heights, states{m}.result(t).h(ixs) + z];
            end
            in_legend = [repmat(true, num_models, 1); false];
            polygraph([heights, z], {style1{1:num_models}, 'k'}, ...
                      {'km','m'}, '', x/1e3, ...
                      [min(z), max(z) + (max(z)-min(z))*0.1], in_legend);
            
            set(gca, 'Ydir', 'reverse');
            set(gca, 'FontSize', 15);
        end
        
        if includes(2)
            shift = tnum * sum(includes(1));
            % Plotting top and interface pressure
            subplot(rows, tnum, shift + t_ix);
            ipress = [];
            tpress = [];
            for m = 1:num_models
                ipress = [ipress, states{m}.result(t).pI(ixs)  ];
                tpress = [tpress, states{m}.result(t).pTop(ixs)];
            end
            in_legend = [repmat(true,  num_models, 1); ...
                         repmat(false, num_models, 1)];
            polygraph([ipress, tpress], ...
                      {style1{1:num_models}, style2{1:num_models}}, ...
                      {'km','MPa'}, '', x/1e3, [], in_legend);
            set(gca, 'FontSize', 15);
            set(get(gca, 'XLabel'), 'FontSize', 15);
        end
        
        if includes(3)
            shift = tnum * sum(includes(1:2));
            % plotting density
            subplot(rows, tnum, shift + t_ix);
            idensity = []; tdensity = [];
            for m = 1:num_models
                idensity = [idensity, double(states{m}.result(t).rhoI(ixs))];
                tdensity = [tdensity, double(states{m}.result(t).rhoTop(ixs))];
            end
            in_legend = [repmat(true,  num_models, 1); ...
                         repmat(false, num_models, 1)];
            polygraph([idensity, tdensity], ...
                      {style1{1:num_models}, style2{1:num_models}}, ...
                      {'km', 'kg/m^3'}, '', x/1e3, [], in_legend);
            set(gca, 'FontSize', 15);
        end

        if includes(4)
            shift = tnum * sum(includes(1:3));
            % plotting mass flux
            subplot(rows, tnum, shift + t_ix);
            mflux = [];
            for m = 1:num_models
                mflux = [mflux, states{m}.result(t).fluxCO2(ixs)];
            end
            polygraph(mflux, style1, {'km','kg/m^2/t'}, '', x(2:end)/1e3, ...
                      []);
            set(gca, 'FontSize', 15);
        end
    end
    
    %% Column labeling
    for c = 1:tnum
        subplot(rows, tnum, c);
        xlim = get(gca, 'Xlim');
        ylim = get(gca, 'Ylim');
        xpos = xlim(1) + diff(xlim)/2;
        ypos = ylim(1) - diff(ylim)/5;
        h = text(xpos, ypos, sprintf('year %d', timesteps(c)), ...
                 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    
    %% Row labeling
    row_labels = {'Depth', 'Pressure', 'Density', 'Mass flux'};
    for r = 1:rows
        if includes(r)
            subplot(rows, tnum, sum(includes(1:(r-1))) * tnum + 1);
            xlim = get(gca, 'Xlim');
            ylim = get(gca, 'Ylim');
            xpos = xlim(1) - diff(xlim)/3;
            ypos = ylim(1) + diff(ylim)/2;
            h = text(xpos, ypos, row_labels{r}, ...
                     'FontSize', 20, ...
                     'HorizontalAlignment', 'center', ...
                     'rotation', 90);
        end
    end
end

%{
* TODO
** DONE axis scaling and units
** DONE invert y-axis for height
** graph markers
** DONE main title
** DONE subtitles
** CANCELED two y-axes for tilted case
%}

