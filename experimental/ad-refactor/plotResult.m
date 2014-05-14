function plotResult(results, timesteps)

    % loading results (presumably from same simulation grid and with same
    % timesteps)
    num_models = numel(results);
    Gt = [];
    states = cell(num_models, 1);
    for i = 1:num_models
        states{i} = load(results{i});
        if isempty(Gt) && isfield(states{i}, 'Gt')
            Gt = states{i}.Gt;
        end
    end
    z = zeros(numel(states{1}.states{1}.h), 1);
    x = [1:numel(z)]';
    if ~isempty(Gt)
        z = Gt.cells.z;
        x = Gt.cells.centroids(:,1);
    end
    
    %% plotting graphs
    figure('Color', [1 1 1]);
    style1 = {'r', 'g', 'b'}
    style2 = {'r--', 'g--', 'b--'};
    tnum = numel(timesteps);
    for t_ix = 1:tnum
        t = timesteps(t_ix);
        
        % Plotting height profiles
        subplot(4, tnum, t_ix);
        heights = [];
        for m = 1:num_models
            heights = [heights, states{m}.states{t}.h + z];
        end
        polygraph([heights, z], {style1{:}, 'k'}, ...
                  {'km','m'}, '', x/1e3, ...
                  [min(z), max(z) + (max(z)-min(z))*0.1]);
        
        set(gca, 'Ydir', 'reverse');
        
        
        % Plotting top and interface pressure
        subplot(4, tnum, tnum + t_ix);
        ipress = [];
        tpress = [];
        for m = 1:num_models
            ipress = [ipress, states{m}.states{t}.pressure];
            tpress = [tpress, states{m}.states{t}.extras.pTop];
        end
        polygraph([ipress, tpress], ...
                  {style1{1:num_models}, style2{1:num_models}}, ...
                  {'km','MPa'}, '', x/1e3, []);
        
        % plotting density
        subplot(4, tnum, 2*tnum + t_ix);
        idensity = []; tdensity = [];
        for m = 1:num_models
            idensity = [idensity, double(states{m}.states{t}.extras.rhoI)];
            tdensity = [tdensity, double(states{m}.states{t}.extras.rhoTop)];
        end
        polygraph([idensity, tdensity], ...
                  {style1{1:num_models}, style2{1:num_models}}, ...
                  {'km', 'kg/m^3'}, '', x/1e3, []);

        % plotting mass flux
        subplot(4, tnum, 3*tnum + t_ix);
        mflux = [];
        for m = 1:num_models
            mflux = [mflux, states{m}.states{t}.extras.fluxCO2];
        end
        polygraph(mflux, style1, {'km','kg/m^2/t'}, '', x(2:end)/1e3, []);
    end
    
    %% Column labeling
    for c = 1:tnum
        subplot(4, tnum, c);
        xlim = get(gca, 'Xlim');
        ylim = get(gca, 'Ylim');
        xpos = xlim(1) + diff(xlim)/2;
        ypos = ylim(1) - diff(ylim)/5;
        h = text(xpos, ypos, sprintf('year %d', timesteps(c)), ...
                 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    
    %% Row labeling
    row_labels = {'Depth', 'Pressure', 'Density', 'Mass flux'};
    for r = 1:4
        subplot(4, tnum, (r-1)* tnum + 1)
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

%{
* TODO
** DONE axis scaling and units
** DONE invert y-axis for height
** graph markers
** DONE main title
** DONE subtitles
** CANCELED two y-axes for tilted case
%}

