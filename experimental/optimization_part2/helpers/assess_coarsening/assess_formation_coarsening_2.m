%% Determine optimal cell size of formations such that the structural
% trapping capacity is within a certain percentage of the un-coarsened
% grid's capacity.

% Formation names to assess:
names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));

ng = numel(names);

% Create color map with distinguishable colors for each grid
assert(exist('colorspace')~=0, 'Ensure colorspace exists and is on path.')
func = @(x) colorspace('RGB->Lab',x);
mymap = distinguishable_colors(ng,'w',func);

% Plotting options:
plotStrap       = true;
plotTrapNum     = false;
combinedPlot    = true;

% Set threshold to determine optimal cell size, (i.e., the strap cap must
% be within this threshold to the un-coarsened grid's strap cap.
threshold       = 0.70;

% Set up figure
if combinedPlot
    k = 0;
    clear names_for_legend
    figure(100); clf; hold on; % for strap cap
    figure(101); clf; hold on; % for num of traps
end

% Loop through each formation and load the pre-computed data (strap cap,
% coarsening level, cell size, number of traps)
for i = 1:numel(names)
    clear res

    file2load = ['strap_cap' names{i} '.mat'];
    try
        % Load data
        load(file2load);
        
        % Prepare data
        cell_size = [res{:,4}]; % will remove empty cells
        strap     = [res{:,2}]; % will remove empty cells
        num_traps = [res{:,6}]; % will remove empty cells
        c_level   = [res{:,3}]; % no empty cells exist...
        % Remove rows of c_level that have empty strap:
        inx2keep = [];
        for r = 1:numel(c_level)
            if ~isempty(res{r,2})
                inx2keep = [inx2keep; r];
            end
        end
        c_level = c_level(inx2keep);
        

        if plotStrap
            if ~combinedPlot
                figure; set(gcf,'name',names{i}); hold on 
            end
            
            % We pick the highest coarsening level (or largest cell_size)
            % that meets the following constraint: the strap capacity must
            % be within 80 percent of the un-coarsened grid's strap
            % capacity
            pts = strap(strap>0);
            cs  = cell_size(strap>0);
            cl  = c_level(strap>0);
            [~,inx] = find((pts./pts(end)) >= threshold, 1);
            
            
            legend_str = [names{i} ' (' num2str(cs(inx)) ')'];
            if ~combinedPlot
                h1 = plot(cell_size(strap>0), strap(strap>0), 'x-', 'Color',mymap(i,:));
                plot(cs(inx), pts(inx), 'o','LineWidth',3,'MarkerSize',12,'Color',mymap(i,:))
                ylim([0 max(strap(strap>0))])
                xlabel('Cell size, meters')
                ylabel('Structural trapping capacity, Gt')
                title(['best cell size: ',num2str(cs(inx))])
                
                legend([h1],legend_str)
            else
                k = k + 1;
                names_for_legend{k} = legend_str;
                names_and_cellsizes{k,1} = names{i};
                names_and_cellsizes{k,2} = cs(inx); % cell size
                names_and_cellsizes{k,3} = cl(inx); % coarsening level
                names_and_cellsizes{k,4} = pts(inx);% strap cap
                
                figure(100); hold on
                h1(k) = plot(cell_size(strap>0), strap(strap>0), 'x-', 'Color',mymap(i,:));
                plot(cs(inx), pts(inx), 'o','LineWidth',3,'MarkerSize',12, 'Color',mymap(i,:)) 
                
            end
            
            fprintf('%s, %5.0f\n', names{i}, cs(inx))
        end
        
        if plotTrapNum
            figure; set(gcf,'name',names{i})
            plot(cell_size(strap>0), num_traps(strap>0), 'x-','Color',mymap(i,:))
            %plot(cs(inx), pts(inx), 'o','LineWidth',3,'MarkerSize',12, 'Color',mymap(i,:)) 
            xlabel('Cell size, meters')
            ylabel('Number of traps')
        end
        
%         % determine best cell size
%         Y               = strap(strap>0)./num_traps(strap>0);
%         [v,inx]         = min(abs(diff(Y))./abs(Y(2:end)));
%         X               = cell_size(strap>0);
%         best_cell_size{i}  = X(inx); 
        
    catch
        
    end

end

if combinedPlot
    xlabel('Cell size, meters')
    ylabel('Structural trapping capacity, Gt')
    legend([h1],names_for_legend)
end


%% Using names_and_cellsizes, create table to compare structural capacity
% of un-coarsened with coarsened grid
rhoCref = 760;

% Get details of un-coarsened grids
for i = 1:size(names_and_cellsizes,1)
    
   % Get formation name:
   fmName   = names_and_cellsizes{i,1};
   assert( strcmpi(fmName, names{i}) )
   
   % Get it's cell size and strap cap:
   [Gt, rock2D] = getFormationTopGrid( fmName, 1 );
   if any(isnan(rock2D.perm))
        rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
        fprintf('\n\nUsing default permeability:\n')
    end
    if any(isnan(rock2D.poro))
        rock2D.poro = 0.25 * ones(Gt.cells.num,1);
        fprintf('Using default porosity:')
    end
    seainfo         = getSeaInfo(fmName, rhoCref);
    fmCapacities{i} = getTrappingInfo(Gt, rock2D, seainfo, 'fmName',fmName);
    full_cell_size{i} = sqrt(mean(Gt.cells.volumes(:)));
 
end

% Print to table
for i = 1:size(names_and_cellsizes,1)
    
    % Print to table:
    x = fmCapacities{i}.breakdown;
    y = full_cell_size{i};
    fprintf('%15s & %8.0f & %8.3f & %8.0f & %8.3f & %8.3f ( %2.0f ) \\\\ \n', ...
        names{i}, ...
        y, x.structural_trapping_capacity, ...
        names_and_cellsizes{i,2}, names_and_cellsizes{i,4}, ...
        x.structural_trapping_capacity - names_and_cellsizes{i,4}, ...
        100 - names_and_cellsizes{i,4}/x.structural_trapping_capacity * 100 );
end






