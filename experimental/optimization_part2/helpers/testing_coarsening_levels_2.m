
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


c_levels = [20:-1:1];
rhoCref  = 760 * kilo*gram/meter;

clear res
res = cell(numel(c_levels),2);

for i = 1:numel(names)
    
    fprintf('--- Formation: %s -------\n', names{i});

    for j = 1:numel(c_levels)

        N = c_levels(j);
        
        try
            strap_cap = [];
            cell_size = [];
            tot_area = [];
            
            [Gt, rock2D] = getFormationTopGrid(names{i},N);

            cell_size = sqrt(Gt.cells.volumes(1)); % m
            tot_area = sum(Gt.cells.volumes(:));   % m2

            if any(isnan(rock2D.perm))
                rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
                fprintf('\n\nUsing default permeability:\n')
            end
            if any(isnan(rock2D.poro))
                rock2D.poro = 0.25 * ones(Gt.cells.num,1);
                fprintf('Using default porosity:')
            end
            seainfo      = getSeaInfo(names{i}, rhoCref);
            fmCapacities = getTrappingInfo(Gt, rock2D, seainfo, 'fmName',names{i});
            strap_cap    = fmCapacities.breakdown.structural_trapping_capacity;
            num_traps    = numel(unique(fmCapacities.ta.trap_regions(fmCapacities.ta.trap_regions>0)));
        catch
         
        end


        % store results:
        res{j,1} = names{i};
        res{j,2} = strap_cap;
        res{j,3} = N;
        res{j,4} = cell_size;
        res{j,5} = tot_area;
        res{j,6} = num_traps;

    end

    % save results to .mat file
    save(['strap_cap' names{i}],'res')

end


%% Plotting of results:

figure;
hold on
k = 0;
for i = 1:numel(names)
    clear res

    file2load = ['strap_cap_' names{i} '.mat'];
    file2load_b = ['strap_cap_b_' names{i} '.mat'];
    try
        load(file2load);   res_a = res;
        load(file2load_b); res_b = res;
        
        cell_size   = [res_b{:,4}, res_a{:,4}]; % empty cells are removed
        strap       = [res_b{:,2}, res_a{:,2}]; % empty cells are removed
        num_traps   = [res_b{:,6}, res_a{:,6}];
        
        % chop off data for cell_size > 4000 meters
        strap     = strap(cell_size<=4000);
        num_traps = num_traps(cell_size<=4000);
        cell_size = cell_size(cell_size<=4000);
        
        plot(cell_size, strap./num_traps, 'x-')
        
        % determine best cell size
        Y               = strap(strap>0)./num_traps(strap>0);
        [v,inx]         = min(abs(diff(Y))./abs(Y(2:end)));
        X               = cell_size(strap>0);
        best_cell_size  = X(inx);
        k = k + 1;
        names_for_legend{k} = [names{i} ' (',num2str(best_cell_size),')'];
    catch
        
    end

end

legend(names_for_legend)
xlabel('Cell size, meters')
ylabel('Structural trapping capacity / num of traps, Gt per trap')

xlim([0 12000])


