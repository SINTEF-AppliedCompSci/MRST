%% Loop through each formation and loop through coarsening levels
% to get coarsening level versus structural trapping capacity

% Formations to assess:
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

% Coarsening levels to loop through
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
run assess_formation_coarsening_2.m



