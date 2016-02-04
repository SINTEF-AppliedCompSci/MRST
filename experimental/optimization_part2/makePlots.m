%% to generate results figs for optimized formation rates.

clear; close all; clc;

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


for i = 1:numel(names)
    
    try  
    fprintf('-----------------FORMATION: %s -----------------\n', names{i})
    pathname = ['opt_results/' names{i} '/MigYrs2000'];
    load([pathname '/' 'Gt.mat'])
    load([pathname '/' 'init.mat'])
    load([pathname '/' 'optim.mat'])
    load([pathname '/' 'other.mat'])
    
    seainfo = getSeaInfo(names{i}, other.opt.refRhoCO2);
    fmCap   = getTrappingInfo(Gt, other.rock, seainfo, 'mapPlotOn',false, 'fmName',names{i});

    strapCap(i) = fmCap.breakdown.structural_trapping_capacity;
    area(i) = convertTo(sum(Gt.cells.volumes),(kilo*meter)^2);
    num_sptr(i) = max( unique(other.traps.trap_regions( other.traps.trap_regions > 0 )) );
    tot_sptr_area(i) = convertTo(sum(Gt.cells.volumes( other.traps.trap_regions > 0 )),(kilo*meter)^2);
    num_wells(i) = numel([optim.schedule.control(1).W.cells]);
    %N_per_cap(i) = numel([optim.schedule.control(1).W.cells]) / fmCap.breakdown.structural_trapping_capacity;
    %N_per_SA(i) = numel([optim.schedule.control(1).W.cells]) / convertTo(sum(Gt.cells.volumes( other.traps.trap_regions > 0 )),(kilo*meter)^2);

    wcinx = [optim.schedule.control(1).W.cells];
    used_trap_regions = other.traps.trap_regions(wcinx);
    used_trap_regions_area = 0;
    used_trap_regions_strapCap = 0;
    for j=1:numel(used_trap_regions)
        
        % area of trap regions used:
        used_trap_regions_area = used_trap_regions_area + ...
            sum(Gt.cells.volumes( other.traps.trap_regions == used_trap_regions(j) ));
        
        % strap capacity of trap regions used:
        tmp = fmCap.trapcap_tot( fmCap.ta.trap_regions == fmCap.ta.trap_regions(wcinx(j)) );
        tmp = unique(tmp(~isnan(tmp))); % kg
        used_trap_regions_strapCap = used_trap_regions_strapCap + tmp;
    end
    tot_sptr_area_used(i) = convertTo(used_trap_regions_area,(kilo*meter)^2);
    strapCap_used(i) = used_trap_regions_strapCap/1e12; % Gt 
    
    % try to compute the used reachable structural capacity of the wells,
    % but don't count trap capacities more than once.
    
    % one approach: compute based on initial schedule well rates
    inj_time_yr = convertTo( sum( init.schedule.step.val(init.schedule.step.control == 1) ), year );
    inj_rates = [init.schedule.control(1).W.val].*other.opt.refRhoCO2.*(1*year); % kg/yr
    inj_masses = inj_rates .* inj_time_yr; % kg
    tot_inj_mass = sum(inj_masses)/1e12; % Gt
    
%     exploreOptWellNCS_postProcess( Gt, init, optim, other, ...
%         'plotWellRates', true, ...
%         'plotPressureChanges', true, ...
%         'savePlots', true, ...
%         'fmName', other.opt.modelname, ...
%         'figDirName', 'OptimizedRates_figs_2000yrsMigration')
    catch
        fprintf('       Incomplete results.\n')
    end
    fprintf('--------------------------------------------------------\n\n')
    close all; clearvars -except i names strapCap area num_sptr tot_sptr_area num_wells N_per_cap N_per_SA
    
end

fprintf('\n Name | Strap Cap (Gt) | Number Spill-trees | Tot spill-tree area (km2) | Number of wells placed (N) | Spill-tree area utilized (km2) | Strap Cap used (Gt) | N/Gt | N/km2 \n');
for i = 1:numel(names)
    
    fprintf(' %16s  &   %6.2f    &   %4.0f   &    %6.2f    &     %4.0f   &    %6.2f    &    %6.4f   \\\\ \n', ...
        names{i}, ...
        strapCap(i), ... %area(i), ...
        num_sptr(i), ...
        tot_sptr_area(i), ...
        num_wells(i), ...
        tot_sptr_area_used(i), ...    
        strapCap_used(i), ...
        num_wells(i) / strapCap_used(i), ...
        num_wells(i) / tot_sptr_area_used(i) );
     
end