
%% Optimize well placement and rates at Utsira
[Gt, optim, init, history, other] = ...
    optimizeFormation('modelname', 'utsirafm', ...
                      'trapfile_name', 'utsira_subtrap_function_3.mat', ...
                      'num_wells', 1, 
                      'leak);

%% Saving optimized results
savedir = 'results/utsira_optimization/';
if ~isdir(savedir)
   mkdir(schedule_savedir);
end
save([savedir 'init'] , 'init');
save([savedir 'optim'], 'optim');
save([savedir 'Gt'], Gt);
save([savedir, 'history'], 'history');

%% Figure 11: Plot map with the numbered wells


%% Figure 12: Compare default and optimized well rates

%% Figure 13: Plot default and optimized trapping distributions

%% Figure 14: Simulation snapshots

%% Figure 15: Overpressure at early injection

