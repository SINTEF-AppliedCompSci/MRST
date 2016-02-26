%% to generate results figs for optimized formation rates.

clear;
close all; clc;

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



%%% for closed bdrys:
fprintf('Formation       | Total inj. (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) | cp \n');
%%% for open bdrys:
%fprintf('Formation       | Total inj. (Gt) | Leaked (Gt,Perc.) | Cap (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) \n');

%fprintf('Formation       | StrapCap(Gt) | RtrapCap(Gt) | DtrapCap(Gt) | Total inj. (Gt) | Leaked (Gt,Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | DtrapAch (Gt,Perc.) | MovePlume(Gt) \n');
%fprintf('Formation       | Vp (km3) | StrapCap(Gt) (Perc) | RtrapCap(Gt) (Perc) | DtrapCap(Gt) (Perc) | Total (Gt) \n');
for i = 1:numel(names)
    
    
    pathname = ['testing/opt_results_one_per_trap_highest_pt_Pressure_plim90_ClosedBdrys/' names{i} '/InjYrs20_MigYrs10_DissOn_0_lowerInitRates1'];
    
    try
      load([pathname '/' 'Gt.mat'])
      load([pathname '/' 'init.mat'])
      load([pathname '/' 'optim.mat'])
      load([pathname '/' 'other.mat'])
    
      exploreOptWellNCS_postProcess( Gt, init, optim, other, ...
        'plotWellRates', true, ...
        'plotInventory', true, ...
        'plotPressureChanges', true, ...
        'savePlots', false, ...
        'fmName', other.opt.modelname, ...
        'figDirName', pathname, ...
        'warningOn', false, ...
        'outputOn', false)
    
      close all;
      openfig([pathname '/' names{i} '_optDetails.fig']);
      close;
    catch
    end
    clearvars -except i names
    
end
