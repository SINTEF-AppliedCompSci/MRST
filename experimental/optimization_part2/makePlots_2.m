%% to generate results figs for optimized formation rates.

%clear;
%close all; clc;
gravity on;

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

% pre-allocation of stored results:
E       = zeros(numel(names),1);
injYr   = zeros(numel(names),1);

%%% for closed bdrys:
%fprintf('Formation       | Total inj. (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) | cp \n');
%%% for open bdrys:
fprintf('Formation       | Total inj. (Gt) | Leaked (Gt,Perc.) | Cap (Gt) | Seff (Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | MovePlume(Gt) | cp \n');

%fprintf('Formation       | StrapCap(Gt) | RtrapCap(Gt) | DtrapCap(Gt) | Total inj. (Gt) | Leaked (Gt,Perc.) | Perc.FracPress | StrapAchieved (Gt,Perc.) | RtrapAchieved (Gt,Perc.) | DtrapAch (Gt,Perc.) | MovePlume(Gt) \n');
%fprintf('Formation       | Vp (km3) | StrapCap(Gt) (Perc) | RtrapCap(Gt) (Perc) | DtrapCap(Gt) (Perc) | Total (Gt) \n');
for i = 1:numel(names)
    
    
    pathname = ['testing2/opt_results_one_per_trap_highest_pt_Pressure_plim90_OpenBdrys_stricterTol/' names{i} '/InjYrs150_MigYrs2000_DissOn_0_adjustInitRates1'];
    
    try
      load([pathname '/' 'Gt.mat'])
      load([pathname '/' 'init.mat'])
      load([pathname '/' 'optim.mat'])
      load([pathname '/' 'other.mat'])
    
      E(i) = exploreOptWellNCS_postProcess( Gt, init, optim, other, ...
        'plotWellRates', false, ...
        'plotInventory', false, ...
        'plotPressureChanges', false, ...
        'plot_co2_height', false, ...
        'plot_fraction_overburden', false, ...
        'savePlots', false, ...
        'fmName', other.opt.modelname, ...
        'figDirName', pathname, ...
        'warningOn', false, ...
        'outputOn', false);

      %fh = openfig([pathname '/' names{i} '_optDetails.fig'],'visible');
      %fh.Position = [4505 22 568 423];
      %fh.Name = names{i};
      close all;

      inj_steps = cumsum((optim.schedule.step.val(optim.schedule.step.control == 1)));
      injYr(i) = convertTo(inj_steps(end), year);
    catch
    end
    clearvars -except i names E* injYr
    
end

% store results for future plotting
% names, E, injYr

% figure
% bar([E_50, E_150]); legend('50 yr','150 yr')
% bar([E_30, E_40, E_50, E_60, E_70, E_150]); legend('30 yr','40 yr', '50 yr', '60 yr', '70 yr', '150 yr');
% ylabel('Storage efficiency', 'FontSize',16)
% set(gca,'xtick',1:numel(names),'xticklabel',names,'xticklabelrotation',45, 'fontsize',16)
