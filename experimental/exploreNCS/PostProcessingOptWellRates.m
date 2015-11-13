%% Post-processing to compare OPTIMIZED rates VS. UN-OPTIMIZED rates

% NB: there appears to be a unit-conversion mis-match between init and
% optim

%% Specify formation name
fmName = 'Stofm';

pathName = ['opt_results/' fmName '/'];
load([pathName,'Gt.mat']);
load([pathName,'init.mat']);
load([pathName,'optim.mat']);
load([pathName,'other.mat']); % largest file


%% Prepare data before passing into plotToolbar()
% set any values below a tolerance to be NaNs, such that they are not
% plotted as color in the figures

% if initial pressure not saved or returned as output, simply re-calculate
% (assuming it was hydrostatic conditions!)
rhoW = other.fluid.rhoW(0);
initPress = Gt.cells.z * norm(gravity) * rhoW;


%% initial (un-optimized) well locations/rates
sts = init.states;
wcinx = [init.schedule.control(1).W.cells];

for i=1:numel(sts)
    
   sts{i}.s( sts{i}.s(:,2) < 0.00001, 2 ) =  nan;
   
   sts{i}.pressDev = sts{i}.pressure - initPress;
   
   % convert pressure deviation to bars
   sts{i}.pressDev = convertTo(sts{i}.pressDev, barsa);
   sts{i}.pressDev( sts{i}.pressDev < 0.5  ) = nan;
   
end

% call to plotting tool
figure; set(gcf,'name',[fmName ' un-optimized'])
plotToolbar(Gt, sts)
plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
mapPlot(gcf, Gt, 'traps', other.traps.traps, 'trapalpha', 0.2, ...
    'rivers', other.traps.cell_lines, 'rivercolor', [1 0 0], ...
    'maplines', 20, 'wellcells', wcinx, ...
    'well_numbering',true);

% no figure needs to be opened first
plotWellSols(init.wellSols)
set(gcf,'name',[fmName ' un-optimized'])


%% optimized well locations/rates
sts = optim.states;
wcinx = [optim.schedule.control(1).W.cells];

for i=1:numel(sts)
    
   sts{i}.s( sts{i}.s(:,2) < 0.01, 2 ) =  nan;
   
   sts{i}.pressDev = sts{i}.pressure - initPress;
   
   % convert pressure deviation to bars
   sts{i}.pressDev = convertTo(sts{i}.pressDev, barsa);
   %sts{i}.pressDev( sts{i}.pressDev < 0.5  ) = nan;
   
end

% call to plotting tool
figure; set(gcf,'name',[fmName ' optimized'])
plotToolbar(Gt, sts)
plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
mapPlot(gcf, Gt, 'traps', other.traps.traps, 'trapalpha', 0.2, ...
    'rivers', other.traps.cell_lines, 'rivercolor', [1 0 0], ...
    'maplines', 20, 'wellcells', wcinx, ...
    'well_numbering',true);

% no figure needs to be opened first
plotWellSols(optim.wellSols)
set(gcf,'name',[fmName ' optimized'])


%% determine how much leaked out of formation

% Init case:
dh = []; % for subscale trapping?
figure; set(gcf,'name',[fmName ' un-optimized'])
plot(1); ax = get(gcf, 'currentaxes');
% NB: {var.initState, states{:}}
reports = makeReports(Gt, {other.initState, init.states{:}}, ...
                         other.rock, other.fluid, init.schedule, ...
                         other.residual, ...
                         other.traps, dh);
% reports contains soln states; could be used for plotting results.
directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
xlabel('Year')
ylabel('Mass (Mt)')
set(gca,'FontSize',14)


% Optimized case:
dh = []; % for subscale trapping?
figure; set(gcf,'name',[fmName ' optimized'])
plot(1); ax = get(gcf, 'currentaxes');
% NB: {var.initState, states{:}}
reports = makeReports(Gt, {other.initState, optim.states{:}}, ...
                         other.rock, other.fluid, optim.schedule, ...
                         other.residual, ...
                         other.traps, dh);
% reports contains soln states; could be used for plotting results.
directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
xlabel('Year')
ylabel('Mass (Mt)')
set(gca,'FontSize',14)
