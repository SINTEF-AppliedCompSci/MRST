function compareTwoSims( sim1str, Gt1, rock2D1, wellSols1, states1, sim_report1, opt1, var1, ...
    sim2str, Gt2, rock2D2, wellSols2, states2, sim_report2, opt2, var2)
% Compare results between two simulations.

% Results are plotted using plotToolbar(Gt, states),
% plotWellSols(wellSols), and directPlotTrappingDistribution()
%
% Implementation could be extended to handle more than two simulation
% results for comparison.

moduleCheck('mrst-gui')

defaultFigPlacement = false;

%% Prepare data before passing into plotToolbar()
% set any values below a tolerance to be NaNs, such that they are not
% plotted as color in the figures



%% sim 1 -
Gt          = Gt1;
initState   = var1.initState;
schedule    = var1.schedule;
sts         = {initState, states1{:}};
ta          = trapAnalysis(Gt,false);
wellSols    = wellSols1;
rock2D      = rock2D1;
fluid       = var1.fluid;
residual    = [var1.fluid.res_water, var1.fluid.res_gas]; % or [opt1.sw, opt1.sr];

wcinx       = [schedule.control(1).W.cells];
maxPressDev = 0;
maxCO2sat   = 0;
timeYr      = [0; convertTo( cumsum(schedule.step.val), year)];

for i=1:numel(sts)
    
   sts{i}.s( sts{i}.s(:,2) < 0.01, 2 ) =  nan;
   
   sts{i}.pressDev = sts{i}.pressure - initState.pressure;
   
   % convert pressure deviation to bars or MPa
   sts{i}.pressDev = convertTo(sts{i}.pressDev * Pascal, mega*Pascal);
   %sts{i}.pressDev( sts{i}.pressDev < 0.5  ) = nan;
   
   % find max pressure dev and when it occurred
   maxPressDev_curr = max(sts{i}.pressDev);
   if maxPressDev_curr > maxPressDev
       maxPressDev          = maxPressDev_curr;
       maxPressDev_timeYr   = timeYr(i);
       maxPressDev_stepNum  = i;
   end
   
   % max co2 saturation and when it occurred
   maxCO2sat_curr = max(sts{i}.s(:,2));
   if maxCO2sat_curr > maxCO2sat
       maxCO2sat            = maxCO2sat_curr;
       maxCO2sat_timeYr     = timeYr(i);
       maxCO2sat_stepNum    = i;
   end
end

% call to plotting tool
figure; set(gcf,'name',sim1str);
if ~defaultFigPlacement; set(gcf,'Position',[3275 421 651 907]); end
plotToolbar(Gt, sts)
plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
    'maplines', 20, 'wellcells', wcinx, ...
    'well_numbering',true);

% no figure needs to be opened first
plotWellSols(wellSols)
set(gcf,'name',sim1str)
if ~defaultFigPlacement; set(gcf,'Position',[2571 936 698 420]); end


fprintf('\n Max press deviation of %d MPa occurred at %d years (step %d).\n', maxPressDev, maxPressDev_timeYr, maxPressDev_stepNum )
fprintf('\n Max CO2 saturation of %d occurred at %d years (step %d).\n', maxCO2sat, maxCO2sat_timeYr, maxCO2sat_stepNum )


% leakage report
dh = []; % for subscale trapping?
figure; set(gcf,'name',sim1str)
if ~defaultFigPlacement; set(gcf,'Position',[2571 422 695 420]); end
plot(1); ax = get(gcf, 'currentaxes');
% NB: {var.initState, states{:}}
reports1 = makeReports(Gt, {initState, states1{:}}, ...
                         rock2D, fluid, schedule, ...
                         residual, ...
                         ta, dh);
% reports contains soln states; could be used for plotting results.
directPlotTrappingDistribution(ax, reports1, 'legend_location', 'northwest');
xlabel('Year')
ylabel('Mass (Mt)')
set(gca,'FontSize',14)

plotFormationPressureChanges( states1, var1, opt1, 'figname', sim1str )

%compareBHPandWellCellPressure( states1, wellSols1, var1.wellCellIndex )


%% sim 2 -
% execute if second sim results were passed in
if exist('Gt2','var')
    
Gt          = Gt2;
initState   = var2.initState;
schedule    = var2.schedule;
sts         = {initState, states2{:}};
ta          = trapAnalysis(Gt,false);
wellSols    = wellSols2;
rock2D      = rock2D2;
fluid       = var2.fluid;
residual    = [var2.fluid.res_water, var2.fluid.res_gas];

wcinx       = [schedule.control(1).W.cells];
maxPressDev = 0;
maxCO2sat   = 0;
timeYr      = [0; convertTo( cumsum(schedule.step.val), year)];

for i=1:numel(sts)
    
   sts{i}.s( sts{i}.s(:,2) < 0.01, 2 ) =  nan;
   
   sts{i}.pressDev = sts{i}.pressure - initState.pressure;
   
   % convert pressure deviation to bars or MPa
   sts{i}.pressDev = convertTo(sts{i}.pressDev * Pascal, mega*Pascal);
   %sts{i}.pressDev( sts{i}.pressDev < 0.5  ) = nan;
   
   % find max pressure dev and when it occurred
   maxPressDev_curr = max(sts{i}.pressDev);
   if maxPressDev_curr > maxPressDev
       maxPressDev          = maxPressDev_curr;
       maxPressDev_timeYr   = timeYr(i);
       maxPressDev_stepNum  = i;
   end
   
   % max co2 saturation and when it occurred
   maxCO2sat_curr = max(sts{i}.s(:,2));
   if maxCO2sat_curr > maxCO2sat
       maxCO2sat            = maxCO2sat_curr;
       maxCO2sat_timeYr     = timeYr(i);
       maxCO2sat_stepNum    = i;
   end
end

% call to plotting tool
figure; set(gcf,'name',sim2str)
if ~defaultFigPlacement; set(gcf,'Position',[4473 422 647 907]); end
plotToolbar(Gt, sts)
plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
mapPlot(gcf, Gt, 'traps', ta.traps, 'trapalpha', 0.2, ...
    'rivers', ta.cell_lines, 'rivercolor', [1 0 0], ...
    'maplines', 20, 'wellcells', wcinx, ...
    'well_numbering',true);

% no figure needs to be opened first
plotWellSols(wellSols)
set(gcf,'name',sim2str,'Position',[3776 933 695 419])


fprintf('\n Max press deviation of %d MPa occurred at %d years (step %d).\n', maxPressDev, maxPressDev_timeYr, maxPressDev_stepNum )
fprintf('\n Max CO2 saturation of %d occurred at %d years (step %d).\n', maxCO2sat, maxCO2sat_timeYr, maxCO2sat_stepNum )


% leakage report
dh = []; % for subscale trapping?
figure; set(gcf,'name',sim2str)
if ~defaultFigPlacement; set(gcf,'Position',[3902 421 560 420]); end
plot(1); ax = get(gcf, 'currentaxes');
% NB: {var.initState, states{:}}
reports2 = makeReports(Gt, {initState, states2{:}}, ...
                         rock2D, fluid, schedule, ...
                         residual, ...
                         ta, dh);
% reports contains soln states; could be used for plotting results.
directPlotTrappingDistribution(ax, reports2, 'legend_location', 'northwest');
xlabel('Year')
ylabel('Mass (Mt)')
set(gca,'FontSize',14)

plotFormationPressureChanges( states2, var2, opt2, 'figname', sim2str )

%compareBHPandWellCellPressure( states2, wellSols2, var2.wellCellIndex )

end


end

function compareBHPandWellCellPressure(states, wellSols, wellCellIndex)
% debug to handle array of wellCellIndex for multiple wells.

    for i=1:numel(states)
       p_bhp(i) = wellSols{1,i}.bhp; % @@ takes well 1
       p_wcl(i) = states{i}.pressure(wellCellIndex(1));
    end
    
    figure
    hold on
    plot(p_bhp,'xr')
    plot(p_wcl,'xb')
    legend(['bottom hole pressure of well cell index ',num2str(wellCellIndex(1))'],'cell pressure') % @@
    ylabel('pressure (Pascals)')

end

