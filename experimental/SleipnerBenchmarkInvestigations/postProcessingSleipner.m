%% Sleipner Benchmark Post-Processing:
% note: some variables or functions might have to be loaded
% specifically, i.e., > load(fileName,'sim_report','wellSols','states','model')
% TODO: check size of function handle model, and consider loading (or
% saving) only parts needed

% manually --> load .mat file containing simulation results


% .mat file containing G, Gt, rock, rock2D is loaded
gridfile2load = [name 'SleipnerGlobalCoords_numRef' num2str(opt.refineLevel) '.mat'];
load(gridfile2load)


% Specify which plots to generate:
plotPanelVE                     = false;
plotModelGrid                   = false;
plotInitialPressure             = false;
plotActualVsSimInjectLocation   = false;
plotInjectRateOverTime          = false;
plotBHPvsTime                   = false;
plotAccumCO2vsTime              = false;
plotEndOfSimResults             = false;
plotCO2simVsCO2obsData          = true; ZoomIntoPlume = false;%true; % if false, entire grid is plotted
plotTrappingInventory           = false;
plotTrapProfiles                = false;
plotTrapAndPlumeCompare         = false;
showTableOfTrapDetails          = false;
plotSideProfileCO2heights       = false;
    

% Trapping analysis method:
% true to use cell-based method, false for node
isCellBasedMethod = false;


% FOR PLOTS:
CO2plumeOutline_SatTol  = (0.01/100); % adjust this value if patch error occurs (which happens when no massCO2 present at specified sat tol)
press_deviation = 0;  % from hydrostatic (percent) --> used for trapping capacity calculation, not simulation


% For plotting of CO2 plumes
% bounds of 2008 plume:
ZoomX1 = 0.4375e6;
ZoomY1 = 6.47e6;
ZoomX2 = 0.4395e6;
ZoomY2 = 6.474e6;


% Coordinate that wellCellIndex corresponds to:
wellCoord_x = Gt.cells.centroids(wellCellIndex,1);
wellCoord_y = Gt.cells.centroids(wellCellIndex,2);
wellCoord_z = 0;


% Call to trap analysis, used for a few plotting functions
if ~exist('isCellBasedMethod','var')
    isCellBasedMethod = false;
end
ta = trapAnalysis(Gt, isCellBasedMethod); % true for cell-based method


%% plotPanelVE
% See also: migrateInjection, plotPanelVE
if plotPanelVE

    Years2plot = [1999; 2001; 2002; 2004; 2006; 2008];

    [ hfig, hax ] = makeSideProfilePlots( Years2plot, inj_year, schedule, ...
        G, Gt, sim_report, states, rock2D, fluid, rhoCref, wellCellIndex, ta );

end


    

%% Plots cooresponding to grid and inital set-up:

if plotModelGrid
    [ hfig, hax ] = plot3DandTopGrids( G, Gt );
end

if plotInitialPressure
    figure;
    plotCellData(Gt, initState.pressure, 'EdgeColor','none')
    title('Initial Pressure','fontSize', 18);
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Pascals', 'fontSize', 18 );
    axis off tight equal
end

if plotActualVsSimInjectLocation
    [ hfig, hax ] = plotRealVsDiscreteInjLoc(Gt,bf,wellXcoord,wellYcoord,wellCoord_x,wellCoord_y);
end

if plotInjectRateOverTime
    [ hfig, hax ] = plotInjectRateVsTime(schedule,inj_year,rhoCref);
end
    


%% Plots cooresponding to VE simulation results:

% BHP VS TIME
if plotBHPvsTime
    time = sim_report.ReservoirTime;
    bhp = zeros(numel(wellSols),1);
    for i = 1:numel(wellSols)
        bhp(i) = wellSols{i}.bhp; % bhp is in Pa=N/m^2
    end
    figure;
    plot(time/365/24/60/60,bhp,'x--')
    xlabel('Reservoir time, years'); ylabel('well bhp, Pascals=10^{-5}bars');
end


% ACCUM CO2 VS TIME (compare this plot against Cavanagh 2013, fig 3)
if plotAccumCO2vsTime
    time = sim_report.ReservoirTime;
    accumCO2sat = zeros(numel(states),1);
    accumCO2mass = zeros(numel(states),1);
    for i = 1:numel(states)
        accumCO2sat(i) = sum( states{i}.s(:,2).*model.G.cells.volumes ); % sat vals are in terms of pore volume
        satCO2          = states{i}.s(:,2);
        densityCO2      = fluid.rhoG(states{i}.pressure); 
        accumCO2mass(i) = sum( model.rock.poro .* model.G.cells.volumes .* model.G.cells.H .* satCO2 .* densityCO2 );
    end
    figure;
    plot(time/365/24/60/60,accumCO2mass/1e9,'o-')
    xlabel('Reservoir time, years'); ylabel('Accumlated CO2 mass, Mt (or 10^9 kg)');
end


% END of SIMULATION PROFILES
% use 'final' or the year
if plotEndOfSimResults
    [ hfig, hax ] = plotProfilesAtGivenTime('final', inj_year, Gt, states, initState, fluid, model, sim_report, caprock_temperature);
end



% INVENTORY (from exploreSimulation.m)
if plotTrappingInventory
    dh = []; % for subscale trapping?
    h2 = figure; plot(1); ax = get(h2, 'currentaxes');
    reports = makeReports(model.G, {initState, states{:}}, model.rock, model.fluid, schedule, ...
                             [model.fluid.res_water, model.fluid.res_gas], ...
                             ta, dh);
    % reports contains soln states; could be used for plotting results.
    directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
    
    %ax = gca;
    %ax.XTickLabel = ax.XTick + inj_year(1)-1;
    % use R2014a and earlier releases syntax to ensure backwards compatibility 
    ax  = get(gca, 'XTick');
    axl = arrayfun(@(a) sprintf('%d', a + inj_year(1)), ax, 'UniformOutput', false);
    set(gca, 'XTickLabel', axl)
    xlabel('Year')
    ylabel('Mass (Mt)')
    set(gca,'FontSize',14)
end


%% Line plots of CO2 migrating plume data
% Note: to run the following function, first downloaded the plume .mat
% files from https://bitbucket.org/mrst/mrst-co2lab/downloads, and place on
% current working directory path
plume = getLayer9CO2plumeOutlines();


% PROFILES AT SELECT TIME
if plotCO2simVsCO2obsData
    
    Years2plot = [1999; 2001; 2002; 2004; 2006; 2008];
    %Years2plot = [2002; 2006; 2008];
    
    [ hfig, hax ] = subplotCO2simVsCO2obsData_basic(Years2plot, inj_year, plume, sim_report, ...
            Gt, states, fluid, model, ...
            wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, ta, ...
            ZoomIntoPlume, ZoomX1, ZoomX2, ZoomY1, ZoomY2, ...
            CO2plumeOutline_SatTol);
        % note the basic function plots in kg, not Mt

end



%% Structural trapping plots

if plotTrapProfiles

    ta_volumes = volumesOfTraps(Gt, ta);
    
    
    % To display analysis method used.
    if isCellBasedMethod
        disp('Trap analysis done using cell-based method.')
    elseif ~isCellBasedMethod
        disp('Trap analysis done using node-based method.')
    end
    
    
    % To display refinement level used, if any.
    if ( exist('myresolution','var') && strcmpi(myresolution,'useRefinedGrid') ) || ( exist('useRefinedGrid','var') && useRefinedGrid )
        fprintf('Refinement level %d:\n', refineLevel);
    elseif ( exist('myresolution','var') && strcmpi(myresolution,'none') ) || ( exist('useRefinedGrid','var') && ~useRefinedGrid )
        disp('No refinement of grid performed.')
    end
    
    
    % Other output.
    fprintf('  Num. global traps: %d\n', numel(ta_volumes));
    fprintf('  Total trap volume: %e m3\n', sum(ta_volumes));
    fprintf('  Avg. global trap size: %e m3\n', mean(ta_volumes));

 
    % PLOT TRAPS COLORED BY CO2 MASS STORAGE CAPACITY
    figure; set(gcf,'Position',[1 1 3000 500])
    hfig = gcf;
    
    %
    subplot(1,5,1); hsub1 = gca; hfsub1 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);

    trapcells = ta.traps~=0;
    cellsTrapVol = zeros(Gt.cells.num,1);
    cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));
    plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0, 'EdgeColor','none')

    set(gca,'DataAspect',[1 1 1/100])
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Trap Volume, km^3', 'fontSize', 18 );
    grid; axis tight; set(gca, 'fontSize', 10);


    % GET TRAPPING BREAKDOWN: structural, residual, dissoluion
    % first, compute theoretical capacity (upper bound):
    [ capacityOutput ] = getTrappingCapacities(Gt, rock2D, ta, ...
        rhoCref, water_density, seafloor_temp, seafloor_depth, ...
        temp_gradient, press_deviation, sr, sw, dis_max);

    % Distributed CO2 mass under structural traps: 
    cellsTrapCO2Mass = zeros(Gt.cells.num,1);
    cellsTrapCO2Mass(trapcells) = capacityOutput.strap_mass_co2(trapcells);

    % Cumulative CO2 mass under structural traps:
    trapcaps = accumarray(ta.traps(trapcells), capacityOutput.strap_mass_co2(trapcells));
    trapcap_tot = zeros(Gt.cells.num,1); %ones(size(ta.traps)) * NaN;
    trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));


    %
    subplot(1,5,2); hsub2 = gca; hfsub2 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, cellsTrapCO2Mass/1e9, cellsTrapCO2Mass~=0, 'EdgeColor','none')

    set(gca,'DataAspect',[1 1 1/100])
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Distributed CO2 Mass under Trap, Mt', 'fontSize', 18 );
    grid; axis tight;
    set(gca, 'fontSize', 10); % check for R2014a

    
    %
    subplot(1,5,3); hsub3 = gca; hfsub3 = gcf;
    hold on
    %plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0, 'EdgeColor','none')

    set(gca,'DataAspect',[1 1 1/100])
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Accumulated CO2 Mass under Trap, Mt', 'fontSize', 18 );
    grid; axis tight; set(gca, 'fontSize', 10);



    % PLOT REACHABLE CAPACITY
    trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, 'removeOverlap', false);
    tvols = [trees.value]; %#ok
    int_tr = find(ta.trap_regions); %#ok ixs of cells spilling into interior trap
    [dummy, reindex] = sort([trees.root], 'ascend'); %#ok

    structural_mass_reached = zeros(Gt.cells.num, 1);
    for i = 1:numel(ta.trap_z) % loop over each trap

        % ix of cells spilling directly into this trap
        cix = find(ta.trap_regions == i);

        % cell indices of all cells of this trap, and its upstream traps
        aix = find(ismember(ta.traps, [trees(reindex(i)).traps]));

        % compute total structural trap capacity (in mass terms) of these
        % cells, and store result
        structural_mass_reached(cix) = sum(capacityOutput.strap_mass_co2(aix)); %#ok

    end

    %
    subplot(1,5,4); hsub4 = gca; hfsub4 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, structural_mass_reached/1e3/1e6, 'EdgeColor','none');
    
    set(gca,'DataAspect',[1 1 1/100])
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Reachable structural capacity, Mt', 'fontSize', 18 );
    grid; axis tight; set(gca, 'fontSize', 10);


    % PLOT SPILL PATHS AND TOPOLOGY
    subplot(1,5,5); hsub5 = gca; hfsub5 = gcf;
    hold on
    mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);

    grid; axis equal tight;

% does not work for R2014a/earlier!
%     % For making plotting adjustments to subplots
%     axesHandles = get(gcf,'children');
%     
%     % Add Injection Location In Each Subplot:
%     for i=1:numel(axesHandles)
%         if strcmpi(axesHandles(i).Type,'axes')
%             
%             subplot(axesHandles(i))
%             % actual location
%             plot(wellXcoord,wellYcoord,'o', ...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10)
%             % simulated location
%             plot(wellCoord_x,wellCoord_y,'x', ...
%                 'LineWidth',3,  ...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','k',...
%                 'MarkerSize',10)
%         end
%     end
    
    hfig = gcf;
    hax  = gca;

end


%% Show the structural trap:
% After the simulation has completed, we are interested in seeing how the
% location of the CO2 plume after a long migration time corresponds to the
% trapping estimates produced by trapAnalysis. This is done by finding the
% trap index of the well injection cell and then plotting the trap along
% with the final CO2 plume.

% Plot the areas with any significant CO2 height along with the trap in red
if plotTrapAndPlumeCompare
 
    % Well in 2D model
    WVE = reports(end).W; % take well of last control since first control might be initial conditions (no well)

    % Generate traps and find the trap corresponding to the well cells
    trap = ta.trap_regions([WVE.cells]);
    
    figure; set(gcf,'Position',[1 1 1600 1000])
    plotCellData(Gt, reports(end).sol.h, reports(end).sol.h > 0.01)
    plotGrid(Gt, ta.traps == trap, 'FaceColor', 'red', 'EdgeColor', 'w')
    plotGrid(Gt, 'FaceColor', 'None', 'EdgeAlpha', .1);

    legend({'CO2 Plume', 'Trap'})
    set(gca,'DataAspect',[1 1 1/10])
    axis tight off
    view(20, 25)
    title('End of simulation CO2 compared to algorithmically determined trap')
    
    % Create textarrow
    annotation(gcf,'textarrow',[0.382421875 0.418396875000001],...
    [0.77 0.857197640117995],'String',{'North'});

end

%% Basic capacity estimates and Show table of Structural trapping details
if showTableOfTrapDetails

    if ~exist('mycase','var')
        if useIEAGHG_model
            mycase = 'IEAGHG';
        elseif useOriginal_model
            mycase = 'GHGT';
        end
    end
    if ~exist('myresolution','var')
        if useRefinedGrid
            myresolution = 'useRefinedGrid';
        else
            myresolution = 'none';
        end
    end

       fprintf('------------------------------------------------\n');
       fprintf('Processing case: %s , %s (numRef=%d) ....\n', mycase, myresolution, refineLevel);

       tan     = trapAnalysis(Gt, false);
       tac     = trapAnalysis(Gt, true);

       %tan_volumes = volumesOfTraps(Gt, tan);
       %tac_volumes = volumesOfTraps(Gt, tac);

       i = 1;
       res{i}.name      = mycase;
       if strcmpi(myresolution,'useRefinedGrid')
           res{i}.refLevel  = refineLevel;
       else
           res{i}.refLevel  = 0;
       end
       res{i}.cells     = Gt.cells.num;
       res{i}.zmin      = min(Gt.cells.z);
       res{i}.zmax      = max(Gt.cells.z);
       res{i}.volume    = sum(G.cells.volumes);
       res{i}.surfarea  = sum(Gt.cells.volumes);
       res{i}.ctrapvols = volumesOfTraps(Gt,tac);
       res{i}.ccapacity = sum(res{i}.ctrapvols);
       res{i}.ntrapvols = volumesOfTraps(Gt,tan);
       res{i}.ncapacity = sum(res{i}.ntrapvols);
       fprintf('done\n');

    % create table:
       fprintf('\n\n------------------------------------------------\n');
       fprintf('%-20s& Refined & Cells  & Min  & Max  & Volume   & Capacity  & Percent &  Capacity & Percent\\\\\n', 'Name');

       fprintf('%-20s&   %2d    & %6d & %4.0f & %4.0f & %4.2e & %4.2e  & %5.2f   & %4.2e  & %5.2f \\\\\n',...
          res{i}.name, res{i}.refLevel, res{i}.cells, res{i}.zmin, res{i}.zmax, res{i}.volume, ...
          res{i}.ncapacity, res{i}.ncapacity/res{i}.volume*100, ...
          res{i}.ccapacity, res{i}.ccapacity/res{i}.volume*100);
       fprintf('------------------------------------------------\n');


      fprintf('\n\n---------------Node-based------------------------\n');
      fprintf('%-20s& Refined & Num. global traps & Tot. trap vol. (m3) & Avg. global trap vol. (m3)\\\\\n', 'Name');

      fprintf('%-20s&   %2d    &     %6d        &     %d    &    %d           \\\\\n',...
          res{i}.name, res{i}.refLevel, ...
          numel(res{i}.ntrapvols), ...
          sum(res{i}.ntrapvols), ...
          mean(res{i}.ntrapvols) );
      fprintf('------------------------------------------------\n');

        fprintf('\n\n---------------Cell-based------------------------\n');
      fprintf('%-20s& Refined & Num. global traps & Tot. trap vol. (m3) & Avg. global trap vol. (m3)\\\\\n', 'Name');

      fprintf('%-20s&   %2d    &     %6d        &     %d    &    %d           \\\\\n',...
          res{i}.name, res{i}.refLevel, ...
          numel(res{i}.ctrapvols), ...
          sum(res{i}.ctrapvols), ...
          mean(res{i}.ctrapvols) );
      fprintf('------------------------------------------------\n');

end

%% Side Vertical Profiles through specified cell, i.e., well cell index
% If all states of simulation are passed in, only last state is plotted. If
% the state of a particular year is to be plotted, pass in that state only.

% sim_report.ReservoirTime contains time (in seconds since start of
% simulation) corresponding to the states{}. To find Year (i.e., 2004)
% corresponding to a state:
%I = 5;
%YearOfStateI = inj_year(1) + sim_report.ReservoirTime(I)/(60*60*24*365) - 1;
%state2plot = { states{I} };

% inj_year (the year) correspondes to states
if plotSideProfileCO2heights

    % To plot a specific injection year:
    YearOfStateToPlot = 2008;
    fprintf('\n Plotting year %d \n', YearOfStateToPlot);
    state2plot = { states{ logical(inj_year==YearOfStateToPlot) } };
    [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, state2plot, fluid, 'SleipnerBounded',true,'legendWithFreeCO2Only',true);
    
    % To plot a specific state (could be a migration year):
    %I = 5;
    %YearOfStateI = inj_year(1) + sim_report.ReservoirTime(I)/(60*60*24*365) - 1;
    %state2plot = { states{I} };
    %fprintf('\n Plotting year %d \n', YearOfStateI);
    %[ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, { states{end} }, fluid, 'SleipnerBounded',true);
    
    % To plot the final state:
    YearOfFinalState = inj_year(1) + sim_report.ReservoirTime(end)/(60*60*24*365) - 1;
    fprintf('\n Plotting final state, year %d \n', YearOfFinalState);
    [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, { states{end} }, fluid, 'SleipnerBounded',true);
    
    % To see all years profiles:
%     for i = 1:numel(states)
%         [ hfig ] = makeSideProfilePlots_CO2heights( G, Gt, wellCellIndex, { states{i} }, fluid, 'SleipnerBounded',true);
%     end
end


% end of post-processing