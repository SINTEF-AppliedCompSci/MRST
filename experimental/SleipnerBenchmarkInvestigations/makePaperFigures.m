%% Figure generation for Sleipner Draft Paper.
% Notes:
%   - export_fig is the routine used to save all figures. This routine
%   can be downloaded from MATLAB's file exchange:
%   <http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig>
%   - use png for mesh (smaller file sizes)
%   - use png or pdf for grids
%   - eps figures are better than png using export_fig. pdf is also okay.
%   - use '-transparent' to ensure figure background is white (some figure
%   viewers will show transparent area as checkboard).
%   - need to fix text that contains super- or sub- scripts.

%%
moduleCheck('co2lab');
mrstVerbose on

% directory for all saved figures
figDirName = 'SleipnerFigs';
mkdir(figDirName)


%% Figure 1
% Observed CO2 outlines (plumes) superimposed onto a grid (Gt) to show mis-match
clearvars -except figDirName
plumes = getLayer9CO2plumeOutlines();
[ ~, Gt, ~, ~ ] = makeSleipnerModelGrid('modelName','IEAGHGmodel', 'refineLevel',1, 'plotsOn',false);
[plumes, topsurface, topfit, hCO2] = makeSurfaceDataAndPlots(plumes, Gt, 'plotsOn',true);



%% Figure 2
% CO2 entry rates into layer 9
clearvars -except figDirName
%[opt, varSPE, model, schedule, initState] = studySleipnerBenchmarkFUN('ratecase','SPE', 'runSimulation',false);
%[opt, varOrig, model, schedule, initState] = studySleipnerBenchmarkFUN('ratecase','original', 'runSimulation',false);
[ hfig1, hfig2 ] = plotRatesForComparison;
%set(hfig1,'Color','white');
export_fig(hfig1, [figDirName '/' 'Rate_volumetric'],'-eps','-transparent')
export_fig(hfig2, [figDirName '/' 'Rate_mass'],'-eps','-transparent')

% OR:
%[ schedule, var ] = getInjRatesAndTimeSteps()


%% Figure 3 - 4
% Fig3: GHGT model and IEAGHG model grids
% Fig4: Top surface comparison and re-centered difference
clearvars -except figDirName
[ hfig1, hfig2, hfigA, hfigB, hfigC ] = inspectSleipnerGridModels('refineLevels',1,'add2008Plume',true,'addlegend',false);
export_fig(hfig1, [figDirName '/' 'GHGTgrid'],'-png','-transparent')           % save as png to avoid a large eps file size
export_fig(hfig2, [figDirName '/' 'IEAGHGgrid'],'-png','-transparent')         % save as png to avoid a large eps file size
export_fig(hfigA, [figDirName '/' 'CompareGridsA'],'-eps','-transparent')
export_fig(hfigB, [figDirName '/' 'CompareGridsB'],'-eps','-transparent')
export_fig(hfigC, [figDirName '/' 'CompareGridsC'],'-png','-transparent')      % save this fig as png to avoid strange coloring behavior that occurs when saving as pdf (and when eps gets converted to pdf in latex document)


%% Figure 5
% Structural traps and their capacity in the GHGT and IEAGHG grids
clearvars -except figDirName




%% Figure 6 - 12: Simulation results (CO2 mass)
clearvars -except figDirName
simsDirName = 'SleipnerSims';

% check to see if results directory exists. If it doesn't exist, create it
% and run all simulation cases. Otherwise, load the simulation results.
if ~exist(simsDirName,'dir')
    fprintf('\n Directory of simulations does not exist. \n')
    fprintf('\n Will make directory and run all simulation cases. \n')
    mkdir(simsDirName)
    variousInjectScenarios
    
else disp('Directory exists')
    % load results of simulation cases
    % first get all mat file names in current directory (where results
    % should have been saved)
    s = what(simsDirName);
    files2load = s.mat;
    
end



% Some specifics required for plot generation:
SimStartYear        = 1999;
Years2plot          = [1999.5; 2001.5; 2002.5; 2004.5; 2006.5; 2008.5];
Year2plot_sideView  = 2008.5;
plumeOutlines       = getLayer9CO2plumeOutlines();
ZoomIntoPlume       = true;
plumeOutline_SatTol = (0.01/100); % adjust this value if patch error occurs
                                  % (which can happen when plotting a year
                                  % with minimal massCO2)
                                  

for j=1:numel(files2load)

    % load .mat results file
    load([simsDirName '/' files2load{j}])
    
%     % check rates and time steps
%     % note: when plotting rates using states.wellSol.val, the rates correspond
%     % to the entire time step that was used to get that particular state
%     % result. So plot these rates with a shift in the x-axis to see the time
%     % these rates actually began to take place.
%     offSet = sim_report.ReservoirTime(1)/convertFrom(1,year);
%     figure('name',files2load{j});
%     subplot(2,1,1)
%     hold on; xlabel('Years since simulation start'); ylabel('inj rate (m^3/yr)')
%     for i=1:numel(sim_report.ReservoirTime)
%         plot(sim_report.ReservoirTime(i)/convertFrom(1,year)-offSet, states{i}.wellSol.val, 'xb')
%     end
%     grid
%     subplot(2,1,2)
%     hold on; xlabel('Year'); ylabel('inj rate (m^3/yr)'); title(['Sim Started Jan1st, ',num2str(SimStartYear)])
%     for i=1:numel(sim_report.ReservoirTime)
%         plot(sim_report.ReservoirTime(i)/convertFrom(1,year)+SimStartYear-offSet, states{i}.wellSol.val, 'xb')
%     end
%     bar(var.inj_year+0.5, var.inj_rates, 'EdgeColor','red','FaceColor','none');
%     grid
% 
%     % check injected (accumulated) mass
%     figure('name',files2load{j}); 
%     sp1 = subplot(2,1,1);
%     hold on; xlabel('Years since simulation start'); ylabel('accumulated mass (kg)')
%     sp2 = subplot(2,1,2);
%     hold on; xlabel('Year'); ylabel('accumulated mass (kg)'); title(['Sim Started Jan1st, ',num2str(SimStartYear)])
%     for i=1:numel(sim_report.ReservoirTime)
%         densityCO2  = model.fluid.rhoG(states{i}.pressure);
%         satCO2      = states{i}.s(:,2);
%         totmassCO2(i)  = sum(model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2);
%     end
%     subplot(sp1); plot(sim_report.ReservoirTime/convertFrom(1,year), totmassCO2, 'xb'); grid;
%     subplot(sp2); plot(sim_report.ReservoirTime/convertFrom(1,year)+SimStartYear, totmassCO2, 'xb'); grid;



    % trap analysis of grid
    % TODO: implement so it is only performed once per grid/grid
    % resolution, rather than re-computing for each iteration
    ta = trapAnalysis(model.G, false);
    

    % --------------- Plot of CO2 mass versus plume outlines --------------
    [ ~, ~ ] = subplotCO2simVsCO2obsData_basic(Years2plot, SimStartYear, ...
                plumeOutlines, sim_report, model.G, states, model.fluid, ...
                model, var.wellXcoord, var.wellYcoord, var.wellCoord_x, ...
                var.wellCoord_y, ta, ZoomIntoPlume, plumeOutline_SatTol, ...
                'figname',files2load{j});
    
    % Prepare the plot name for saving:
    % (first remove any '.' from fileName to make it latex-friendly, then
    % remove the part of fileName which contains clock info and file
    % extension, and add a specific plot title)
    plotname = fileName;
    plotname(plotname=='.') = [];
    plotname = [plotname(1:end-18) '_CO2vsObsPlume'];
    
    % Save figure
    %set(gcf,'PaperPositionMode','auto'); %print(plotname,'-dpng','-r0');
    set(gcf,'Color','white');
    export_fig(gcf, [figDirName '/' plotname],'-png')
    % there was an issue exporting this fig using '-transparent', so
    % background color set to white manually
    % ---------------------------------------------------------------------
    
    
%     % for side profile
%     stateYears = sim_report.ReservoirTime/convertFrom(1,year)+SimStartYear;
%     %state2plot = { states{ logical(var.inj_year == Year2plot_sideView) } };
%     state2plot = { states{ logical(stateYears == Year2plot_sideView) } };
%     [ ~ ] = makeSideProfilePlots_CO2heights( model.G.parent, model.G, model.wellmodel.W.cells, ...
%         state2plot, model.fluid, 'SleipnerBounded',true,'legendWithFreeCO2Only',true, ...
%         'figname',files2load{j});
    

    % ------------ Plot of CO2 mass in 2008 with side profiles ------------
    % for combined top view and side profiles, in one figure:
    % specify XY coordinate(s) you want to draw EW-slice(s) through:
    XYcoord1 = [var.wellXcoord, var.wellYcoord + 500];
    %XYcoord2 = [var.wellXcoord, var.wellYcoord + 900];
    %XYcoord3 = [var.wellXcoord, var.wellYcoord + 1200];
    XYcoords = [XYcoord1];%; XYcoord2];% XYcoord3];
    cellIndex = zeros(numel(XYcoords(:,1)),1);
    % get index of closest cell to your specified XY coordinate(s):
    for p = 1:numel(XYcoords(:,1))
        dv              = bsxfun(@minus, model.G.cells.centroids(:,1:2), XYcoords(p,:));
        [v, cinx]       = min(sum(dv.^2, 2));   
        cellIndex(p,:)  = cinx;
    end
    
    inxPts = [model.wellmodel.W.cells; cellIndex];
%     
%     [ ~, ~ ] = subplotCO2simVsCO2obsData_withSideProfiles(2008.5, SimStartYear, plumeOutlines, sim_report, ...
%                 model.G, states, model.fluid, model, ...
%                 var.wellXcoord, var.wellYcoord, var.wellCoord_x, var.wellCoord_y, ta, ...
%                 plumeOutline_SatTol, ...
%                 'sliceCellIndex',inxPts, ...
%                 'SleipnerBounded',true, ...
%                 'figname',files2load{j});
            
    [ ~, ~ ] = subplotCO2simVsCO2obsData_withSideProfiles_basic(2008.5, ...
                SimStartYear, plumeOutlines, sim_report, model.G, states, ...
                model.fluid, model, var.wellXcoord, var.wellYcoord, ...
                var.wellCoord_x, var.wellCoord_y, ta, plumeOutline_SatTol, ...
                'sliceCellIndex',inxPts, 'SleipnerBounded',true, ...
                'figname',files2load{j});
            
    % Prepare the plot name for saving:
    % (first remove any '.' from fileName to make it latex-friendly, then
    % remove the part of fileName which contains clock info and file
    % extension, and add a specific plot title)
    plotname = fileName;
    plotname(plotname=='.') = [];
    plotname = [plotname(1:end-18) '_CO2vsObsPlume_2008_withSideProfiles'];
    
    % Save figure
    %set(gcf,'PaperPositionMode','auto'); %print(plotname,'-dpng','-r0');
    set(gcf,'Color','white');
    export_fig(gcf, [figDirName '/' plotname],'-png')
    % there was an issue exporting this fig using '-transparent', so
    % background color set to white manually
    % ---------------------------------------------------------------------
    
end






%% Figure X-X
% Sensitivities:
clearvars -except figDirName

% 1. run a simulation using smodel, and get the CO2 height data.
% (Alternatively, you could load any previously run simulation results,
% which used smodel.)
[opt, var, smodel, schedule, initState, wellSols, states, sim_report] = ...
    studySleipnerBenchmarkFUN(  'gridname',         'ORIGINALmodel', ...
                                'refineLevel',      1, ... %'num_years',        11, ...
                                'useSensModel',     true, ...
                                'extractSubgrid',   true, ...
                                'modifyPerm',       true, ...
                                'modifyPoro',       true, ...
                                'modifyRhoCO2',     true );
states = addHeightData(states, smodel.G, smodel.fluid);

% Keep only states for the years of plume outline data
summary = zeros(1:numel(var.inj_year),4);
summary(:,1) = var.inj_year;
summary(:,2) = var.inj_rates;
summary(:,3) = sim_report.ReservoirTime;
summary(:,4) = sim_report.ReservoirTime/(60*60*24*365);

% plot CO2 height under top-surface of grid (final year)
% TODO: confirm plotting of correct year with correct outline.
figure; plotCellData(smodel.G, states{end}.h, 'EdgeColor','none')
title({'CO2 height (meters)';['year ',num2str(var.inj_year(1)-1 + sim_report.ReservoirTime(end)/(60*60*24*365))]})
axis equal tight off; colorbar

% 2. assess match between the 'simulated' states.h just obtained and the
% 'observed' states.h (i.e., from another simulation or using observed
% plume height data wrt a given grid)

% a) get plume outlines
plumes_base = getLayer9CO2plumeOutlines();

% b) compute CO2 heights under grid for these plume outlines, and make
% plots to assess how plume outline fits to the grid top-surface. Also plot
% difference between simulated and observed CO2 heights under topsurface
% grid. Plumes_base_comb are the plume polygons combined into the same year
% (i.e., 2006 has two polygons).
[ plumes_base, topsurface, ~, ~ ] = makeSurfaceDataAndPlots(plumes_base, smodel.G, 'plotsOn',true);

[ plumes_base_comb ] = plotDiff( smodel.G, plumes_base, sim_report, states, topsurface );

% c) pass into sensitivity function, to compare 'states.h' with 'plumes_base.h'
dobj_dz = studySleipnerSensitivitiesFUN( initState, smodel, schedule, wellSols, states, 'plumes_base',plumes_base );

% plot dobj_dz
figure; plotCellData(smodel.G, dobj_dz, 'EdgeColor','none');
axis equal tight off; colorbar
ny = numel(plumes_base);
lp1 = line(plumes_base{ny}.outline(:,1), plumes_base{ny}.outline(:,2), ...
    topsurface(plumes_base{ny}.outline)-3, 'LineWidth',3, 'Color','r');
legend([lp1],{['year ',num2str(plumes_base{ny}.year),' outline']})
title('d (obj) / dz')






