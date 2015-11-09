% Explore storage capacity and injection scenearios for formations found in
% Norwegian Sea

% DESCRIPTION:
%   The Bat Group is made up of Are, Tilje, Ror (sealing) and Tofte
%   (missing data?) formations. Are and Tilje are treated as one aquifer at
%   the regional scale due to lack of regional sealing shales (Compiled
%   Atlas, Chp 5, pg 93). Modeling in Are formation done (pg 100-101).
%
%   The Fangst Group is made up of Ile, Not (sealing), and Garn formations.
%   In NPD's modeling, Ile and Garn were treated as one aquifer (Compiled
%   Atlas, Chp 5, pg 94). An estimated 400 million tonnes CO2 can be
%   stored in Garn/Ile (8 million per year, 50 years, using 4 wells) (pg
%   96).
%
%   Trondelag Platform contains a portion of the Fangst Group (Ile, Not,
%   Garn). It is one of the main structural elements of this sea, and is
%   considered to have the best storage potential (pg. 78)
%
%   Froan Basin is a sub-element of the Trondelag Platform (pg. 96), and is
%   located in the southern part of the Garn/Ile aquifer (pg. 105). CO2
%   migration modeling performed (pg. 97 in Atlas, and in Riis and Halland,
%   2014 pg. 5262). **The simulation grid in the Froan Basin (FBS) is
%   depicted in the map on pg. 95.**
%
%   See map on pg. 93 to see relative location of Froan and Trondelag.
%
%   Prospects have been identified by NPD (pg. 104), and storage estimates
%   tabulated.
%
% *** NB: Need to get heterogeneous rock properties from Andreas at NPD.***

moduleCheck('co2lab', 'opm_gridprocessing');
mrstVerbose on
gravity on

% directory for saving figures
figDirName = 'NorwegianSeaFigs';
mkdir(figDirName)


%% Load formation grids and rock properties
% NB: un-coarsened grids take about 1 minute each to process. Thus grids
% are saved to .mat file to avoid re-generation.
coarsening  = 5;
fileName    = ['NorwegianSeaGrids_ref', num2str(coarsening),'.mat'];
try
    fprintf('\n Loading grids if file exists.\n')
    load(fileName);
    
catch
    fprintf('\nNo such file exists. Processing grids to be saved.\n')
    fmNames     = {'Arefm';'Tiljefm';'Rorfm';'Ilefm';'Notfm';'Garnfm'};
    ng          = numel(fmNames);

    [Gts, rock2Ds, tans] = deal(cell(ng,1));

    for i = 1:ng
        [Gts{i}, rock2Ds{i}] = getFormationTopGrid(fmNames{i}, coarsening);
        tans{i}              = trapAnalysis(Gts{i}, false);
    end

    save(fileName,'fmNames','Gts','rock2Ds','tans','-v7.3');
    
end

%% Try to find any trap(s) that match NPD's identified prospects
%% ------------------------------------------------------------------------
% 1. Assess capacity with GUI:
% NB: keep coarse resolution as these formations are large
% NB: need to confirm sea floor data @@
info = getSeaInfo('NorwegianSea');
exploreCapacity('default_formation',   'Garnfm',             ...
                'grid_coarsening',     5,                    ...
                'seafloor_depth',      info.seafloor_depth,  ...
                'seafloor_temp',       info.seafloor_temp,   ...
                'temp_gradient',       info.temp_gradient     );
            
%% ------------------------------------------------------------------------
% 2. Launch interactiveTrapping:
% NB: use formation name in order to load avg rock data via getAtlasGrid()
N    = 5;
name = 'Garnfm';
interactiveTrapping(name, 'coarsening',N) 
% then use icons in figure to toggle various plotting options.

% re-size figure for clarity
set(gcf,'Position',[1 1 822 839])
view(2)
hfig = gcf;
set(findobj(hfig.Children,'Type','axes'),'FontSize',14)
set(findobj(hfig.Children,'Type','Legend'),'FontSize',14,'Orientation','vertical')
haxs = hfig.Children; % pie axes should be haxs(5)?
set(findobj(haxs(4).Children,'Type','Text'),'FontSize',14) % pie chart might be haxes(4)

export_fig(gcf,[figDirName '/' 'interactiveTrapping' name '_ref',num2str(N)], '-png','-transparent')
close gcf

%% ------------------------------------------------------------------------
% 3. Using results from trapAnalysis() and getTrappingCapacities():
garnInd = find(strcmpi(fmNames,'Garnfm'));
ileInd = find(strcmpi(fmNames,'Ilefm'));
tiljeInd = find(strcmpi(fmNames,'Tiljefm'));
areInd = find(strcmpi(fmNames,'Arefm'));

fmind = [garnInd, ileInd, tiljeInd, areInd]; % no trapAnalysis required for sealing formations
for i = 1:numel(fmind)
    
    Gt      = Gts{fmind(i)};
    ta      = tans{fmind(i)};
    rock    = rock2Ds{fmind(i)};

    % useful plots:
    % - caprock topology (structural traps, spill paths, contours)
    % - structural trap capacity (Mt)
    % - caprock CO2 density (gas vs liquid)
    fprintf('  Formation name: %s\n', fmNames{fmind(i)});
    [ capOutput, hfig, hax ] = getTrappingPlots(Gt, ta, rock, 'NorwegianSea');
    set(hfig,'name', fmNames{fmind(i)});
    
    % save fig
    export_fig(gcf,[figDirName '/' fmNames{fmind(i)} 'Trapping_ref' num2str(coarsening)], '-png','-transparent')
    close gcf

end



% Plot Garn/Not/Ile and Tilje/Are aquifers:

Gt_gn = Gts{ logical(strcmpi(fmNames,'Garnfm')) };
Gt_nt = Gts{ logical(strcmpi(fmNames,'Notfm')) };
Gt_il = Gts{ logical(strcmpi(fmNames,'Ilefm')) };
Gt_tl = Gts{ logical(strcmpi(fmNames,'Tiljefm')) };
Gt_ar = Gts{ logical(strcmpi(fmNames,'Arefm')) };
    
% Garn/Not/Ile
figure; set(gcf,'Position',[1 1 1000 800])
hp1 = plotGrid(Gt_il.parent,'FaceColor','c','EdgeColor','none');
hp2 = plotGrid(Gt_nt.parent,'FaceColor','m','EdgeColor','none');
hp3 = plotGrid(Gt_gn.parent,'FaceColor',[1 .9 .9],'EdgeColor','none');
light('Position',[-1 -1 1],'Style','infinite'); lighting phong
view([-85 45]); axis equal tight off
hl = legend([hp3,hp2,hp1],{'Garn','Not','Ile'},'Location','East');
set(hl,'Fontsize',20)
set(gca,'DataAspect',[1 1 0.02])
set(gca,'Fontsize',20)
export_fig(gcf,[figDirName '/' 'GarnIle3D_ref' num2str(coarsening)], '-png','-transparent')

% Tilje/Are
figure; set(gcf,'Position',[1 1 1000 800])
hp1 = plotGrid(Gt_ar.parent,'FaceColor','m','EdgeColor','none');
hp2 = plotGrid(Gt_tl.parent,'FaceColor',[1 .9 .9],'EdgeColor','none');
light('Position',[-1 -1 1],'Style','infinite'); lighting phong
view([-85 45]); axis equal tight off
hl = legend([hp2,hp1],{'Tilje',[char(197),'re']},'Location','East');
set(hl,'Fontsize',20)
set(gca,'DataAspect',[1 1 0.02])
set(gca,'Fontsize',20)
export_fig(gcf,[figDirName '/' 'TiljeAre3D_ref' num2str(coarsening)], '-png','-transparent')
    

%% ------------------------------------------------------------------------
% 4. Tabulate results using capOutput
% Use sub-functions to retrieve NPD data for comparison (TODO)



%% Set up injection scenario(s):

%% ------------------------------------------------------------------------
% 1. Simulate injection with GUI:
% NB: keep coarse resolution as these formations are large
% NB: need to confirm sea floor data @@
info = getSeaInfo('NorwegianSea');
exploreSimulation('default_formation',   'Garnfm',           ...
                'grid_coarsening',     5,                    ...
                'seafloor_depth',      info.seafloor_depth,  ...
                'seafloor_temp',       info.seafloor_temp,   ...
                'temp_gradient',       info.temp_gradient,   ...
                'inj_time',            100 * year   );
            
%% ------------------------------------------------------------------------
% 2. Set-up injection scenario manually:
Gt_gn = Gts{ logical(strcmpi(fmNames,'Garnfm')) };
rock2D_gn = rock2Ds{ logical(strcmpi(fmNames,'Garnfm')) };

Gt      = Gt_gn;
rock2D  = rock2D_gn;
ta      = trapAnalysis(Gt,'false');

[ capOutput, hfig, hax ] = getTrappingPlots(Gt, ta, rock2D, 'NorwegianSea');

seainfo = getSeaInfo('NorwegianSea');
wellinfo = getWellInfo(Gt, capOutput);
export_fig(gcf,[figDirName '/' 'WellsGarn_ref' num2str(coarsening)], '-png','-transparent')

% assess impact of presence of producer wells:
% - first, run case with producer wells
[ wellSols, states, sim_report, opt, var ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo );

% - second, run another case where producer wells are turned off
wellinfo_off                 = wellinfo;
wellinfo_off.wellCoords_prod = [];
[ wellSols_off, states_off, sim_report_off, opt_off, var_off ] = ...
    runGenericInjectionScenario( Gt, rock2D, seainfo, wellinfo_off );

% options to inspect results:
figure;
plotToolbar(Gt, states_off)
clf
plotToolbar(Gt, states{end}, states{end}.s(:,2)>0.01)
plotGrid(Gt, 'FaceColor','none')

figure;
plotWellSols(wellSols_off)

% a) using 4 wells in Ile/Not/Garn to match NPD's estimated 400 million
% tonne (0.4 Gt) capacity (Riis and Halland, 2014, pg 5262). This capacity
% is reportedly low (Atlas, chp 5, pg 105). **The simulation grid in the
% Froan Basin (FBS) is depicted in the map on pg. 95 of the Atlas.**



% b) Also, make comparison to Lothe et al 2014 which studied injection into
% Garn.


%% Post-processing
% initState, schedule added to var structure

    dh = []; % for subscale trapping?
    figure; plot(1); ax = get(gcf, 'currentaxes');
    % NB: {var.initState, states{:}}
    reports = makeReports(var.model.G, {var.initState, states{:}}, ...
                             var.model.rock, var.model.fluid, var.schedule, ...
                             [var.model.fluid.res_water, var.model.fluid.res_gas], ...
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



