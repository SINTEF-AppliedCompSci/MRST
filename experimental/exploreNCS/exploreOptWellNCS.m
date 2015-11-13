%% Get optimized well rates in NCS formations

fmNames = {'Arefm';'Tiljefm';'Ilefm';'Garnfm';...
           'Bjarmelandfm';...
           'Stofm';'Tubaenfm'};
       
for i = 1:numel(fmNames);
    
    fn = fmNames(i);
    fn = fn{:};     % to make fn a string
    
    %% Call to optimize
    [Gt, optim, init, history, other] = ...
    optimizeFormation(  'modelname',        fn, ...
                        'coarse_level',     5, ...
                        'num_wells',        4, ...
                        'trapfile_name',    'utsira_subtrap_function_3.mat' );
                    
                    

    %% Saving optimized results
    savedir = fullfile('opt_results/', fn,'/');
    if ~isdir(savedir)
       mkdir(savedir);
    end
    save([savedir  'init'] , 'init');
    save([savedir  'optim'], 'optim');
    save([savedir  'Gt'], 'Gt');
    save([savedir, 'history'], 'history');
    save([savedir, 'other'], 'other');
    
    clearvars init optim Gt history other

end


%% Get optimized well rates for an array of wells in trap region

% directory for saving figures
figDirName = 'NorwegianSeaFigs';
mkdir(figDirName)


% 1. Load formation grids and rock properties
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

% 2. Set-up well injection sites:
fn                  = 'Garnfm';

Gt                  = Gts{ logical(strcmpi(fmNames,fn)) };
rock2D              = rock2Ds{ logical(strcmpi(fmNames,fn)) };
ta                  = trapAnalysis(Gt,'false');
[ capOutput, ~, ~ ] = getTrappingPlots(Gt, ta, rock2D, 'NorwegianSea');
seainfo             = getSeaInfo('NorwegianSea');
wellinfo            = getWellInfo(Gt, capOutput, 'prod',false, 'setInjRates',true);
wcells  = wellinfo.cinx_inj;
qtot    = wellinfo.vols_inj / 10000;  %@@
isteps = 10;
msteps = 31;
itime  = 50   * year;
mtime  = 3000 * year;

% put wells, rates, time steps into schedule
schedule = setSchedule(Gt, rock2D, wcells, qtot, ...
                                isteps, itime, msteps, mtime, true, ...
                                'minval', sqrt(eps));
% pass schedule into optimizeFormation(), instead of using the default
% schedule inside of that function. NB: bdry conditions will be added to
% schedule, however no option for bdry type yet.

% 3. Call to optimize
% NB: dryrun = true gives init case, doesn't proceed to optim
[Gt, optim, init, history, other] = ...
      optimizeFormation2('modelname',       fn, ...
                        'coarse_level',     coarsening, ...
                        'num_wells',        4, ...
                        'trapfile_name',    'utsira_subtrap_function_3.mat', ... %); %, ...
                        'schedule',         schedule, ...
                        'dryrun',           true);
                    
% 4. Saving optimized results
savedir = fullfile('opt_results/', fn, '/Arrays/NTrapRegions/');
if ~isdir(savedir)
   mkdir(savedir);
end
save([savedir  'init'] , 'init');
save([savedir  'optim'], 'optim');
save([savedir  'Gt'], 'Gt');
save([savedir, 'history'], 'history');
save([savedir, 'other'], 'other');




