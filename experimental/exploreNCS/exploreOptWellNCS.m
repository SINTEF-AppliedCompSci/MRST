%% Get optimized well rates in NCS formations

fmNames = {'Arefm';'Tiljefm';'Ilefm';'Garnfm';...
           'Bjarmelandfm';...
           'Stofm';'Tubaenfm'};
       
useWellArrays = true;
coarsening = 5;
       
for i = 1:numel(fmNames);
    
    fn = fmNames(i);
    fn = fn{:};     % to make fn a string
    
    if useWellArrays
        
        % 1. Load formation grids and rock properties
        % NB: implement st checks if Gt.mat, rock2D.mat exists for fn
        [Gt, rock2D] = getFormationTopGrid(fn, coarsening);
        
        
        % 2. Set-up well injection sites:
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

        % 3. Call to optimize
        % pass schedule into optimizeFormation(), instead of using the default
        % schedule inside of that function. NB: bdry conditions will be added to
        % schedule, however no option for bdry type yet.
        [Gt, optim, init, history, other] = ...
              optimizeFormation2('modelname',       fn, ...
                                'coarse_level',     coarsening, ...
                                'num_wells',        4, ...
                                'trapfile_name',    'utsira_subtrap_function_3.mat', ... %); %, ...
                                'schedule',         schedule);

        % 4. Dir for saving results
        savedir = fullfile('opt_results/', fn, '/Arrays/NTrapRegions/');
        
        
    else
    
        % 1. Call to optimize
        [Gt, optim, init, history, other] = ...
        optimizeFormation(  'modelname',        fn, ...
                            'coarse_level',     coarsening, ...
                            'num_wells',        4, ...
                            'trapfile_name',    'utsira_subtrap_function_3.mat' );



        % 2. Dir for saving results
        savedir = fullfile('opt_results/', fn,'/');

    end
    
    %% Saving optimized results
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