%% Get optimized well rates in NCS formations

% fmNames = {'Arefm';'Tiljefm';'Ilefm';'Garnfm';...
%            'Bjarmelandfm';...
%            'Stofm';'Tubaenfm'};
fmNames = {'Stofm'};
       
useWellArrays = true;
coarsening = 5;
       
for i = 1:numel(fmNames);
    
    fn = fmNames(i);
    fn = fn{:};     % to make fn a string
    
    if useWellArrays
        % Using trap structure of formation, an array of wells are placed
        % according to varargin passed in (i.e., to cover entire formation,
        % set 'limits' = 'none', or to cover only best N catchment area,
        % set 'limits' = 'bestReachCap' and 'numTopReachTrap' = N).
        
        
        % 1. Load formation grids and rock properties
        % NB: implement st checks if Gt.mat, rock2D.mat exists for fn
        [Gt, rock2D] = getFormationTopGrid(fn, coarsening);
        
        % 2. Set-up well injection sites:
        ta                  = trapAnalysis(Gt,'false');
        [ capOutput, ~, ~ ] = getTrappingPlots(Gt, ta, rock2D, 'NorwegianSea');
        seainfo             = getSeaInfo('NorwegianSea');
        wellinfo            = getWellInfo(Gt, capOutput, 'prod',false, 'setInjRates',true, 'limits','none');
        wcells  = wellinfo.cinx_inj;
        qtot    = wellinfo.vols_inj / 100;  %@@
        isteps = 10;
        msteps = 31;
        itime  = 50   * year;
        mtime  = 3000 * year;

        % put wells, rates, time steps into schedule
        schedule = setSchedule(Gt, rock2D, wcells, qtot/seainfo.co2_density, ...
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
                                'schedule',         schedule, ...
                                'leakPenalty',      10, ...
                                'dryrun',           false );

        % 4. Dir for saving results
        savedir = fullfile('opt_results/', fn, '/Arrays/NTrapRegions/leakPen10/');
        
        
    else
        % Use default options that are specified inside optimizeFormation()
        % i.e., a catchment area contains at most 1 well, and num_wells may
        % be specified.
        
        % 1. Call to optimize
        [Gt, optim, init, history, other] = ...
        optimizeFormation2(  'modelname',        fn, ...
                            'coarse_level',     coarsening, ...
                            'num_wells',        4 );

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