%% Get optimized well rates in NCS formations

%fmNames = {'Arefm';'Tiljefm';'Ilefm';'Garnfm';...
%           'Bjarmelandfm';...
%           'Stofm';'Tubaenfm'};
fmNames = {'Utsirafm';'Stofm';'Garnfm'};
       
useWellArrays = false; %true; %true;
coarsenings = [5, 5, 10];
       
for i = 1:numel(fmNames);
    
    fn = fmNames(i);
    fn = fn{:};     % to make fn a string
    coarsening = coarsenings(i);
    
    if useWellArrays
        % Using trap structure of formation, an array of wells are placed
        % according to varargin passed in (i.e., to cover entire formation,
        % set 'limits' = 'none', or to cover only best N catchment area,
        % set 'limits' = 'bestReachCap' and 'numTopReachTrap' = N).
        
        
        % 1. Load formation grids and rock properties
        % NB: implement st checks if Gt.mat, rock2D.mat exists for fn
        [Gt, rock2D] = getFormationTopGrid(fn, coarsening);
        
        % 2. Set-up well injection sites:
        ta                  = trapAnalysis(Gt,false); % @@ syntax must be false, not 'false'
        trapCapacities      = getTrappingInfo(fn, coarsening, 'plotsOn',false);
        seainfo             = getSeaInfo(fn, 760);
        wellinfo            = getWellInfo(Gt, trapCapacities, 'prod',false, ...
            'setInjRates',true, 'limits','none', 'plotsOn',false, ...
            'DX',10*5000, 'DY',10*5000);
        wcells  = wellinfo.cinx_inj;
        %qtot    = wellinfo.vols_inj / 100;  %@@
        qtot    = wellinfo.vols_inj; % already in m^3 so do not need to divide by co2_density
        isteps = 10;
        msteps = 30;
        itime  = 50   * year;
        mtime  = 3000 * year;

        % put wells, rates, time steps into schedule
        % (4th varargin should be in m^3)
        schedule = setSchedule(Gt, rock2D, wcells, qtot, ...
                                        isteps, itime, msteps, mtime, true, ...
                                        'minval', sqrt(eps));
        

        % 3. Call to optimize
        % pass schedule into optimizeFormation(), instead of using the default
        % schedule inside of that function. NB: bdry conditions will be added to
        % schedule, however no option for bdry type yet.
        [Gt, optim, init, history, other] = ...
              optimizeFormation3('modelname',       fn, ...
                                'coarse_level',     coarsening, ...
                                'schedule',         schedule, ...
                                'leakPenalty',      10, ...
                                'dryrun',           true );

        % 4. Dir for saving results
        savedir = fullfile('opt_results/', fn, ['/ref' num2str(coarsening)], '/Arrays/', ['leakPenalty' num2str(other.leakPenalty) '/']);
        
        
    else
        % Use default options that are specified inside optimizeFormation()
        % i.e., a catchment area contains at most 1 well, and num_wells may
        % be specified.
        
        % 1. Call to optimize: use optForm3() since fluid was created using
        % 'smooth' top surface. Otherwise subscale trap file must be passed
        % in as varargin to avoid error from occurring.
        [Gt, optim, init, history, other] = ...
        optimizeFormation3(  'modelname',        fn, ...
                            'coarse_level',     coarsening, ...
                            'num_wells',        10, ...
                            'maximise_boundary_distance' , true);

        % 2. Dir for saving results
        savedir = fullfile('opt_results/',fn,['/ref' num2str(coarsening)],'/OneWellPerCatch/',['leakPenalty' num2str(other.leakPenalty) '/']);

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