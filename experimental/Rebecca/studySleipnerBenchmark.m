%% Sleipner simulation
% The following studies the Sleipner benchmark data set which comes from
% Singh et al 2010 (SPE paper).

% There are several options available to run this benchmark, which are
% explained below:

    % The first option is related to which Sleipner grid is loaded.
    % Currently, this script can handle the IEAGHG model (when
    % 'useIEAGHG_model' is set to true) and the ORIGINAL model (when
    % 'useOriginal_model' is set to true). Only the ORIGINAL model has the
    % additional option of performing model refinement (i.e., grid, perm,
    % poro, etc).

    % The second option is which annual injection rates to use. The
    % available rates are either those which come from Singh et al 2010
    % (SPE 134891) (i.e., rates between 1999 and 2009) or those which came
    % with the original Sleipner benchmark files (i.e., rates between 1999
    % and 2030).

    % The third option is whether to modify the rock and/or fluid
    % properties which have been loaded from the model file. Currently,
    % permeability and porosity are modified when mod_rock=true, and CO2
    % reference density is modified when mod_rhoCO2=true. The modification
    % factors must be defined.

    % Other options available in this script coorespond to the typical
    % input parameters that must be defined such as sea floor temperature,
    % time step sizes, migration period, etc.
    
    
% Notes about the files needed to load the specified grid:

    % To use the IEAGHG model (i.e., when 'useIEAGHG_model=true'), the
    % necessary grdecl files are:
        % M9X1.grdecl, M9X1_perm_X_mD_.inc, M9X1_perm_Y_mD_.inc,
        % M9X1_poro___.inc, SLEIPNER.DATA
    % The above files should be downloaded and placed in:
        % co2lab/data/sleipner/

        
    % To use the ORIGINAL model (i.e., when 'useORIGINAL_model=true'), the
    % necessary grdecl files are:
        % sleipner_prep.data
    % The above file(s) should be placed in:
        % co2lab/data/sleipner/original/
    

% Notes about other files needed for this script to run:

    % CO2 plume outlines - layer9_polygons_XXXX.mat (where XXXX is the
    % year, such as 1999) files should be created and placed in current
    % working directory.
    
    % Injection rates - SleipnerOriginalInjectionRates.mat should be
    % created and placed in current working directory.


% Note that the following script will read the appropriate grdecl files and
% make the MRST-type grids (G and Gt). The first time the grids and rock
% data are generated, they are saved to /co2lab/data/mat/ to avoid
% re-generation every time this script is run. Since the grids are in the
% physical coordinate system, the injection location is specified in this
% same coordinate system.

%%


moduleCheck('ad-core');
mrstVerbose on
gravity on;


% ******************** START OF USER OPTIONS ******************************


% selection of what will be plotted:
plotModelGrid                   = true;
plotInitialPressure             = true;
plotActualVsSimInjectLocation   = true;
plotInjectRateOverTime          = true;
plotBHPvsTime                   = true;
plotAccumCO2vsTime              = true;
plotEndOfSimResults             = true;
plotCO2simVsCO2obsData          = true; ZoomIntoPlume = true; % if false, entire grid is plotted
plotTrappingInventory           = true;
plotTrapProfiles                = true;
plotTrapAndPlumeCompare         = true;


% Trapping analysis method
isCellBasedMethod = false; % true to use cell-based method, false to use node-based method


% FOR PLOTS:
CO2plumeOutline_SatTol  = (0.01/100); % adjust this value if patch error occurs (which happens when no massCO2 present at specified sat tol)
press_deviation = 0;  % from hydrostatic (percent) --> used for trapping capacity calculation, not simulation
N = 1; % coarsening level of grid. (no option yet to re-fine or coarsen Sleipner grid). ==> in progress


% For plotting of CO2 plumes
% bounds of 2008 plume:
ZoomX1 = 0.4375e6;
ZoomY1 = 6.47e6;
ZoomX2 = 0.4395e6;
ZoomY2 = 6.474e6;


% OPTION - Select the grid model to load/use:
useIEAGHG_model     = false;
useOriginal_model   = true; useRefinedGrid = true;
if useRefinedGrid
    refineLevel = 4;
else
    refineLevel = []; % or 1?
end


% Physical coordinate of injection well (Singh et al. 2010)
wellXcoord      = 4.38516e5;
wellYcoord      = 6.47121e6;


% OPTION - Well injection rate:
useRatesFromSPE134891 = true; extrapolateRates = false; % TODO - implement for true (i.e., extrapolation)
useUserDefinedRates = false; % TODO - implement
useSleipnerOriginalInjectionRates = false; 

if useRatesFromSPE134891
    % Note: variable annual injection rates are given in Singh et al 2010,
    % however it is likely their paper contains a typo. The injection rates
    % were reported as in units of meter^3/day, however a more realistic value
    % implies the units are meter^3/year.
    inj_year   = [1999; 2000; 2001; 2002; 2003; 2004; 2005; 2006; 2007; 2008; 2009];
    inj_rates  = [2.91e4; 5.92e4; 6.35e4; 8.0e4; 1.09e5; 1.5e5; 2.03e5; 2.69e5; 3.47e5; 4.37e5; 5.4e5] .* meter^3/year;
    % inj_rates is in meter^3/s
    
elseif useUserDefinedRates
    %error('TODO: implement user defined rates.')
    
elseif useSleipnerOriginalInjectionRates
    
    [ CumSurfVol_m3, Mass_kg, ReservoirVol_m3, year ] = getSleipnerOriginalInjectionRates();
    % OLD: load SleipnerOriginalInjectionRates.mat
    
    % We only take the years between 1999 and 2030 as injection years since
    % ReservoirVol_m3 is the cumulative value as of Jan 1st of each year
    % 1999 is taken as first injection year since ReservoirVol_m3 is
    % non-zero starting in 2000 (Jan 1st) which implies there was an
    % injection rate for 1999. The ReservoirVol_m3 amount corresponding to
    % 2031 is used to determine the injection rate in the previous year
    % (2030), thus we assume no injection occurs in 2031 (or a value could
    % be extrapolated).
    inj_year    = year(2:end-1);              clear year
    inj_rates   = zeros(1,numel(inj_year));
    
    ReservoirVol_m3 = ReservoirVol_m3(2:end);
    for i = 1:numel(ReservoirVol_m3)-1
        inj_rates(i) = ( ReservoirVol_m3(i+1) - ReservoirVol_m3(i) ); % annual rate
    end
    inj_rates = inj_rates * meter^3/year; % in meter^3/second
    
else
    error('The injection rate option was either unvalid or not selected.')

end


% Specify and compute time step size for injection period.
% ***Note: inj_time and inj_steps are applied to each inj_rate given***
inj_time    = 1 * year; % DEFAULT. CAN ONLY ADJUST NUMBER OF STEPS.
inj_steps   = 5;
dTi         = inj_time / inj_steps; % timestep size in seconds

% Specify and compute time step size for migration period. 
mig_time    = 1 * year; % CAN ADJUST.
mig_steps   = 1;        % CAN ADJUST.
dTm         = mig_time / mig_steps; % timestep size in seconds


% Specify fluid properties:
[rho, mu, sr, sw]   = getValuesSPE134891();
water_density       = rho(1) * kilogram / meter ^3;
rhoCref             = rho(2) * kilogram / meter ^3;

seafloor_temp       = 7; % Celsius
seafloor_depth      = 100; % meters
temp_gradient       = 35.6; % Celsius / km
water_compr_val     = 0; %4.3e-5/barsa; % will convert to compr/Pa
pvMult              = 0; %1e-5/barsa;
isDissOn            = false;
dis_max             = (53 * kilogram / meter^3) / rhoCref; % from CO2store

% kwm? 0.75, 0.54 in Appendix of Singh et al 2010.


% OPTION - Select which parameters to modify from original data:
mod_rock        = true;    
mod_rhoCO2      = true;

% Then, set parameter modifier factors:
por_mod     = 0.6;
perm_mod    = 3;
rhoCO2_mod  = 2/3;



% ************************ END OF USER OPTIONS ****************************


%% 1. Load formation
% loads or generates file containing G, Gt, rock, rock2D

if useIEAGHG_model

    try

        disp(' -> Reading SleipnerGlobalCoords.mat');
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,'SleipnerGlobalCoords')); clear datadir
        %return;

    catch %#ok<*CTCH>
        disp('    SleipnerGlobalCoords.mat has not yet been created.');
        disp('    Building G, Gt, and rock2D from grdecl files...')

        % First loading of Sleipner Eclipse grid (to get PERMX, PERMZ,
        % PORO)
        sdir    = fullfile('data', 'sleipner');
        disp([' -> Reading data from: ' sdir]);
        grdecl  = readGRDECL(fullfile(mrstPath('co2lab'), sdir, 'SLEIPNER.DATA'));
        % this grdecl contains: cartDims, COORD, ZCORN, ACTNUM, PERMX,
        % PERMZ, PORO
        clear sdir

        % Reshaping
        lines = reshape(grdecl.COORD,6,[]);
        grdecl.COORD = lines(:); clear lines

        % Then, we remove the bottom and top layers that contain shale
        grdecl.ACTNUM(grdecl.PERMX<200) = 0;


        % Second loading of Sleipner Eclispe grid, to get MAPAXES
        moduleCheck('deckformat', 'mex');
        sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'M9X1.grdecl'); % IEAGHG 
        fn      = fopen(sl_file);
        gr  = readGRID(fn, fileparts(sl_file), initializeDeck());
        % this grdecl contains: GRID, and others. grdecl.GRID contains
        % MAPUNITS, MAPAXES, cartDims, COORD, ZCORN, ACTNUM
        fclose(fn);


        % Add data loaded from first loading of Sleipner Eclispe grid
        grdecl.MAPAXES = gr.GRID.MAPAXES;
        clear gr sl_file


        % Recompute X and Y coordinates in terms of the provided axes
        % (depths, Z, do not require any recomputation)
        coords        = reshape(grdecl.COORD,3,[])';
        coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
        coords        = coords';
        grdecl.COORD  = coords(:); clear coords


        % Next, we process the grid and compute geometry
        mrstModule add libgeometry opm_gridprocessing
        G = processgrid(grdecl);
        G = mcomputeGeometry(G);

        % Adding tags needed by topSurfaceGrid
        G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];

        % Construct petrophysical model
        rock = grdecl2Rock(grdecl, G.cells.indexMap);
        rock.perm = convertFrom(rock.perm, milli*darcy);
        clear grdecl

        % Construct top-surface grid
        disp(' -> Constructing top-surface grid');
        [Gt, G] = topSurfaceGrid(G);
        rock2D  = averageRock(rock, Gt);


        % Store data
        disp(' ')
        disp(' -> Writing SleipnerGlobalCoords.mat')
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save(fullfile(datadir,'SleipnerGlobalCoords'), 'G', 'Gt', 'rock', 'rock2D');
        clear datadir

    end

elseif useOriginal_model && useRefinedGrid
    
    % there should be a directory named "sleipner/original/" which contains
    % files such as sleipner_prep.data, injection rates, etc.
    try

        disp([' -> Reading OriginalSleipnerGlobalCoords_numRef', num2str(refineLevel), '.mat']);
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,['OriginalSleipnerGlobalCoords_numRef', num2str(refineLevel)])); clear datadir
        %return;

    catch %#ok<*CTCH>
        disp(['    OriginalSleipnerGlobalCoords_numRef', num2str(refineLevel), '.mat has not yet been created.']);
        disp('    Building REFINED G, Gt, and rock2D from grdecl files...')

        % Open and read original/sleipner_prep.data
        moduleCheck('deckformat', 'mex');
        sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'original', 'sleipner_prep.data'); % IEAGHG original
        fn      = fopen(sl_file);
        grdecl  = readGRID(fn, fileparts(sl_file), initializeDeck());
        % this grdecl contains: GRID, and others. grdecl.GRID contains
        % cartDims, MAPUNITS, MAPAXES, COORD, ZCORN, ACTNUM as well as
        % PERMX, PERMY, PERMZ, and PORO
        fclose(fn);
        
        % Rename
        grdecl = grdecl.GRID;

        % Reshaping
        lines = reshape(grdecl.COORD,6,[]);
        grdecl.COORD = lines(:); clear lines

        % Then, we remove the bottom and top layers that contain shale
        grdecl.ACTNUM(grdecl.PERMX<200) = 0;

        
        % Then we perform a few steps for GRID REFINEMENT:
        % - first cut the grdecl by removing the top and bottom layers (i.e.,
        % keep layers 2 to 6, thus remove layers 1 and 7). There will be 5
        % layers remaining.
        lowerBound = [1 1 2];
        upperBound = [grdecl.cartDims(1) grdecl.cartDims(2) 6];
        ind = [lowerBound; upperBound]';
        grdecl_cut = cutGrdecl(grdecl, ind);
        % grdecl_new contains: cartDims, COORD, ZCORN, PERMX, PERMY, PERMZ,
        % PORO, ACTNUM. These all coorespond to the 5 remaining layers.
        
        % - second coarsen the grdecl in the z direction since we don't
        % need to keep the resolution of 5 layers (i.e., we treat 5 layers
        % as 1 homogeneous layer for VE). In order to ensure PERMX, PERMY,
        % PERMZ, and PORO are coarsened, 'only_grid' must be set to false.
        dim = [1 1 5];
        grdecl_coarsened = coarseGrdecl(grdecl_cut, dim, 'only_grid',false);
        % grdecl_new now contains: cartDims, COORD, ZCORN, ACTNUM, PERMX,
        % PERMY, PERMZ, PORO which correspond to the single layer.
        
        % - third refine the grdecl in the X and Y direction (Z already has
        % a resolution of 1 cell) by the refinement level specified
        dim = [refineLevel; refineLevel; 1];
        grdecl_refined = refineGrdecl(grdecl_coarsened, dim);
        grdecl_refined.MAPAXES = grdecl.MAPAXES;
        grdecl_refined.MAPUNITS = grdecl.MAPUNITS;
        
        
        % Recompute X and Y coordinates in terms of the provided axes
        % (depths, Z, do not require any recomputation)
        coords        = reshape(grdecl_refined.COORD,3,[])';
        coords(:,1:2) = mapAxes(coords(:,1:2), grdecl_refined.MAPAXES);
        coords        = coords';
        grdecl_refined.COORD  = coords(:); clear coords


        % Next, we process the grid and compute geometry
        mrstModule add libgeometry opm_gridprocessing
        G = processGRDECL(grdecl_refined); % note: processgrid() didn't work
        G = mcomputeGeometry(G);

        % Adding tags needed by topSurfaceGrid
        G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];

        % Construct petrophysical model
        rock = grdecl2Rock(grdecl_refined, G.cells.indexMap);
        rock.perm = convertFrom(rock.perm, milli*darcy);
        clear grdecl grdecl_cut grdecl_coarsened grdecl_refined

        % Construct top-surface grid
        disp(' -> Constructing top-surface grid');
        [Gt, G] = topSurfaceGrid(G);
        rock2D  = averageRock(rock, Gt);


        % Store data
        disp(' ')
        disp([' -> Writing OriginalSleipnerGlobalCoords_numRef', num2str(refineLevel), '.mat'])
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save( fullfile(datadir,['OriginalSleipnerGlobalCoords_numRef', num2str(refineLevel), '.mat']), ...
            'G', 'Gt', 'rock', 'rock2D', ...
            '-v7.3'); % use flag to save as v7.3 to avoid issues with compression of a lot of data
        clear datadir
        
    end
        
    
elseif useOriginal_model
    
    % there should be a directory named "sleipner/original/" which contains
    % files such as sleipner_prep.data, injection rates, etc.
    try

        disp(' -> Reading OriginalSleipnerGlobalCoords.mat');
        datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
        load(fullfile(datadir,'OriginalSleipnerGlobalCoords')); clear datadir
        %return;

    catch %#ok<*CTCH>
        disp('    OriginalSleipnerGlobalCoords.mat has not yet been created.');
        disp('    Building G, Gt, and rock2D from grdecl files...')

        % Open and read original/sleipner_prep.data
        moduleCheck('deckformat', 'mex');
        sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'original', 'sleipner_prep.data'); % IEAGHG original
        fn      = fopen(sl_file);
        grdecl  = readGRID(fn, fileparts(sl_file), initializeDeck());
        % this grdecl contains: GRID, and others. grdecl.GRID contains
        % cartDims, MAPUNITS, MAPAXES, COORD, ZCORN, ACTNUM as well as
        % PERMX, PERMY, PERMZ, and PORO
        fclose(fn);
        
        % Rename
        grdecl = grdecl.GRID;

        % Reshaping
        lines = reshape(grdecl.COORD,6,[]);
        grdecl.COORD = lines(:); clear lines

        % Then, we remove the bottom and top layers that contain shale
        grdecl.ACTNUM(grdecl.PERMX<200) = 0;


        % Recompute X and Y coordinates in terms of the provided axes
        % (depths, Z, do not require any recomputation)
        coords        = reshape(grdecl.COORD,3,[])';
        coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
        coords        = coords';
        grdecl.COORD  = coords(:); clear coords


        % Next, we process the grid and compute geometry
        mrstModule add libgeometry opm_gridprocessing
        G = processGRDECL(grdecl); % note: processgrid() didn't work
        G = mcomputeGeometry(G);

        % Adding tags needed by topSurfaceGrid
        G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];

        % Construct petrophysical model
        rock = grdecl2Rock(grdecl, G.cells.indexMap);
        rock.perm = convertFrom(rock.perm, milli*darcy);
        clear grdecl

        % Construct top-surface grid
        disp(' -> Constructing top-surface grid');
        [Gt, G] = topSurfaceGrid(G);
        rock2D  = averageRock(rock, Gt);


        % Store data
        disp(' ')
        disp(' -> Writing OriginalSleipnerGlobalCoords.mat')
        if ~isdir(datadir)
           mkdir(datadir);
        end
        save(fullfile(datadir,'OriginalSleipnerGlobalCoords'), 'G', 'Gt', 'rock', 'rock2D');
        clear datadir
         
    end
    
    
else
    error('A model to read or generate has not been specified.')

end

    
%% Modify original parameters (optional) and visualize model grids
if mod_rock
    disp('Original rock parameters are being modified ...')
    
    % modify parameters
    rock.poro   = rock.poro .* por_mod;
    rock2D.poro = rock2D.poro .* por_mod;
    rock.perm   = rock.perm .* perm_mod;
    rock2D.perm = rock2D.perm .* perm_mod;
end
    
if plotModelGrid
    [ hfig, hax ] = plot3DandTopGrids( G, Gt );
    
end


% Get boundary faces of formation (or grid region)
bf = boundaryFaces(Gt);


%% 2. Basic routine to perform VE simulation, using simulateScheduleAD().
% _________________________________________________________________________
% a) set up initial state, OR get literature data:


if mod_rhoCO2
    disp('Original CO2 density value is being modified ...')
    rhoCref = rhoCref * rhoCO2_mod; 
end

initState.pressure  = Gt.cells.z * norm(gravity) * water_density;   % hydrostatic pressure, in Pa=N/m^2
initState.s         = repmat([1 0], Gt.cells.num, 1);               % sat of water is 1, sat of CO2 is 0
initState.sGmax     = initState.s(:,2);                             % max sat of CO2 is initially 0
initState.rs        = 0 * initState.sGmax;                          % initially 0

if plotInitialPressure
    figure;
    plotCellData(Gt, initState.pressure)
    title('Initial Pressure (hydrostatic)'); axis off tight equal
    hcb = colorbar;
    hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)
end


% _________________________________________________________________________
% b) set up schedule (wells, bc, etc.).

% WELLS:

dv = bsxfun(@minus, Gt.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
[v,i] = min(sum(dv.^2, 2));

wellCellIndex = i; % or Gt.cells.indexMap(i);

[i, j] = ind2sub(Gt.cartDims, wellCellIndex);


% Check coordinate that wellCellIndex corresponds to:
wellCoord_x = Gt.cells.centroids(wellCellIndex,1);
wellCoord_y = Gt.cells.centroids(wellCellIndex,2);
wellCoord_z = 0;

% Compare actual against simulated injection location
if plotActualVsSimInjectLocation
    [ hfig, hax ] = plotRealVsDiscreteInjLoc(Gt,bf,wellXcoord,wellYcoord,wellCoord_x,wellCoord_y);
end

inj_rates_MtPerYr = inj_rates.*(rhoCref/1e9*365*24*60*60); % Mt/year


% Put into schedule fields --> [inj period 1; inj period 2; etc...; migration period]
for i = 1:numel(inj_rates)
    schedule.control(i).W = addWell([], Gt.parent, rock2D, wellCellIndex, ...
        'name', sprintf('W%i', i), 'Type', 'rate', 'Val', inj_rates(i), 'comp_i', [0 1]);
end
schedule.control(end+1).W       = schedule.control(1).W;
schedule.control(end).W.name    = 'W_off';
schedule.control(end).W.val     = 0;



% BOUNDARY CONDITIONS: (TODO - put in function and give options for
% different bdry condition types)
% First get the faces of the boundaries. face.neighbors are the indices of
% the cells on either side of the faces, i.e., face.neighbor(100,1) and
% face.neighbor(100,2) give the index of the cells on either side of face
% with index 100. Any 0 cell index means there is no cell, i.e., the face
% is along an external boundary of the domain. Thus bdryFaces may be
% obtained by finding all the face indices that contain a 0 cell index on
% either side.
bdryFaces = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );

bdryType = 'pressure';
bdryVal  = Gt.faces.z(bdryFaces) * water_density * norm(gravity);
% Then use function bc = addBC(bc, faces, type, value, varargin)
bc = addBC( [], bdryFaces, bdryType, bdryVal, 'sat', [1 0] );


% Put into schedule fields --> [injection period; migration period]
for i = 1:numel(schedule.control)
    schedule.control(i).bc = bc;
end
             

% TIME STEP:

% For simulation schedule
istepvec = repmat( ones(inj_steps, 1) * dTi , [numel(inj_rates) 1] );
mstepvec = ones(mig_steps, 1) * dTm;

% schedule.step.val and schedule.step.control are same size arrays:
% schedule.step.val is the timestep (size) used for that control step.
schedule.step.val       = [istepvec; mstepvec];
% schedule.step.control is a index (1,2,...) indicating which control (i.e., schedule.control) is to be used for the timestep.
schedule.step.control = [];
for i = 1:numel(schedule.control)
    
    if schedule.control(i).W.val ~= 0
        % an injection period
        schedule.step.control = [schedule.step.control; ones(inj_steps, 1) * i];

    elseif schedule.control(i).W.val == 0
        % a migration period
        schedule.step.control = [schedule.step.control; ones(mig_steps, 1) * i];
        
    end

end


% confirm inj_rate, inj_year, and time steps
if plotInjectRateOverTime
    [ hfig, hax ] = plotInjectRateVsTime(schedule,inj_year,rhoCref);
end


% _________________________________________________________________________
% c) set up model (grid, rock and fluid properties).

caprock_temperature = 273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient; % Kelvin

% pressure at which water density equals the reference density:
ref_p           = mean(initState.pressure); % use mean pressure as ref for linear compressibilities


fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                              'fixedT'      , caprock_temperature, ...
                              'wat_mu_ref'  , mu(1), ...
                              'co2_mu_ref'  , mu(2), ...
                              'wat_rho_ref' , water_density, ...
                              'co2_rho_ref' , rhoCref, ...
                              'wat_rho_pvt' , [water_compr_val, ref_p], ...
                              'pvMult_p_ref', ref_p, ...
                              'pvMult_fac'  , pvMult, ...
                              'residual'    , [sw, sr] , ...
                              'dissolution' , isDissOn, 'dis_max', dis_max);
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);


% Prepare plotting (from runSleipner.m)
% We will make a composite plot that consists of several parts: a 3D plot
% of the plume, a pie chart of trapped versus free volume, a plane view of
% the plume from above, and two cross-sections in the x/y directions
% through the well
%wellIn = % in 3D
%opts = {'slice', wellIx, 'Saxis', [0 1-fluidVE.sw], ...
%   'maxH', 5, 'Wadd', 10, 'view', [130 50]};
%plotPanelVE(G, Gt, W, sol, 0.0, zeros(1,6), opts{:});


% _________________________________________________________________________
% d) call to simulateScheduleAD().
[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);

% TODO: write output file
% to save: G, Gt, rock, rock2D, all input parameters, all outputs (results)



% _________________________________________________________________________
%% Post-Processing

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
if plotEndOfSimResults
    [ hfig, hax ] = plotProfilesAtGivenTime('final', Gt, states, initState, fluid, model, caprock_temperature);
end


% Call to trap analysis:
ta = trapAnalysis(Gt, isCellBasedMethod); % true for cell-based method


% INVENTORY (from exploreSimulation.m)
if plotTrappingInventory
    dh = []; % for subscale trapping?
    h2 = figure; plot(1); ax = get(h2, 'currentaxes');
    reports = makeReports(model.G, {initState, states{:}}, model.rock, model.fluid, schedule, ...
                             [model.fluid.res_water, model.fluid.res_gas], ...
                             ta, dh);
    % reports contains soln states; could be used for plotting results.
    directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
    %xlabel('Years since simulation start (1999)')
    ax.XTickLabel = ax.XTick + inj_year(1);
end


%% Line plots of CO2 migrating plume data
% Note: to run the following function, first downloaded the plume .mat
% files from https://bitbucket.org/mrst/mrst-co2lab/downloads, and place on
% current working directory path
plume = getLayer9CO2plumeOutlines();


% PROFILES AT SELECT TIME
if plotCO2simVsCO2obsData
    
    Years2plot = [1999; 2001; 2002; 2004; 2006; 2008];
    
    [ hfig, hax ] = subplotCO2simVsCO2obsData(Years2plot, inj_year, plume, sim_report, ...
            Gt, states, fluid, model, ...
            wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, ta, ...
            ZoomIntoPlume, ZoomX1, ZoomX2, ZoomY1, ZoomY2, ...
            CO2plumeOutline_SatTol);

end



%% Structural trapping plots

if plotTrapProfiles

    ta_volumes = volumesOfTraps(Gt, ta);
    
    fprintf('Coarsening level %d:\n', N);
    fprintf('  Num. global traps: %d\n', numel(ta_volumes));
    fprintf('  Total trap volume: %e m3\n', sum(ta_volumes));
    fprintf('  Avg. global trap size: %e m3\n', mean(ta_volumes));
    
    
%     % PLOT FACES OF STRUCTURAL TRAPS (redundant with following plots)
%     figure;
%     
%     plotCellData(Gt, ta.traps, ta.traps~=0, 'FaceColor', 'r', 'EdgeColor','none')
%     bf = boundaryFaces(Gt, ta.traps~=0);
%     plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
% 
%     title({['Sleipner, coarsening level=' num2str(N)]; ...
%         [num2str(numel(ta_volumes)) ' structural traps shown in red']}, 'FontSize', 14)
%     axis equal tight

 
    % PLOT TRAPS COLORED BY CO2 MASS STORAGE CAPACITY
    figure; set(gcf,'Position',[1000 1000 3000 500])
    
    %
    subplot(1,5,1)
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);

    trapcells = ta.traps~=0;
    cellsTrapVol = zeros(Gt.cells.num,1);
    cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));
    plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0)

    set(gca,'DataAspect',[1 1 1/100])
    hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
    hcb.Label.String = 'Trap Volume, km^3';
    grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);


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
    subplot(1,5,2)
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, cellsTrapCO2Mass/1e9, cellsTrapCO2Mass~=0)

    set(gca,'DataAspect',[1 1 1/100])
    hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
    hcb.Label.String = 'Distributed CO2 Mass under Trap, Mt';
    grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);

    
    %
    subplot(1,5,3)
    hold on
    %plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0)

    set(gca,'DataAspect',[1 1 1/100])
    hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
    hcb.Label.String = 'Accumulated CO2 Mass under Trap, Mt';
    grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);



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
    subplot(1,5,4)
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, structural_mass_reached/1e3/1e6, 'EdgeColor','none');

    set(gca,'DataAspect',[1 1 1/100])
    hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
    hcb.Label.String = 'Reachable structural capacity, Mt';
    grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);


    % PLOT SPILL PATHS AND TOPOLOGY
    subplot(1,5,5)
    hold on
    mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);

    grid; axis equal tight;
    
    % For making plotting adjustments to subplots
    axesHandles = get(gcf,'children');
    
    % Add Injection Location In Each Subplot:
    for i=1:numel(axesHandles)
        if strcmpi(axesHandles(i).Type,'axes')
            
            subplot(axesHandles(i))
            % actual location
            plot(wellXcoord,wellYcoord,'o', ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
            % simulated location
            plot(wellCoord_x,wellCoord_y,'x', ...
                'LineWidth',3,  ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',10)
        end
    end
    
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
    
    figure; set(gcf,'Position',[1000 1000 1600 1000])
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


