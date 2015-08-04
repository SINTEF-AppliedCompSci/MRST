%% Sleipner simulation
% The following studies the Sleipner benchmark data set which comes from
% Singh et al 2010 (SPE paper).

% The M9Z1.grdecl file and other necessary files should be downloaded and
% placed in co2lab/data/sleipner/. The following script will read those
% files to make the appropriate MRST-type grids (G and Gt). Also 'rock2D'
% is required. The first time the grids and rock data are generated, they
% are saved to /co2lab/data/mat/SleipnerGlobalCoords.mat to avoid
% re-generation every time this script is run. Since the grids are in the
% physical coordinate system, the injection location is specified in this
% same coordinate system.

moduleCheck('ad-core');
mrstVerbose on

% selection of what will be plotted:
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
% Default parameters/inputs:

% select the grid model to load/use:
useIEAGHG_model     = false;
useOriginal_model   = true;



% Merge user-specified inputs with other defaults:

% Select which parameters to modify from original data:
mod_rock        = true;    
mod_rhoCO2      = true;

% Then, set parameter modifier factors:
por_mod     = 0.6;
perm_mod    = 3;
rhoCO2_mod  = 2/3;



%


%% 1. Load formation
% need G, Gt, rock, rock2D

% Load and use Sleipner Benchmark model directly:

% If using Sleipner.mat, which has been created using
% makeSleipnerVEmodel():
%load Sleipner.mat
%OriginX = 436000; % used to line up plots between actual coordinates and Sleipner benchmark model centroids
%OriginY = 6469000;
%wellCellIndex = sub2ind(Gt.cartDims, 36,78);

% If using SleipnerGlobalCoords.mat, first check to see
% SleipnerGlobalCoords.mat exists. Otherwise generate it. (note: code taken
% from resTiltUtsira.m)
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
    

if mod_rock
    disp('Original rock parameters are being modified ...')
    
    % modify parameters
    rock.poro   = rock.poro * por_mod;
    rock2D.poro = rock2D.poro * por_mod;
    rock.perm   = rock.perm * perm_mod;
    rock2D.perm = rock2D.perm * perm_mod;
end
    
%     % to inspect set-up of 3D grid and 2D top-surface grid
%     figure;
%     subplot(1,3,1)
%     plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
%     plotCellData(G, G.cells.centroids(:,3))
%     
%     subplot(1,3,2)
%     plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
%     plotCellData(Gt, Gt.cells.z)
%     
%     subplot(1,3,3)
%     plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
%     plotCellData(G, G.cells.centroids(:,3))
%     plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(3)
    
    % Inspect orientation of grid (compare with orientation shown in Singh
    % et al 2010, Cavanagh 2013, etc.)
    figure; set(gcf,'Position',[1000 1000 1000 1000])
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(3)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x'); ylabel('y'); zlabel('z');
    % Create textarrow
    view(44,22)
    annotation(gcf,'textarrow',[0.4734 0.5391],[0.7825 0.81],'String',{'North'});
    set(gca,'FontSize',14)
    
    
    % Inspect orientation of grid: x-sectional views
    figure; set(gcf,'Position',[1000 1000 1100 500])
    % facing 
    subplot(2,2,1)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(270,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing East')
    
    % facing 
    subplot(2,2,3)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(90,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing West')
    
    % facing 
    subplot(2,2,2)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(180,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing South')
    
    % facing 
    subplot(2,2,4)
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8 .8 .8]); view(3)
    plotCellData(G, G.cells.centroids(:,3))
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'r'); view(0,0)
    set(gca,'DataAspect',[1 1 0.02]); grid
    xlabel('x, meters'); ylabel('y, meters'); zlabel('z, meters');
    title('Facing North')
    



% Plot of Sleipner Benchmark Region:
% bounds to zoom into:
zoomX1 = 0.436e6; % from Fig 2 in Singh et al 2010.
zoomY1 = 6.469e6;
zoomX2 = 0.441e6;
zoomY2 = 6.476e6;

% bounds of 2008 plume:
ZoomX1 = 0.4375e6;
ZoomY1 = 6.47e6;
ZoomX2 = 0.4395e6;
ZoomY2 = 6.474e6;

% Get boundary faces of formation (or grid region)
bf = boundaryFaces(Gt);


%% 2. Basic routine to perform VE simulation, using simulateScheduleAD().
% _________________________________________________________________________
% 1) set up initial state, OR get literature data:
[rho, mu, sr, sw] = getValuesSPE134891();

% rename variables:
water_density = rho(1) * kilogram / meter ^3;
rhoCref = rho(2) * kilogram / meter ^3;

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
% 2) set up schedule (wells, bc, etc.).

% WELLS:

%if useOptionA
    
    % Well location at Sleipner: (x,y) = (4.38516e5, 6.47121e6) (Singh et al.
    % 2010)
    wellXcoord = 4.38516e5;
    wellYcoord = 6.47121e6;

%     % get cell index of specified coordinate location:
%     wellIndexI_candidates = find(Gt.cells.centroids(:,1)<wellXcoord);
%     wellIndexJ_candidates = find(Gt.cells.centroids(:,2)<wellYcoord);
% 
%     PossibleIndexI = zeros(Gt.cells.num,1);
%     PossibleIndexI(wellIndexI_candidates) = Gt.cells.centroids(wellIndexI_candidates,1);
% 
%     PossibleIndexJ = zeros(Gt.cells.num,1);
%     PossibleIndexJ(wellIndexJ_candidates) = Gt.cells.centroids(wellIndexJ_candidates,2);
% 
%     % get the matching PossibleIndexI and PossibleIndexJ that are non-zero, and
%     % get the cell index of the maximum matching value
%     matchingIndex = PossibleIndexI.*PossibleIndexJ;
%     maxMatchingIndex = find(matchingIndex,1,'last');
% 
%     wellCellIndex = maxMatchingIndex;


    % OR: compute distances between all cell centroids to specified well
    % location, and use the cell centroid which has the minimum squared norm.
    % TODO.
    dv = bsxfun(@minus, Gt.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
    [v,i] = min(sum(dv.^2, 2));
    
    wellCellIndex = i; % or Gt.cells.indexMap(i);

    [i, j] = ind2sub(Gt.cartDims, wellCellIndex);

%end

% Check coordinate that wellCellIndex corresponds to:
wellCoord_x = Gt.cells.centroids(wellCellIndex,1);
wellCoord_y = Gt.cells.centroids(wellCellIndex,2);
wellCoord_z = 0;

% Compare actual against simulated injection location
if plotActualVsSimInjectLocation
figure; hold on

% actual location
plot(wellXcoord,wellYcoord,'ok', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
% simulated location
plot(wellCoord_x,wellCoord_y,'xk', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
legend('Actual injection location','Simulation injection location')

rectangle('Position',[zoomX1 zoomY1 zoomX2-zoomX1 zoomY2-zoomY1], 'EdgeColor','r', 'LineWidth',3,...
          'LineStyle','--')
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
axis equal %tight

end


% Well injection rate:
%annual_inj_rate = 1; % Mt/year
%inj_rate = annual_inj_rate * 1e9 / rhoCref / (365*24*60*60); % m^3/s

% Note: variable annual injection rates are given in Singh et al 2010,
% however it is likely their paper contains a typo. The injection rates
% were reported as in units of meter^3/day, however a more realistic value
% implies the units are meter^3/year.
inj_year = [1999; 2000; 2001; 2002; 2003; ...
         2004; 2005; 2006; 2007; 2008; 2009];
inj_rates  = [2.91e4; 5.92e4; 6.35e4; 8.0e4; 1.09e5; ...
         1.5e5; 2.03e5; 2.69e5; 3.47e5; 4.37e5; 5.4e5] * meter^3/year;
     % meter^3/s, automatically converted by specifying units.
inj_rates_MtPerYr = inj_rates.*(rhoCref/1e9*365*24*60*60); % Mt/year
% TODO - use sleipner_rates.txt:
% 1.9980000e+03   1.9990000e+03   2.0000000e+03   2.0010000e+03   2.0020000e+03   2.0030000e+03   2.0040000e+03   2.0050000e+03   2.0060000e+03   2.0070000e+03   2.0080000e+03   2.0090000e+03   2.0100000e+03   2.0110000e+03   2.0120000e+03   2.0130000e+03   2.0140000e+03   2.0150000e+03   2.0160000e+03   2.0170000e+03   2.0180000e+03   2.0190000e+03   2.0200000e+03   2.0210000e+03   2.0220000e+03   2.0230000e+03   2.0240000e+03   2.0250000e+03   2.0260000e+03   2.0270000e+03   2.0280000e+03   2.0290000e+03   2.0300000e+03   2.0310000e+03
% 0.0000000e+00   0.0000000e+00   9.7268000e+06   2.9532000e+07   5.0765000e+07   7.7537600e+07   1.1396160e+08   1.6414880e+08   2.3221100e+08   3.2226000e+08   4.3840760e+08   5.8476560e+08   7.6544580e+08   9.8456000e+08   1.2462200e+09   1.5545376e+09   1.9136246e+09   2.3275928e+09   2.8005540e+09   3.3366200e+09   3.9399026e+09   4.6145136e+09   5.3645648e+09   6.1941680e+09   7.1074350e+09   8.1084776e+09   9.2014076e+09   1.0390337e+10   1.1679377e+10   1.3072640e+10   1.4574238e+10   1.6188282e+10   1.7918884e+10   1.9770156e+10
% 0.0000000e+00   0.0000000e+00   1.8189116e+07   5.5224840e+07   9.4930550e+07   1.4499531e+08   2.1310819e+08   3.0695826e+08   4.3423457e+08   6.0262620e+08   8.1982221e+08   1.0935117e+09   1.4313836e+09   1.8411272e+09   2.3304314e+09   2.9069853e+09   3.5784780e+09   4.3525985e+09   5.2370360e+09   6.2394794e+09   7.3676179e+09   8.6291404e+09   1.0031736e+10   1.1583094e+10   1.3290903e+10   1.5162853e+10   1.7206632e+10   1.9429930e+10   2.1840435e+10   2.4445837e+10   2.7253824e+10   3.0272087e+10   3.3508313e+10   3.6970192e+10
% 0.0000000e+00   0.0000000e+00   2.6137918e+04   7.9358577e+04   1.3641603e+05   2.0835953e+05   3.0623833e+05   4.4110169e+05   6.2399886e+05   8.6597911e+05   1.1780917e+06   1.5713858e+06   2.0569108e+06   2.6457158e+06   3.3488502e+06   4.1773632e+06   5.1423041e+06   6.2547221e+06   7.5256664e+06   8.9661863e+06   1.0587331e+07   1.2400150e+07   1.4415692e+07   1.6645007e+07   1.9099144e+07   2.1789152e+07   2.4726081e+07   2.7920979e+07   3.1384896e+07   3.5128881e+07   3.9163983e+07   4.3501252e+07   4.8151738e+07   5.3126488e+07

% Put into schedule fields --> [inj period 1; inj period 2; etc...; migration period]
for i = 1:numel(inj_rates)
    schedule.control(i).W = addWell([], Gt.parent, rock2D, wellCellIndex, ...
        'name', sprintf('W%i', i), 'Type', 'rate', 'Val', inj_rates(i), 'comp_i', [0 1]);
end
schedule.control(end+1).W       = schedule.control(1).W;
schedule.control(end).W.name    = 'W_off';
schedule.control(end).W.val     = 0;



% BOUNDARY CONDITIONS:
% First get the faces of the boundaries. face.neighbors are the indices of
% the cells on either side of the faces, i.e., face.neighbor(100,1) and
% face.neighbor(100,2) give the index of the cells on either side of face
% with index 100. Any 0 cell index means there is no cell, i.e., the face
% is along an external boundary of the domain. Thus bdryFaces may be
% obtained by finding all the face indices that contain a 0 cell index on
% either side.
bdryFaces = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );
% Then use function bc = addBC(bc, faces, type, value, varargin)
bc = addBC( [], bdryFaces, ...
    'pressure', Gt.faces.z(bdryFaces) * water_density * norm(gravity), ...
    'sat', [1 0] );
% Put into schedule fields --> [injection period; migration period]
for i = 1:numel(schedule.control)
    schedule.control(i).bc = bc;
end
             

% TIME STEP:
% Specify and compute time step size for injection period.
% ***Note: applied to each inj_rate***
inj_time = 1 * year; % gets converted to seconds
inj_steps = 1;
dTi = inj_time / inj_steps; % timestep size in seconds

% Specify and compute time step size for migration period. 
mig_time = 1 * year;
mig_steps = 1;
dTm = mig_time / mig_steps; % timestep size in seconds

% For simulation schedule
istepvec = repmat( ones(inj_steps, 1) * dTi , [numel(inj_rates) 1] );
mstepvec = ones(mig_steps, 1) * dTm;

% schedule.step.val and schedule.step.control are same size arrays:
% schedule.step.control is a index (1,2,...) indicating which control (i.e., schedule.control) is to be used for the timestep.
schedule.step.control   = [ones(inj_steps, 1); ...
    ones(inj_steps, 1) * 2; ...
    ones(inj_steps, 1) * 3; ... 
    ones(inj_steps, 1) * 4; ...
    ones(inj_steps, 1) * 5; ...
    ones(inj_steps, 1) * 6; ...
    ones(inj_steps, 1) * 7; ...
    ones(inj_steps, 1) * 8; ...
    ones(inj_steps, 1) * 9; ...
    ones(inj_steps, 1) * 10; ...
    ones(inj_steps, 1) * 11; ...
    ones(mig_steps, 1) * 12];
% schedule.step.val is the timestep (size) used for that control step.
schedule.step.val       = [istepvec; mstepvec];


% confirm inj_rate, inj_year, and time steps
if plotInjectRateOverTime
figure; set(gcf,'Position',[1000 1000 1500 400])
subplot(1,3,1)
timeSinceInjLast = 0;
for i = 1:numel(schedule.step.val)
    timeSinceInj(i)  = schedule.step.val(i)/365/24/60/60 + timeSinceInjLast; % years
    rateNow(i)       = schedule.control(schedule.step.control(i)).W.val; % m^3/s
    timeSinceInjLast = timeSinceInj(i);
end
plot(timeSinceInj, rateNow, 'o--k')
xlabel('time, yr','FontSize',14)
ylabel('injection rate, m^3/s','FontSize',14)
hold off
ax = gca;
ax.XTickLabel = ax.XTick + inj_year(1)-1;

subplot(1,3,2)
plot(timeSinceInj, rateNow*(rhoCref/1e9*365*24*60*60), 'o--k') % Mt/yr
xlabel('time, yr','FontSize',14')
ylabel('injection rate, Mt/year','FontSize',14)
hold off
ax = gca;
ax.XTickLabel = ax.XTick + inj_year(1)-1;

subplot(1,3,3) %(compare this plot against Cavanagh 2013, fig 3)
massNowLast = 0;
%X = 1; % need to use timestep size here!! % need units in years
for i = 1:numel(schedule.step.val)
    timeStepNow(i) = schedule.step.val(i)/(365*24*60*60); % yr
    massNow(i) = rateNow(i)*(rhoCref/1e9*365*24*60*60)*timeStepNow(i) + massNowLast; % Mt/yr
    massNowLast = massNow(i);
end
plot(timeSinceInj, massNow, 'o--k')
xlabel('time, yr','FontSize',14)
ylabel('Accumulated mass injected, Mt','FontSize',14)
hold off
ax = gca;
ax.XTickLabel = ax.XTick + inj_year(1)-1;
end


% _________________________________________________________________________
% 3) set up model (grid, rock and fluid properties).
seafloor_temp = 7; % Celsius
seafloor_depth = 100; % meters
temp_gradient = 35.6; % Celsius / km
caprock_temperature =  273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient; % Kelvin

% water density, pvt info:
water_compr_val = 0; %4.3e-5/barsa; % will convert to compr/Pa

% pressure at which water density equals the reference density:
ref_p           = mean(initState.pressure); % use mean pressure as ref for linear compressibilities

% pore volume multiplier:
pvMult          = 0; %1e-5/barsa;

% dissolution:
dis_max         = (53 * kilogram / meter^3) / rhoCref; % from CO2store

% kwm? 0.75, 0.54 in Appendix of Singh et al 2010.


fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
                              'fixedT'      , caprock_temperature, ...
                              'wat_mu_ref'  , mu(1), ...
                              'co2_mu_ref'  , mu(2), ...
                              'wat_rho_ref' , water_density, ...
                              'co2_rho_ref' , rhoCref, ...
                              ... %'wat_rho_pvt' , [water_compr_val, ref_p], ...
                              'pvMult_p_ref', ref_p, ...
                              'pvMult_fac'  , pvMult, ...
                              'residual'    , [sw, sr] , ...
                              'dissolution' , false); %, 'dis_max', dis_max);
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
% 4) call to simulateScheduleAD().
[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);


% _________________________________________________________________________
% 5) Look at results:

% BHP VS TIME
time = sim_report.ReservoirTime;
bhp = zeros(numel(wellSols),1);
for i = 1:numel(wellSols)
    bhp(i) = wellSols{i}.bhp; % bhp is in Pa=N/m^2
end
if plotBHPvsTime
figure;
plot(time/365/24/60/60,bhp,'x--')
xlabel('Reservoir time, years'); ylabel('well bhp, Pascals=10^{-5}bars');
end


% ACCUM CO2 VS TIME (compare this plot against Cavanagh 2013, fig 3)
accumCO2sat = zeros(numel(states),1);
accumCO2mass = zeros(numel(states),1);
for i = 1:numel(states)
    accumCO2sat(i) = sum( states{i}.s(:,2).*model.G.cells.volumes ); % sat vals are in terms of pore volume
    
    satCO2          = states{i}.s(:,2);
    densityCO2      = fluid.rhoG(states{i}.pressure); 
    accumCO2mass(i) = sum( model.rock.poro .* model.G.cells.volumes .* model.G.cells.H .* satCO2 .* densityCO2 );
end
if plotAccumCO2vsTime
figure;
plot(time/365/24/60/60,accumCO2mass/1e9,'o-')
xlabel('Reservoir time, years'); ylabel('Accumlated CO2 mass, Mt (or 10^9 kg)');
end


% END of SIMULATION PROFILES
if plotEndOfSimResults
% meaningful profiles
press       = states{end}.pressure;
pressDiffFromHydrostatic = press - initState.pressure;
densityCO2  = fluid.rhoG(states{end}.pressure);  % fluid.rhoG is function handle to get CO2 density
satCO2      = states{end}.s(:,2);
massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg

bf = boundaryFaces(Gt);

h = figure;
subplot(1,6,1)
plotCellData(Gt, caprock_temperature, 'EdgeColor','none')
title('Caprock Temperature'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Kelvin'; set(hcb, 'fontSize', 18)

subplot(1,6,2)
plotCellData(Gt, press, 'EdgeColor','none')
title('Pressure'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,3)
plotCellData(Gt, pressDiffFromHydrostatic, 'EdgeColor','none')
title({'Pressure diff. from hydrostatic';'i.e., the initial condition'}); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pascals'; set(hcb, 'fontSize', 18)

subplot(1,6,4)
plotCellData(Gt, densityCO2, 'EdgeColor','none')
title('CO2 density at Caprock'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'kg/m^3'; set(hcb, 'fontSize', 18)

subplot(1,6,5)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, satCO2, satCO2>(0.1/100), 'EdgeColor','none') %satCO2~=0)
title('Saturation of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Pore Volume CO2 Saturation'; set(hcb, 'fontSize', 18)

subplot(1,6,6)
%plotGrid(Gt, 'FaceColor', 'white', 'EdgeAlpha', 0)
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, massCO2/1e9, satCO2>(0.1/100), 'EdgeColor','none') % only plot plume that has sat > 0.1 percent
title('Mass of CO2'); axis off equal tight
hcb = colorbar;
hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)

end


% INVENTORY (from exploreSimulation.m)
if plotTrappingInventory
trapstruct = trapAnalysis(Gt, false); % true for cell-based method
dh = []; % for subscale trapping?
h2 = figure; plot(1); ax = get(h2, 'currentaxes');
reports = makeReports(model.G, {initState, states{:}}, model.rock, model.fluid, schedule, ...
                         [model.fluid.res_water, model.fluid.res_gas], ...
                         trapstruct, dh);
% reports contains soln states; could be used for plotting results.
directPlotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
%xlabel('Years since simulation start (1999)')
ax.XTickLabel = ax.XTick + inj_year(1);
end


%% Line plots of CO2 migrating plume data
% TODO: implement code to download the following .mat files from a repo
% (and put .mat files onto repo)
load layer9_polygons_1999.mat;
plume{1}.outline = CO2plumeOutline;
plume{1}.year    = 1999;

load layer9_polygons_2001.mat;
plume{2}.outline = CO2plumeOutline;
plume{2}.year    = 2001;

load layer9_polygons_2002.mat;
plume{3}.outline = CO2plumeOutline;
plume{3}.year    = 2002;

load layer9_polygons_2004.mat;
plume{4}.outline = CO2plumeOutline;
plume{4}.year    = 2004;

load layer9_polygons_2006.mat;
plume{5}.outline = CO2plumeOutline;
plume{5}.year    = 2006;

%load layer9_polygons_2006a.mat;
%plume{6}.outline = CO2plumeOutline;
%plume{6}.year    = 2006.3;

%load layer9_polygons_2006b.mat;
%plume{7}.outline = CO2plumeOutline;
%plume{7}.year    = 2006.6;

load layer9_polygons_2008.mat;
plume{6}.outline = CO2plumeOutline;
plume{6}.year    = 2008;


% PROFILES AT SELECT TIME
Years2plot          = [1999; 2001; 2002; 2004; 2006; 2008];
ReservoirTime2plot  = (Years2plot - inj_year(1)+1 ).*(365*24*60*60); % seconds
CO2plumeOutline_SatTol = (0.01/100); % adjust this value if patch error occurs (which happens when no massCO2 present at specified sat tol)
if plotCO2simVsCO2obsData
h = figure; set(gcf, 'Position', [1000 1000 1500 1100])
hold on

for i = 1:numel(ReservoirTime2plot)
    
    % get reservoir time index
    [rti,J] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));

    % meaningful profiles
    press       = states{rti}.pressure;
    pressDiffFromHydrostatic = press - initState.pressure;
    densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{rti}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
    maxMassCO2(i)= max(massCO2);
    
    subplot(2,numel(ReservoirTime2plot)/2,i)
    hold on

    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, massCO2/1e9, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none') % only plot plume that has sat > tolerance specified 
    title({'Mass of CO2 at';['year ', num2str(Years2plot(i))]}, 'fontSize', 18); axis equal
    hcb = colorbar; hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)

    % add Sleipner region box:
    %rectangle('Position',[zoomX1 zoomY1 zoomX2-zoomX1 zoomY2-zoomY1], 'EdgeColor','r', 'LineWidth',3,...
    %      'LineStyle','--')
    %axis tight
    %xlim([zoomX1 zoomX2]); ylim([zoomY1 zoomY2])
    
    % add CO2 plume outline (check matching year):
    if plume{i}.year == Years2plot(i)
        disp('Plotting Observed CO2 plume outline...')
        %if useOptionA
            line(plume{i}.outline(:,1), plume{i}.outline(:,2), ...
            'LineWidth',3, 'Color','r')
        %else
%             disp('But first converting observed CO2 plume outline coords...')
%             line(plume{i}.outline(:,1)-OriginX, plume{i}.outline(:,2)-OriginY, ...
%             'LineWidth',4, 'Color','r')
        %end
    end
    
    % add injection point:
    %if exist('wellXcoord')
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    %end
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)
    
    axis tight
    
    % We visualize the spill paths between structural traps
    mapPlot(gcf, Gt, 'traps', trapstruct.traps, 'rivers', trapstruct.cell_lines, ...
        'maplines',20); % default maplines is 40

end

% make plotting adjustments to subplots
axesHandles = get(gcf,'children');

% set caxis to be between 0 and the max CO2 mass value plotted in any of
% the subplots
cmax = max(maxMassCO2./1e9);
for i=1:numel(axesHandles)
    if strcmpi(axesHandles(i).Type,'axes')
        axesHandles(i).CLim = [0 cmax];
        
        if ZoomIntoPlume
           axesHandles(i).XLim = [ZoomX1 ZoomX2];
           axesHandles(i).YLim = [ZoomY1 ZoomY2];
        end
    end
end

end


% % PLOT OF CO2 PLUME DATA ON ITS OWN (NO SIMULATION RESULTS)
% figure
% for i=1:numel(plume)
%     
%     subplot(1,numel(plume),i)
%     hold on
%     line(plume{i}.outline(:,1), plume{i}.outline(:,2), ...
%         'LineWidth',2, 'Color','k')
%     title(['Year ', num2str(plume{i}.year)])
%     
%     % add injection point:
%     % actual location
%     plot(wellXcoord,wellYcoord,'ok', ...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','g',...
%         'MarkerSize',10)
%     % simulated location
% %     plot(wellCoord_x,wellCoord_y,'xk', ...
% %         'MarkerEdgeColor','k',...
% %         'MarkerFaceColor','g',...
% %         'MarkerSize',10)
%     hold off
%     
%     rectangle('Position',[zoomX1 zoomY1 zoomX2-zoomX1 zoomY2-zoomY1], 'EdgeColor','r', 'LineWidth',3,...
%           'LineStyle','--')
%     axis equal tight
%     %ylim([6.425e6 6.55e6])
%     %xlim([0.41e6 0.5e6])
% end






%% 3. Analyze traps with spill-point analysis, and plot structural traps.
%ta = trapAnalysis(Gt, false); % true for cell-based method
if plotTrapProfiles
ta = trapstruct;

figure;
plotCellData(Gt, ta.traps, ta.traps~=0, 'FaceColor', 'r', 'EdgeColor','none')
bf = boundaryFaces(Gt, ta.traps~=0);
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);

ta_volumes = volumesOfTraps(Gt, ta);

N = 1;
fprintf('Coarsening level %d:\n', N);
fprintf('  Num. global traps: %d\n', numel(ta_volumes));
fprintf('  Total trap volume: %e m3\n', sum(ta_volumes));
fprintf('  Avg. global trap size: %e m3\n', mean(ta_volumes));

% plot of structural traps in red:
title({['Sleipner, coarsening level=' num2str(N)]; ...
    [num2str(numel(ta_volumes)) ' structural traps shown in red']}, 'FontSize', 14)

axis equal tight

% ****************
% We can use plotCellData and color the faces according to the trap volume
% in km^3. The cell of Gt is associated with a particular trap by ta.traps,
% and is 0 otherwise. The volume of a particular trap is given by
% ta_volumes.
figure;
subplot(1,5,1)
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')

trapcells = ta.traps~=0;
cellsTrapVol = zeros(Gt.cells.num,1);
cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));
plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Trap Volume, km^3';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);

% Note that it may be more useful to color the faces according to the
% trapped mass of CO2. This requires using the CO2 density data which can
% vary with depth (since pressure and temperature gradients exist in the
% formation), to determine the cooresponding CO2 mass trapped by the traps.
% Porosity might also be required, in order to compute both trapping
% volumes and CO2 mass trapping.
% See functions which produce Utsira plots of traps, catchment area, etc:
%   co2 = CO2props();
%   utsiraCapacity('rhoCfun', @co2.rho);
% The above functions also compute breakdown of CO2 trapping types
% (structural, residual, and dissolution) without VE simulation.






%% 3. b) get trapping breakdown (structural, residual, dissoluion)
% The following section computes the breakdown of structural, residual, and
% dissolution trapping in the formation, assuming CO2 is injected and fills
% the formation completely, i.e., all spill paths are utilized. The
% contents of this section closely follows the implementation found in
% exploreCapacity(), which is a GUI. Note: some fluid parameters are fixed
% to those used in exploreCapacity.

% first, use default values to compute theoretical capacity (upper bound):
out_default = exploreParameterRanges(Gt, rock2D, ta);

% ****************
% We can visualize the structural traps in terms of CO2 mass, not volume.

% Distributed CO2 mass under structural traps: 
cellsTrapCO2Mass = zeros(Gt.cells.num,1);
cellsTrapCO2Mass(trapcells) = out_default.strap_mass_co2(trapcells);

% Cumulative CO2 mass under structural traps:
trapcaps = accumarray(ta.traps(trapcells), out_default.strap_mass_co2(trapcells));
trapcap_tot = zeros(Gt.cells.num,1); %ones(size(ta.traps)) * NaN;
trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));

% Boundary of formation
bf = boundaryFaces(Gt);

subplot(1,5,2)
hold on
%plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, cellsTrapCO2Mass/1e9, cellsTrapCO2Mass~=0)

    % add injection point:
    %if exist('wellXcoord')
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    %end
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Distributed CO2 Mass under Trap, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);

subplot(1,5,3)
hold on
%plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0)

    % add injection point:
    %if exist('wellXcoord')
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    %end
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Accumulated CO2 Mass under Trap, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);



% ****************
% We can also visualize the reachable structural capacity, which is based
% on the output of trapAnalysis: ta. Then we compute cumulative structural
% trap volume reached when injecting at a given cell
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
    structural_mass_reached(cix) = sum(out_default.strap_mass_co2(aix)); %#ok

end

subplot(1,5,4)
hold on
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
plotCellData(Gt, structural_mass_reached/1e3/1e6, 'EdgeColor','none');

    % add injection point:
    %if exist('wellXcoord')
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    %end
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar; %hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Reachable structural capacity, Mt';
grid; axis tight; set(hcb, 'fontSize', 18); set(gca, 'fontSize', 10);


% ****************
% We visualize the spill paths between structural traps
subplot(1,5,5)
hold on
mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);

    % add injection point:
    %if exist('wellXcoord')
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    %end
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)
    
grid; axis equal tight;

end


%% Show the structural trap
% After the simulation has completed, we are interested in seeing how the
% location of the CO2 plume after a long migration time corresponds to the
% trapping estimates produced by trapAnalysis. This is done by finding the
% trap index of the well injection cell and then plotting the trap along
% with the final CO2 plume.

% Well in 2D model
%WVE = convertwellsVE(W, Gt.parent, Gt, rock2D);
WVE = reports(end).W; % take well of last control since first control might be initial conditions (no well)

% Generate traps and find the trap corresponding to the well cells
trap = trapstruct.trap_regions([WVE.cells]);

% Plot the areas with any significant CO2 height along with the trap in red
if plotTrapAndPlumeCompare
figure;
plotCellData(Gt, reports(end).sol.h, reports(end).sol.h > 0.01)
plotGrid(Gt, trapstruct.traps == trap, 'FaceColor', 'red', 'EdgeColor', 'w')
plotGrid(Gt, 'FaceColor', 'None', 'EdgeAlpha', .1);

legend({'CO2 Plume', 'Trap'})
axis tight off
view(-20, 75)
title('End of simulation CO2 compared to algorithmically determined trap')
end


