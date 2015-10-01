function studySleipnerBenchmarkFUN( varargin )
% The following studies the Sleipner benchmark data set which comes from
% Singh et al 2010 (SPE paper).

% SYNOPSIS:
%
% studySleipnerBenchmarkFUN()
% studySleipnerBenchmarkFUN('mycase','useOriginal_model')
% studySleipnerBenchmarkFUN('refineLevel', 2)
% studySleipnerBenchmarkFUN('myInjRates', 'useSleipnerOriginalInjectionRates')
% studySleipnerBenchmarkFUN('mod_rock_perm', true)
% studySleipnerBenchmarkFUN('mod_rock_perm', true, 'perm_fac', 2.75)


% OPTIONS:


% SEE ALSO:
%   runSleipner, analyseSleipner (co2lab/experimental/project/tests/)

%%
mrstModule add co2lab
moduleCheck('ad-core','opm_gridprocessing','mex','deckformat', ...
    'coarsegrid','upscaling','incomp','mrst-experimental');
mrstVerbose on
gravity on;


%%
opt = struct(   'mycase',               'useIEAGHG_model', ...
                'refineLevel',          1, ...
                'myInjRates',           'useRatesFromSPE134891', ...
                'mod_rock_perm',        false, ...
                'mod_rock_poro',        false, ...
                'mod_rhoCO2',           false, ...
                'perm_fac',             [], ...
                'poro_fac',             [], ...
                'rhoCO2_fac',           [], ...
                'pauseBeforeSolver',    false               );
opt = merge_options(opt, varargin{:});


% Default parameter modifier factors:
if opt.mod_rock_perm
    if isempty(opt.perm_fac)
        opt.perm_fac    = 3;
    end
end
if opt.mod_rock_poro
    if isempty(opt.poro_fac)
        opt.poro_fac     = 0.6;
    end
end  
if opt.mod_rhoCO2
    if isempty(opt.rhoCO2_fac)
        opt.rhoCO2_fac  = 2/3;
    end  
end


%%
% selection of what will be plotted before simulation starts:
plotModelGrid                   = false;
plotInitialPressure             = false;
plotActualVsSimInjectLocation   = false;
plotInjectRateOverTime          = false;


% Trapping analysis method (used for Post-processing, not simulation).
isCellBasedMethod = false; % true to use cell-based method, false to use node-based method


% FOR PLOTS:
CO2plumeOutline_SatTol  = (0.01/100); % adjust this value if patch error occurs (which happens when no massCO2 present at specified sat tol)
press_deviation = 0;  % from hydrostatic (percent) --> used for trapping capacity calculation, not simulation


% For plotting of CO2 plumes
% bounds of 2008 plume:
ZoomX1 = 0.4375e6;
ZoomY1 = 6.47e6;
ZoomX2 = 0.4395e6;
ZoomY2 = 6.474e6;


% Physical coordinate of injection well (Singh et al. 2010)
wellXcoord      = 4.38516e5;
wellYcoord      = 6.47121e6;


% OPTION - Well injection rate:
switch opt.myInjRates
    
    % Note, inj_rates are in terms of reservoir rates (i.e., the volumetric
    % rate of CO2 entering layer 9, not the volumetric surface rate).
    % Seismic imaging provided estimates of how much CO2 accumlated in the
    % pore space of layer 9. These volumes were likely converted into a
    % mass using an infered CO2 density, and then into a surface rate using
    % the CO2 density at the surface. Specifying the inj_rates in terms of
    % reservoir volume instead of surface volume allows one to test other
    % CO2 densities without the need to modify a surface volume injection
    % rate.
    
    case 'useRatesFromSPE134891'
        % See Singh et al 2010 for more info about how they determined
        % these rates. Note: the injection rates were reported as surface
        % rates. Both volume and mass were given, thus surface density can
        % be calculated (=1.87 kg/m3). The CO2 density at reservoir
        % conditions was reported as 760 kg/m3.
        inj_year   = [1999; 2000; 2001; 2002; 2003; 2004; 2005; 2006; 2007; 2008; 2009];
        inj_rates  = [2.91e4; 5.92e4; 6.35e4; 8.0e4; 1.09e5; 1.5e5; 2.03e5; 2.69e5; 3.47e5; 4.37e5; 5.4e5] .* meter^3/day;
        % inj_rates is in meter^3/s
        % Convert to rate at reservoir conditions
        inj_rates  = inj_rates.*(1.87/760);
        
    case 'useSleipnerOriginalInjectionRates'
        % See "Injection rates Layer 9.xls" under
        % co2lab/data/sleipner/original for more info about these rates.
        % Note: the CO2 density at reservoir conditions was reported as
        % 695.89 kg/m3, and the surface density was 1.87 kg/m3.
        [ inj_year, inj_rates ] = getSleipnerOriginalInjectionRates();
        % inj_rates is in meter^3/s
        
    otherwise
        error('The injection rate option was either invalid or not selected.')
        
end

% Plot inj_rates over inj_year. Plotting later occurs using schedule fields
figure
plot(inj_year, inj_rates, 'o')
ylabel('Reservoir rates, m^3/year')
xlabel('Year')


% Specify and compute time step size for injection period.
% ***Note: inj_time and inj_steps are applied to each inj_rate given***
inj_time    = 1 * year; % DEFAULT. CAN ONLY ADJUST NUMBER OF STEPS.
inj_steps   = 1;
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




%% 1. Load formation
% makeSleipnerModelGrid() looks for file or generates it from grdecl files
% and writes .mat file. Output is the variables G, Gt, rock, rock2D.

% get case info:
switch opt.mycase
    
    case 'useIEAGHG_model'
        modelname   = 'IEAGHGmodel';

    case 'useOriginal_model'
        modelname   = 'ORIGINALmodel';
        
    case 'useInhouse_model'
        modelname   = 'INHOUSEmodel';
        
    otherwise
        error('No such case')
end


% make grid model:
fprintf(['\nYour case is set to ' opt.mycase '.\n'])
fprintf(['You have chosen to refine the model grid ',num2str(opt.refineLevel),' times.\n'])
fprintf('\nGetting grid...\n\n')
[ G, Gt, rock, rock2D ] = makeSleipnerModelGrid('modelName', modelname, 'refineLevel',opt.refineLevel);
fprintf('\n\nGrid obtained.\n')

    
%% 2. Modify original parameters (optional) and visualize model grids
if opt.mod_rock_perm
    disp('Original rock permeabilities are being modified ...')
    rock.perm   = rock.perm .* opt.perm_fac;
    rock2D.perm = rock2D.perm .* opt.perm_fac;
end

if opt.mod_rock_poro
    disp('Original rock porosities are being modified ...')
    rock.poro   = rock.poro .* opt.poro_fac;
    rock2D.poro = rock2D.poro .* opt.poro_fac;
end
    
if plotModelGrid
    [ hfig, hax ] = plot3DandTopGrids( G, Gt );
    
end


% Get boundary faces of formation (or grid region)
bf = boundaryFaces(Gt);


%% 3. Basic routine to perform VE simulation, using simulateScheduleAD().
if opt.mod_rhoCO2
    disp('Original CO2 density value is being modified ...')
    rhoCref = rhoCref * opt.rhoCO2_fac; 
end

initState.pressure  = Gt.cells.z * norm(gravity) * water_density;   % hydrostatic pressure, in Pa=N/m^2
initState.s         = repmat([1 0], Gt.cells.num, 1);               % sat of water is 1, sat of CO2 is 0
initState.sGmax     = initState.s(:,2);                             % max sat of CO2 is initially 0
initState.rs        = 0 * initState.sGmax;                          % initially 0

if plotInitialPressure
    figure;
    plotCellData(Gt, initState.pressure, 'EdgeColor','none')
    title('Initial Pressure','fontSize', 18);
    % setColorbarHandle() is able to deal with handles of class 'double'
    % (pre-R2014) and graphic objects (post-R2014)
    [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Pascals', 'fontSize', 18 );
    axis off tight equal
end


% WELLS:
dv = bsxfun(@minus, Gt.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
[v,i] = min(sum(dv.^2, 2));
wellCellIndex = i; % or Gt.cells.indexMap(i);
[i, j] = ind2sub(Gt.cartDims, wellCellIndex);


% Compare actual against simulated injection location
if plotActualVsSimInjectLocation
    [ hfig, hax ] = plotRealVsDiscreteInjLoc(Gt,bf,wellXcoord,wellYcoord,wellCoord_x,wellCoord_y);
end



% Put into schedule fields --> [inj period 1; inj period 2; etc...; migration period]
for i = 1:numel(inj_rates)
    schedule.control(i).W = addWell([], Gt.parent, rock2D, wellCellIndex, ...
        'name', sprintf('W%i', i), 'Type', 'rate', 'Val', inj_rates(i), 'comp_i', [0 1]); % inj_rate should be mass rate / fluid.rhoGS 
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
% either side. (But will this include 'top' and 'bottom' faces?)
bdryFaces = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );

bdryVal  = Gt.faces.z(bdryFaces) * water_density * norm(gravity);
% Then use function bc = addBC(bc, faces, type, value, varargin)
bc = addBC( [], bdryFaces, 'pressure', bdryVal, 'sat', [1 0] );


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
% schedule.step.control is a index (1,2,...) indicating which control
% (i.e., schedule.control) is to be used for the timestep.
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
    [ hfig, hax, timeSinceInj, massNow ] = plotInjectRateVsTime(schedule,inj_year,rhoCref);
end


caprock_temperature = 273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient; % Kelvin

% pressure at which water density equals the reference density:
ref_p           = mean(initState.pressure); % use mean pressure as ref for linear compressibilities


%% 5. fluid and model set-up, call to solver
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

if opt.pauseBeforeSolver
    disp('do you wish to proceed to solver?')
    pause
end

[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);


%% 6. Save variables in workspace
% first, close all figures
close all
if strcmpi(opt.mycase(1:3), 'use') && strcmpi(opt.mycase(end-5:end), '_model')
    name = opt.mycase(4:end-6);
end
if strcmpi(opt.myInjRates(1:3), 'use')
    rateName = opt.myInjRates(4:end);
end
[modstr1, modstr2, modstr3] = deal('_');
if opt.mod_rock_perm
    modstr1 = ['_PermFac', num2str(opt.perm_fac), '_'];
end
if opt.mod_rock_poro
    modstr2 = ['PoroFac', num2str(opt.poro_fac), '_'];
end
if opt.mod_rhoCO2
    modstr3 = ['RhoFac', num2str(opt.rhoCO2_fac), '_'];
end
modstrs = [modstr1 modstr2 modstr3];
fileName = [name 'refNum' num2str(opt.refineLevel) '_' rateName modstrs datestr(clock,30) '.mat'];
save(fileName);
% TODO: only save certain variables to .mat file




end

