function setup = fivespot_geothermal(varargin)
%Conceptual model for geothermal heat storage in five-spot pattern

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    % Description
    description = ['Conceptual model for geothermal heat storage in ', ...
                   'five-spot pattern'];
    % Optional input arguments
    month = 365*day/12;
    K0    = 273.15*Kelvin;
    options = struct( ...
        'cartDims'     , [31, 31, 32]                , ... % Grid dimensions
        'physDim'      , [10, 10, 20]*meter          , ... % Physical dimensions
        'permMat'      , 1*milli*darcy               , ... % Matrix permeability
        'poroMat'      , 0.05                        , ... % Matrix porosity
        'numFrac'      , 3                           , ... % Number of fractures
        'dfm'          , false                       , ... % Discrete fractures
        'poroFrac'     , 0.5                         , ... % Matrix porosity
        'apertureFrac' , 1e-3                        , ... % Fracture aperture
        'lambdaR'      , 2*watt/(meter*Kelvin)       , ... % Thermal conductivity
        'rhoR'         , 2700*kilogram/meter^3       , ... % Density
        'CpR'          , 1000*joule/(kilogram*Kelvin), ... % Specific heat capacity
        'lambdaF'      , 0.6*watt/(meter*Kelvin)     , ... % Thermal conductivity
        'Cp'           , 4.2*joule/(Kelvin*gram)     , ... % Heat capacity
        'cT'           , 207e-6/Kelvin               , ... % Thermal expansion
        'TRef'         , (273.15 + 20)*Kelvin        , ... % Density reference temperature
        'useEOS'       , false                       , ... % Use EOS for fluid props
        'rateCharge'   , 0.5*litre/second            , ... % Charge rate
        'rateDischarge', 0.5*litre/second            , ... % Discharge rate
        'bhp'          , 1*atm                       , ... % Produced BHP
        'tempCharge'   , K0 + 30*Kelvin              , ... % Charge temperature
        'tempDischarge', K0 + 10*Kelvin              , ... % Discharge temperature
        'timeCharge'   , 6*month                     , ... % Charging time
        'timeDischarge', 1*month                     , ... % Discharging time
        'dtCharge'     , 7*day                       , ... % Timesteps during charging
        'dtDischarge'  , 1*day                       , ... % Timestep during discharging
        'numCycles'    , 1                           , ... % Number of cycles (charge, discharge)
        'chargeOnly'   , false                       , ... % Simulate single charge period
        'dischargeOnly', false                       , ... % Simulate single discharge period
        'useGroupCtrl' , false                       , ... % Use group control
        'initPres'     , 1*atm                       , ... % Initial pressure
        'initTemp'     , 273.15 + 10                   ... % Initial temperature
    );
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    options = checkOptions(options);
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Define module dependencies
    require ad-core ad-props compositional geothermal
    % Enable gravity
    gravity reset on;
    % Make Cartesian grid
    G = computeGeometry(cartGrid(options.cartDims, options.physDim));
    % Make rock
    if ~options.dfm
        % Compute effective perm/poro if not DFM
        options.permMat = effectivePerm(options);
        options.poroMat = effectivePoro(options);
    end
    rock = makeRock(G, options.permMat, options.poroMat);
    rock = addThermalRockProps(rock, 'lambdaR', options.lambdaR, ...
                                     'rhoR'   , options.rhoR   , ...
                                     'CpR'    , options.CpR    );
    G0 = G;
    if options.dfm
        % Add fractures if DFM
        [G, rock] = addFractures(G, rock, options);
    end
                              
    % Make fluid
    fluid = initSimpleADIFluid('phases', 'W'                  , ...
                               'n'     , 1                    , ... 
                               'rho'   , 1000*kilogram/meter^3, ...
                               'mu'    , 1*centi*poise        );
    fluid = addThermalFluidProps(fluid, 'lambdaF', options.lambdaF, ...
                                        'Cp'     , options.Cp     , ...
                                        'cT'     , options.cT     , ...
                                        'useEOS' , options.useEOS );
    % Make model
    model = GeothermalModel(G, rock, fluid);
    model = model.validateModel();
    model.FacilityModel.implicitConnectionDP = true;
    % Make schedule
    schedule = setUpSchedule(G, G0, rock, fluid, options);
    % Set up initial state
    state0 = setUpInitialState(model, schedule.control(1).W, options);

    % Set plotting specs
    plotOptions = {'Size'              , [400, 600], ...
                   'PlotBoxAspectRatio', [1,1,2]   , ...
                   'View'              , [-15, 22] , ...
                   'Box'               , true      };
    %---------------------------------------------------------------------%
    
     % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
    %---------------------------------------------------------------------%

end

%-------------------------------------------------------------------------%
function options = checkOptions(options)
    
    assert(~(options.chargeOnly && options.dischargeOnly), ...
        'Cannot simulate only charge and only discharge at the same time');
    
end

%-------------------------------------------------------------------------%
function permEff = effectivePerm(options)
    
    permFrac = options.apertureFrac^2/12;
    f = options.numFrac*options.apertureFrac./options.physDim(3);
    permEff = f.*permFrac + (1-f).*options.permMat;
    permEff = [permEff, permEff, options.permMat];

end

%-------------------------------------------------------------------------%
function poroEff = effectivePoro(options)

    f = options.numFrac*options.apertureFrac./options.physDim(3);
    poroEff = f.*options.poroFrac + (1-f).*options.poroMat;

end

%-------------------------------------------------------------------------%
function [G, rock] = addFractures(G, rock, options)

    % Find fracture faces
    N   = G.faces.neighbors + 1;
    ncl = prod(G.cartDims(1:2));
    nl  = ceil(G.cartDims(3)/(options.numFrac+1));
    ix  = reshape(repmat(1:options.numFrac+1, ncl*nl, 1), [], 1);
    ix  = ix(1:G.cells.num);
    ix  = [inf; ix];
    is_frac = diff(ix(N),1,2) == 1;
    % Set faces to fractures
    [G, rock] = setFacesToFractures(G, rock, is_frac, ...
                                  'aperture', options.apertureFrac, ...
                                  'poro'    , options.poroFrac    , ...
                                  'lambdaR' , options.lambdaR     , ...
                                  'CpR'     , options.CpR         , ...
                                  'rhoR'    , options.rhoR        );

end

%-------------------------------------------------------------------------%
function W = setUpWells(G, rock, fluid, options)

    W = [];
    % Add one cold well in each corner
    addColdWell = @(W, ii, jj, name) ...
        verticalWell(W, G, rock, ii, jj, [], ...
                     'Name'  , name         , ...
                     'Radius', 5*centi*meter, ...
                     'type'  , 'bhp'        , ...
                     'val'   , options.bhp  );
    W = addColdWell(W,             1,             1, 'Cold-1-1');
    W = addColdWell(W, G.cartDims(1),             1, 'Cold-2-1');
    W = addColdWell(W, G.cartDims(1), G.cartDims(2), 'Cold-2-2');
    W = addColdWell(W,             1, G.cartDims(2), 'Cold-1-2');
    % Add a hot well in the center
    ii = floor(G.cartDims(1:2)/2) + 1;
    W  = verticalWell(W, G, rock, ii(1), ii(2), [], ...
                    'Name', 'Hot', ...
                    'Radius', 5*centi*meter, ...
                    'type', 'rate', ...
                    'val', options.rateCharge);
    % Set injection temperature
    W = addThermalWellProps(W, G, rock, fluid, 'T', options.tempCharge);
    % Set groups if we use group control
    if options.useGroupCtrl
        [W.group] = deal({'cold', 'cold', 'cold', 'cold', 'hot'});
    end
    
end

%-------------------------------------------------------------------------%
function schedule = setUpSchedule(G, G0, rock, fluid, options)

    W = setUpWells(G0, rock, fluid, options);
    
    if options.dfm
        W = addHybridCellsToWell(W, G, rock);
    end
    
    [W(1:4).type] = deal('bhp');
    [W(1:4).val ] = deal(options.bhp);
    [W(1:4).sign] = deal(-1);
    W(5).type     = 'rate';
    W(5).val      = options.rateCharge;
    W(5).T        = options.tempCharge;
    W(5).sign     = 1;
    
    dtCharge       = rampupTimesteps(options.timeCharge, options.dtCharge);
    scheduleCharge = simpleSchedule(dtCharge, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end
    
    [W(1:4).type] = deal('rate');
    [W(1:4).val ] = deal(options.rateDischarge);
    [W(1:4).T   ] = deal(options.tempDischarge);
    [W(1:4).sign] = deal(1);
    W(5).type     = 'bhp';
    W(5).val      = options.bhp;
    W(5).sign     = -1;
    
    dtDischarge       = rampupTimesteps(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end
    
    if options.chargeOnly
        schedule = scheduleCharge;
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    else
        schedule = combineSchedules(scheduleCharge, scheduleDischarge, 'makeConsistent', false);
        schedule = repmat({schedule}, 1, options.numCycles);
        schedule = combineSchedules(schedule{:}, 'makeConsistent', false);
    end

end

%-------------------------------------------------------------------------%
function state0 = setUpInitialState(model, W, options)
    
    [p,T] = initializeGeothermalEquilibrium(model, ...
                            'datumPressure'   , options.initPres, ...
                            'datumTemperature', options.initTemp);
    state0   = initResSol(model.G, 1, 1);
    state0.pressure = p;
    state0.T        = T;
    % Initialize well solutions
    wellSol = initWellSolAD(W, model, state0);
    [wellSol.bhp] = deal(1*atm);
    state0.wellSol = wellSol;
     
end
%-------------------------------------------------------------------------%

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}