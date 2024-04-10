function setup = convection_geothermal(varargin)
% Setup function demonstrating simulation of geothermal convection
%
% SYNOPSIS:
%   setup = magmatic_intrusion_geothermal('pn1', pv1, ...)
%   setup = magmatic_intrusion_geothermal(fullSetup, 'pn1', pv1, ...)

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

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    % Description
    description = 'Geothermal conduction in a 2000 x 4000 m 2D reservoir';
    
    % Optional input arguments
    options = struct( ...
        'massInjection', false, ...
        'timeHeat'     , inf    ...
    );
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
     %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Define module dependencies
    require ad-core ad-props geothermal compositional

    gravity reset on
    
    % PEBI grid
    L = 1000*meter;
    boundary = [-L, -L; L, -L; L, 3*L; -L, 3*L];
    dx = L/15;
    
    r = 500*meter;
    n = r*2*pi/dx;
    dtheta = 2*pi/n;
    theta = (dtheta/2:dtheta:2*pi)' - pi/2;
    x = r.*[cos(theta), sin(theta)];
    
    rng(20230321);
    G = pebiGrid2D(dx, [L,L], 'polyBdr', boundary, 'faceConstraints', {x});
    G = computeGeometry(G);
    srcCells = find(sqrt(sum(G.cells.centroids.^2,2)) <= r);
    
    n0 = G.cells.num;
    G = makeLayeredGrid(G, 1);
    G.nodes.coords(n0+1:end,3) = G.nodes.coords(n0+1:end,3)*0.1*L;
    G.nodes.coords(:, 2:3) = fliplr(G.nodes.coords(:,2:3));
    G.nodes.coords(:, 3) = -G.nodes.coords(:, 3);
    G.nodes.coords(:, 3) = G.nodes.coords(:, 3) - min(G.nodes.coords(:,3));
    G = computeGeometry(G);
    
    rock = makeRock(G, 100*milli*darcy, 0.4);
    % Add thermal properties
    rock = addThermalRockProps(rock, 'lambdaR', 2*watt/(meter*Kelvin)       , ...
                                     'rhoR'   , 2700*kilogram/meter^3       , ...
                                     'CpR'    , 1000*joule/(kilogram*Kelvin));
    K0   = 273.15*Kelvin;
    % Make fluid
    fluid = initSimpleADIFluid(                     ...
                   'phases', 'W'                  , ... % Water only
                   'mu'    , 0.5*centi*poise      , ... % Viscosity
                   'rho'   , 1000*kilogram/meter^3, ... % Reference density
                   'c'     , 4.4e-10/Pascal       , ... % Compressibility
                   'pRef'  , 1*atm                );    % Ref pressure
    % Add thermal properties
    fluid = addThermalFluidProps(fluid            , ...
                'Cp'     , 4.2*joule/(gram*Kelvin), ... % Heat capacity
                'lambdaF', 0.6*watt/(meter*Kelvin), ... % Thermal cond
                'cT'     , 207e-6/Kelvin          , ... % Thermal expansion
                'TRef'   , K0 + 10*Kelvin         );    % Reference temp
    
    
    model = GeothermalModel(G, rock, fluid);
    model.radiogenicHeatFluxDensity = 2.79*micro*watt/meter^3;
    
    top = find(G.faces.centroids(:,3) == 0);
    bc = addBC([], top, 'pressure', 1*atm, 'sat', 1);
    bc = addThermalBCProps(bc, 'T', convertFromCelcius(20));
    
    timeTot = 3000*year;
    [timeHeat, dt] = deal(min(options.timeHeat, timeTot), 100*year);
    timeRest = timeTot - timeHeat;
    rate = 0;
    if options.massInjection
        rate = 0.1*sum(model.operators.pv(srcCells))/timeHeat;
    end
    src = addSource([], srcCells, rate, 'sat', 1);
    src = addThermalSourceProps(src, 'T', convertFromCelcius(500));
    
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(timeHeat, dt), 'src', src, 'bc', bc);
    if timeRest > 0
        dtRest = rampupTimesteps(timeRest, dt, 0);
        nRest  = numel(dtRest);
        schedule.control(2) = schedule.control;
        schedule.control(2).src = [];
        schedule.step.val     = [schedule.step.val; dtRest];
        schedule.step.control = [schedule.step.control; repmat(2, nRest, 1)];
    end
    
    % Initial state
    state0 = initResSol(G, 1, 1);
    [state0.pressure, state0.T] = initializeGeothermalEquilibrium(model);
    
    % Plotting options
    plotOptions = {
        'PlotBoxAspectRatio', [2, 0.1, 4]   , ...
        'Projection'        , 'orthographic', ...
        'View'              , [0,0]         , ...
        'Size'              , [450, 800]      ...
    };
    
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename, ...
        'description', description, ...
        'options'    , options    , ...
        'state0'     , state0     , ...
        'model'      , model      , ...
        'schedule'   , schedule   , ...
        'plotOptions', plotOptions  ...
    );
    %---------------------------------------------------------------------%

end