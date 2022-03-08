function setup = benchmark_1d_geothermal(varargin)
%Test case describing four 1D geothermal benchmarks from Weis et al (2014)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    description = ['Test case describing four 1D geothermal '       , ...
        'benchmarks from Weis et al. "Hydrothermal, multiphase '    , ...
        'convection of H2O-NaCl fluids from ambient to magmatic '   , ...
        'temperatures: A new numerical scheme and benchmarks for '  , ...
        'code comparison." Geofluids (2014). DOI: 10.1111/gfl.12080'];
        
    options = struct('case'   , 'a'       , ...
                     'length' , 2000*meter, ...
                     'dx'     , 10*meter  , ...
                     'dt'     , 5*year    , ...
                     'gravity', false    );
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                    options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    require geothermal compositional
    if options.gravity
        gravity reset x on;
        g = gravity(); gravity(-g);
    else
        gravity reset off;
    end
    % Define initial state, model and schedule
    G = computeGeometry(cartGrid([round(options.length/options.dx), 1], ...
                                 [options.length, options.dx]         ));
    fluid = initSimpleADIFluid('phases', 'WG', 'n', [1,1], 'mu', [1,1], 'rho', [1,1]);
    fluid.krW = coreyPhaseRelpermAD(1, 0.3, 1, 0.3);
    fluid.krG = coreyPhaseRelpermAD(1, 0  , 1, 0.3);
    fluid = addThermalFluidProps(fluid);
    fluid = setIAPWSIF97Properties(fluid);
    
    fluid.lambdaF = fluid.lambdaF(5*mega*Pascal, 273.15 + 200);
    fluid.CpW = 4.2*joule/(Kelvin*gram);
    fluid.CpG = 4.2*joule/(Kelvin*gram);
    fluid.Cp = 4.2*joule/(Kelvin*gram);

    rock = makeRock(G, 1e-15, 0.1);
    Watt = joule/second;
    rock = addThermalRockProps(rock                                  , ...
                               'CpR'    , 880*joule/(kilogram*Kelvin), ...
                               'lambdaR', 2*Watt/(meter*Kelvin)      , ...
                               'rhoR'   , 2700*kilogram/meter^3      );
    model = GeothermalModel(G, rock, fluid);
    model.thermalFormulation = 'enthalpy';
    
    K0 = 273.15*Kelvin;
    switch options.case
        case 'a'
            p0  = 25*mega*Pascal;
            pbc = 50*mega*Pascal;
            T0  = K0 + 150*Kelvin;
            Tbc = K0 + 350*Kelvin;
            if options.gravity, t = 750*year; else, t = 250*year; end
        case 'b'
            p0  = 20*mega*Pascal;
            pbc = 40*mega*Pascal;
            T0  = K0 + 300*Kelvin;
            Tbc = K0 + 450*Kelvin;
            if options.gravity, t = 350*year; else, t = 120*year; end
        case 'c'
            p0  = 1*mega*Pascal;
            pbc = 15*mega*Pascal;
            T0  = K0 + 350*Kelvin;
            Tbc = K0 + 500*Kelvin;
            t   = 1500*year;
        case 'd'
            p0  = 1*mega*Pascal;
            pbc = 20*mega*Pascal;
            T0  = K0 + 150*Kelvin;
            Tbc = K0 + 400*Kelvin;
            if options.gravity, t = 1000*year; else, t = 200*year; end
        case 'e'
            assert(~options.gravity);
            p0  = 1*mega*Pascal;
            pbc = 4*mega*Pascal;
            T0  = K0 + 150*Kelvin;
            Tbc = K0 + 300*Kelvin;
            t   = 2000*year;
    end
    
    state0 = initResSol(G, p0, [0.5, 0.5]);
    state0.T = repmat(T0, G.cells.num, 1);
    state0.components = ones(G.cells.num, 1);
        
    bc = [];
    bc = pside(bc, G, 'left' , pbc, 'sat', [0.5, 0.5]);
    bc = pside(bc, G, 'right', p0, 'sat', [0.5, 0.5]);
    bc = addThermalBCProps(bc, 'T', [Tbc; T0]);
    
    time = 1500*year;
    dt = rampupTimesteps(time, options.dt);
    
    schedule = simpleSchedule(dt, 'bc', bc);
    
    options.tplot = t;
    
    plotOptions = {'plot1d', true};
    
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