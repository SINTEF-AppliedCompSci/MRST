function [description, options, state0, model, schedule, plotOptions] = benchmark_2d_geothermal(varargin)

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
    description = '';
    options = struct('case'    , 'a'              , ...
                     'physDims', [9, 3]*kilo*meter, ...
                     'cartDims', [20,20]        , ...
                     'dt'      , 250*year    , ...
                     'gravity' , false    );
                     
    options = merge_options(options, varargin{:});
    cartDims = options.cartDims.*options.physDims/(kilo*meter);
    if nargout <= 2, return; end
    require geothermal compositional
    gravity reset on y;
    % Define initial state, model and schedule
    G = computeGeometry(cartGrid(cartDims, options.physDims));

    fluid = initSimpleADIFluid('phases', 'WG', 'n', [1,1], 'mu', [1,1], 'rho', [1,1]);
    fluid.krW = coreyPhaseRelpermAD(1.0, 0.3, 1.0, 0.3);
    fluid.krG = coreyPhaseRelpermAD(1.0, 0.0, 1.0, 0.3);

    fluid = addThermalFluidProps(fluid);
    fluid = setIAPWSIF97Properties(fluid, 'numX', 1000, 'numY', 1000);
    
    fluid.lambdaF = 0;%fluid.lambdaF(5*mega*Pascal, 273.15 + 200);
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

    model = model.validateModel();
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
    
    gpd = model.FlowDiscretization.getStateFunction('GravityPotentialDifference');
    gpd.saturationWeighting = true;
    model.FlowDiscretization = model.FlowDiscretization.setStateFunction('GravityPotentialDifference', gpd);
    
    pte = model.FlowPropertyFunctions.getStateFunction('PhaseThermalEnergy');
    pte.useEnthalpy = true;
    model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('PhaseThermalEnergy', pte);
    
    K0 = 273.15*Kelvin;
    switch options.case
        case 'a'
            p0  = 1*atm;
            pbc = 50*mega*Pascal;
            T0  = K0 + 10*Kelvin;
            Tbc = K0 + 10*Kelvin;
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
    

    p = @(z) 1*atm + fluid.rhoWS.*norm(gravity())*z;
    state0 = initResSol(G, p(G.cells.centroids(:,2)), [0.5, 0.5]);
    state0.T = repmat(T0, G.cells.num, 1);
    state0.components = ones(G.cells.num, 1);
    
    bc = [];
    bc = fluxside(bc, G, 'yMax', 0, 'sat', [0.5, 0.5]);
    bc = pside(bc, G, 'yMin', p0, 'sat', [0.5, 0.5]);

    Watt = joule/second;
    K0 = 273.15*Kelvin;


    hf = nan(numel(bc.face),1);
    T  = nan(numel(bc.face),1);

    bottom = G.faces.centroids(bc.face,2) == options.physDims(2);
    
    T(~bottom) = T0;

    hf(bottom) = 0.05*Watt/(meter^2).*G.faces.areas(bc.face(bottom));

    mid = bottom & abs(G.faces.centroids(bc.face,1) - options.physDims(1)/2) <= 0.5*kilo*meter;    

    hf(mid) = 5*Watt/(meter^2).*G.faces.areas(bc.face(mid));


    bc = addThermalBCProps(bc, 'Hflux', hf, 'T', T);


    
    time = 50000*year;
    dt = rampupTimesteps(time, options.dt);
    
    schedule = simpleSchedule(dt, 'bc', bc);
    
    
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {'YDir', 'reverse'};
end