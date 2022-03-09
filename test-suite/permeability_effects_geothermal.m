function setup = permeability_effects_geothermal(varargin)
%Test case illustrating temperature-dependent permeability

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
    description = 'Test case illustrating temperature-dependent permeability';
    options = struct('ncells', 100);
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%

    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    K0 = 273.15*Kelvin;
    rockSPE10 = getSPE10rock(40);
    rockSPE10.poro = max(rockSPE10.poro, 0.01);
    
    G = pebiGrid2D(1/options.ncells, [0.45, 1]);
    
    G.nodes.coords = G.nodes.coords*1000;
    G = computeGeometry(G);
    
    perm = reshape(rockSPE10.perm(:,1), [60,220]);
    perm = sampleFromBox(G, perm);
    
    poro = reshape(rockSPE10.poro(:,1), [60,220]);
    poro = sampleFromBox(G, poro);
    poro = max(poro, 0.01);
   
    rock = makeRock(G, perm, poro);
    permMult = 1e-6;
    
    Tmin = K0 + 100*Kelvin;
    Tmax = K0 + 200*Kelvin;
%     tau = @(T) (min(max(T, Tmin), Tmax) - Tmin)./(Tmax - Tmin);
%     rock.perm0 = rock.perm;
%     rock.perm = @(p,T) rock.perm.*((1-tau(T)) + tau(T)*permMult);
    
    rock = addThermalRockProps(rock);
    
    fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'rho', 1, 'mu', 1);
    fluid = addThermalFluidProps(fluid, 'useEOS', true);
    
    model = GeothermalModel(G, rock, fluid);
    model.maximumTemperature = K0 + 275;
    
    bc = [];
    faces = boundaryFaces(G);
    bottom = G.faces.centroids(faces,2) == min(G.faces.centroids(faces,1));
    top    = G.faces.centroids(faces,2) == max(G.faces.centroids(faces,2));
    
    time = 100*year;
    rate = 0.25*sum(model.operators.pv)/time;
    
    bc = addBC(bc, faces(bottom), 'pressure', 200*barsa);
    bc = addBC(bc, faces(top), 'pressure', 10*barsa);
    Tres = Tmin - 50;
    Tbc = [repmat(Tmax, nnz(bottom), 1);
           repmat(Tres     , nnz(top), 1)];
    bc = addThermalBCProps(bc, 'T', Tbc);

    schedule = simpleSchedule(rampupTimesteps(time, time/100), 'bc', bc);
    state0 = initResSol(G, 100*barsa, 1);
    state0.T = repmat(Tres, G.cells.num, 1);
    
    plotOptions = {'PlotBoxAspectRatio', [1,2,1]   , ...
                   'Size'              , [430, 720]};
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,          ...
                      'description', description, ...
                      'options'    , options    , ...
                      'state0'     , state0     , ...
                      'model'      , model      , ...
                      'schedule'   , schedule   , ...
                      'plotOptions', plotOptions);
    %---------------------------------------------------------------------%
    
end