function [description, options, state0, model, schedule, plotOptions] = saigup_wo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

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
    % One-line description
    description = ['SAIGUP Model (Manzocchi et al, 2008, doi: 10.1144/1354-079307-790) ', ...
                   'with water injection into reservoir filled with oil and water'      ];
    % Optional input arguments
    options = struct();
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Read grdecl
    gravity reset on
    grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
    grdecl = readGRDECL(grdecl);
    usys   = getUnitSystem('METRIC');
    grdecl = convertInputUnits(grdecl, usys);
    % Get grid
    G = processGRDECL(grdecl);
    G = computeGeometry(G);
    % Get rock
    rock = grdecl2Rock(grdecl, G.cells.indexMap);
    % Make fluid
    fluid = initSimpleADIFluid('phases', 'WO'                         , ...
                                'mu'   ,  [1,5]*centi*poise           , ...
                                'rho'  ,  [1014, 859]*kilogram/meter^3, ...
                                'n'    ,  [2, 2]                      , ...
                                'c'    ,  [1e-6, 1e-5]/barsa          );
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);    
    % Eight vertical injectors around the perimeter of the model
    nz = G.cartDims(3);
    I = [ 3, 20,  3, 25,  3, 30,  5, 29]; % I-indices
    J = [ 4,  3, 35, 35, 70, 70,113,113]; % J-indices
    R = [ 1,  3,  3,  3,  2,  4,  2,  3]*500*meter^3/day; % Rates
    W = [];
    refDepth = min(G.nodes.coords(:,3));
    for i = 1 : numel(I)
        W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
            'Val', R(i), 'Radius', .1*meter, 'Comp_i', [1 0], ...
            'name', ['I', int2str(i)], 'refDepth', refDepth);
    end
    % Set six vertical producers
    I = [15, 12, 25, 21, 29, 12]; % I-indices
    J = [25, 51, 51, 60, 95, 90]; % J-indices
    bhp = 200*barsa(); % Bottom-hole pressure
    for i = 1 : numel(I)
        W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
            'Val', bhp, 'Radius', .1*meter, ...
            'name', ['P', int2str(i)], 'Comp_i',[0 1], 'refDepth', refDepth, 'sign', -1);
    end
    % Make schedule
    dt = rampupTimesteps(30*year, 100*day);
    schedule = simpleSchedule(dt, 'W', W);
    % Initial state: oil in the shallow regions, water in the rest
    xmax = max(model.G.nodes.coords);
    xmin = min(model.G.nodes.coords);
    dz = xmax(3) - xmin(3);
    ix = G.cells.centroids(:,3) < xmin(3) + 0.5*dz;
    state0 = initResSol(G, 350*barsa, [1,0]);
    state0.s(ix,:) = state0.s(ix,:) + [-1,1];
    % Plotting
    plotOptions  = {'PlotBoxAspectRatio', [1,3,0.5] , ...
                    'View'              , [-70, 30] , ...
                    'Size'              , [800, 400], ...
                    'field'             , 's:1'     };
end