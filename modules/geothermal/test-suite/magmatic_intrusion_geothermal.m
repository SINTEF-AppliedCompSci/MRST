function setup = magmatic_intrusion_geothermal(varargin)
% Setup function for a quarter five-spot water-oil model
%
% SYNOPSIS:
%   setup = magmatic_intrusion_geothermal('pn1', pv1, ...)
%   setup = magmatic_intrusion_geothermal(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%
% RETURNS:
%
% SEE ALSO:

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
    % One-line description
    description = '';
    
    % Optional input arguments
    options = struct( ...
    );
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil

    % Cartesian grid covering a 1000 x 1000 m domain
    
    L = 1000*meter;
    boundary = [-L, -L; L, -L; L, 3*L; -L, 3*L];
    dx = L/15;
    
    r = 500*meter;
    n = r*2*pi/dx;
    dtheta = 2*pi/n;
    theta = (dtheta/2:dtheta:2*pi)';
    x = r.*[cos(theta), sin(theta)];
    x = [x; x(1,:)];
    
    G = pebiGrid2D(dx, [L,L], 'polyBdr', boundary, 'faceConstraints', {x});
    
    G     = computeGeometry(cartGrid([1,1].*options.ncells, [1000, 1000]*meter));
    rock  = makeRock(G, 100*milli*darcy, 0.4); % Rock with 100 md perm and 0.4 poro
    fluid = initSimpleADIFluid('phases', 'WO'             , ... % Water and oil
                               'n'     , [1,1]*options.nkr, ... % Relperm exponents
                               'mu'    , [1,1]*centi*poise, ... % Viscosity
                               'rho'   , [1,1]            );    % Density (no gravity)

    % Construct two-phase model with oil and water
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    
    % Wells
    rate = options.pvi*sum(poreVolume(G, rock))/options.time;
    bhp  = 50*barsa;
    W = [];
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'comp_i', [1,0], 'sign', 1);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', bhp, 'comp_i', [1,0], 'sign', -1);
    
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(options.time, options.dt), 'W', W);
    
    % Initial state
    state0 = initResSol(G, 100*barsa, [0,1]);
    
    % Pack setup
    setup = packTestCaseSetup(mfilename, ...
        'description', description, ...
        'options'    , options    , ...
        'state0'     , state0     , ...
        'model'      , model      , ...
        'schedule'   , schedule     ...
    );

end