function setup = ensemble_base_problem_3d_reservoir(varargin)
% Creates a 3D model with two-phase flow between two injectors and two
% producers.
%
% SYNOPSIS:
%   setup = ensemble_base_problem_3d_reservoir('pn1', pv1, ...)
%   setup = ensemble_base_problem_3d_reservoir(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   
%   Very simple example that generates a 3D reservoir with a two-phase flow
%   problem. It has two injectors and two producers, of which one of the
%   injectors are not in a corner. May be used as a stand-alone example 
%   definition, or to construct an instance of `TestCase` as 
%   example = TestCase('ensemble_base_problem_3D_reservoir');
%   At the time of writing, the main purpose is to use this example for a
%   base in an example ensemble simulation.
%
% OPTIONAL PARAMETERS:
%   This example currently does not take any optional inputs
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%           options, state0, model, schedule, and plotOptions.
%           If the optional input fullSetup (see synopsis) is false, the
%           returned setup only contains name, description, and options.
%
% SEE ALSO:
%   `TestCase`, `testcase_template`, `testSuiteTutorial`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    % Each example must start with the description and options, followed by
    % an nargout check that returns if we only asked for the description
    % and options
    description = 'two-phase flow in fixed 3D reservoir between two injectors and two producers';
    
    % Each example can have any number of optional input arguments, and
    % must return a (possibly empy) options struct
    options = struct('longWells', false, ...
                     'randomSchedule', false);
 
                 
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end    
   
    
    % Define any module dependencies for the example. The following are
    % typically always needed
    require ad-core ad-props ad-blackoil incomp
    
    %% Define model
    length = 275; 
    G = cartGrid([10, 6, 4], [length, length/5, length/5]*meter^3);
    G = computeGeometry(G);
    
    % Define rock and fluid
    porosity = gaussianField(G.cartDims, [0.2 0.4]); 
    permeability = porosity.^3.*(1e-5)^2./(0.81*72*(1-porosity).^2);
    rock = makeRock(G, permeability(:), porosity(:));
    
    % Two-phase model
    fluid = initSimpleADIFluid('phases', 'wo', ...
                               'mu',     [1, 5]*centi*poise, ...
                               'rho',    [1000, 700]*kilogram/meter^3, ...
                               'n',      [2, 2], ...
                               'cR',     1e-8/barsa, ...
                               'sMin',   [.2, .2]);

    % Generic black oil model
    gravity off
    model = GenericBlackOilModel(G, rock, fluid);
    model.gas = false;
    model.OutputStateFunctions = {};
        
    %% Schedule
    % Add wells
    wellDepths = {1, 2, 3, 4};
    if options.longWells
        wellDepths = {1, (2:3), (1:3), (1:4)};
    end
    W = verticalWell([], G, rock, 3, 3, wellDepths{1}, ...
        'Type', 'rate', 'Val', 300*meter^3/day, ...
        'Radius', 0.1, 'Name', 'I1', 'Comp_i', [1, 0], 'sign', 1);

    W = verticalWell(W, G, rock, 1, G.cartDims(2), wellDepths{2}, ...
        'Type', 'rate', 'Val', 300*meter^3/day, ...
        'Radius', 0.1, 'Name', 'I2', 'Comp_i', [1, 0], 'sign', 1);

    W = verticalWell(W, G, rock, G.cartDims(1), 1, wellDepths{3}, ...
        'Type', 'bhp', 'Val', 150*barsa, ...
        'Radius', 0.1, 'Name', 'P1', 'Comp_i', [0, 1], 'sign', -1);

    W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), wellDepths{4}, ...
        'Type', 'bhp', 'Val', 150*barsa, ...
        'Radius', 0.1, 'Name', 'P2', 'Comp_i', [0, 1], 'sign', -1);

    
    %% Standard or random schedule
    ts = ones(20, 1)*30*day;
    schedule = simpleSchedule(ts, 'W', W);

    if options.randomSchedule
        schedule.step.control= ceil(0.25*(1:numel(schedule.step.val))');
        for n=2:max(schedule.step.control)
            schedule.control(n) = schedule.control(1);
            for i=1:numel(schedule.control(n).W)
                W = schedule.control(n).W(i);
                switch W.type
                    case 'rate'
                        W.val = (.75 + .5*rand)*W.val;
                    case 'bhp'
                        %if rand < 0.2
                        %    W.status = false;
                        % else
                            W.val = (.95 + 0.1*rand)*W.val;
                        %end
                end
                schedule.control(n).W(i) = W;
            end
        end
    end
    
    %% Initial state (state0)
    initPressure = 200*barsa;
    initSaturation = [0, 1];

    state0 = initState(G, W , initPressure, initSaturation);


    %% Plot options
    % plotOptions are only by TestCase. In case of empty plotOptions,
    % TestCase will attempt to set reasonable defaults
    plotOptions = {};
    
     %% Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);

end
