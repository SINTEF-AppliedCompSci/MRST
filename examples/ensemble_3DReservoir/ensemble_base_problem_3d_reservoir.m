function [description, options, state0, model, schedule, plotOptions] = ensemble_base_problem_3d_reservoir(varargin)
% Creates a 3D model with two-phase flow between two injectors and two
% producers.
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions] = base_function_1d_reservoir('pn1', pv1, ...)
%
% DESCRIPTION:
%   
%   Very simple example that generates a 3D reservoir with a two-phase flow
%   problem. It has two injectors and two producers, of which one of the
%   injectors are not in a corner. May be used as a stand-alone example 
%   definition, or to construct an instance of `MRSTExample` as 
%   example = MRSTExample('ensemble_base_problem_3D_reservoir');
%   At the time of writing, the main purpose is to use this example for a
%   base in an example ensemble simulation.
%
% OPTIONAL PARAMETERS:
%   This example currently does not take any optional inputs
%
% RETURNS:
%   description - One-line example description, displayed in list-examples,
%                 and the only input argument if the function is called as
%                 description = my_example_wog()
%
%   options     - A struct of the optional input parameters, with defaults
%                 for all arguments that were not passed as optional
%                 parameters. Returned for convenient access to the example
%                 configuration.
%
%   state0, model, schedule - Initial state, model, and simulation schedule
%                             that can be passed to `simulateScheduleAD`
%
%   plotOptions - Cell array on the form {'pn1', pv1, ...} with arguments
%                 that can be used in any of the following ways
%                   - set(myAxis, 'pn1, vn1, ...)
%                   - figure('pn1', vn1, ...)
%                   - plotToolbar(G, state, 'pn1', vn1, ...)
%                 In addition to the standard optional parameters of
%                 `figure`, {'Size', [width, height]} can also be provided,
%                 which `MRSTExample` interprets as
%                 [pos(1:2), [width, height]], where
%                 pos = get(0, 'DefaultFigurePosition')
%
% SEE ALSO:
%   `MRSTExample`, `listExamples`, `exampleSuiteTutorial`

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
    options = struct();
    if nargout <= 2, return; end
    
    
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
    W = verticalWell([], G, rock, 3, 3, [1], ...
        'Type', 'rate', 'Val', 300*meter^3/day, ...
        'Radius', 0.1, 'Name', 'I1', 'Comp_i', [1, 0], 'sign', 1);

    W = verticalWell(W, G, rock, 1, G.cartDims(2), [2], ...
        'Type', 'rate', 'Val', 300*meter^3/day, ...
        'Radius', 0.1, 'Name', 'I2', 'Comp_i', [1, 0], 'sign', 1);

    W = verticalWell(W, G, rock, G.cartDims(1), 1, [3], ...
        'Type', 'bhp', 'Val', 150*barsa, ...
        'Radius', 0.1, 'Name', 'P1', 'Comp_i', [0, 1], 'sign', -1);

    W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [4], ...
        'Type', 'bhp', 'Val', 150*barsa, ...
        'Radius', 0.1, 'Name', 'P2', 'Comp_i', [0, 1], 'sign', -1);

    ts = ones(15, 1)*30*day;
    schedule = simpleSchedule(ts, 'W', W);
    
    %% Initial state (state0)
    initPressure = 200*barsa;
    initSaturation = [0, 1];

    state0 = initState(G, W , initPressure, initSaturation);


    %% Plot options
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
end
