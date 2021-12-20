function [description, options, state0, model, schedule, plotOptions] = pinn_da_case_two_producers(varargin)
% A simple 2D square water-oil reservoir with injector and producer in
% opposing corners 
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions] = ...
% pinn_da_case_const_perm('pn1', pv1, ...)
%
% DESCRIPTION:
% A simple 2D square water-oil reservoir with injector and producer in
% opposing corners. With injector at [-1 -1] and producer at [1 1].
%
%
% OPTIONAL PARAMETERS:
%   The example function may take in any number of optional input
%   arguments to facilitate different configurations of the same example,
%   e.g., number of cells, number of timesteps, and so on:
%   Example: my_example_wog('ncells', 2020, 'nsteps', 5);
%   Suggested standard names:
%   
%   'ncells' - Number of cells in each axial direction
%   'nsteps' - Number of timesteps
%   'tend'   - Total simulation time
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
    % Each example must start with the description and options, followed by
    % an nargout check that returns if we only asked for the description
    % and options
    description = 'PINN DA case - two producers';
    % Each example can have any number of optional input arguments, and
    % must return a (possibly empy) options struct
    options = struct('ncells', 50, ...
                     'nsteps', 100, ...
                     'tend', 0.7);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define any module dependencies for the example. The following are
    % typically always needed
    require ad-core ad-props ad-blackoil
    % Define initial state, model and schedule
    
    G = cartGrid([options.ncells, options.ncells], [2 2]);
    G.nodes.coords = G.nodes.coords - 1;
    G = computeGeometry(G);
    
    p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);
    K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
    %K = gaussianField(G.cartDims, [0.1 0.8], [9 9], 2.5);
    
    rock = makeRock(G, K(:), p(:));
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1 1], 'mu', [1 1], 'rho', [1 1]);
      
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    
    dt = repmat(options.tend/options.nsteps, options.nsteps, 1);
    
    W = verticalWell([], G, rock, floor(G.cartDims(1)/2), 1, 1, ...
        'Type', 'rate', 'Val', 1, ...
        'Radius', 2/(options.ncells*10), ...
        'Name', 'I1', 'Comp_i', [1, 0], 'sign', 1);

    % Producer (with pressure control)
    W = verticalWell(W, G, rock, 1, G.cartDims(2), 1, ...
        'Type', 'bhp', 'Val', 860, ...
        'Radius', 2/(options.ncells*10), ...
        'Name', 'P1', 'Comp_i', [0, 1], 'sign', -1);
    W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), 1, ...
        'Type', 'bhp', 'Val', 860, ...
        'Radius', 2/(options.ncells*10), ...
        'Name', 'P2', 'Comp_i', [0, 1], 'sign', -1);


    schedule = simpleSchedule(dt, 'W', W);
    
    state0 = initResSol(G, 1, [0 1]);
    
    for i = 1:numel(schedule.control.W)
        schedule.control.W(i).WI = schedule.control.W(i).WI; % * 100; 
    end
    
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
end