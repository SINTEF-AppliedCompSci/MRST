function [description, options, state0, model, schedule, plotOptions] = injector_to_producer_network(varargin)
% Injector-producer networks created from a full reservoir model. See
% description below.
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions] = ...
%       injector_to_producer_network('referenceModel', referenceModel, ...
%                                    'cellsPerConnection', cellsPerConnection)
%
% DESCRIPTION:
%   This function maps a full reservoir model into a redced GPSNet-type 
%   network model with a *single* communication path between all pairs of
%   injector and producer wells. The full reference model is also required 
%   to be a MRSTExample object.
%
% REQUIRED PARAMETERS:
%   'referenceModel' - A full reservoir model as a MRSTExample object.
% 
% OPTIONAL PARAMETERS:   
%   'cellsPerConnection' - Number of cells in each connection
%   'p0'   - Initial network pressure
%   's0'   - Initial network saturation
%   'perm' - Permeability in each connection
%   'poro' - Porosity in each connection
%   'plotNetwork' - true or false
%
% RETURNS:
%   description - One-line example description, displayed in list-examples,
%                 and the only input argument if the function is called as
%                 description = injector_to_producer_network(referenceModel);
%
%   options     - A struct of the optional input parameters, with defaults
%                 for all arguments that were not passed as optional
%                 parameters. Returned for convenient access to the example
%                 configuration.
%                 Also includes the NetworkModel object itself
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
    
    % Optional input parameters with default values
    options = struct('referenceExample', {{}}, ...
                     'cellsPerConnection', 10,  ...
                     'p0', 400*barsa, ...
                     's0', [0.2, 0.8], ...
                     'perm', 200*milli*darcy, ...
                     'poro', 0.1, ...
                     'plotNetwork',  false);
    options = merge_options(options, varargin{:});
    
    if isempty(options.referenceExample)
        error("'referenceExample' is a required input to this example but is currently missing");
    end
    if iscell(options.referenceExample)
        options.referenceExample = options.referenceExample{1};
    end
    
    if ~isa(options.referenceExample, 'MRSTExample')
        error("'referenceExample' is required to be a MRSTExample object");
    end
    
    % One line description
    description = ['Creates a injector_to_producer network of ', options.referenceExample.name];

    if nargout <= 2, return; end
    
    % Module dependencies
    require ad-core ad-props ad-blackoil network-models
    
    %% Create network
    % To create the network, we only need a single perforation from each well,
    % which we somewhat arbitrarily pick from the top.
    Wnetwork = options.referenceExample.schedule.control.W;
    for i = 1:numel(Wnetwork)
        Wnetwork(i).cells = Wnetwork(i).cells(end);
    end

    % Creating the network by defining which of the wells are injectors and
    % which are producers
    networkType = 'injectors_to_producers';
    injectors = find([Wnetwork.sign] ==   1);
    producers = find([Wnetwork.sign] ==  -1);
    network = Network(Wnetwork, options.referenceExample.model.G, ...
                      'type', networkType, ...
                      'injectors', injectors, 'producers', producers);

    % Plot the network
    if options.plotNetwork
        figure; network.plotNetwork()
    end
    
    %% Build network model
    % We split each flow path into ten cells and define a computational 2D grid
    % in which each flow path is represented with a single row. We use the same
    % fluid model as in the full reference model.
    % 
    % The grid is set to have an aspect ratio of [5 1 1] and a volum that
    % matches the bulk volume of the reservoir. 
    L = nthroot(sum(options.referenceExample.model.operators.pv./options.referenceExample.model.rock.poro)*25, 3);
    G = cartGrid([options.cellsPerConnection, 1, numedges(network.network)], ...
                 [L, L/5, L/5]*meter^3);
    G = computeGeometry(G);
    rock = makeRock(G, options.perm, options.poro);

    % Fluid model
    gravity off
    fluid = options.referenceExample.model.fluid;
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    model.OutputStateFunctions = {};

    % Map the network onto the MRST model
    networkModel = NetworkModel(model, options.cellsPerConnection, network.network, Wnetwork);
    model = networkModel.model;
    model = model.validateModel();
    W = networkModel.W;
    state0 = initState(G, W, options.p0, options.s0);
    schedule = simpleSchedule(options.referenceExample.schedule.step.val(:), 'W', W);
    
    % Return the NetworkModel through the options struct
    options.networkModel = networkModel;

    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
end