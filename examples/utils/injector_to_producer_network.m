function setup = injector_to_producer_network(varargin)
% Injector-producer networks created from a full reservoir model. See
% description below.
%
% SYNOPSIS:
%   setup = injector_to_producer_network('referenceCase', referenceCase)
%   setup = injector_to_producer_network(fullSetup, 'referenceCase', referenceCase)
%
% DESCRIPTION:
%   This function maps a full reservoir model into a reduced GPSNet-type 
%   network model with a *single* communication path between all pairs of
%   injector and producer wells. The full reference case is also required 
%   to be a TestCase object.
%
% REQUIRED PARAMETERS:
%   'referenceCase' - A full reservoir model as a TestCase object.
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
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
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
    
    % Optional input parameters with default values
    options = struct('referenceCase', {{}}, ...
                     'cellsPerConnection', 10,  ...
                     'p0', 400*barsa, ...
                     's0', [0.2, 0.8], ...
                     'perm', 200*milli*darcy, ...
                     'poro', 0.1, ...
                     'plotNetwork',  false, ...
                     'plottype', 'default');
    description = ['Creates an injector_to_producer network of model X'];
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    
    if isempty(options.referenceCase)
        error("'referenceCase' is a required input to this example but is currently missing");
    end
    if iscell(options.referenceCase)
        options.referenceCase = options.referenceCase{1};
    end
    
    if ~isa(options.referenceCase, 'TestCase')
        error("'referenceCase' is required to be a TestCase object");
    end
    
    % One line description
    setup.description = ['Creates a injector_to_producer network of ', options.referenceCase.name];

    if ~fullSetup, return; end
    
    % Module dependencies
    require ad-core ad-props ad-blackoil network-models
    
    %% Create network
    % To create the network, we only need a single perforation from each well,
    % which we somewhat arbitrarily pick from the top.
    Wnetwork = options.referenceCase.schedule.control.W;
    for i = 1:numel(Wnetwork)
        Wnetwork(i).cells = Wnetwork(i).cells(end);
    end

    % Creating the network by defining which of the wells are injectors and
    % which are producers
    networkType = 'injectors_to_producers';
    injectors = find([Wnetwork.sign] ==   1);
    producers = find([Wnetwork.sign] ==  -1);
    network = Network(Wnetwork, options.referenceCase.model.G, ...
                      'type', networkType, ...
                      'injectors', injectors, 'producers', producers);

    % Plot the network
    if options.plotNetwork
        figure; network.plotNetwork(options.plottype)
    end
    
    %% Build network model
    % We split each flow path into ten cells and define a computational 2D grid
    % in which each flow path is represented with a single row. We use the same
    % fluid model as in the full reference model.
    
    
    % Creat the network model
    networkModel = GPSNet(options.referenceCase.model, network, ...
                          options.referenceCase.schedule.control.W, ...
                          options.cellsPerConnection, ...
                          'p0', options.p0, ...
                          'S0', options.s0);
    
    % Create schedule for the network model
    networkSetup = gpsNetSimulationSetup(networkModel, options.referenceCase.schedule);
        
    % Return the GPSNet network model through the options struct
    options.networkModel = networkModel;

    % plotOptions are only by TestCase. In case of empty plotOptions,
    % TestCase will attempt to set reasonable defaults
    plotOptions = {};
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', setup.description, ...
                              'options'    , options, ...
                              'state0'     , networkSetup.state0, ...
                              'model'      , networkSetup.model, ...
                              'schedule'   , networkSetup.schedule, ...
                              'plotOptions', plotOptions);
end