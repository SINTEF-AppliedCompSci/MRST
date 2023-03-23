
function [description, options, state0, model, schedule, plotOptions, connectionIndices] ...
    = network_model_template(varargin)
% Creates a GPSNET model with two-phase flow between the injectors and
% producers found in the Egg model.
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions, connectionIndices] ...
%        = network_model_template('pn1', pv1, ...)
%
% DESCRIPTION:
%   
%   Network model connecting all injectors to all producors of the Egg 
%   model for two-phase flow. May be used as a stand-alone example 
%   definition, or to construct an instance of `MRSTExample` as 
%   example = MRSTExample('full_egg_network');
%
% OPTIONAL PARAMETERS:
%   'cellsPerConnection'
%   'gpsnetPerm'
%   'gpsnetPoro'
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
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
    description = 'two-phase flow in a fixed (full) network reduced version of the Egg model';
    
    % Each example can have any number of optional input arguments, and
    % must return a (possibly empy) options struct
    options = struct('networkGraph', graph(), ...
                     'ModelGrid',[],...
                     'deleteOldResults', false, ...
                     'cellsPerConnection', 10, ...
                     'gpsnetPoro', 0.2, ...
                     'gpsnetPerm', 200*milli*darcy, ...
                     'plotNetwork', false, ...
                     'realization', 0, ...
                     'fullSchedule', true, ...
                     'fullExampleName', []);
    options = merge_options(options, varargin{:});
    
    
    if nargout <= 2, return; end
    
    
    % Define any module dependencies for the example. The following are
    % typically always needed
    require ad-core ad-props ad-blackoil incomp network-models ensemble diagnostics
    
    %% Define model
    % Create full 3D reservoir based on a full 3D example reservoir 
    fullExample = MRSTExample(options.fullExampleName, ...
                              'realization', options.realization);
    
    % Use a simpler schedule 
    if ~options.fullSchedule
        fullExample.schedule = simpleSchedule(fullExample.schedule.step.val(1:15), 'W', fullExample.schedule.control.W);
    end
                           
    
    % Defining a well structure with only once well cell per well, to impose a
    % Network that has only once well-cell per well.
    W_ref_V2 =fullExample.schedule.control.W;
    for i = 1:numel(W_ref_V2) % TODO maybe this part should be done more systematically.
        W_ref_V2(i).cells = W_ref_V2(i).cells(1);
        num_cells = numel(W_ref_V2(i).cells);
       % W_ref_V2(i).cells=W_ref_V2(i).cells([1 num_cells]);
    end    

    ntwkr =  Network(W_ref_V2,fullExample.model.G,...
                                 'type','injectors_to_producers',...
                                 'injectors',[1:7],...
                                 'producers',[8:18]);   
    
    %% Creating data driven model
    L= nthroot(sum(fullExample.model.operators.pv./fullExample.model.rock.poro)*25,3)  ;                            

    G = cartGrid([options.cellsPerConnection, 1, numedges(ntwkr.network)], ...
                 [L, L/5 ,L/5]*meter^3);
    G = computeGeometry(G);

    fluid =  fullExample.model.fluid;
    rock = makeRock(G, options.gpsnetPerm, options.gpsnetPoro);

    gravity off
    tmp_model = GenericBlackOilModel(G, rock, fluid);
    tmp_model.gas = false;
    tmp_model.OutputStateFunctions = {};

  
    
    %evalc("[model, W, connectionIndices] = createDDmodel_1(tmp_model, options.cellsPerConnection, options.networkGraph, fullExample.schedule.control.W)");
    evalc("obj = NetworkModel(tmp_model,options.cellsPerConnection, ntwkr.network,fullExample.schedule.control.W)");

    model = obj.model;
    W     = obj.W;
    
    connectionIndices.faces = obj.Graph.Edges.Face_Indices;
    connectionIndices.cells = obj.Graph.Edges.Cell_Indices;
    
    model = model.validateModel();
    state0 = initState(model.G, W , 215*barsa,[0.05, 0.95]); 
    ts = fullExample.schedule.step.val;
    schedule = simpleSchedule(ts, 'W', W);

    problem = packSimulationProblem(state0, model, schedule, 'network');
    
    if options.deleteOldResults
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    

    %% Plot options
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
    
    %% Output connectionIndices ++
    %  To be able to run optimization algorithms or history matching on
    %  this model, we need to output the connectionIndices (lists of grid
    %  indices to identify each connection in G).
    %  Since this function is called from inside `MRSTExample`, we can't
    %  add new output variables, but has to attach it to some of those we
    %  have.
    %  Here, we will add the connectionIndices to the options output
    %  parameter.
    
    options.connectionIndices = connectionIndices;
    
end
