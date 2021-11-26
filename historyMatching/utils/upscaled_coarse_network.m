function [description, options, state0, model, schedule, plotOptions] = upscaled_coarse_network(varargin)
% A rich network model represented by a very coarse grid obtained through
% upscaling of a full reservoir model. See description below.
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions] = ...
%       upscaled_coarse_network('referenceModel', referenceModel)
%
% DESCRIPTION:
%   This function maps a full reservoir model into a rich network model
%   represented by a very coarse computational grid. The full reference model is also required 
%   to be a MRSTExample object.
%
% REQUIRED PARAMETERS:
%   'referenceModel' - A full reservoir model as a MRSTExample object.
% 
% OPTIONAL PARAMETERS:   
%   'plotCoarseModel' - true or false
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
                     'partition', [], ...
                     'wellUpscaleMethod', 'mean', ...
                     'plotCoarseModel',  false);
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
    description = ['Creates a coarse model of ', options.referenceExample.name];

    if nargout <= 2, return; end
    
    % Module dependencies
    require ad-core ad-props ad-blackoil upscaling coarsegrid
    
    %% Create coarse model
    if isempty(options.partition)
        % Just specify an arbitrary partition based on the reference model
        options.partition = max(floor(options.referenceExample.model.G.cartDims./10), [1 1 1]);
    end
    
    % Make coarse grid based on the partition
    blockIx = partitionUI(options.referenceExample.model.G, options.partition);
    blockIx = processPartition(options.referenceExample.model.G, blockIx);
    blockIx = compressPartition(blockIx);
    
   % Perform a simple upscaling to obtain a coarse model
   model = upscaleModelTPFA(options.referenceExample.model, blockIx);
   model.AutoDiffBackend = AutoDiffBackend();
   
    % We want to be able to include rel-perm scaling as tunabale 
    % parameters, so includethese for the coarse model. These parameters 
    % have no effect for theinitial coarse model (they are set equal to the 
    % ones given by the rel-perm curves).
    pts = model.fluid.krPts;
    scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
               'SOWCR', pts.ow(2), 'KRW',  pts.w(4), 'KRO', pts.ow(4)};  
    model = imposeRelpermScaling(model, scaling{:});
    
    % Set tighter tolerance to improve gradient accuracy
    model.toleranceCNV = 1e-6; 
    
    %% Plot comparison of grids from reference model and coarse model
    if options.plotCoarseModel
        figure('position',[100 100 1000 400])
        axes('position',[.02 .05 .48 .9]);
        plotGrid(options.referenceExample.model.G, 'EdgeAlpha',.2); 
        view(174,60);
        title('Fine-scale grid (18553 cells)')
        plotWell(options.referenceExample.model.G, ...
                 options.referenceExample.schedule.control(1).W, ...
                 'Color', 'k', 'FontSize', 10); 
        axis off tight
        camlight headlight

        axes('position',[.5 .05 .48 .9]);
        %plotCellData(cModel.G, cModel.rock.poro, 'EdgeColor', 'none');
        plotGrid(model.G, 'EdgeAlpha',.8);
        title('Coarse-scale grid (33 cells)')
        view(174,60); 
        plotWell(options.referenceExample.model.G, ...
                 options.referenceExample.schedule.control(1).W, ...
                 'Color', 'k', 'FontSize', 10); 
        axis off tight
        camlight headlight
    end
    
    %% Compute state and schedule using upscaling
    state0 = upscaleState(model, ...
                          options.referenceExample.model, ...
                          options.referenceExample.state0);
    schedule = upscaleSchedule(model, ...
                               options.referenceExample.schedule, ...
                               'wellUpscaleMethod', options.wellUpscaleMethod);
    
    % Temporarily store all alternatives for upscaling:
    wellUpscalingMethods = {'sum', 'mean', 'harmonic'};
    for i=1:numel(wellUpscalingMethods)
        tmpSchedule = upscaleSchedule(model, ...
                               options.referenceExample.schedule, ...
                               'wellUpscaleMethod', wellUpscalingMethods{i});
        options.(wellUpscalingMethods{i}) = [tmpSchedule.control(1).W.WI];
    end

    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
end