function setup = upscaled_coarse_network(varargin)
% A rich network model represented by a very coarse grid obtained through
% upscaling of a full reservoir model. See description below.
%
% SYNOPSIS:
%   setup = upscaled_coarse_network('referenceCase', referenceCase)
%   setup = upscaled_coarse_network(fullSetup, 'referenceCase', referenceCase)
%
% DESCRIPTION:
%   This function maps a full reservoir model into a rich network model
%   represented by a very coarse computational grid. The full reference 
%   model is also required to be a TestCase object.
%
% REQUIRED PARAMETERS:
%   'referenceCase' - A full reservoir model as a TestCase object.
% 
% OPTIONAL PARAMETERS:   
%   'plotCoarseModel' - true or false
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
%      Note that setup.options also includes the NetworkModel object itself
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

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
                     'partition', [], ...
                     'wellUpscaleMethod', 'mean', ...
                     'plotCoarseModel',  false);

    description = ['Creates a coarse network model of model X'];
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
    setup.description = ['Creates a coarse network model of ', options.referenceCase.name];

    if ~fullSetup, return; end
    
    % Module dependencies
    require ad-core ad-props ad-blackoil upscaling coarsegrid
    
    %% Create coarse model
    if isempty(options.partition)
        % Just specify an arbitrary partition based on the reference model
        options.partition = max(floor(options.referenceCase.model.G.cartDims./10), [1 1 1]);
    end
    
    % Make coarse grid based on the partition
    blockIx = partitionUI(options.referenceCase.model.G, options.partition);
    blockIx = processPartition(options.referenceCase.model.G, blockIx);
    blockIx = compressPartition(blockIx);
    
   % Perform a simple upscaling to obtain a coarse model
   model = upscaleModelTPFA(options.referenceCase.model, blockIx);
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
        plotGrid(options.referenceCase.model.G, 'EdgeAlpha',.2); 
        view(174,60);
        title(['Fine-scale grid (', num2str(options.referenceCase.model.G.cells.num), ' cells)'])
        plotWell(options.referenceCase.model.G, ...
                 options.referenceCase.schedule.control(1).W, ...
                 'Color', 'k', 'FontSize', 10); 
        axis off tight
        camlight headlight

        axes('position',[.5 .05 .48 .9]);
        %plotCellData(cModel.G, cModel.rock.poro, 'EdgeColor', 'none');
        plotGrid(model.G, 'EdgeAlpha',.8);
        title(['Coarse-scale grid (', num2str(model.G.cells.num), ' cells)'])
        view(174,60); 
        plotWell(options.referenceCase.model.G, ...
                 options.referenceCase.schedule.control(1).W, ...
                 'Color', 'k', 'FontSize', 10); 
        axis off tight
        camlight headlight
    end
    
    %% Compute state and schedule using upscaling
    state0 = upscaleState(model, ...
                          options.referenceCase.model, ...
                          options.referenceCase.state0);
    schedule = upscaleSchedule(model, ...
                               options.referenceCase.schedule, ...
                               'wellUpscaleMethod', options.wellUpscaleMethod);
    
    % Temporarily store all alternatives for upscaling:
    wellUpscalingMethods = {'sum', 'mean', 'harmonic'};
    for i=1:numel(wellUpscalingMethods)
        tmpSchedule = upscaleSchedule(model, ...
                               options.referenceCase.schedule, ...
                               'wellUpscaleMethod', wellUpscalingMethods{i});
        options.(wellUpscalingMethods{i}) = [tmpSchedule.control(1).W.WI];
    end

    % plotOptions are only by TestCase. In case of empty plotOptions,
    % TestCase will attempt to set reasonable defaults
    plotOptions = {};
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', setup.description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end