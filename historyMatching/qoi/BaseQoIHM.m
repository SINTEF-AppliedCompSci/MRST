classdef BaseQoIHM 
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    %
    % NOTE:
    %   Not intended for direct use.
    %
    % DESCRIPTION:
    %   This class (and its super classes) is used within a MRSTEnsemble
    %   to extract, store, and work with quantities of interest. 
    %   This base class defines the main API for interacting with QoI's.
    %
    % SEE ALSO:
    %   `WellQoI`, `ReservoirStateQoI`, `MRSTExample`, `BaseSamples`
    
    properties
        
        % Properties related to history matching
        observationResultHandler % Result handler that holds observations.
        truthResultHandler       % Result handler that holds the true data. 
                          % The data read by these result handler must be 
                          % on a format that matches the given QoI class.
        observationCov    % Observation error covariance matrix, several  
                          % forms might be valid depending on the relevant
                          % QoI implementation. This can be either
                          % - scalar: uncorrelated observation with the
                          % same variance
                          % - vector: Uncorrelated observations, but
                          % different values for different parts of the QoI
                          % - matrix: The full error covariance matrix
        
    end
    
    methods
        

        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi)
            % Validate the history matching part of the quantity of 
            % interest. Checks that it is used in a subclass that also
            % inherits from BaseQoI, and validate the properties related to
            % history matching.
            % 
            % This function is called from the relevant QoIHM's validateQoI
            % implementation, after its call to its standard QoI's
            % validateQoI, like this:
            % qoi = validateQoI@BaseQoIHM(qoi);
            %
            % SYNOPSIS:
            %   qoi = qoi.validateQoI()
            %
            % PARAMETERS:
            %   problem - An mrst problem of the same nature as what this
            %             QoI class will be used for.
            
            % Validate that qoi has also inherited BaseQoI by checking for
            % one of the properties defined there
            assert(isprop(qoi, 'ResultHandler'), ...
                'qoi does not seem to be subclass of BaseQoI');
                
            % Validation related to history matching
            if ~isempty(qoi.observationResultHandler)
                assert(~isempty(qoi.observationCov), ...
                    'Observation error covariance matrix must be provided along with observationResultHandler');
            
            end
        end
        
        function qoi = setObservations(qoi, referenceProblem, varargin)
            % Takes a reference problem as input and extracts the QoI from
            % it as observations for history matching. 
            opt = struct('perturb', true);
            [opt, extra] = merge_options(opt, varargin{:});
            
            disp('Simulating the reference problem');
            simulatePackedProblem(referenceProblem);
            
            referenceObservations = qoi.getQoI(referenceProblem);
            
            qoi.truthResultHandler = ResultHandler('dataPrefix', 'trueQoI', ...
                                   'writeToDisk', qoi.ResultHandler.writeToDisk,...
                                   'dataDirectory', qoi.ResultHandler.dataDirectory, ...
                                   'dataFolder', qoi.ResultHandler.dataFolder, ...
                                   'cleardir', false);
            qoi.truthResultHandler{1} = {referenceObservations};
            
            perturbedObservations = referenceObservations;
            if opt.perturb
                perturbedObservations = qoi.perturbQoI(referenceObservations);
            end
            
            qoi.observationResultHandler = ResultHandler('dataPrefix', 'observedQoI', ...
                                   'writeToDisk', qoi.ResultHandler.writeToDisk,...
                                   'dataDirectory', qoi.ResultHandler.dataDirectory, ...
                                   'dataFolder', qoi.ResultHandler.dataFolder, ...
                                   'cleardir', false);
            qoi.observationResultHandler{1} = {perturbedObservations};
            
            
        end
        
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function obs = getObservationVector(qoi, varargin)
            
            % Check that the observation is valid
            assert(~isempty(qoi.observationResultHandler), ...
                'qoi.observationResultHandler is missing');
            assert(numel(qoi.observationResultHandler.getValidIds) > 0, ...
                'No available data in the observationResultHandler');
            
            obs = qoi.observationResultHandler{1};
            obs = qoi.qoi2vector(obs, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function truth = getTrueObservation(qoi, varargin)
            
            % Check that the observation is valid
            assert(~isempty(qoi.truthResultHandler), ...
                'qoi.observationResultHandler is missing');
            assert(numel(qoi.truthResultHandler.getValidIds) > 0, ...
                'No available data in the observationResultHandler');
            
            truth = qoi.truthResultHandler{1};
            truth = qoi.qoi2vector(truth, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function [obs, scaling] = getObservationAndScaling(qoi, varargin)
            % Returns the observation and values that shows how these
            % observations should be scaled for best results in history
            % matching. Optional parameters are:
            %  - 'vectorize': If true, the outputs will be one dimensional 
            %         vectors, or they will be structured as the given QoI.

            error('Template class not meant for direct use!');
        end
            
            
        %-----------------------------------------------------------------%
        function R = getObservationErrorCov(qoi)
            % Returns the full observation error covariance matrix
            error('Template class not meant for direct use!');
        end
            
        
        
        %-----------------------------------------------------------------%
        % Functions related to plotting
        % TODO: How to reuse from BaseQoI here???
        %       For now, reimplemented and defined as sealed so that the
        %       correct version is inherited for *HM classes
        %-----------------------------------------------------------------%
    end
    methods (Sealed = true)
        
        %-----------------------------------------------------------------%
        function h = plotEnsembleQoI(qoi, ensemble, h, varargin)
            % Create a meaningful plot of the ensemble based on the
            % relevant QoI
            %
            % SYNOPSIS:
            %   h = qoi.plotEnsembleQoI(ensemble, h);
            %
            % OPTIONAL PARAMETERS:
            %   ensemble - ensemble of which this QoI object of.
            %   h        - Figure handle
            %   'range'  - Subset of ensemble member IDs, if only parts of
            %              the ensemble is to be plotted.
            %   Extra parameters might depending on the actual QoI and
            %   others acceptable for `plot`.
            
            opt = struct('range'     , inf         , ...
                         'subplots'  , false       , ...
                         'subplotDir', 'horizontal', ...
                         'clearFigure', true       , ...
                         'plotObservation', true   , ...
                         'plotTruth', false        , ...
                         'legend'     , {{}}       , ...
                         'plotWells'  , []         , ...
                         'alreadyOpenFigures', 0   , ...
                         'Position',  get(0, 'DefaultFigurePosition'), ...
                         'savefig', false          , ...
                         'savefolder', ''          );
            [opt, extra] = merge_options(opt, varargin{:});
            [u_mean, u_var, u]  = qoi.getQoIMean(opt.range);
            numQoIs      = numel(u_mean);
            numFieldsThatAreNotFields = 0;
            if isfield(u_mean, 'name')
                numFieldsThatAreNotFields = numFieldsThatAreNotFields + 1;
            end
            if isfield(u_mean, 'cost')
                numFieldsThatAreNotFields = numFieldsThatAreNotFields + 1;
            end
            numSubQoIs   = numel(fieldnames(u_mean)) - numFieldsThatAreNotFields;
            
            % The plotWells parameter is a boolean structure to show which
            % wells and which fields should be plotted. Should follow the
            % same structure as the qoi. If empty, we create an all true
            % structure:
            if isempty(opt.plotWells)
                for w = 1:numQoIs
                    for f = 1:numel(qoi.names)
                        opt.plotWells(w).(qoi.names{f}) = true;
                    end
                end
            end
            % Make similar structure to keep track of figure IDs
            figIds = [];
            figIdCounter = opt.alreadyOpenFigures + 1; 
            for w = 1:numQoIs
                for f = 1:numel(qoi.names)
                    if opt.plotWells(w).(qoi.names{f})
                        figIds(w).(qoi.names{f}) = figIdCounter;
                        figIdCounter = figIdCounter + 1;
                    else
                        figIds(w).(qoi.names{f}) = 0;
                    end
                end
            end
            
            if opt.savefig 
                if strcmp(opt.savefolder, '')
                    opt.savefolder = fullfile(ensemble.mainDirectory(), ...
                                              strcat('history_matching_figures_', ...
                                                     datestr(now,'yyyy_mm_dd_HHMMSS')));
                end
                if ~exist(opt.savefolder, 'dir')
                    mkdir(opt.savefolder);
                end
                fprintf('Saving figures to %s\n', opt.savefolder);    
            end
                        
            plotQoI = @(u, i, k, varargin) qoi.plotQoI(ensemble, u(i).(qoi.names{k}), ...
                    'cellNo', i, 'subCellNo', k, varargin{:});
            
            numSamples = numel(u);
            if nargin < 2, ensemble = []; end
            if nargin < 3 || isempty(h)
                if opt.subplots
                    h = nan(numSubQoIs,1);
                else
                    h = nan(numQoIs*numSubQoIs,1);
                end
            end
            if opt.subplots
                switch opt.subplotDir
                    case 'vertical'
                        nr = numQoIs; nc = 1;
                    case 'horizontal'
                        nr = 1; nc = numQoIs;
                end
            end
            for i = 1:numQoIs
                for k = 1:numSubQoIs
                    if opt.subplots
                        figureId = k;
                    else
                        if ~opt.plotWells(i).(qoi.names{k})
                            continue;
                        end
                        figureId = figIds(i).(qoi.names{k});
                    end
                    if isnan(h(figureId))
                        h(figureId) = qoi.figure(ensemble, 'Position', opt.Position);
                    else
                        set(0, 'CurrentFigure', h(figureId));
                        if ~opt.subplots && opt.clearFigure
                            clf(h(figureId));
                        end
                    end
                    if opt.subplots
                        subplot(nr, nc, i);
                    end
                    hold on
                    for j = 1:numSamples
                        plotQoI(u{j}, i, k, 'isMean', false, extra{:});
                    end
                    plotQoI(u_mean, i, k, extra{:}, 'tag', 'mean');
                    
                    if opt.plotObservation && ~isempty(qoi.observationResultHandler)
                        obs = qoi.getObservationVector('vectorize', false);
                        plotQoI(obs, i, k, 'isObservation', true, 'tag', 'obs', extra{:});
                    end
                    
                    if opt.plotTruth && ~isempty(qoi.truthResultHandler)
                        truth = qoi.getTrueObservation('vectorize', false);
                        plotQoI(truth, i, k, 'isTruth', true, 'tag', 'truth', extra{:}); 
                    end
                    
                    hold off
                    
                    % Stack the lines so that the mean(s) come on top
                    qoi.organizePlots(opt.legend);
                    if opt.savefig
                        qoi.saveFigure(h(figureId), opt.savefolder, ...
                                       qoi.wellNames{i}, qoi.names{k});
                    end
                end
            end
        end
        
         
        %-----------------------------------------------------------------%
        function organizePlots(qoi, legendText) %#ok
            % Ensure that the line representing the mean comes on top, but
            % still under the observations, if any.
            lines = get(gca, 'Children');
            meansID = [];
            obsID = [];
            truthID = [];
            for line=1:numel(lines)
                if strcmp(lines(line).Tag, 'mean')
                    meansID = [meansID, line]; 
                elseif strcmp(lines(line).Tag, 'obs')
                    obsID = [obsID, line];
                elseif strcmp(lines(line).Tag, 'truth')
                    truthID = [truthID, line];
                end
            end
            uistack(lines(meansID), 'top');
            if ~isempty(truthID)
                uistack(lines(truthID), 'top');
            end
            if ~isempty(obsID)
                uistack(lines(obsID), 'top');
            end
            
            
            % Legend is added as
            % legend([chi(obsID(1)); chi(meansID)], {'Obs', 'It N', ..., 'It 2', 'It 1'})

            if ~isempty(legendText)
                if numel(legendText) == numel(meansID)+1 && ~isempty(obsID)
                    legend([lines(obsID(1)); lines(meansID)], ...
                           legendText, 'Location', 'Best');
                elseif numel(legendText) == numel(meansID)+1 && ~isempty(truthID)
                    legend([lines(truthID(1)); lines(meansID)], ...
                           legendText, 'Location', 'Best');
                elseif numel(legendText) == numel(meansID)+2 && ~isempty(obsID) && ~isempty(truthID)
                    legend([lines(obsID(1)); lines(truthID(1)); lines(meansID)], ...
                           legendText, 'Location', 'Best');
                elseif numel(legendText) == numel(meansID)
                    legend(lines(meansID), legendText, 'Location', 'Best');
                else
                    warning('mismatch between number of legendText and elements to name');
                end
            end    
        end
    
        function n = norm(qoi, u)
            % Compute norm n of the quantity of interest u.
            % 
            % SYNOPSIS
            %   n = qoi.norm(u)
            %
            if isstruct(u)
                % We got a full QoI struct, compute norm for each well and
                % each field by calling qoi.norm for each of them
                n = u;
                for i = 1:numel(u)
                    for fn = qoi.names
                        n(i).(fn{1}) = qoi.norm(u(i).(fn{1}));
                    end
                end
                return;
            else
                n = sqrt(sum(u.*u));
            end
        end
        
        
        
    end
    
    methods (Access = protected, Sealed = true)

        
        function saveFigure(qoi, figNr, folderName, wellName, field)
            
            %figHandles = findobj('Type', 'figure')
            %figHandle = figHandles(figNr);
            figHandle = gcf;
            
            fileExt = {'.fig', '.png', '.eps'}; %, '.pdf'};
            formats = {'fig', 'png', 'epsc'}; %, 'pdf'};
            
            filenameBase = strcat(wellName, '_', field);

            for showLegend = 0:1
                legendText = '';
                if showLegend
                    figHandle.CurrentAxes.Legend.Visible = 1;
                    legendText = '_legend';
                else
                    figHandle.CurrentAxes.Legend.Visible = 0;
                end

                for fe = 1:numel(fileExt)
                    savefilename = fullfile(folderName, strcat(filenameBase, legendText, fileExt{fe}))
                    saveas(figHandle, savefilename, formats{fe});
                end

            end      
        end

    end
    
    methods (Access = protected) 
        function perturbedU = perturbQoI(qoi, u)
            error('Template class not meant for direct use!');
        end
    end
end
    
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
