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
                         'legend'     , {{}} );
            [opt, extra] = merge_options(opt, varargin{:});
            [u_mean, u]  = qoi.computeMean(opt.range);
            numQoIs      = numel(u_mean);
            numSubQoIs   = 1;
                        
            plotQoI = @(u, i, k, varargin) qoi.plotQoI(ensemble, u{i}, ...
                'cellNo', i, varargin{:});
            
            if iscell(u_mean{1})
                numSubQoIs = numel(u_mean{1});
                plotQoI = @(u, i, k, varargin) qoi.plotQoI(ensemble, u{i}{k}, ...
                    'cellNo', i, 'subCellNo', k, varargin{:});
            end
            
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
                        figureId = (i-1)*numSubQoIs + k;
                    end
                    if isnan(h(figureId))
                        h(figureId) = qoi.figure(ensemble);
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
            
        
        
        
    end
end
    
%{
#COPYRIGHT#
%}