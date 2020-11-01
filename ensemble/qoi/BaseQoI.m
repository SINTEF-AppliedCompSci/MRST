classdef BaseQoI
    % Template class for extracting a quantity of interest from a simulated
    % problem.
    %
    % NOTE:
    %   Not intended for direct use.
    %
    % DESCRIPTION:
    %   This class (and its super classes) is used within a MRSTEnsemble
    %   to extract, store, and work with quantities of interest. 
    %   This base class defines functions 
    %
    % SEE ALSO:
    %   `WellQoI`, `ReservoirStateQoI`, `MRSTExample`, `BaseSamples`
    
    properties
        ResultHandler % Handler for writing/reading QoIs to/from file
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = BaseQoI()
            % Constructor is intentionally empty
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem, varargin)
            % Validate the quantity of interest. BaseQoI sets up an
            % appropriate ResultHandler, so that all subclass
            % implementations of this function should start with
            % qoi = validateQoI@BaseQoI(qoi, problem);
            %
            % SYNOPSIS:
            %   qoi = qoi.validateQoI(problem)
            %
            % PARAMETERS:
            %   problem - An mrst problem of the same nature as what this
            %             QoI class will be used for.
            opt = struct('qoiNo', []);
            opt = merge_options(opt, varargin{:});
            if isempty(qoi.ResultHandler)
                prefix = ['qoi', num2str(opt.qoiNo), '_'];
                % Set up ResultHandler. By default, data is stored in the
                % the dataDirectory of the problem output handler
                dataDir = problem.OutputHandlers.states.dataDirectory;
                qoi.ResultHandler = ResultHandler('dataDirectory', dataDir, ...
                                                  'dataFolder'   , ''     , ...
                                                  'dataPrefix'   ,  prefix); %#ok
            end
            % Check that output is stored with the correct name
            assert(strcmp(qoi.ResultHandler.dataPrefix(1:3), 'qoi'), ...
                   'ResultHandler data prefix must begin with ''qoi''.');
        end
        
        %-----------------------------------------------------------------%
        function u = getQoI(qoi, problem)
            % Get quantity of interest (QoI) for a given problem. The 
            % structure of the relevant QoI is given by the implementation
            % of qoi.computeQoI(problem).
            %
            % SYNOPSIS:
            %   u = qoi.getQoI(problem)
            %
            % PARAMETERS:
            %   problem - The specific problem for which we will extract
            %             relevant quantity of interest. 
            %
            % RETURNS:
            %   u - quantity of interest for the given problem
            %
            % See also:
            %    `qoi.computeQoI(problem)`
            seed = str2double(problem.OutputHandlers.states.dataFolder);
            if qoi.isComputed(seed)
                % QoI already computed - read from file
                u = qoi.ResultHandler{seed};
            else
                % Compute QoI and store to file
                u  = qoi.computeQoI(problem);

                % TODO: Check that problem has been simulated successfully, and
                % issue a warning if it is not
                
                us = u; % Handle special case when u is a cell array
                if iscell(us), us = {us}; end 
                qoi.ResultHandler{seed} = us;
            end 
        end
        
        %-----------------------------------------------------------------%
        function ok = isComputed(qoi, seed)
            % Check if qoi for a given seed it computed
            %
            % SYNOPSIS:
            %   ok = qoi.isComputed(seed)
            ids = qoi.ResultHandler.getValidIds();
            ok  = any(ids == seed);
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem) %#ok
            % Compute quantity of interest for a given problem
            %
            % SYNOPSIS:
            %   u = qoi.computeQoI(problem)
            error('Template class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u) %#ok
            % Compute norm n of the quantity of interest u.
            % 
            % SYNOPSIS
            %   n = qoi.norm(u)
            %
            n = abs(u);
        end
        
        %-----------------------------------------------------------------%
        function [u_mean, u] = computeMean(qoi, range)
            % Computes the mean according to the ensemble given by the
            % range of ensemble member IDs (if any).
            %
            % SYNOPSIS:
            %   u_mean = qoi.computeMean(range)
            %   [u_mean, u] = qoi.computeMean(range)
            %
            % OPTINAL PARAMETERS:
            %   range - A range of ensemble IDs enabling the possibility to
            %           compute the mean only for a part of the ensemble.
            %           If not provided, all available QoI's will be used
            %           to compute the mean.
            %
            % RETURNS:
            %   u_mean - Mean values of the quantity of interest
            %   u      - Cell array of the QoIs for the individual ensemble
            %            members.

            
            % Check if we have computed any QoIs so far
            assert(~isempty(qoi.ResultHandler), 'No QoIs computed yet!')
            ids = qoi.ResultHandler.getValidIds();
            if nargin < 2 || any(isinf(range))
                % No range give, compute mean for all QoIs
                range = ids;
            else
                % Range given, check for computed QoIs
                keep = ismember(range, ids);
                if ~all(keep)
                    warning(['Only a subset of the QoIs in range have '         , ...
                             'been computed (%d %%). Computing mean for subset'], ...
                             nnz(keep)/numel(range)*100                        )
                end
                range = range(keep);
            end
            % Compute mean of all QoIs in range
            u_mean = qoi.ResultHandler{range(1)};
            % Return all QoIs if requested
            if nargout > 1
                u    = cell(numel(range),1);
                u{1} = u_mean;
            end
            for i = 2:numel(range)
                u_tmp = qoi.ResultHandler{range(i)};
                for j = 1:numel(u_tmp)
                    u_mean{j} = computeMean(u_mean{j}, u_tmp{j}, i-1, 1);
                end
                if nargout > 1
                    u{i} = u_tmp;
                end
            end
        end
        
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
            
            opt = struct('range', inf);
            [opt, extra] = merge_options(opt, varargin{:});
            [u_mean, u] = qoi.computeMean(opt.range);
            numQoIs    = numel(u_mean);
            numSamples = numel(u);
            if nargin < 2,               ensemble = [];      end
            if nargin < 3 || isempty(h), h = nan(numQoIs,1); end
            for i = 1:numQoIs
                if isnan(h(i))
                    if ~isempty(ensemble)
                        h(i) = ensemble.setup.figure();
                    else
                        h(i) = figure();
                    end
                else
                    set(0, 'CurrentFigure', h(i)), clf(h(i))
                end
                hold on
                for j = 1:numSamples
                    qoi.plotQoI(ensemble, u{j}{i}, 'isMean', false, extra{:});
                end
                qoi.plotQoI(ensemble, u_mean{i}, extra{:});
                hold off
            end
        end
        
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin)
            % Plot a single QoI u in current figure. This function is meant
            % to be implemented on a per-class-basis to generate suitable
            % plots for the QoI in question.
            warning('BaseQoI:notImplemented', ['Method plotQoI is not '    , ...
                    'implemented for this QoI - using simple 1d plotting. ', ...
                    'This warning message will be turned off for '         , ...
                    'subsequent calls.']                                   );
            warning('off', 'BaseQoI:notImplemented');
            % Optional input arguments. Can be used to pass arguments
            % directly to e.g., plot or plotCellData
            opt = struct('isMean', true);
            [opt, extra] = merge_options(opt, varargin{:});
            color = [1,1,1]*0.8*(1-opt.isMean); % Plot mean in distinct color
            plot(u, 'lineWidth', 2, 'color', color, extra{:});
        end
        
        %-----------------------------------------------------------------%
        function h = plotQoIHistogram(qoi, edges, varargin)
            % Plots the distribution of the QoI of the ensemble in the form
            % of a histogram. If the QoI is nonscalar, norm(QoI) is used.
            %
            % SYNOPSIS:
            %   h = plotQoIHistogram(egdes, ...)
            %
            if nargin < 2, edges = 10; end
            opt = struct('range'      , inf  , ...
                         'log10'      , false, ...
                         'includeMean', false, ...
                         'includeRMSE', false);
            [opt, extra] = merge_options(opt, varargin{:});
            % Get QoIs and mean
            [u_mean, u] = qoi.computeMean(opt.range);
            numQoIs     = numel(u_mean);
            % Compute norm
            n_mean = qoi.norm(u_mean);
            n      = cell2mat(cellfun(@(u) qoi.norm(u), u, 'UniformOutput', false));
            if opt.log10
                % Logarithmic transformation
                n_mean = log10(abs(n_mean));
                n      = log10(abs(n));
            end
            for i = 1:numQoIs
                % Plot each QoI in separate figure
                h = histogram(n(:,i), edges, extra{:});
                if opt.includeMean
                    % Plot mean as vertical, dashed line
                    hold on
                    plot([n_mean(i), n_mean(i)], h.Parent.YLim, '--k', 'lineWidth', 1);
                    hold off
                end
                if opt.includeRMSE
                    error('Not implemented yet')
                end
            end
        end
            
    end
end
    
%{
#COPYRIGHT#
%}