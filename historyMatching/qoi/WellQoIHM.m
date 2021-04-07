classdef WellQoIHM < BaseQoIHM & WellQoI
    % Class for extracting well related quantity of interests from
    % (an ensemble of) simulated problems.
    %
    % DESCRIPTION:
    %   This class is used within a MRSTEnsemble to extract, store, and 
    %   work with quantities of interest related to well solutions.
    %
    % SYNOPSIS
    %   qoi = WellQoI('wellNames', {'P1', 'P2'}, ...);
    %   qoi = WellQoI('wellIndices', [3 4], ...);
    %
    % PARAMETERS
    %   'wellNames' - Names of the wells that we are interested in
    %   'wellIndices' - Indices of the wells that we are interested in
    % 
    %    Either 'wellNames' or 'wellIndices' must be provided. If both are
    %    provided, 'wellIndices' is ignored.
    %
    % Properties related to history matching
    %   'observationResultHandler' - Result handler that holds 
    %                                observations.
    %   'truthResultHandler'       - Result handler that holds the true 
    %                                data. 
    %       The data read by these result handler must be on a format that 
    %       matches the given QoI class.
    %
    %    observationCov - Observation error covariance matrix, several  
    %                     forms might be valid depending on the relevant
    %                     QoI implementation. This can be either
    %                     - scalar: uncorrelated observation with the
    %                       same variance
    %                     - vector: Uncorrelated observations with each
    %                     element refering to each fieldname, fieldname
    %                     x well, or fieldname x well x timestep, depending
    %                     on the size of the vector.
    %                     - matrix: The full error covariance matrix
    %
    % OPTIONAL PARAMETERS
    %   'fldname' - cell array of the well output fields that we are
    %               interested in. Valid values are well solution field
    %               names. Default: {'qOs'}
    %   'cumulative' - Quantity of interest is cumulative production data
    %                  per well. Default: false
    %   'total' - Quantity of interest is total production data (scalar per
    %             field and well). Default: false
    %   'combined' - Quantity of interest is combined production data over
    %                the given wells. Default: false
    %
    % 
    % SEE ALSO:
    %   `BaseQoI`, `ReservoirStateQoI`, `MRSTExample`, `BaseSamples`
    
    properties

        historyMatchDtRange
                
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = WellQoIHM(varargin)
            qoi = qoi@WellQoI(varargin{:});
            %qoi = merge_options(qoi, varargin{:});
            % Check input
            %assert(xor(isempty(qoi.wellIndices), isempty(qoi.wellNames)), ...
            %    'Please provide either well indices or well names');
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Check that the configs that are inserted to the constructor
            % makes sense for the base problem for the ensemble, and
            % updates remaining fields.
            %
            % SYNOPSIS:
            %   qoi = qoi.validateQoI(problem)
            %
            % PARAMETERS:
            %   problem - MRST problem that is representative for the
            %             (ensemble member) problems that will be used with
            %             this QoI instance. Especially, the well setups
            %             should be representative, as well as the
            %             production schedule.
            %
            % NOTE:
            %   When used in an MRSTEnsemble, this function is called by
            %   the ensemble constructor.
            
            qoi = validateQoI@WellQoI(qoi, problem);
            
            qoi = validateQoI@BaseQoIHM(qoi);
            
            % Get timesteps
            qoi.dt = problem.SimulatorSetup.schedule.step.val;
            if isempty(qoi.historyMatchDtRange)
                qoi.historyMatchDtRange = 1:numel(qoi.dt);
            end
        end
        
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin) %#ok
            % Plot a single well QoI u in current figure.
            opt = struct('color'         , [0,0,0], ...
                         'lineWidth'     , 2      , ...
                         'alpha'         , 0.8    , ...
                         'isMean'        , true   , ...
                         'timescale'     , day    , ...
                         'labels'        , true   , ...
                         'title'         , true   , ...
                         'cellNo'        , 1      , ...
                         'subCellNo'     , 1      , ...
                         'isObservation' , false  , ...
                         'isTruth'       , false  , ...
                         'observationIndices', [] );
            [opt, extra] = merge_options(opt, varargin{:});
            
            color = opt.color; % Plot mean in distinct color
            if ~opt.isMean
                color = color.*(1-opt.alpha) + opt.alpha;
            end
            
            % Invert rate so that production is plotted as positive values.
            plotScale = 1;
            if strcmpi(qoi.fldname{opt.subCellNo}, 'qOs') || ...
               strcmpi(qoi.fldname{opt.subCellNo}, 'qWs')
                plotScale = -1;
            end
            
            is_timeseries = true;
            if is_timeseries
                
                time = cumsum(qoi.dt)./opt.timescale;
                if opt.isObservation
                    if isempty(opt.observationIndices)
                        plot(time, u*opt.scale, ...,
                             'x', 'color', [0 0 0], extra{:});
                    else
                        plot(time(opt.observationIndices), u(opt.observationIndices)*plotScale, ... ...
                            'x', 'color', [0 0 0], extra{:});
                        unobservedIndices = setdiff(1:numel(u),opt.observationIndices);
                        if ~isempty(unobservedIndices)
                            plot(time(unobservedIndices), u(unobservedIndices)*plotScale, ...
                                 'o', 'color', [0 0 0], extra{:});
                        end
                    end
                elseif opt.isTruth
                    plot(time, u*plotScale, ...
                         'o', 'color', [1 1 1].*0.5, extra{:});
                else
                    plot(time(1:numel(u)), u*plotScale, ...
                         'color'    , color, ...
                         'lineWidth', opt.lineWidth, ...
                         extra{:}         );
                end
                
                xlim([time(1), time(end)]);
                box on, grid on
                if opt.title
                    if qoi.combined
                        title(sprintf('Combined produced %s', qoi.fldname{opt.subCellNo}));
                    else
                        title(sprintf('%s for well %s', qoi.fldname{opt.subCellNo}, qoi.wellNames{opt.cellNo}));
                    end
                end
                if opt.labels
                    xlabel(sprintf('Time (%s)', formatTime(opt.timescale)));
                    ylabel(sprintf('%s', qoi.fldname{opt.subCellNo}));
                end
            end
        end
        
       
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        
        function u = getQoIVector(qoi, in, varargin)
            % Returns the qoi as a single vector.
            % Input can either be a seed or a problem.
            
            opt = struct('dtIndices', qoi.historyMatchDtRange);
            [opt, extra] = merge_options(opt, varargin{:});
            u = getQoIVector@WellQoI(qoi, in, 'dtIndices', opt.dtIndices, extra{:});
        end
        
        %-----------------------------------------------------------------%
        function [obs, scaling] = getObservationAndScaling(qoi, varargin) 
            
            opt = struct('vectorize', true );
            [opt, extra] = merge_options(opt, varargin{:});
            
            obs = qoi.getObservationVector('vectorize', false);
            
            % Apply one scaling per field name
            fieldScales = cell(1, numel(qoi.fldname));
            for f = 1:numel(qoi.fldname)
                fieldScales{f} = 0.0;
            end
            
            % Find scale values based on maximums
            for w = 1:numel(qoi.wellNames)
                for f = 1:numel(qoi.fldname)
                    fieldScales{f} = max(fieldScales{f}, max(abs(obs{w}{f}(:))));
                end
            end
            
            for w = 1:numel(qoi.wellNames)
                for f = 1:numel(qoi.fldname)
                    scaling{w}{f} = ones(size(obs{w}{f}(:)))*fieldScales{f};
                end
            end
            
            % The call to getObservationVector already extracts the correct
            % time indices, so we must avoid doing the same again here
            if opt.vectorize
                obs = qoi.qoi2vector(obs, 'dtIndices', []);
                scaling = qoi.qoi2vector(scaling, 'dtIndices', []);
            end
    
        end
        
        
        %-----------------------------------------------------------------%
        function R = getObservationErrorCov(qoi)
            
            assert(~isempty(qoi.observationCov), ...
                'qoi.observationCov is missing');
            
            % Check how observationCov matches the observation
            u = qoi.getObservationVector('vectorize', false);
            u_vec = qoi.qoi2vector(u, 'dtIndices', []);
            numObs = size(u_vec,1);
            
            if isscalar(qoi.observationCov)
                R = speye(numObs)*qoi.observationCov;
            
            elseif isvector(qoi.observationCov)
                rdiag = [];
            
                if numel(qoi.observationCov) == numel(qoi.fldname)*numel(qoi.wellNames)
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldname)
                            covindex = w*(numel(qoi.fldname)-1) + f;
                            rdiag = [rdiag ; repmat(qoi.observationCov(covindex), [numel(u{w}{f}), 1])];
                        end
                    end  
                elseif numel(qoi.observationCov) == numel(qoi.fldname)
                    % Assume that there are different covariances for each
                    % fieldName
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldname)
                            rdiag = [rdiag ; repmat(qoi.observationCov(f), [numel(u{w}{f}), 1])];
                        end
                    end    
                elseif numel(qoi.observationCov) == numObs
                    rdiag = qoi.observationCov;
                end
                assert(numel(rdiag) == numObs, ...
                    'Wrong number of diagonal elements after mapping the vector in qoi.observationCov');
                R = sparse([1:numObs], [1:numObs], rdiag, numObs, numObs); 
                
            else % observationCov is matrix
                assert(all(size(qoi.observationCov) == [numObs, numObs]), ...
                    'The matrix in qoi.observationCov is of wrong size');
                R = qoi.observationCov;
            end
        end
        
        function u = qoi2vector(qoi, u, varargin)
            opt = struct('dtIndices', qoi.historyMatchDtRange);
            [opt, extra] = merge_options(opt, varargin{:});
            u = qoi2vector@WellQoI(qoi, u, 'dtIndices', opt.dtIndices, extra{:});
        end
    
        
        
        
    end % methods
    
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------%
        function wellOutputOut = interpolateWellOutput(qoi, dtProblem, wellOutput)
            % If wellOutput is given with time intervals dtProblem, this
            % function can be used to give wellOutputs at the timesteps 
            % found in qoi.dt
            
            % Note: WellOutput interpolation is not done if the QoI is used
            % to history match part of the total timesteps only.
            % We therefore simply return the input if the history matching
            % dt range does not match the entire dt list.
            if numel(qoi.dt) ~= numel(qoi.historyMatchDtRange)
                wellOutputOut = wellOutput;
                return
            end
            
            wellOutputOut = interpolateWellOutput@WellQoI(qoi, dtProblem, wellOutput);
        end

        
    end
    
            

    
end

%-------------------------------------------------------------------------%
% Helpers
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function dt = getTimestepsFromProblem(problem)
    reports = reshape(problem.OutputHandlers.reports(:), [], 1);
    dt = zeros(numel(reports),1);
    for i = 1:numel(reports)
        dtstep = 0;
        stepReports = reports{i}.StepReports;
        for j = 1:numel(stepReports)
            dtstep = dtstep + stepReports{j}.Timestep;
        end
        dt(i) = dtstep;
    end     
end

%-------------------------------------------------------------------------%
function str = formatTime(timescale)
    t   = [second, minute, day, year];
    str = {'Second', 'Minute', 'Day', 'Year'};
    [d, ix] = min(abs(t - timescale));
    if d == 0
        str = str{ix};
    else
        warning('Non-standard timescale, using ''Second''');
        str = str{1};
    end
end

%{
#COPYRIGHT#
%}