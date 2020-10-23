classdef WellQoI < BaseQoI
    
    properties
        fldname      = {'qOs'}
        wellIndices
        wellNames    
                               
        numTimesteps
        cumulative   = false % Cumulative production
        total        = false % Total production
        combined     = false % Combine all wells, or give per well values.
        timesteps
        
        historyMatchFromTimestep = 1 
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = WellQoI(varargin)
            qoi = qoi@BaseQoI();
            qoi = merge_options(qoi, varargin{:});
            % Check input
            assert(xor(isempty(qoi.wellIndices), isempty(qoi.wellNames)), ...
                'Please provide either well indices or well names');
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Check that the configs that are inserted to the constructor
            % makes sense for the base problem for the ensemble, and
            % updates remaining fields.
            qoi = validateQoI@BaseQoI(qoi, problem);
            % Either, we provide the indices of the chosen wells, or their
            % names. In either case, we find the other.
            names = {problem.SimulatorSetup.schedule.control(1).W.name};
            if ~isempty(qoi.wellNames)
                % Find well indices
                ix = cellfun(@(wn) find(strcmpi(wn, names)), ...
                                    qoi.wellNames, 'UniformOutput', false);
                % Verify that all wells exist in schedule
                found = ~cellfun(@isempty, ix);
                assert(all(found), 'Did not find well %s in schedule', ...
                                                qoi.wellNames{~found});
                % Assign well indices
                qoi.wellIndices = cell2mat(ix);
            elseif ~isempty(qoi.wellIndices)
                % Verify well indices
                assert(min(qoi.wellIndices) >= 1 && max(qoi.wellIndices) <= numel(names), ...
                    'Well indices must be in the range [1, numel(W)] = [1, %d]', numel(names));
                % Find well names
                qoi.wellNames = names(qoi.wellIndices);
            else
                % We should never get here
                error('Please provide either well names or well indices');
            end
            % Get timesteps
            qoi.timesteps = problem.SimulatorSetup.schedule.step.val;
            if ~isempty(qoi.numTimesteps)
                qoi.timesteps = qoi.timesteps(1:qoi.numTimesteps);
            else
                qoi.numTimesteps = numel(qoi.timesteps);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Reads the well solutions from the given problem and extract
            % the relevant data, which is then saved to file.
            
            % Read well solutions
            wellSols = reshape(problem.OutputHandlers.wellSols(:), [], 1);

            % Organize as matrix with dimensions:
            % (numTimesteps, numWells, numFields)
            wellOutputs = getWellOutput(wellSols, qoi.fldname, qoi.wellNames);

            dt = getTimestepsFromProblem(problem);
            
            numAvailableTimesteps = size(dt, 1);
            if qoi.numTimesteps < numAvailableTimesteps
                wellOutputs = wellOutputs(1:qoi.numTimesteps,:,:);
                dt = dt(1:qoi.numTimesteps);
            end
            
            % Compute total or cumulative if requested
            if qoi.total
                wellOutputs = sum(wellOutputs.*dt,1);
            elseif qoi.cumulative
                wellOutputs = cumsum(wellOutputs.*dt,1);
            end
            
            if qoi.combined && numel(qoi.wellNames) > 1
                wellOutputs = sum(wellOutputs, 2);
            end
            
            u = {wellOutputs};
            
        end

        %-----------------------------------------------------------------%
        function mean_u = getEnsembleMean(qoi, ensemble)
            mean_u = zeros(numel(qoi.timesteps), numel(qoi.wellNames), numel(qoi.fldname));
            for i = 1:ensemble.num
                mean_u = mean_u + ensemble.qoi.ResultHandler{i}{1};
            end
            mean_u = mean_u / ensemble.num;
        end
        
        %-----------------------------------------------------------------%
        function plotEnsemble(qoi, ensemble)
            % Plots well properties for the ensemble and ensemble mean.
            % Creates one figure per field
            % Organizes results from each well in a subplot
            
            num_fields = numel(qoi.fldname);
            num_wells  = numel(qoi.wellNames);
            
            num_wells_horizontal = ceil(sqrt(num_wells));
            num_wells_vertical   = ceil(num_wells/num_wells_horizontal);
            
            mean_qoi = qoi.getEnsembleMean(ensemble);
            
            t = cumsum(qoi.timesteps);
            for fld = 1:num_fields
                figure
                for w = 1:num_wells
                    subplot(num_wells_horizontal, num_wells_vertical, w)
                    hold on
                    for i = 1:ensemble.num
                        plot(t, ensemble.qoi.ResultHandler{i}{1}(1:qoi.numTimesteps,w,fld),   'color', [1 1 1]*0.6);
                    end
                    plot(t, mean_qoi(1:qoi.numTimesteps, w, fld), 'color', 'red')
                    title(strcat(qoi.fldname{fld}, " for well ", qoi.wellNames{w}));
                end
            end
            
        end
        
        function n = norm(qoi, u)
            if ~qoi.total
                n = cellfun(@(u) sum(u.*qoi.timesteps), u);
            else
                n = norm@BaseQoI(qoi, u);
            end
        end
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function u = getObservationVector(qoi, varargin)
            
            opt = struct('vectorize', true);
            [opt, extra] = merge_options(opt, varargin{:});
            
            % Check that the observation is valid
            assert(~isempty(qoi.observationResultHandler), ...
                'qoi.observationResultHandler is missing');
            assert(numel(qoi.observationResultHandler.getValidIds) > 0, ...
                'No available data in the observationResultHandler');
            
            u = qoi.observationResultHandler{1};
            
            assert(size(u,1) <= qoi.numTimesteps, ...
                'The observation has too few timesteps to match the QoI class');
            assert(size(u,2) == numel(qoi.wellNames), ...
                'The observation does not match the number of wells in QoI');
            assert(size(u,3) == numel(qoi.fldname), ...
                'The observation does not match the number of fldnames in QoI');
            
            u = u(qoi.historyMatchFromTimestep:qoi.numTimesteps, :, :);
            if opt.vectorize
                u = qoi2vector(u);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = getQoIVector(qoi, seed)
            u = qoi.ResultHandler{seed}{1};
            u = u(qoi.historyMatchFromTimestep:qoi.numTimesteps, :, :);
            u = qoi2vector(u);
        end
            
        %-----------------------------------------------------------------%
        function R = getObservationErrorCov(qoi)
            
            assert(~isempty(qoi.observationCov), ...
                'qoi.observationCov is missing');
            
            % Check how observationCov matches the observation
            u = qoi.getObservationVector('vectorize', false);
            u_vec = qoi2vector(u);
            numObs = size(u,1);
            
            if isscalar(qoi.observationCov)
                R = speye(numObs)*qoi.observationCov;
            
            elseif isvector(qoi.observationCov)
                rdiag = [];
            
                if numel(qoi.observationCov) == numel(qoi.fldnames)*numel(qoi.wellNames)
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldnames)
                            covindex = w*(numel(qoi.fldnames-1)) + f;
                            rdiag = [rdiag ; repmat(qoi.observationCov(covindex), [size(u,1), 1])];
                        end
                    end  
                elseif numel(qoi.observationCov) == numel(qoi.fldnames)
                    % Assume that there are different covariances for each
                    % fieldName
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldnames)
                            rdiag = [rdiag ; repmat(qoi.observationCov(f), [size(u,1), 1])];
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
        
    end % methods
    

end

%-------------------------------------------------------------------------%
% Helpers
%-------------------------------------------------------------------------%
function u = qoi2vector(u)
    % For multiple fields and wells, this vectorization will result in
    % u = [ (well1, field1), (well1, field2), (well2, field1), (well2,
    % field2) ...]
    
    u = u(:);
end

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
