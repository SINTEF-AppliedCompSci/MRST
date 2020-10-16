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
                        plot(t, ensemble.qoi.ResultHandler{i}{1}(:,w,fld),   'color', [1 1 1]*0.6);
                    end
                    plot(t, mean_qoi(:, w, fld), 'color', 'red')
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
        
    end
end

%-------------------------------------------------------------------------%
% Helpers
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
