classdef WellQoI < BaseQoI
    
    properties
        fldname      = 'qOs'
        wells    = 1 % This can either be a cell array of well 
                               % names or indices.
        numTimesteps
        cumulative   = false % Cumulative production
        total        = false % Total production
        combined     = false % Combine all wells, or give per well values.
        timesteps
    end
    
    methods
        
        function qoi = WellQoI(varargin)
            qoi = qoi@BaseQoI();
            qoi = merge_options(qoi, varargin{:});
        end
        
        function qoi = validateQoI(qoi, problem)
            % Check that the configs that are inserted to the constructor 
            % makes sense for the base problem for the ensemble, and
            % updates reminding fields.
            
            qoi = validateQoI@BaseQoI(qoi, problem);
            
            qoi.timesteps = problem.SimulatorSetup.schedule.step.val;
            if ~isempty(qoi.numTimesteps)
                qoi.timesteps = qoi.timesteps(1:qoi.numTimesteps);
            else
                qoi.numTimesteps = numel(qoi.timesteps);
            end
        end
        
        function u = computeQoI(qoi, problem)
            % Reads the well solutions from the given problem and extract
            % the relevant data, which is then saved to file.
            
            % Read well solutions
            wellSols = reshape(problem.OutputHandlers.wellSols(:), [], 1);

            % Organize as matrix with dimensions:
            % (numTimesteps, numWells, numFields)
            wellOutputs = getWellOutput(wellSols, qoi.fldname, qoi.wells);

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
            
            if qoi.combined && numel(qoi.wells) > 1
                wellOutputs = sum(wellOutputs, 2);
            end
            
            u = wellOutputs;
            
        end
        
        
        
    end
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
