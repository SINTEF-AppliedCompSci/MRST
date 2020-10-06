classdef WellQoI < BaseQoI
    
    properties
        fldname      = 'qOs'
        wellIndices  = 1
        wellNames    
                               
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
            
            % Either, we provide the indices of the chosen wells, or their
            % names. If wellNames are provided, any wellIndices input is
            % overwritten
            if ~isempty(qoi.wellNames)
                numWellsTotal = numel(problem.SimulatorSetup.schedule.control.W);
                qoi.wellIndices = zeros(1, numel(qoi.wellNames));
                for wi = 1:numel(qoi.wellNames)
                    found = false;
                    for j = 1:numWellsTotal
                        if strcmp(qoi.wellNames{wi}, problem.SimulatorSetup.schedule.control.W(j).name);
                            qoi.wellIndices(wi) = j;
                            found = true;
                            break
                        end
                    end
                    assert(found, 'Did not find given well name');
                end
            else
                qoi.wellNames = cell(numel(qoi.wellIndices, 1));
                for wi = 1:numel(qoi.wellIndices)
                    qoi.wellNames{wi} = problem.SimulatorSetup.schedule.control.W(wi).name;
                end
            end % well names and indices
            
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
