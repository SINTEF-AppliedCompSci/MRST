classdef WellProductionQuantityOfInterest < QuantityOfInterestBase
   
    
    properties
        fldname      = 'qOs'
        wellIndices          % Not the same as wellProductionIndex
        wellNames    = {'P1'}
        numTimesteps
        cumulative   = false % Cumulative production
        total        = false % Total production 
        combined     = false % Combine all wells, or give per well values.
        timesteps  
        
        
    end
    
    methods
        function wellQoI = WellProductionQuantityOfInterest(simulatorSetup, varargin)
            wellQoI = wellQoI@QuantityOfInterestBase();
            wellQoI = merge_options(wellQoI, varargin{:});
            
            % Either, we provide the indices of the chosen wells,
            % or the names. If wellIndices is provided, any wellNames input
            % is overwritten
            if ~isempty(wellQoI.wellIndices)
                wellQoI.wellNames = cell(numel(wellQoI.wellIndices, 1));
                for wi = 1:numel(wellQoI.wellIndices)
                    wellQoI.wellNames{wi} = simulatorSetup.schedule.control.W(wi).name;
                end
            else
                numWellsTotal = numel(simulatorSetup.schedule.control.W);
                wellQoI.wellIndices = zeros(1, numel(wellQoI.wellNames));
                for wi = 1:numel(wellQoI.wellNames)
                    found = false;
                    for j = 1:numWellsTotal
                        if strcmp(wellQoI.wellNames{wi}, simulatorSetup.schedule.control.W(j).name);
                            wellQoI.wellIndices(wi) = j;
                            found = true;
                            break
                        end
                    end
                    assert(found, 'Did not find given well name');
                end
            end % well names and indices
            
            wellQoI.timesteps = simulatorSetup.schedule.step.val;
            if ~isempty(wellQoI.numTimesteps)
                wellQoI.timesteps = wellQoI.timesteps(1:wellQoI.numTimesteps);
            else
                wellQoI.numTimesteps = numel(wellQoI.timesteps);
            end
                
        end
        
        
        
        function wellOutputs = getQoI(wellQoI, problem)
            % getQoI reads the result files of the given problem and
            % extracts the quantity of interest, which is returned.
            % Get well output from the wells specified in the constructor
            
            wellSols = reshape(problem.OutputHandlers.wellSols(:), [], 1);
            
            % Organize as matrix with dimensions:
            % (numTimesteps, numWells, numFields)
            wellOutputs = getWellOutput(wellSols, wellQoI.fldname, wellQoI.wellIndices);
            
            dt = getTimestepsFromProblem(problem);
            
            numAvailableTimesteps = size(dt, 1);
            if wellQoI.numTimesteps < numAvailableTimesteps
                wellOutputs = wellOutputs(1:wellQoI.numTimesteps,:,:);
                dt = dt(1:wellQoI.numTimesteps);
            end
            
            % Compute total or cumulative if requested
            if wellQoI.total
                wellOutputs = sum(wellOutputs.*dt,1);
            elseif wellQoI.cumulative
                wellOutputs = cumsum(wellOutputs.*dt,1);
            end
            
            if wellQoI.combined && numel(wellQoI.wellIndices) > 1
                wellOutputs = sum(wellOutputs, 2);
            end
            
            
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
