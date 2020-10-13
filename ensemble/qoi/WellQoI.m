classdef WellQoI < BaseQoI
    
    properties
        fldname      = {'qOs'}
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
                        if strcmp(qoi.wellNames{wi}, problem.SimulatorSetup.schedule.control.W(j).name)
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
        
        function mean_u = getEnsembleMean(qoi, ensemble)
            mean_u = zeros(numel(qoi.timesteps), numel(qoi.wellNames), numel(qoi.fldname));
            for i = 1:ensemble.num
                mean_u = mean_u + ensemble.qoi.ResultHandler{i}{1};
            end
            mean_u = mean_u / ensemble.num;
        end
        
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
