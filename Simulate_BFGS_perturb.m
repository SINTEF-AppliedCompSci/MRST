function [misfitVal,varargout] = Simulate_BFGS_perturb(p,parameters,model,schedule,state0, wellSols_ref,weighting,objScaling,varargin)
    opt = struct('Verbose',           mrstVerbose());

    opt = merge_options(opt, varargin{:});     
 

             [model,schedule,state0] = control2problem(p,model,schedule,state0, parameters);     
            
             lsolve = AMGCL_CPRSolverAD('maxIterations', 200, 'tolerance', 1e-3);
                lsolve.setSRelaxation('ilu0');
                % We can set coarsening and solver options as well
                lsolve.setCoarsening('aggregation');
                lsolve.setSolver('bicgstab');
             solver = NonLinearSolver('LinearSolver',AMGCL_CPRSolverAD(),'continueOnFailure',true);
             
            [ wellSols,states] = simulateScheduleAD(state0, model, schedule);
        
% compute misfit function value (first each summand corresonding to each time-step)
    misfitVals = matchObservedOW(model.G, wellSols, schedule, wellSols_ref, weighting{:});
    
    objh = @(tstep)matchObservedOW(model.G, wellSols, schedule, wellSols_ref, 'computePartials', true, 'tstep', tstep, weighting{:});

% sum values to obtiain scalar objective 
    misfitVal = (objScaling - sum(vertcat(misfitVals{:}))) / objScaling ;
    
    if nargout > 1 % then the gradient is required
        np = numel(parameters);
        initStateSensitivity =  false;
        for k = 1:np
            params{k}     = parameters{k}.name;
            paramTypes{k} = parameters{k}.type; 
            boxLims{k}  = parameters{k}.boxLims;
            paramsDist{k} = parameters{k}.distribution;
            if  strcmp(parameters{k}.name , 'initSw')
                initStateSensitivity =  true;
            end
        end
        %Lets Assume param{4} is for state0
        grad = cell(1, numel(parameters));
        reel = 1;
        for k = 1:np% Numel Parameters


            grad{k} = zeros(size(parameters{k}.Indx,1), 1);
            valk    = zeros(size(parameters{k}.Indx,1), 1);
            Gradost = zeros(size(parameters{k}.Indx,1), 1);
            dispif(opt.Verbose, 'Solving for parameter  %d of %s', k, parameters{k}.name);

            parfor i = 1:size(parameters{k}.Indx,1) % Numel in each parameter

                dispif(opt.Verbose, 'Processing parameter %d of %s\n', i, parameters{k}.name);

                %% TODO: Define epsilon for the corresponding parameter 
                     e = 1e-7;       
                valk(i) = run_parallel(p,(reel+i-1),e,model,schedule,state0, parameters,wellSols_ref, weighting,objScaling)

                %% DONE: Calculate Gradient
                Gradost(i) = (valk(i)  -misfitVal)/e;
            end
            grad{k}=Gradost;
            reel = reel + numel(parameters{k}.Indx);
        end

        % Concatenate distributed gradient
        varargout{1} =   vertcat(grad{:}) ;  % Adjoinf functions gives the negative of the Gradient                                   
        if nargout > 2
            varargout{2} = wellSols;
            if nargout > 3
                varargout{3} = states;
            end
         end
    end

end

