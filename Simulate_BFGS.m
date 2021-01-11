function [misfitVal,varargout] = Simulate_BFGS(p,parameters,model,schedule,state0, states_ref,objScaling, obj, varargin)
    opt = struct('Verbose',           mrstVerbose(),...
                 'Gradient',   'AdjointAD');

    opt = merge_options(opt, varargin{:});     
 

             [model,schedule,state0] = control2problem(p,model,schedule,state0, parameters);     
            
             lsolve = AMGCL_CPRSolverAD('maxIterations', 200, 'tolerance', 1e-3);
                lsolve.setSRelaxation('ilu0');
                % We can set coarsening and solver options as well
                lsolve.setCoarsening('aggregation');
                lsolve.setSolver('bicgstab');
             solver = NonLinearSolver('LinearSolver',AMGCL_CPRSolverAD(),'continueOnFailure',true);
             
            [ wellSols,states, schedulereport, model] = simulateScheduleAD(state0, model, schedule);
        
% compute misfit function value (first each summand corresonding to each time-step)
    misfitVals = obj(model, states, schedule, states_ref, false, []);
    
    objh = @(tstep) obj(model, states, schedule, states_ref, true, tstep);%'computePartials', true, 'tstep', tstep, weighting{:});

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

        switch opt.Gradient
                case 'AdjointAD'
                gradient = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                             'Parameters'    , {params{1:3}}, ...
                                             'ParameterTypes', {paramTypes{1:3}},...
                                             'initStateSensitivity',initStateSensitivity);
                case  'PerturbationAD'
                gradient = computeGradientPerturbationAD(state0, states, model, schedule, objh, ...
                                             'Parameters'    , {params{1:3}}, ...
                                             'ParameterTypes', {paramTypes{1:3}},...
                                             'initStateSensitivity',initStateSensitivity);
            otherwise   
                warning('Greadient method %s is not implemented',opt.Gradient)
        end
                                         

        
        reel = 1;
        % Distribute gradient for each parameter params{k} and scale it
        for k = 1:np
            switch paramsDist{k}
                case 'cell' %parameter disstribution per cell
                    switch params{k}
                            case 'porevolume'
                               I_pv = parameters{k}.Indx;
                               dBox   = boxLims{k}(2) - boxLims{k}(1);
                               scaled_gradient{k} = (dBox/objScaling)*gradient.porevolume(I_pv);
                            case 'initSw'
                               I_sw = parameters{k}.Indx;
                               dBox   = boxLims{k}(2) - boxLims{k}(1);
                               scaled_gradient{k} = (dBox/objScaling)*gradient.init.sW(I_sw);   
                            otherwise
                               warning('Parameter %s is not implemented',params{k})
                    end
               case  'connection'     
                   switch params{k}
                            case 'transmissibility'
                               for i =  1 : numel(parameters{k}.Indx)
                                   I_Tr = parameters{k}.Indx{i};
                                   dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                   scaled_grad_tr(i,1) = (dBox/objScaling)*sum(gradient.transmissibility(I_Tr));
                               end                               
                               scaled_gradient{k} = scaled_grad_tr;
                            case 'porevolume'
                               for i =  1 : numel(parameters{k}.Indx)
                                   I_pv = parameters{k}.Indx{i};
                                   dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                   scaled_grad_pv(i,1) = (dBox/objScaling)*sum(gradient.porevolume(I_pv));
                               end                               
                               scaled_gradient{k} = scaled_grad_pv; 
                            case 'permeability'
                               for i =  1 : numel(parameters{k}.Indx)
                                   I_perm = parameters{k}.Indx{i};
                                   if (parameters{k}.log == true)
                                        dBox_0 = boxLims{k}(i,2) - boxLims{k}(i,1);
                                        dBox   = (10^(  (p(i)*dBox_0)   + boxLims{k}(i,1)))*...
                                                     log(10)*dBox_0;
                                   else                                      
                                        dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                   end
                                   scaled_grad_perm(i,1) = (dBox/objScaling)*sum(gradient.permx(I_perm));
                               end                               
                               scaled_gradient{k} = scaled_grad_perm;                                   
                            otherwise
                               warning('Parameter %s is not implemented',params{k})
                   end
                case 'general'
                   switch params{k}
                            case 'transmissibility'
                                I_Tr =  parameters{k}.Indx;
                                 dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.transmissibility(I_Tr));
                            case 'porevolume'
                                I_pv =  parameters{k}.Indx;
                                dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.porevolume(I_pv));
                            case 'permeability'
                                I_perm =  parameters{k}.Indx;
                                dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.permx(I_perm));
                            case 'conntrans'
                                I_wi =  parameters{k}.Indx;
                                for i =  1 : size(parameters{k}.Indx,1)
                                    dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                    scaled_grad_well(i,1) = (dBox/objScaling)*gradient.conntrans(i);
                                end
                                scaled_gradient{k} = scaled_grad_well;
                            otherwise
                               warning('Parameter %s is not implemented',params{k})
                   end
              otherwise
            warning('Parameter distribution %s is not implemented',paramDist{k})
           end
        end
        % Concatenate distributed gradient
        varargout{1} =   vertcat(scaled_gradient{:}) ;  % Adjoinf functions gives the negative of the Gradient                                   
         if nargout > 2
            varargout{2} = wellSols;
            if nargout > 3
                varargout{3} = states;
            end
         end
    end

end

