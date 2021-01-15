function [misfitVal,varargout] = Simulate_BFGS(p,parameters,model_org,schedule_org,state0_org, states_ref,objScaling, obj, varargin)
    opt = struct('Verbose',           mrstVerbose(),...
                 'Gradient',   'AdjointAD',...
                 'NonlinearSolver', [],...
                 'pertub',1e-4);

    opt = merge_options(opt, varargin{:});     
 

             [model,schedule,state0] = control2problem(p,model_org,schedule_org,state0_org, parameters);     
            
             %lsolve = AMGCL_CPRSolverAD('maxIterations', 200, 'tolerance', 1e-3);
             %   lsolve.setSRelaxation('ilu0');
             % We can set coarsening and solver options as well
             %   lsolve.setCoarsening('aggregation');
             %    lsolve.setSolver('bicgstab');
             %solver = NonLinearSolver('LinearSolver',AMGCL_CPRSolverAD(),'continueOnFailure',true);
             
            [ wellSols,states, schedulereport, model] = simulateScheduleAD(state0, model, schedule,'NonLinearSolver',opt.NonlinearSolver);
        
% compute misfit function value (first each summand corresonding to each time-step)
    misfitVals = obj(model, states, schedule, states_ref, false, []);
    
    objh = @(tstep) obj(model, states, schedule, states_ref, true, tstep);%'computePartials', true, 'tstep', tstep, weighting{:});

% sum values to obtiain scalar objective 
    misfitVal = (objScaling - sum(vertcat(misfitVals{:}))) / objScaling ;
    %misfitVal = - sum(vertcat(misfitVals{:})) / objScaling ;
    
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
                                             'Parameters'    , {params{:}}, ...
                                             'ParameterTypes', {paramTypes{:}},...
                                             'initStateSensitivity',initStateSensitivity);
                 case  'PerturbationAD'
                       gradient = computeGradientPerturbationAD(state0, states, model, schedule, objh, ...
                                                     'Parameters'    , {params{:}}, ...
                                                     'ParameterTypes', {paramTypes{:}},...
                                                     'initStateSensitivity',initStateSensitivity);                               
                case  'PerturbationADNUM'
                    % do manual pertubuation of the defiend control variabels
                    eps_scale=1e-2;
                    p_org=p;
                    val=nan(size(p));
                    dp=nan(size(p));
                    for i=1:numel(p)
                        p_pert=p_org;
                        dp(i) = max(p_pert(i)*eps_scale,eps_scale);
                        p_pert(i) = p_pert(i) + dp(i);
                        val(i) = Simulate_BFGS(p_pert,parameters,model_org,schedule_org,state0_org, states_ref,objScaling, obj, 'Gradient', 'none' );
                    end
                    gradient= (val-misfitVal)./dp;
                    varargout{1} =   gradient;  % Adjoinf functions gives the negative of the Gradient                                   
                    if nargout > 2
                        varargout{2} = wellSols;
                        if nargout > 3
                            varargout{3} = states;
                        end
                    end
                    
                    return
                case 'none'
                    % No gradient caluculations to used for numerical
                    % differentiation
                    varargout{1} =   NaN ;  % Adjoinf functions gives the negative of the Gradient                                   
                    if nargout > 2
                        varargout{2} = wellSols;
                    if nargout > 3
                        varargout{3} = states;
                    end
                    return
                end

            otherwise   
                error('Greadient method %s is not implemented',opt.Gradient)
        end
                                         

        
        reel = 1;
        % Distribute gradient for each parameter params{k} and scale it
        for k = 1:np
            switch paramsDist{k}
                case 'cell' %parameter disstribution per cell
                    dBox   = boxLims{k}(2) - boxLims{k}(1);
                    Indx = parameters{k}.Indx;
                    switch params{k}
                            case {'porevolume'}
                               scaled_gradient{k} = (dBox/objScaling)*gradient.porevolume(Indx);
                            case 'initSw'
                               scaled_gradient{k} = (dBox/objScaling)*gradient.init.sW(Indx);   
                            otherwise
                               error('Parameter %s is not implemented',params{k})
                    end
                    reel = reel +  length(parameters{k}.Indx);
               case  'connection'
                   reel_old=reel;
                   switch params{k}
                            case 'transmissibility'
                               for i =  1 : numel(parameters{k}.Indx)
                                   I_Tr = parameters{k}.Indx{i};
                                   dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                   scaled_grad_tr(i,1) = (dBox/objScaling)*sum(gradient.transmissibility(I_Tr));
                               end                               
                               scaled_gradient{k} = scaled_grad_tr;
                               reel = reel + numel(parameters{k}.Indx) ;
                            case 'porevolume'
                               for i =  1 : numel(parameters{k}.Indx)
                                   I_pv = parameters{k}.Indx{i};
                                   dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                   scaled_grad_pv(i,1) = (dBox/objScaling)*sum(gradient.porevolume(I_pv));                                                      
                               end                               
                               scaled_gradient{k} = scaled_grad_pv;
                               reel = reel + numel(parameters{k}.Indx) ;
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
                               reel = reel + numel(parameters{k}.Indx) ;
                            case 'conntrans'                              
                                for i =  1 : numel(parameters{k}.Indx)
                                    indx=parameters{k}.Indx{i};
                                    dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                    scaled_grad_well(i,1) = (dBox/objScaling)*sum(gradient.conntrans(indx(:,3)));
                                end
                                scaled_gradient{k} = scaled_grad_well;
                                reel = reel + numel(parameters{k}.Indx)
                            otherwise
                               error('Parameter %s is not implemented',params{k})
                   end
                   
                   for i =  1 : (reel-reel_old)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       val =  p(reel_old + i-1)*(umax-umin)+umin; 
                       switch parameters{k}.type
                           case 'value'
                                
                           case 'multiplier'
                               scaled_gradient{k}(i,1) =   scaled_gradient{k}(i,1)/val;
                           otherwise
                               error('Type not valid: %s ', parameters{k}.type);
                       end           
                   end
                   
                case 'general'
                   switch params{k}
                            case 'transmissibility'
                                I_Tr =  parameters{k}.Indx;
                                 dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.transmissibility(I_Tr));
                                reel = reel +1;
                            case 'porevolume'
                                I_pv =  parameters{k}.Indx;
                                dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.porevolume(I_pv));
                                reel = reel +1;
                            case 'permeability'
                                I_perm =  parameters{k}.Indx;
                                dBox   = boxLims{k}(2) - boxLims{k}(1);
                                scaled_gradient{k} = (dBox/objScaling)*sum(gradient.permx(I_perm));
                                reel = reel +1;
                            case 'conntrans'
                                I_wi =  parameters{k}.Indx;
                                for i =  1 : size(parameters{k}.Indx,1)
                                    dBox   = boxLims{k}(i,2) - boxLims{k}(i,1);
                                    scaled_grad_well(i,1) = (dBox/objScaling)*gradient.conntrans(i);
                                end
                                scaled_gradient{k} = scaled_grad_well;
                                reel = reel + size(parameters{k}.Indx,1);
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
                if nargout > 4
                    varargout{4} = model;
                end
            end
         end
    end

end

