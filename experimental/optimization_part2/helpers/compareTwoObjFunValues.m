function compareTwoObjFunValues( Gt, other, optim1, type1, ...
    optim2, type2, optim3, type3 )
% Compare the obj values at all time steps of two different simulation runs

% for example: simulation A is run for 2000 yrs migration with penalize
% leakage at each time step, while simulation B is run for 100 yrs
% migration with penalize future leakage (at each time step). The obj
% values are plotted to determine at what time they intersect.

    
    % Reconstruct the fluid structure which wasn't saved to avoid large
    % .mat files (~400 MB)
    if ~isfield(other,'fluid') %isempty(other.fluid)
        other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );
    end
    
    % items which should be consistent between all sim 1, 2, 3...
    model.fluid = other.fluid;
    model.G = Gt;
    model.rock = other.rock;
    
    figure; hold on
    plot(1);
    ylabel('J(t), (Gt CO2)')
    xlabel('time step')
    
    
    
    %% Simulation A:
    %%% recompute obj val at each step:
    if strcmpi(type1,'leak_penalizer')
    	obj_val_steps_1 = leak_penalizer_Rerun(model, optim1.wellSols, ...
                    optim1.states, optim1.schedule, other.opt.leak_penalty);
    
    elseif strcmpi(type1,'leak_penalizer_at_infinity')
        obj_val_steps_1 = leak_penalizer_at_infinity_Rerun(model, ...
                    optim1.wellSols, optim1.states, optim1.schedule, ...
                    other.opt.leak_penalty, other.opt.surface_pressure, ...
                    other.opt.rhoW);
    end
    
    plot([obj_val_steps_1{:}], 'x');
    legend('')
    
    %% Simulation B:
    if strcmpi(type2,'leak_penalizer')
    	obj_val_steps_2 = leak_penalizer_Rerun(model, optim2.wellSols, ...
                    optim2.states, optim2.schedule, other.opt.leak_penalty);
    
    elseif strcmpi(type2,'leak_penalizer_at_infinity')
        obj_val_steps_2 = leak_penalizer_at_infinity_Rerun(model, ...
                    optim2.wellSols, optim2.states, optim2.schedule, ...
                    other.opt.leak_penalty, other.opt.surface_pressure, ...
                    other.opt.rhoW);
    end
    
    plot([obj_val_steps_2{:}], 'x');
    legend('')
    
    %% Simulation C:
    if strcmpi(type3,'leak_penalizer')
    	obj_val_steps_3 = leak_penalizer_Rerun(model, optim3.wellSols, ...
                    optim3.states, optim3.schedule, other.opt.leak_penalty);
    
    elseif strcmpi(type3,'leak_penalizer_at_infinity')
        obj_val_steps_3 = leak_penalizer_at_infinity_Rerun(model, ...
                    optim3.wellSols, optim3.states, optim3.schedule, ...
                    other.opt.leak_penalty, other.opt.surface_pressure, ...
                    other.opt.rhoW);
    end
    
    plot([obj_val_steps_3{:}], 'x');
    legend('')
    

%     %% plotting:
%     figure; hold on
%     
%     time_yr = convertTo(cumsum(optim.schedule.step.val), year);
%     plot(time_yr, [optim.obj_val_steps_recompute{:}], 'x')
%     
%     time_yrB = convertTo(cumsum(optimB.schedule.step.val), year);
%     plot(time_yrB, [optimB.obj_val_steps_recompute{:}], 'o')
%     
%     legend('sim A', 'sim B')
%     ylabel('J(t)'); xlabel('time (yr)'); title('optimal simulation')
    
% 
%         figure; hold on
%     time_yr = convertTo(cumsum(init.schedule.step.val), year);
%     plot(time_yr, init.obj_val_steps, 'x')
%     plot(time_yr, [init.obj_val_steps_recompute{:}], 'o')
%         legend(['J(1:end-1) = Mi * (1-C), J(end) = Mi * (1-C) + C * Ma; total=',num2str(init.obj_val_total)], ...
%         ['J = Mi * (1-C) + C * Ma; total=',num2str(init.obj_val_total_recompute)])
%     ylabel('J(t)'); xlabel('time (yr)'); title('initial simulation')
%     
%     figure; hold on
%     plot(time_yr, [optim.obj_val_steps{:}], 'x')
%     plot(time_yr, [optim.obj_val_steps_recompute{:}], 'o')
%     legend(['J(1:end-1) = Mi * (1-C), J(end) = Mi * (1-C) + C * Ma; total=',num2str(optim.obj_val_total)], ...
%         ['J = Mi * (1-C) + C * Ma; total=',num2str(optim.obj_val_total_recompute)])
%     ylabel('J(t)'); xlabel('time (yr)'); title('optimal simulation')


end

