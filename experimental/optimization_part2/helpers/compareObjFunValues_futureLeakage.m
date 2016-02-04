%% compare the obj values against obj values wrt:
% (a) future leakage penalized, or
% (b) pressure penalized

plotPenalizeFutureLeakageObj = true;
plotPenalizePressureObj = false;

% First load a simulation result (Gt, init, optim, other), 
% then run the following:

% Reconstruct the fluid structure which wasn't saved to avoid large
% .mat files (~400 MB)
if ~isfield(other,'fluid') %isempty(other.fluid)
    other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );
end
% Renaming
model.fluid = other.fluid;
model.G = Gt;
model.rock = other.rock;

if plotPenalizeFutureLeakageObj
    
    % Obtain the obj values: (Xyrs is the number of migration years)
    obj_val_steps_Xyrs = leak_penalizer_Rerun(model, optim.wellSols, ...
                        optim.states, optim.schedule, other.opt.leak_penalty); 

    [obj_val_steps_Xyrs_future, vol, vol_inf] = leak_penalizer_at_infinity_Rerun(model, ...
                        optim.wellSols, optim.states, optim.schedule, ...
                        other.opt.leak_penalty, other.opt.surface_pressure, ...
                        other.opt.rhoW);

    % Plot:
    figure; hold on
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs{:}], '+')
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs_future{:}], 'o')
    hl = legend('X yrs migration, with leakage penalized',...
        'X yrs migration, with future leakage penalized');
    set(hl,'Location','best')
    ylabel('J(t), optimal, (Gt CO2)'); xlabel('time (years)')
    title(other.opt.modelname)
    box; grid;
    set(gca,'FontSize',16)

    
    
elseif plotPenalizePressureObj
    
    % compute pressure limit (unless exists in 'other')
    [P_over, p_limit] = computeOverburdenPressure(Gt, other.rock, ...
                                    other.opt.ref_depth, model.fluid.rhoWS);
    %P_limit = mean(P_over); % @@ try using an array for P_limit
    P_limit = P_over; % an array
    
    % set pressure penalty (unless exists in 'other.opt')
    pressure_penalty = 1e-9; % @@ adjustable
    
    % Obtain the obj values: (Xyrs is the number of migration years)
    obj_val_steps_Xyrs = leak_penalizer_Rerun(model, optim.wellSols, ...
                        optim.states, optim.schedule, other.opt.leak_penalty); 

    obj_val_steps_Xyrs_pressure = pressure_penalizer_Rerun(model, ...
                        optim.states, optim.schedule, ...
                        pressure_penalty, P_limit);

    obj_val_steps_combined = cell2mat(obj_val_steps_Xyrs) + cell2mat(obj_val_steps_Xyrs_pressure);
    % Plot:
    figure; hold on
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs{:}], '+')
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs_pressure{:}], 'o')
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        obj_val_steps_combined, 'x')
    hl = legend('X yrs migration, with leakage penalized',...
        'X yrs migration, with pressure penalized', ...
        'summation of obj fun values');
    set(hl,'Location','best')
    ylabel('J(t), optimal, (Gt CO2)'); xlabel('time (years)')
    title(other.opt.modelname)
    box; grid;
    set(gca,'FontSize',16)
    

    
    
end