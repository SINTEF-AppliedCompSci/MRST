%% compare the obj values against obj values wrt:
% (a) future leakage penalized, or
% (b) pressure penalized

plotPenalizeFutureLeakageObj = true;
plotPenalizePressureObj = false;

% First load a simulation result (Gt, init, optim, other), 
% then run the following:

plotOptim = true; % otherwise plots Init

% Reconstruct the fluid structure which wasn't saved to avoid large
% .mat files (~400 MB)
if ~isfield(other,'fluid') %isempty(other.fluid)
    other.fluid = reconstructFluidStructure( Gt, other.rock, other.opt, other.dh );
end
% Renaming
model.fluid = other.fluid;
model.G = Gt;
model.rock = other.rock;

if plotOptim
    wellSols = optim.wellSols;
    states = optim.states;
    schedule = optim.schedule;
else
    wellSols = init.wellSols;
    states = init.states;
    schedule = init.schedule;  
end

if plotPenalizeFutureLeakageObj
    
    % Obtain the obj values: (Xyrs is the number of migration years)
    % Also, mass or volume inventories can be shown
    [obj_val_steps_Xyrs, Mi, Mi_tot, Ma] = leak_penalizer_Rerun(model, wellSols, ...
                        states, schedule, other.opt.leak_penalty); 
    [obj_val_steps_Xyrs_noPenalty] = leak_penalizer_Rerun(model, wellSols, ...
                        states, schedule, 1);

    [obj_val_steps_Xyrs_future, Mi_tot_inf, Ma_inf] = leak_penalizer_at_infinity_Rerun(model, ...
                        wellSols, states, schedule, ...
                        other.opt.leak_penalty, ...
                        other.opt.surface_pressure, ...
                        other.opt.rhoW, other.traps);
    [obj_val_steps_Xyrs_future_noPenalty, ~, ~] = leak_penalizer_at_infinity_Rerun(model, ...
                        wellSols, states, schedule, ...
                        1, ...
                        other.opt.surface_pressure, ...
                        other.opt.rhoW, other.traps, 'plotsOn',false);
    % NB: Mi_tot and Mi_tot_inf should be the same given the same injection schedule 
    assert(all(Mi_tot == Mi_tot_inf))
    
    % Compare obj values:
    figure; hold on; set(gcf,'Position',[2749 166 1120 411])
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs{:}], 'o') % Gt
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs_future{:}], 'o')
    hl = legend(['leakage penalized (c_p=',num2str(other.opt.leak_penalty),')'], ...
        ['future leakage penalized (c_p=',num2str(other.opt.leak_penalty),')']);
    set(hl,'Location','best')
    ylabel('J(t), (Gt CO2)'); xlabel('time (years)')
    %title([other.opt.modelname ', C=',num2str(other.opt.leak_penalty)])
    title(other.opt.modelname)
    box; grid;
    set(gca,'FontSize',16)
    
    % Compare obj values with c=1 (which corresponds to mass inventory):
    figure; hold on; set(gcf,'Position',[3873 165 1119 412])
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs_noPenalty{:}], 'o') % Gt
    plot(convertTo(cumsum(init.schedule.step.val), year), ...
        [obj_val_steps_Xyrs_future_noPenalty{:}], 'o')
    plot(convertTo(cumsum(init.schedule.step.val), year), Mi_tot , 'x')
    %plot(convertTo(cumsum(init.schedule.step.val), year), Ma , '+')
    %plot(convertTo(cumsum(init.schedule.step.val), year), Ma_inf , '+')
    assert(all([obj_val_steps_Xyrs_noPenalty{:}]' == Ma), 'obj fun values should match with mass remaining')
    assert(all([obj_val_steps_Xyrs_future_noPenalty{:}]' == Ma_inf), 'obj fun values for future leakage should match with future mass remaining')
    hl = legend('leakage penalized (c_p=1), mass remaining', ...
        'future leakage penalized (c_p=1), future mass remaining', ...
        'tot Mi(t), mass injected');
        %'Ma(t) (mass remaining)', ...
        %'future Ma(t) (future mass remaining)');
    set(hl,'Location','best')
    ylabel('J(t), (Gt CO2)'); xlabel('time (years)')
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