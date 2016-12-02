function saveResults(dir, Gt, states, start_at, schedule, traps, rock, fluid, ...
                     sr, sw, dh, varargin)
    
    opt.initState = [];
    opt = merge_options(opt, varargin{:});
    
    if ~isdir(dir)
        mkdir(dir)
    end
    
    %% saving grid

    % Removing (inessential) functions whose implementations will in any case
    % be inaccessible upon load.
    Gt = rmfield(Gt, {'grav_pressure', 'primitives'});

    % Saving
    save(sprintf('%s/simulation_grid', dir), 'Gt');
    
    %% saving states

    W = [schedule.control(:).W]; 
    
    times   = cumsum(schedule.step.val);
    tnum    = numel(states);
    
    tot_inj = 0;
    for t = 1:(start_at - 1)
        cnum = schedule.step.control(t);
        dt = schedule.step.val(t);
        tot_inj = tot_inj + (sum([W(:,cnum).val]) * fluid.rhoGS * dt);
    end

    no_dissol = ~isfield(states{1}, 'sGmax'); 
    
    for t = 1:tnum
        t_global = t + (start_at - 1);
        cnum = schedule.step.control(t_global);
        cur_sol = states{t};
        
        if no_dissol
            cur_smax = cur_sol.smax(:,2);
        else
            cur_smax = cur_sol.sGmax;
        end
        rs = 0;
        if isfield(cur_sol, 'rs');
           rs = cur_sol.rs;
        end
        
        [cur_sol.h, cur_sol.h_max] = ...
            upscaledSat2height(cur_sol.s(:,2), cur_smax, Gt, 'resSat', [sw, sr]);
        
        mass_dist = ...
            massTrappingDistributionVEADI(Gt, cur_sol.pressure, cur_sol.s(:,2), ...
                                          cur_sol.s(:,1), cur_sol.h, cur_sol.h_max, ...
                                          rock, fluid, traps, dh, 'rs', rs); 

        dt      = schedule.step.val(t_global);
        tot_inj = tot_inj + (sum([W(:,cnum).val]) * fluid.rhoGS * dt);
        
        report.t      = times(t_global);
        report.sol    = cur_sol;
        report.W      = W(:,cnum);
        report.masses = [mass_dist, tot_inj - sum(mass_dist)];
        
        save(sprintf('%s/report_%i', dir, t_global), 'report');
    end
    
    % Saving initial state, if provided
    if ~isempty(opt.initState)
        report.t = 0;
        report.sol = opt.initState;
        report.W = W(:,1);
        save(sprintf('%s/report_0', dir),  'report');
    end
end
