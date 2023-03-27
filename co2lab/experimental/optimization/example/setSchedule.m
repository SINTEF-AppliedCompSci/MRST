function schedule = setSchedule(Gt, rock, wcells, qtot, ...
                                isteps, itime, msteps, mtime, single_control, ...
                                varargin)
                                    
   opt.minval = 0;
   opt = merge_options(opt, varargin{:});
   
    assert(isteps>0);
    if msteps == 1
        msteps = 2;
        warning(['If migration is happening, we need at least two steps, in ' ...
                 'order to have a transition step.  Migration step has been ' ...
                 'increased to two.']);
    end
    
    % constructing wells
    W = [];
    wellradius = 0.3;
    dummy = 0;
    
    % computing fixed rates
    rate = qtot / itime;  
    
    for i = 1:numel(wcells)
        W = addWell(W, Gt, rock, wcells(i), ...
                    'Type', 'rate', ...
                    'Val', rate(i), ...
                    'Radius',  wellradius, ...
                    'Comp_i', [0, 1], ...
                    'name', ['I', num2str(i)]);
    end

    % constructing schedule
    cpos = 1;
    ctrls = [];
    if single_control && isteps > 0 
        schedule.control(cpos).W = W;
        cpos = cpos+1;
        ctrls = ones(isteps,1);
    else
        for i = 1:isteps
            schedule.control(cpos).W = W;
            cpos = cpos+1;
        end
        ctrls = [1:isteps]';
    end    
    
    if msteps > 0
        schedule.control(cpos).W = W;
        for i = 1:numel(schedule.control(cpos).W)
            schedule.control(cpos).W(i).val = opt.minval;
        end
        cpos = cpos+1;
    end

    mig_ctrl = numel(schedule.control); % control step for migration, only
                                        % relevant (and correct) if msteps > 0

    dTi = itime/isteps;
    dTm = mtime/(msteps-1);
    
    istepvec = ones(isteps, 1) * dTi; 
    mstepvec = [dTi; ones(msteps-1,1) * dTm];
    mstepvec(2) = mstepvec(2) - dTi;

    schedule.step.val = [istepvec; mstepvec];
    schedule.step.control = [ctrls; mig_ctrl * ones(msteps,1)];
end
