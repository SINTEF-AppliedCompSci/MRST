function [model, schedule] = assignCoarseTransmissibility(p, model, schedule, T_struct, varargin)
    opt = struct('timeMethod',    'max', ...
                 'timeFunction',  [], ...
                 'phaseMethod',   'max', ...
                 'phaseFunction', []);
    opt = merge_options(opt, varargin{:});
    
    
    r = T_struct.reservoir;
    T_r = processConnection(schedule.step.val, r.T, r.T_ph, opt);
    
    model = upscaleModelTPFA(model, p, 'transCoarse', T_r, 'neighborship', r.N);
    schedule = upscaleSchedule(model, schedule);
    if nargout > 1
        w = T_struct.wells;
        for ctrl = 1:numel(schedule.control)
            WI = processConnection(schedule.step.val, w.WI, w.WI_ph, opt);
            for wNo = 1:numel(schedule.control(ctrl).W)
                schedule.control(ctrl).W(wNo).WI = WI(w.perf2well == wNo);
            end
        end
    end
end

function T = fixTrans(T)
    T(~isfinite(T)) = 0;
    if any(T > 0)
        T(T<0) = 1e-3*min(T(T>0));
    else
        T(T<0) = 0;
        warning('All transmissibilities are negative!');
    end
end

function T = processConnection(dt, T, T_ph, opt)
    nph = numel(T_ph);
    
    T = fixTrans(T);
    for i = 1:nph
        T_ph{i} = fixTrans(T_ph{i});
    end
    if isempty(opt.phaseFunction)
        switch lower(opt.phaseMethod)
            case 'max'
                T = -inf;
                for i = 1:nph
                    T = max(T, T_ph{i});
                end
            case 'min'
                T = inf;
                for i = 1:nph
                    T = min(T, T_ph{i});
                end
            case 'mean'
                T = 0;
                for i = 1:nph
                    T = T + T_ph{i};
                end
                T = T./nph;
            case 'none'
                % Just use total trans
            otherwise
                error(['Unknown matching method "', opt.phaseMethod, '".']);
        end
    else
        T = opt.phaseFunction(T, T_ph, dt);
    end
    
    if isempty(opt.timeFunction)
        switch lower(opt.timeMethod)
            case 'max'
                T = max(T, [], 2);
            case 'min'
                T = min(T, [], 2);
            case 'mean'
                T = bsxfun(@times, T, dt');
                T = sum(T, 2);
                T = T./sum(dt);
            otherwise
                error(['Unknown matching method "', opt.timeMethod, '".']);
        end
    else
        T = opt.timeMethod(T, T_ph, dt);
    end
end