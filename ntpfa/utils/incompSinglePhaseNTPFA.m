function state = incompSinglePhaseNTPFA(model, varargin)
    if mod(numel(varargin), 2) > 0
        state0 = varargin{1};
        varargin = varargin(2:end);
    else
        state0 = initResSol(model.G, 0, [1, 0]);
    end
    opt = struct('Wells', [], 'src', [], 'bc', []);
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.src)
        for i = 1:numel(opt.src)
            opt.src(i).sat = repmat([1, 0], size(opt.src(i).cell));
        end
    end
    
    if ~isempty(opt.bc)
        opt.bc.sat = repmat([1, 0], size(opt.bc.face));
    end
    solver = NonLinearSolver('maxIterations', 1000);

    schedule = simpleSchedule(1, 'W', opt.Wells, 'bc', opt.bc, 'src', opt.src);
    
    
%     solver.minIterations = 100;
    [~, states] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', solver);
    state = states{1};
    state.flux = sum(state.flux, 2);
    if isfield(state, 'wellSol')
        for i = 1:numel(state.wellSol)
            state.wellSol(i).flux = sum(state.wellSol(i).flux, 2);
            state.wellSol(i).pressure = state.wellSol(i).bhp;
        end
    end
end