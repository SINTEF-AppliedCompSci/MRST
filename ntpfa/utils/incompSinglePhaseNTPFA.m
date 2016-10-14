function state = incompSinglePhaseNTPFA(model, varargin)
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
    state0 = initResSol(model.G, 0, [1, 0]);
    
%     solver.minIterations = 100;
    [~, states] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', solver);
    state = states{1};
end