function [convergence, values] = checkWellConvergence(model, problem)
    isperf = problem.indexOfType('perf');
    iswell = problem.indexOfType('well');
    
    
    ratevals = problem.equations(isperf);
    bhpvals = problem.equations(iswell);

    values = [cellfun(@(x) norm(double(x), inf), ratevals), ...
              cellfun(@(x) norm(double(x), inf), bhpvals)];
    
    tmp = find(isperf | iswell);
    isperf = isperf(tmp);
    iswell = iswell(tmp);
    
    convergence = false(size(tmp));
    convergence(isperf) = values(isperf) < model.toleranceWellRate;
    convergence(iswell) = values(iswell) < model.toleranceWellBHP;
end
