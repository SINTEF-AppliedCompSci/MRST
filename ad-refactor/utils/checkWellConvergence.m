function [convergence, values] = checkWellConvergence(model, problem)
    ratevals = problem.equations(problem.indexOfType('perf'));
    
    bhpvals = problem.equations(problem.indexOfType('well'));
    
    convergence = [cellfun(@(x) all(double(x) < model.toleranceWellRate), ratevals), ...
                   cellfun(@(x) all(double(x) < model.toleranceWellBHP), bhpvals)];
               
    values = [cellfun(@(x) norm(double(x), inf), ratevals), ...
            cellfun(@(x) norm(double(x), inf), bhpvals)];
end
