% Parallel Computing Toolbox 
function model = ParEvalTimer(model, maxTime)
    % Invoke your function on a worker
    fut = parfeval(@Run, 1, model);
    % Block for up to maxTime seconds waiting for a result
    didFinish = wait(fut, 'finished', maxTime);
    if ~didFinish
        % Execution didn't finish in time, cancel this iteration
        cancel(fut);
        model.calculatedp = ...
            zeros(1,height(model.experiment.schedule.procedure))';
    else
        % Did complete, retrieve results
        model = fetchOutputs(fut);
    end
end