function [state, dx, fail] = linesearchADI(state0, dx0, system, getEqs, updateState, isBO)

    ni = system.nonlinear.lineIter;
    target = system.nonlinear.lineRelTol;

    fail = true;
    i = 0;
    alph = 0;

    e = @(eqs) cellfun(@(x) norm(x, 'inf'), {eqs{system.cellwise}});

    if isBO
        [eqs, history, explTrms] = getEqs(state0);
    else
        eqs = getEqs(state0);
    end
    err0 = e(eqs);


    dx = dx0;

    target = target*norm(err0);
    while fail && (i < ni)
        dx = cellfun(@(x) pow2(x, alph), dx0, 'UniformOutput', false);
        if isBO
            state = updateState(dx, explTrms);
            [eqs, history, explTrms] = getEqs(state);
        else
            state = updateState(dx);
            eqs = getEqs(state);
        end
        err = e(eqs);
        alph = alph - 1;
        i    = i + 1;

        fail = ~(norm(err) < target);

    end
    if fail
        state = state0;
        dx = dx0;
    end
end
