function [state, dx, fail] = linesearchADI(state0, dx0, system, getEqs, updateState, isBO)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
