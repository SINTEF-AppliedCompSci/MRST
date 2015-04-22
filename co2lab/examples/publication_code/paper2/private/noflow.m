function state = noflow(state)
% Trivial helper function. Sets flux to zero and pressure to empty.
% This ensures that any functions which use the pressure later on will not
% run, as it is not defined.

    state.flux     = state.flux*0.0;
    state.pressure = [];
end
