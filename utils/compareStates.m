function stateDiff = compareStates(state1, state2)

    stateDiff = state1;
    flds = fieldnames(state1);
    for fNo = 1:numel(flds)
        v = flds{fNo};
        if isfield(state2, v)    && ...
           isnumeric(state2.(v)) && ...
           all(size(state2.(v)) == size(state1.(v)))
           stateDiff.(v) = abs(state1.(v) - state2.(v));
        end
    end
       
end