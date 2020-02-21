function stateDiff = compareStates(state1, state2, varargin)

    opt = struct('relative', false);
    opt = merge_options(opt, varargin{:});

    stateDiff = state1;
    flds = fieldnames(state1);
    for fNo = 1:numel(flds)
        v = flds{fNo};
        if isfield(state2, v)    && ...
           isnumeric(state2.(v)) && ...
           all(size(state2.(v)) == size(state1.(v)))
           stateDiff.(v) = abs(state1.(v) - state2.(v));
           if opt.relative
               stateDiff.(v) = stateDiff.(v)./state1.(v);
           end
        else
            stateDiff = rmfield(stateDiff, v);
        end
    end
       
end