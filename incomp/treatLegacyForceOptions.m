function opt = treatLegacyForceOptions(opt)
    % Internal function for ensuring that both W and Wells are supported as
    % keyword arguments to incompressible solvers.
    hasWells = isfield(opt, 'wells');
    hasW = isfield(opt, 'W');
    if hasWells
        Wells = opt.wells;
    else
        Wells = [];
    end
    if hasW
        W = opt.W;
    else
        W = [];
    end
    
    if ~isempty(W)
        assert(isempty(Wells), ...
            'Solver cannot have both legacy ''wells'' argument and ''W'' argument at the same time');
        opt.wells = W;
    end
    if ~isempty(Wells)
        assert(isempty(W), ...
            'Solver cannot have both legacy ''wells'' argument and ''W'' argument at the same time');
        opt.W = Wells;
    end
end