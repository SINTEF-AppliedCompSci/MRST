function scheduleMRST = convertDeckScheduleToMRST(G, rock, scheduleDeck, varargin)
    opt = struct('StepLimit', inf);
    opt = merge_options(opt, varargin{:});


    if isfield(scheduleDeck, 'RUNSPEC') &&...
       isfield(scheduleDeck, 'SCHEDULE')
       % Support passing deck directly
       scheduleDeck = scheduleDeck.SCHEDULE;
    end
    scheduleMRST = struct('step', scheduleDeck.step);
    
    nc = numel(scheduleDeck.control);
    
    tmp = cell(nc, 1);
    scheduleMRST.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    for i = 1:nc
        scheduleMRST.control(i).W = ...
            processWellsLocal(G, rock, scheduleDeck.control(i));
    end
    
    if ~isinf(opt.StepLimit)
        scheduleMRST.step.val     = scheduleMRST.step.val(1:opt.StepLimit);
        scheduleMRST.step.control = scheduleMRST.step.control(1:opt.StepLimit);
    end
end