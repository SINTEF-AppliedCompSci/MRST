function scheduleMRST = convertDeckScheduleToMRST(G, rock, scheduleDeck)
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
end