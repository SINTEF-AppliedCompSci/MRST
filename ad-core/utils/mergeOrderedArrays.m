function ind = mergeOrderedArrays(old, new)
    % Merge two sets of cells that are similar in that they may contain
    % elements from the same superset in the same order, but each set may
    % be missing one or more elements that the other has.
    %
    % This is done by having two simple pointers that are incremented as we
    % go along trying to merge the two sets.
    N = numel(unique([old; new]));
    if isempty(old)
        ind = new;
        return
    elseif isempty(new)
        ind = old;
        return
    end
    if isempty(setdiff(new, old))
        if ~all(new == old)
            warning('Index sets are permuted versions of each other. Taking old ordering uncritically.');
        end
        ind = old;
        return
    end
    ind = nan(N, 1);
    
    iOld = 1;
    iNew = 1;
    for i = 1:N
        nv = new(iNew);
        no = old(iOld);
        if nv == no
            ind(i) = nv;
            iNew = iNew + 1;
            iOld = iOld + 1;
        else
            oOld = searchForward(old(iOld+1:end), nv);
            oNew = searchForward(new(iNew+1:end), no);
            
            if isinf(oNew)
                ind(i) = no;
                iOld = iOld + 1;
            elseif isinf(iOld)
                ind(i) = nv;
                iNew = iNew + 1;
            else
                error('Unable to correctly reassign indices based on given data. Consider sorting them first.')
            end
            assert(~(isinf(oOld) && isinf(oNew)));
        end
    end
end

function offset = searchForward(d, i)
    offset = find(d == i);
    if isempty(offset)
        offset = inf;
    end
end
