function ind = mergeOrderedArrays(old, new)
    % Merge two sets of cells that are similar in that they may contain
    % elements from the same superset in the same order, but each set may
    % be missing one or more elements that the other has.
    %
    % This is done by having two simple pointers that are incremented as we
    % go along trying to merge the two sets.

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

    combined = unique([old; new]);
    N = numel(combined);
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
    no = numel(old);
    nn = numel(new);
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
            elseif isinf(oOld)
                ind(i) = nv;
                iNew = iNew + 1;
            else
                warning('Unable to correctly reassign indices based on given data. Consider sorting them first.')
                ind = unique([old; new]);
                return
            end
            assert(~(isinf(oOld) && isinf(oNew)));
        end
        % If one of the arrays stops, we just finish by taking the rest of
        % the other one
        if iNew > nn && iOld < no
            ind(i+1:end) = old(iOld:end);
            break
        end
        if iNew < nn && iOld > no
            ind(i+1:end) = new(iNew:end);
            break
        end
    end
end

function offset = searchForward(d, i)
    offset = find(d == i);
    if isempty(offset)
        offset = inf;
    end
end
