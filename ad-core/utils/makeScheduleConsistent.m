function schedule = makeScheduleConsistent(schedule, varargin)
%Ensure that a schedule is consistent in terms of well counts/perforations
%
% SYNOPSIS:
%   schedule = makeScheduleConsistent(schedule)
%
% DESCRIPTION:
%   For a given schedule with varying amount of wells and perforated cells
%   per well, this schedule makes the schedule internally consistent so
%   that all wells are defined at each control step. Some wells will be
%   disabled at different points, but they are always present and thus the
%   simulator output will be normalized and easier to work with.
%
% PARAMETERS:
%   schedule - Schedule with possibly inconsistent numbers of wells and
%              perforations.
%
% RETURNS:
%   schedule - Equivialent schedule that is consistent in the well and cell
%              numberings.
%
% SEE ALSO:
%   convertDeckScheduleToMRST

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
    opt = struct('perforationFields', {{'WI', 'dZ'}}, ...
                 'fixSign',           true);

    opt = merge_options(opt, varargin{:});
    
    % Find superset of well names
    W_all = [];
    % Names of all wells in the whole schedule
    names = [];
    % Flag indicating if a well has changing wells 
    cellsChangedFlag = [];
    
    perffields = opt.perforationFields;
    
    % First, we loop over all controls, adding any wells we haven't seen
    % before to the superset of all wells. If there are mismatches in the
    % perforated cells, we try to reconcile the differences without
    % changing the ordering. Other properties are ignored: These will be
    % replicated afterwards.
    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        currentNames = {W.name};
        
        [newNames, subs] = setdiff(currentNames, names);
        
        % Wells we have seen before
        other = true(size(W));
        other(subs) = false;
        other = find(other);

        for j = 1:numel(other)
            ind_all = find(strcmp(names, currentNames{other(j)}));
            c_all = W_all(ind_all).cells;
            c = W(other(j)).cells;
            
            % The number of cells / the actual cells have changed. We call
            % the remap function and hope for the best.
            if ~(numel(c) == numel(c_all) && all(c == c_all))
                % The perforations have changed
                cellsChangedFlag(ind_all) = true; %#ok
                W_all(ind_all).cells = remapIndices(c_all, c); %#ok
            end
        end
        
        % Wells that are new to us
        if ~isempty(newNames)
            names = [names, newNames]; %#ok
            W_all = [W_all; W(subs)]; %#ok
            cellsChangedFlag = [cellsChangedFlag; false(numel(subs), 1)]; %#ok
        end
    end
    
    % Create disabled/closed wells where neither the well nor the
    % perforations are flagged as active, but with the otherwise correct
    % dimensions.
    for i = 1:numel(W_all)
        c = W_all(i).cells;
        % Assign zero perforation values that will be filled in as we go
        for j = 1:numel(perffields)
            pf  = perffields{j};
            W_all(i).(pf) = zeros(numel(c),  size(W_all(i).(pf), 2)); %#ok
        end
        W_all(i).cstatus = false(numel(c), 1);
        W_all(i).status = false; %#ok
    end
    W_closed = W_all;
    
    % The schedule can now be updated from the original schedule, using the
    % superset wells that contains all wells from both the present, past
    % and future.
    passed = false(numel(W_all), 1);
    for i = 1:numel(schedule.control)       
        W = schedule.control(i).W;
        active = false(numel(W_all), 1);
        
        % Restfields are all fields that are not explicitly handled by the
        % this script and are not defined as perforation variables.
        restfields = setdiff(fieldnames(W), perffields);
        restfields = setdiff(restfields, {'cells', 'cstatus'});
        
        for j = 1:numel(W_all)
            % Find where we fit in the global well list
            sub = find(strcmp({W.name}, W_all(j).name));
            active(sub) = true;
            w = W(sub);
            % The completion active status can be defined from our already
            % computed well cell list.
            if ~isempty(sub)
                if cellsChangedFlag(j)
                    % Only some perforations are actually active.
                    isActivePerf = ismember(W_all(j).cells, w.cells(w.cstatus));
                else
                    % This well does not change the number of completions
                    isActivePerf = true(size(W_all(j).cells));
                end
                
                W_all(j).cstatus = ismember(W_all(j).cells, w.cells(w.cstatus)); %#ok
                
                for k = 1:numel(perffields)
                    pf = perffields{k};
                    % Take the values from the active perforations
                    W_all(j).(pf)(isActivePerf, :) = W(sub).(pf); %#ok
                end
                % Treat rest of the fields, whatever they may be
                for k = 1:numel(restfields)
                    fn = restfields{k};
                    W_all(j).(fn) = W(sub).(fn); %#ok
                end
                
                if ~passed(j) && W_all(j).status
                    % Grab the first active value for the closed set
                    W_closed(j).val  = w.val;
                    W_closed(j).type = w.type;
                    W_closed(j).sign = w.sign;
                    
                    passed(j) = true;
                end
            end
        end
        W = W_all;
        W(~active) = W_closed(~active);
        schedule.control(i).W = W;
    end
    
    % At this point, the schedule contains all the controls. We now want to
    % ensure that the disabled wells contain the controls from the first
    % time they are activated (so that any initialized well solutions are
    % reasonable for when they appear). Do another pass through, and ensure
    % that we also have the correct signs for all wells.
    for i = 1:numel(schedule.control)
        active = vertcat(schedule.control(i).W.status);
        W = schedule.control(i).W;
        W(~active) = W_closed(~active);
        schedule.control(i).W = setWellSign(W);
    end
end

function ind = remapIndices(old, new)
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
            elseif isinf(oOld)
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
