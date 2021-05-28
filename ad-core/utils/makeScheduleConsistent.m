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
    opt = struct('perforationFields', {{'WI', 'dZ', 'dir', 'r', 'rR'}}, ...
                 'DepthReorder',      false, ...
                 'ReorderStrategy',   {{}}, ...
                 'G',                 [], ...
                 'fixSign',           true);

    opt = merge_options(opt, varargin{:});
    if isempty(schedule.control)
        return
    end

    % First, we loop over all controls, adding any wells we haven't seen
    % before to the superset of all wells. If there are mismatches in the
    % perforated cells, we try to reconcile the differences without
    % changing the ordering. Other properties are ignored: These will be
    % replicated afterwards.
    
    % Account for positional order in controls
    ctrl_order = getControlOrdering(schedule);
    % Calculate the superset of all wells
    [W_all, cellsChangedFlag] = getWellSuperset(schedule, ctrl_order, opt);
    % Update the schedule with superset of wells, setting inactive/active
    % wells and perforations
    schedule = updateSchedule(schedule, ctrl_order, W_all, cellsChangedFlag, opt);
    % Perform alternate ordering of well cells
    schedule = reorderWellsPerforations(schedule, opt);
    if ~isempty(W_all)
        % Well may have been updated, we need a representative well to work
        % with when setting reference depths.
        W_all = schedule.control(1).W;
    end
    % Fix defaulted reference depths
    schedule = setReferenceDepths(schedule, W_all, opt);
end

function schedule = setReferenceDepths(schedule, W_all, opt)
    for wNo = 1:numel(W_all)
        w = W_all(wNo);
        if (isfield(w.defaulted, 'refDepth') && w.defaulted.refDepth) || isnan(w.refDepth)
            assert(~isempty(opt.G), 'Grid must be provided when refDepth is defaulted');
            z = opt.G.cells.centroids(:, 3);
            c = w.cells;
            refDepth = z(c(1));
            dZ = z(w.cells) - refDepth;
            for i = 1:numel(schedule.control)
                schedule.control(i).W(wNo).refDepth = refDepth;
                schedule.control(i).W(wNo).dZ = dZ;
            end
        end
    end
end

function ctrl_order = getControlOrdering(schedule)
    nctrl = numel(schedule.control);
    found = false(nctrl, 1);
    ctrl_order = nan(nctrl, 1);
    pos = 1;
    for i = 1:numel(schedule.step.val)
        ctrl = schedule.step.control(i);
        if ~found(ctrl)
            ctrl_order(pos) = ctrl;
            found(ctrl) = true;
            pos = pos + 1;
        end
        if pos > nctrl
            break
        end
    end
    % NaN means entries were not found, we assign them increasing values
    not_found = isnan(ctrl_order);
    last_found = max(max(ctrl_order), 0);
    ctrl_order(not_found) = last_found + (1:sum(not_found))';
end

function [W_all, cellsChangedFlag] = getWellSuperset(schedule, ctrl_order, opt)
    perffields = opt.perforationFields;
    % Find superset of well names
    W_all = [];
    % Names of all wells in the whole schedule
    names = [];
    % Flag indicating if a well has changing wells
    cellsChangedFlag = [];
    
    for i = 1:numel(ctrl_order)
        order = ctrl_order(i);
        W = schedule.control(order).W;
        if isempty(W)
            continue
        end
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
            
            % The number of cells / the actual cells have changed.
            if numel(c) == numel(c_all) && ...
               (all(c == c_all) || all(sort(c) == sort(c_all)))
                % The two arrays are (possibly permuted) versions of
                % each other. No new cells are encountered. We continue.
                flag = ~all(W(other(j)).cstatus);
            else
                new_cells = setdiff(c, c_all, 'stable');
                flag = true;
                % Expand cells and cell_origin
                W_all(ind_all).cells = [W_all(ind_all).cells; new_cells];
                W_all(ind_all).cell_origin = [W_all(ind_all).cell_origin; ...
                                             repmat(i, size(new_cells))]; %#ok
            end
            cellsChangedFlag(ind_all) = cellsChangedFlag(ind_all) | flag;
        end
        
        % Wells that are new to us
        if ~isempty(newNames)
            names = [names, newNames]; %#ok
            W_new = W(subs);
            for j = 1:numel(W_new)
                W_new(j).cell_origin = repmat(i, numel(W_new(j).cells), 1);
            end
            W_all = [W_all; W_new]; %#ok
            changed_new = arrayfun(@(x) ~all(x.cstatus), W_new);
            cellsChangedFlag = [cellsChangedFlag; changed_new]; %#ok
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
            sample = W_all(i).(pf);
            d2 = size(W_all(i).(pf), 2);
            if isnumeric(sample)
                fn = @zeros;
            elseif islogical(sample)
                fn = @false;
            elseif iscell(sample)
                fn = @cell;
            elseif ischar(sample)
                fn = @(d1, d2) repmat(' ', d1, d2);
            else
                error('Unknown type %s', class(sample));
            end
            new_fld = fn(numel(c), d2);
            W_all(i).(pf) = new_fld;
        end
        W_all(i).cstatus = false(numel(c), 1);
        W_all(i).status = false; 
    end
end

function schedule = updateSchedule(schedule, ctrl_order, W_all, cellsChangedFlag, opt)
    % The schedule can now be updated from the original schedule, using the
    % superset wells that contains all wells from both the present, past
    % and future.
    perffields = opt.perforationFields;
    W_closed = W_all;
    
    passed = false(numel(W_all), 1);
    for i = 1:numel(ctrl_order)
        ctrl = ctrl_order(i);
        W = schedule.control(ctrl).W;
        if isempty(W)
            continue;
        end
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
                nwc = numel(W_all(j).cells);
                if cellsChangedFlag(j)
                    % Only some perforations are actually active.
                    [isActivePerf, order] = ismember(W_all(j).cells, w.cells(w.cstatus));
                    % Positions of currently active perforations in global
                    % cell list
                    order = order(order > 0);
                else
                    % This well has the same number of perforations across
                    % all time-steps. We do not really need to do anything.
                    isActivePerf = true(nwc, 1);
                    order = (1:nwc)';
                end
                W_all(j).cstatus = isActivePerf;
                for k = 1:numel(perffields)
                    pf = perffields{k};
                    % Take the values from the active perforations
                    if any(isActivePerf)
                        if not(isempty(W(sub).(pf)))
                            tmp = W(sub).(pf)(W(sub).cstatus, :);
                            W_all(j).(pf)(isActivePerf, :) = tmp(order, :);
                        end
                    end
                end
                % Treat rest of the fields, whatever they may be
                for k = 1:numel(restfields)
                    fn = restfields{k};
                    W_all(j).(fn) = W(sub).(fn);
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
        schedule.control(ctrl).W = W;
    end
    % At this point, the schedule contains all the controls. We now want to
    % ensure that the disabled wells contain the controls from the first
    % time they are activated (so that any initialized well solutions are
    % reasonable for when they appear). Do another pass through, and ensure
    % that we also have the correct signs for all wells.
    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        if ~isempty(W)
            active = vertcat(W.status);
            W(~active) = W_closed(~active);
            schedule.control(i).W = setWellSign(W);
        end
    end
    
end

function schedule = reorderWellsPerforations(schedule, opt)
    order = opt.ReorderStrategy;
    for wNo = 1:numel(schedule.control(1).W)
        schedule = setUniformDirection(schedule, wNo);
    end
    % One strategy per well is supported
    W = schedule.control(1).W;
    nw = numel(W);
    if isempty(order)
        % Defaulted
        if opt.DepthReorder
            % Backwards compatibility
            order = {'depth'};
        else
            % Default varies based on type of well
            order = cell(nw, 1);
            for i = 1:nw
               w = W(i);
               if all(w.dir == 'Z')
                   order{i} = 'depth';
               else
                   order{i} = 'direction-ijz-heel';
               end
            end
        end
    end
    
    if numel(order) == 1
        val = order{1};
        order = cell(nw, 1);
        [order{:}] = deal(val);
    end
    assert(numel(order) == nw);
    
    for wNo = 1:nw
        dispif(mrstVerbose(), 'Ordering well %d (%s) with strategy "%s".\n', ...
                                wNo, schedule.control(1).W(wNo).name, order{wNo})
        strategy = lower(order{wNo});
        if ~isempty(regexp(strategy, 'direction', 'once'))
            schedule = directionReorder(schedule, wNo, opt, strategy);
        else
            switch strategy
                case 'origin'
                    schedule = originReorder(schedule, wNo, opt);
                case 'depth'
                    schedule = depthReorder(schedule, wNo, opt);
                case 'none'
                    % We are leaving everything to chance!
                otherwise
                    error('Unknown ordering strategy %s', order{wNo});
            end
        end
    end
end

function schedule = setUniformDirection(schedule, wellNo)
    nc = numel(schedule.control(1).W(wellNo).cells);
    dir = repmat(' ', nc, 1);
    % Find all defined dZ values
    for i = 1:numel(schedule.control)
        w = schedule.control(i).W(wellNo);
        ok_dir = w.dir ~= ' ';
        dir(ok_dir) = w.dir(ok_dir);
    end
    dir(dir == ' ') = 'Z';
    % Set uniform direction
    for i = 1:numel(schedule.control)
        schedule.control(i).W(wellNo).dir = dir;
    end
end

function schedule = reorderCellFields(schedule, wellNo, opt, sortIx)
    flds = [opt.perforationFields, 'cells', 'cstatus'];
    for i = 1:numel(schedule.control)
        w = schedule.control(i).W(wellNo);
        for j = 1:numel(flds)
            f = flds{j};
            if isfield(w, f) && not(isempty(w.(f)))
                w.(f) = w.(f)(sortIx);
            end
        end
        schedule.control(i).W(wellNo) = w;
    end
end

function schedule = originReorder(schedule, wellNo, opt)
    W = schedule.control(1).W(wellNo);
    [~, sortIx] = sort(W.cell_origin);
    schedule = reorderCellFields(schedule, wellNo, opt, sortIx);
end

function schedule = depthReorder(schedule, wellNo, opt)
    w = schedule.control(1).W(wellNo);
    z = opt.G.cells.centroids(w.cells, 3);
    if issorted(z)
        return
    end
    [~, sortIx] = sort(z);
    schedule = reorderCellFields(schedule, wellNo, opt, sortIx);
end

function schedule = directionReorder(schedule, wellNo, opt, strategy)
    G = opt.G;
    w = schedule.control(1).W(wellNo);
    options = strsplit(strategy, '-');
    nopt = numel(options);
    % Get IJK and centroids for processing
    ijk = value(gridLogicalIndices(G));
    ijk = ijk(w.cells, :);
    centroids = G.cells.centroids(w.cells, :);

    if nopt > 1
        direction = options{2};
    else
        direction = 'xyz';
    end
    
    if nopt > 2
        switch options{3}
            case 'first'
                index = [];
            case 'heel'
                assert(isfield(w, 'heel') && numel(w.heel) > 1, ...
                    'Heel I, J index must be provided.');
                % Pick closest to heel
                ij_dist = sum(abs(bsxfun(@minus, ijk(:, 1:2), w.heel(1:2))), 2);
                mdist = min(ij_dist);
                candidates = find(ij_dist == mdist);
                % Pick top value if multiple cells are equally close
                [~, mz] = min(centroids(candidates, 3));
                index = candidates(mz);
            otherwise
                index = num2str(options{3});
        end
    else
        index = 1;
    end
    assert(not(isempty(G)), 'Grid must be provided for directional ordering');
    nc = numel(w.cells);
    sortIx = nan(nc, 1);
    % Always start with first perforation
    mO = min(w.cell_origin);
    sortIx(1) = find(w.cell_origin == mO, 1, 'first');
    dir = lower(w.dir);
    if numel(dir) == 1
        dir = repmat(dir, nc, 1);
    end
    convention = 'xyz';
    assert(all(ismember(dir, convention)));
    
    pts = zeros(size(centroids));
    ijk_names = 'ijk';
    unitscale = @(x) (x - min(x))./(max(x) - min(x));
    for i = 1:G.griddim
        if strcmpi(direction(i), ijk_names(i))
            pts(:, i) = unitscale(ijk(:, i));
        else
            pts(:, i) = unitscale(centroids(:, i));
        end
    end

    if isempty(index)
        [~, index] = min(pts(:, 3));
    end
    pts = swap(pts, index, 1);
    current_pt = pts(1, :);
    pts(1, :) = inf;
    % Loop over all segments. Assume that direction in previous cell is
    % used to define direction. Take closest point in xyz direction.
    for i = 2:nc
        coord = convention == dir(i-1);
        dist = sum((current_pt(:, coord) - pts(:, coord)).^2, 2);
        [tmp, pos] = min(dist);
        sortIx(i) = pos;
        current_pt = pts(pos, :);
        pts(pos, :) = inf;
    end
    sortIx = swap(sortIx, sortIx == index, sortIx == 1);
    schedule = reorderCellFields(schedule, wellNo, opt, sortIx);
end

function array = swap(array, in, out)
    tmp = array(out, :);
    array(out, :) = array(in, :);
    array(in, :) = tmp;
end
