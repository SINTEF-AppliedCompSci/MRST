function scheduleMRST = convertDeckScheduleToMRST(model, scheduleDeck, varargin)
% Convert deck-type schedule to MRST style schedule
%
% SYNOPSIS:
%   schedule = convertDeckScheduleToMRST(model, deck.SCHEDULE)
%
% DESCRIPTION:
%   Take a schedule in deck-style (from for example the output of
%   readEclipseDeck), parse all wells and create a new schedule suitable
%   for 'simulateScheduleAD'.
%
% REQUIRED PARAMETERS:
%   G       - Valid grid (likely from initEclipseGrid). Must be the same
%           grid as the wells in the schedule are defined for.
% 
%   rock    - Valid rock used to compute the well indices. Typically from
%             initEclipseRock.
%
%   scheduleDeck - Either a deck struct from readEclipseDeck or the
%                  schedule (typically deck.SCHEDULE).
%
% 
% OPTIONAL PARAMETERS:
%
%   'StepLimit' - Only parse the first n control steps.
%
% RETURNS:
%   scheduleMRST - Schedule ready for simulation in 'simulateScheduleAD'.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('StepLimit',        inf, ...
                 'useCpGeometry',    size(model.G.cells.faces, 2) > 1, ...
                 'DepthReorder',     false, ...
                 'skipWarnings',     {{'RefDepth:BelowTopConnection', 'RepRad:FaceGeomMissing'}}, ...
                 'ReorderStrategy',  {{}}, ...
                 'EnsureConsistent', true);
    [opt, wellArg] = merge_options(opt, varargin{:});

    cellDims = [];
    if isfield(scheduleDeck, 'RUNSPEC') &&...
       isfield(scheduleDeck, 'SCHEDULE')
       if opt.useCpGeometry &&...
          isfield(scheduleDeck, 'GRID') &&...
          isfield(scheduleDeck.GRID, 'COORD')
           [~, ~, cellDims] = computeCpGeometry(model.G, scheduleDeck.GRID);
           cellDims = cellDims(model.G.cells.indexMap,:);
       end
       % Support passing deck directly
       scheduleDeck = scheduleDeck.SCHEDULE;
    end
    scheduleMRST = struct('step', scheduleDeck.step);
    isCompostional = isa(model, 'ThreePhaseCompositionalModel');
    if isCompostional
        n_hc = model.EOSModel.getNumberOfComponents();
    end
    nc = numel(scheduleDeck.control);
    tmp = cell(nc, 1);
    scheduleMRST.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    
    % Massage phases in compi to match active components in model
    satVarNames = model.getSaturationVarNames();
    ncomp = numel(satVarNames);
    map = zeros(ncomp, 1);
    for i = 1:ncomp
        switch lower(satVarNames{i})
            case 'sw'
                map(i) = 1;
            case 'so'   
                map(i) = 2;
            case 'sg'
                map(i) = 3;
            case 'ss'
                map(i) = 4;
            otherwise
                map(i) = 4;
                warning('Unknown phase, translation directly from deck difficult, setting zero');
        end
    end
    nwarn = numel(opt.skipWarnings);
    warnStatus = cell(1, nwarn);
    for i = 1:nwarn
        w = opt.skipWarnings{i};
        warnStatus{i} = warning('query', w);
        warning('off', w);
    end
    for i = 1:nc
        % Parse well
        W = processWells(model.G, model.rock, scheduleDeck.control(i), 'OutputDefaulted', true, ...
            'DepthReorder', opt.DepthReorder, wellArg{:}, 'cellDims', cellDims);
        
        for j = 1:numel(W)
            c = [W(j).compi, 0];
            if isfield(W(j), 'solventFrac')
                gasIx = 3;
                solIx = 4;
                c(solIx) = c(gasIx)*W(j).solventFrac;
                c(gasIx) = c(gasIx)*(1-W(j).solventFrac);
            end
            ci = c(map);
            W(j).compi = ci/sum(ci);
            if isCompostional
                W(j).components = ones(1, n_hc)/n_hc;
            end
        end
        
        if isfield(W, 'solventFrac')
            W = rmfield(W, 'solventFrac');
        end
        scheduleMRST.control(i).W = W;
    end
    for i = 1:nwarn
        warning(warnStatus{i}.state, opt.skipWarnings{i});
    end
    
    if ~isinf(opt.StepLimit)
        scheduleMRST.step.val     = scheduleMRST.step.val(1:opt.StepLimit);
        scheduleMRST.step.control = scheduleMRST.step.control(1:opt.StepLimit);
    end

    if opt.EnsureConsistent
        scheduleMRST = makeScheduleConsistent(scheduleMRST, ...
                            'DepthReorder', opt.DepthReorder, ...
                            'ReorderStrategy', opt.ReorderStrategy, ...
                            'G', model.G);
        % Apply productivity modifiers post consistency. Algorithm assumes
        % that the wells have been made consistent.
        ectrls = scheduleDeck.control;
        if isfield(ectrls(1), 'WPIMULT')
            has_wpi = any(arrayfun(@(x) ~isempty(x.WPIMULT), ectrls));
            if has_wpi
                IJK = gridLogicalIndices(model.G);
                W_prev = scheduleMRST.control(1).W;
                WI_prev = applyFunction(@(x) x.WI, W_prev);
                WI_raw = WI_prev;
                for ctrl = 1:nc
                    wpi = ectrls(ctrl).WPIMULT;
                    W = scheduleMRST.control(ctrl).W;
                    [W, WI_raw, WI_prev] = apply_wpimult(W, IJK, wpi, WI_raw, WI_prev);
                    scheduleMRST.control(ctrl).W = W;
                    % [WI_raw{3}, WI_prev{3}, scheduleMRST.control(ctrl).W(3).WI]
                end
            end
        end
    end
end

function [W, WI_raw, WI_now] = apply_wpimult(W, IJK, WPIMULT, WI_raw, WI_now)
    wnames = {W.name};
    % Need to keep track of two things: Last defined unmodified WI and
    % current WI with multipliers. Unmodified is needed to know when we
    % need to reset WI since it has changed from the original definition.
    for well = 1:numel(W)
        WI_r = WI_raw{well};
        WI = WI_now{well};
        % Look at current WI relative to old WI (without mults applied)
        WI_current = W(well).WI;
        % If these are different, and the new one is non-zero, update it,
        % since it means that we found a reset value (or a new one).
        new_perf = WI_r ~= WI_current & WI_current > 0;
        new_vals = WI_current(new_perf);
        WI_r(new_perf) = new_vals;
        WI(new_perf) = new_vals;
        % Update storage
        WI_raw{well} = WI_r;
        WI_now{well} = WI;
    end
    % Loop over and apply WI mults
    for i_w = 1:numel(WPIMULT)
        wpi = WPIMULT{i_w};
        for i = 1:size(wpi, 1)
            [name, val, I, J, K, start, stop] = deal(wpi{i, :});
            well = find(strcmp(wnames, name));
            assert(numel(well) == 1);
            assert(val > 1e-10);
            assert(isfinite(val));
            w = W(well);
            N = numel(w.WI);
            wc = w.cells;
            if start < 1
                start = 1;
            end
            if stop < 1
                stop = N;
            end
            rng = (1:N)';
            active = rng >= start & rng <= stop;
            if I > 0
                active = active & IJK{1}(wc) == I;
            end
            if J > 0
                active = active & IJK{2}(wc) == J;
            end
            if K > 0
                active = active & IJK{3}(wc) == K;
            end
            WI = WI_now{well};
            % Apply mult after masking
            WI(active) = WI(active)*val;
            WI_now{well} = WI;
        end
    end
    % Finally make sure that all wells have the updated WIs
    for well = 1:numel(W)
        w = W(well);
        % Apply cstatus again, just to be sure.
        w.WI = WI_now{well}.*w.cstatus;
        W(well) = w;
    end
end
