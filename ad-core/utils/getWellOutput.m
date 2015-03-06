function [welldata, wellnames, fldnames] = getWellOutput(wellsols, fldnames, wells)
    assert(iscell(wellsols))
    nsteps = numel(wellsols);
    if nsteps == 0
        % No data, exit early
        welldata = [];
        wellnames = {};
        return
    end
        
    ws0 = wellsols{1};
    nw = numel(ws0);
    wnames = arrayfun(@(x) x.name, ws0, 'UniformOutput', false);
    
    % Valid fields are fields with exactly one value per well
    validfields = fieldnames(ws0);
    singledata = cellfun(@(x) numel(ws0(1).(x)) == 1, validfields);
    validfields = validfields(singledata);
    if nargin < 2
        fldnames = validfields;
    end
    
    % If user did not supply third argument, set subset to include all
    % wells by default
    if nargin < 3
        wells = 1:nw;
    end
    
    % Well subset can either be a list of indices or a set of names. We
    % will treat both cases as a set of numeric indices.
    if isnumeric(wells)
        subs = wells;
        if islogical(subs)
            subs = find(subs);
        end
        assert(min(subs) > 0 && max(subs) <= nw, 'Indices for wells outside of valid range');
    elseif ischar(wells) || iscell(wells)
        if ischar(wells)
            wells = {wells};
        end
        subs = cellfun(@(x) find(strcmpi(wnames, x)), wells);
    else
        error('Unknown format for well subset, provide either a list of names or indices');
    end
    wellnames = wnames(subs);
    
    if ~iscell(fldnames)
        fldnames = {fldnames};
    end
    welldata = nan(nsteps, numel(subs), numel(fldnames));
    
    for i = 1:numel(fldnames)
        fld = fldnames{i};
        assert(any(strcmpi(fld, validfields)),...
            'Only well fields with a single output value can be extracted using getWellOutput');
        data = cellfun(@(x) [x(subs).(fld)], wellsols, 'UniformOutput', false);
        welldata(:, :, i) = vertcat(data{:});
    end
end
