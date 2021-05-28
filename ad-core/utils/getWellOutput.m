function [welldata, wellnames, fldnames] = getWellOutput(wellsols, fldnames, wells, index)
%Extract values from wellsols.
%
% SYNOPSIS:
%   [wd, wn, flds]= getWellOutput(wellSols, 'bhp')
%   [wd, wn, flds]= getWellOutput(wellSols, {'bhp', 'qWs'}, 'W1')
%
% DESCRIPTION:
%   Given a cell array of well solution structures representing multiple
%   timesteps, this routine extracts requested values in matrix form ready
%   for plotting / inspection.
%
% REQUIRED PARAMETERS:
%   wellSols - Cell array of NSTEP length, each containing a struct with
%              NWELL entries. All wells must exist at all timesteps.
%                
%   fldnames - (OPTIONAL) Either a single string, or a cell array of
%              desired fields for output. Bottom hole pressures, rates, ...
%              Defaults to all fields.
%
%   wells    - (OPTIONAL) Either a single string, or a cell array of well
%              names for which output is desired. Defaults to all wells.
% RETURNS:
%   welldata - A NSTEP by NWELL by NFIELDS matrix. For instance, for
%   calling
%
%   D = getWellOutput(wellsols, {'bhp', 'qWs'}, {'Injector', 'Producer'})
%   
%   will give bottom hole pressures for the producer in D(:, 1, 2);
%             
%
% SEE ALSO:
%   plotWellSols

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
    assert(iscell(wellsols))
    nsteps = numel(wellsols);
    if nargin < 4
        index = 1;
    end
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
    datacounts = cellfun(@(x) numel(ws0(1).(x)), validfields);
    acceptible = datacounts >= index;
    validfields = validfields(acceptible);
    datacounts  = datacounts(acceptible);
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
    elseif islogical(wells)
        assert(numel(wells) == nw);
        subs = find(wells);
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
        fieldNo = strcmpi(validfields, fld);
        assert(any(strcmpi(fld, validfields)),...
            'Dimensions of %s does not match the given index %d', fld, index);
        if datacounts(fieldNo) == 0
            % No data, return NaN
            continue
        end
        data = cellfun(@(x) [arrayfun(@(y) y.(fld)(index), x(subs))], wellsols, 'UniformOutput', false);
        welldata(:, :, i) = vertcat(data{:});
    end
end
