function output = readEclipseRestartUnFmt(prefix, spec, steps)
%Read unformatted (binary) ECLIPSE unified or multiple restart file(s)
%
% SYNOPSIS:
%   output = readEclipseRestartUnFmt(prefix, spec, steps)
%
% PARAMETERS:
%   prefix - Name (string) of file(s) (omitting extensions) containing 
%            unformatted ECLIPSE results. 
%
%   spec   - Specification (struct) of restart data as obtained by
%            processEclipseRestartSpec. If empty or not given, 
%            processEclipseRestartSpec will be called with default options.
%
%   steps  - List (indices) of steps to read. If empty or not given, all 
%            steps are read.            
%
% RETURNS:
%   output - Data structure containing restart field data for requested 
%            fields/steps. 
%
% SEE ALSO:
%   `readEclipseOutputFileFmt`, `readEclipseRestartSpec`.

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

if nargin < 2 || isempty(spec)
    try
        spec = processEclipseRestartSpec(prefix);
    catch
        warning('No restart spesification found, attempting to use potentially inefficient fallback routine')
        output = readEclipseRestartUnFmt_fallback(prefix);
        fn = fieldnames(output);
        if nargin < 3 || isempty(steps)
            steps = ':';
        end
        for k = 1:numel(fn)
            if isfield(output.(fn{k}), 'values')
                output.(fn{k}) = output.(fn{k}).values(steps);
            else
                output.(fn{k}) = output.(fn{k})(steps);
            end
        end
        return;
    end
end

if nargin < 3 || isempty(steps)
    steps = 1:numel(spec.time);
end

% Check compatibility of steps
isOK = and(steps >= 1, steps <= numel(spec.time));
if any(~isOK)
    warning('Ignoring out-of-range entries of input-vector ''steps''.')
    steps = steps(isOK);
    if isempty(steps)
        output = [];
        return;
    end
end
    
% Preallocate
[maxnum, ix] = max(cellfun(@numel, spec.keywords(steps)));   
nms = cellfun(@fixVarName, spec.keywords{steps(ix)}, 'UniformOutput', false);  
v   = repmat({cell(1, numel(steps))}, maxnum, 1);
output = cell2struct(v, nms);

% Read
if strcmp(spec.type, 'unified')
    fname = [prefix, '.UNRST'];
    [fid, msg] = fopen(fname, 'rb', 'ieee-be');
    if fid < 0
       error('Open:Failure', ...
             'Failed to open Restart File ''%s'': %s', fname, msg);
    end
end
dispif(mrstVerbose, 'Reading restart:     ')

for ks = 1:numel(steps)
    dispif(mrstVerbose, '\b\b\b\b%3d%%', round(100*ks/numel(steps)));
    step = steps(ks);
    if strcmp(spec.type, 'multiple')
        [fid, msg] = fopen(spec.fnames{step}, 'rb',  'ieee-be');
        if fid < 0,
           error('Open:Failure', ...
                 'Failed to open Restart File ''%s'': %s', ...
                 spec.fnames{step}, msg);
        end
    end
    nms  = cellfun(@fixVarName, spec.keywords{step}, 'UniformOutput', false);
    p    = spec.pointers{step};
    n    = spec.num{step};
    prec = spec.prec{step};
    for kf = 1:numel(nms)
        curp = ftell(fid);
        fseek(fid, p(kf)+28-curp, 0);
        if n(kf)>0
        if ~strcmp(prec{kf}, '840*uchar=>char')
            output.(nms{kf}){ks} = fread(fid, n(kf), prec{kf}, 8);
        else
            v = fread(fid, 8*n(kf), prec{kf}, 8);
            output.(nms{kf}){ks} = cellstr(reshape(v , 8, n(kf))');
        end
        end
    end
    if strcmp(spec.type, 'multiple')
        fid = fclose(fid);
    end
end

if strcmp(spec.type, 'unified')
    fclose(fid);
end
dispif(mrstVerbose, ',  done\n')
end
        
%--------------------------------------------------------------------------

function name = fixVarName(name)
    if ~isvarname(name)
        name = regexprep(name, {'+', '-'}, {'p', 'n'});
        name = genvarname(regexprep(name, '\W', '_'));
    end
end
