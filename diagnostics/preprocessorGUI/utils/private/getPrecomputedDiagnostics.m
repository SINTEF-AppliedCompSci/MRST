function precomp = getPrecomputedDiagnostics(casenm, steps, pdir)
% Function accepts either casenm for ECLIPSE input or problem structure for
% MRST input.

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

% if not given, assume precompted diagnositics lies in caseDir/mrst_diagnostics
if ischar(casenm) % ECLIPSE Input
    [caseDir, prefix] = fileparts(casenm);
elseif isstruct(casenm) % MRST Input
    caseDir = casenm.OutputHandlers.states.getDataPath;
    prefix = casenm.BaseName;
end

if nargin < 3
    pdir = fullfile(caseDir, 'mrst_diagnostics');
end

if ~isempty(dir(pdir))
    fn = @(n)fullfile(pdir, [prefix, sprintf('_diagn%0.4d.mat', n)]);
    fail = false;
    fprintf('Loading precomputed diagnostics %3.0d%%', 0)
    precomp = cell(numel(steps), 1);
    for k = 1:numel(steps)
        try
            precomp{k} = load(fn(steps(k)));
        catch
            fail = true;
            break
        end
        fprintf('\b\b\b\b%3.0d%%', round(100*k/numel(steps)))
    end
    fprintf(', done\n')
    if fail
        warning('Not able to load precomputed diagnostics')
        precomp = [];
    end
else
    precomp = [];
end
end
