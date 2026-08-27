function rts = getAllRegressionTests()
%Collect regression tests from every registered MRST module.
%
% DESCRIPTION:
%   Calls ``getModuleRegressionTests()`` in every registered MRST module
%   that provides such a function.  Returns all collected
%   MRSTRegressionTest objects in a single cell array.
%
%   Module authors should expose regression tests by placing a function
%   ``getModuleRegressionTests.m`` at the root of their module.  Legacy
%   ``reg-tests/getModuleRegressionTests.m`` layouts are also supported.
%   The function must return a cell array of MRSTRegressionTest objects::
%
%       function rts = getModuleRegressionTests()
%           rts = { MRSTRegressionTest(@myCase), ... };
%       end
%
% RETURNS:
%   rts - Cell array of MRSTRegressionTest objects.
%
% SEE ALSO:
%   MRSTRegressionTest, runPreReleaseTests

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

    rts  = {};
    mods = mrstPath();

    for i = 1:numel(mods)
        modPath = mrstPath('query', mods{i});
        if isempty(modPath), continue; end

        fnpaths = {
            fullfile(modPath, 'getModuleRegressionTests.m')
            fullfile(modPath, 'reg-tests', 'getModuleRegressionTests.m')
        };
        idx = find(cellfun(@(f) exist(f, 'file') > 0, fnpaths), 1);
        if isempty(idx), continue; end
        fnpath = fnpaths{idx};

        wasOnPath = contains(path(), modPath);
        if ~wasOnPath
            addpath(modPath);
        end
        fnDir = fileparts(fnpath);
        addedFnDir = false;
        if ~strcmp(fnDir, modPath) && ~contains(path(), fnDir)
            addpath(fnDir);
            addedFnDir = true;
        end
        try
            modRTs = getModuleRegressionTests();
            if isempty(modRTs)
                modRTs = {};
            elseif iscell(modRTs)
                % keep as-is
            elseif numel(modRTs) > 1
                modRTs = num2cell(modRTs);
            else
                modRTs = {modRTs};
            end
            rts = [rts, modRTs(:)']; %#ok<AGROW>
        catch ME
            warning('mrst:regressionTestDiscovery', ...
                'Could not load regression tests from module ''%s'': %s', ...
                mods{i}, ME.message);
        end
        if addedFnDir
            rmpath(fnDir);
        end
        if ~wasOnPath
            rmpath(modPath);
        end
    end
end
