function suite = getUnitTestSuiteMRST(varargin)
%Build a unit test suite from mrst-testing and all registered MRST modules.
%
% SYNOPSIS:
%   suite = getUnitTestSuiteMRST()
%   suite = getUnitTestSuiteMRST('modules', {'mod1', 'mod2'})
%
% DESCRIPTION:
%   Collects matlab.unittest.TestCase subclasses from:
%     1. The fixed folder in the mrst-testing module (tests/unitTests).
%        This is always included regardless of the 'modules' filter.
%     2. Any ``tests/unitTests/`` subdirectory in every registered
%        MRST module, provided the directory contains files whose names
%        match the convention  ``Test*.m``  or  ``*Test.m``.
%
% OPTIONAL PARAMETERS:
%   'modules' - Cell array of module names.  When provided, only tests from
%               those modules are included (step 2 above).  The fixed
%               mrst-testing folder (step 1) is always included.
%               Default: {} (include all registered modules).
%
% RETURNS:
%   suite - matlab.unittest.TestSuite aggregating all discovered tests.
%
% NOTE:
%   Discovery silently skips directories that cause errors during scan
%   (e.g. missing dependencies) to keep the runner robust.
%
% SEE ALSO:
%   runPreReleaseTests, getIntegrationTestSuiteMRST

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

    mrstModule add mrst-testing

    opt = struct('modules', {{}});
    opt = merge_options(opt, varargin{:});

    filterMods = opt.modules;
    if ~isempty(filterMods) && ischar(filterMods)
        filterMods = {filterMods};
    end

    % Collect sub-suites in a cell array to avoid vertcat issues with
    % matlab.unittest.TestSuite across MATLAB releases.
    parts = {};

    % --- 1. Fixed mrst-testing folders (always included) -----------------
    mtpath      = mrstPath('query', 'mrst-testing');
    unitfolders = {'tests/unitTests'};
    for i = 1:numel(unitfolders)
        p = fullfile(mtpath, unitfolders{i});
        if isfolder(p)
            try
                s = matlab.unittest.TestSuite.fromFolder(p);
                parts{end+1} = s; %#ok<AGROW>
            catch ME
                warning('mrst:unitTestDiscovery', ...
                    'Skipping folder %s: %s', p, ME.message);
            end
        end
    end

    % --- 2. Auto-discover from every registered module ------------------
    if isempty(filterMods)
        mods = mrstPath();
    else
        mods = filterMods;
        % Always test ad-unittest itself when a module filter is active
        if ~any(strcmpi(mods, 'ad-unittest'))
            mods = ['ad-unittest', mods];
        end
    end
    testDirs = {'tests/unitTests', 'UnitTests'};

    for im = 1:numel(mods)
        modPath = mrstPath('query', mods{im});
        if isempty(modPath), continue; end

        for id = 1:numel(testDirs)
            candidate = fullfile(modPath, testDirs{id});
            if ~isfolder(candidate), continue; end

            % Only include files that follow Test*.m or *Test.m convention
            files = [dir(fullfile(candidate, 'Test*.m')); ...
                     dir(fullfile(candidate, '*Test.m'))];
            if isempty(files), continue; end

            % Add directory to path temporarily so fromFile can resolve it
            wasOnPath = contains(path(), candidate);
            if ~wasOnPath
                addpath(candidate);
            end
            for fi = 1:numel(files)
                fpath = fullfile(candidate, files(fi).name);
                % Skip plain scripts — only load files that define a
                % matlab.unittest.TestCase subclass.  Plain scripts
                % (tutorial/demo files) would otherwise be picked up as
                % script-based tests whose %%-sections run in isolation,
                % causing failures due to missing cross-section state.
                if ~isTestCaseClass(fpath)
                    continue
                end
                try
                    s = matlab.unittest.TestSuite.fromFile(fpath);
                    parts{end+1} = s; %#ok<AGROW>
                catch ME
                    warning('mrst:unitTestDiscovery', ...
                        'Skipping %s: %s', fpath, ME.message);
                end
            end
            if ~wasOnPath
                rmpath(candidate);
            end
        end
    end

    % Combine all collected sub-suites
    if isempty(parts)
        suite = matlab.unittest.TestSuite.empty();
    else
        suite = horzcat(parts{:});
    end
end

% -------------------------------------------------------------------------

function tf = isTestCaseClass(fpath)
%Return true only if fpath defines a classdef that inherits from
%matlab.unittest.TestCase.  Plain scripts return false.
    tf = false;
    try
        txt = fileread(fpath);
    catch
        return
    end
    % Look for:  classdef <Name> < matlab.unittest.TestCase
    % (possibly with extra whitespace / multiple inheritance)
    tf = ~isempty(regexp(txt, ...
        '^\s*classdef\s+\w+\s*<[^%\n]*matlab\.unittest\.TestCase', ...
        'once', 'lineanchors'));
end
