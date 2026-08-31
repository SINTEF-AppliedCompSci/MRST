function suite = getIntegrationTestSuiteMRST()
%Build an integration test suite from all registered MRST modules.
%
% SYNOPSIS:
%   suite = getIntegrationTestSuiteMRST()
%
% DESCRIPTION:
%   Collects matlab.unittest.TestCase subclasses from:
%   Any ``tests/integrationTests/`` subdirectory in every registered
%   MRST module, provided the directory contains files whose names
%   match the convention  ``Test*.m``  or  ``*Test.m``.
%
% RETURNS:
%   suite - matlab.unittest.TestSuite aggregating all discovered tests.
%
% NOTE:
%   Discovery silently skips directories that cause errors during scan
%   (e.g. missing dependencies) to keep the runner robust.
%
% SEE ALSO:
%   runTests, getUnitTestSuiteMRST

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
    import matlab.unittest.TestSuite;

    % Collect sub-suites in a cell array to avoid vertcat issues with
    % matlab.unittest.TestSuite across MATLAB releases.
    parts = {};

    % Auto-discover from every registered module
    mods = mrstPath();
    testDirs = {'tests/integrationTests'};

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
                % matlab.unittest.TestCase subclass.
                if ~isTestCaseClass(fpath)
                    continue
                end
                try
                    s = matlab.unittest.TestSuite.fromFile(fpath);
                    parts{end+1} = s; %#ok<AGROW>
                catch ME
                    warning('mrst:integrationTestDiscovery', ...
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

