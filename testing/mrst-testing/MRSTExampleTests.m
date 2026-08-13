classdef MRSTExampleTests < matlab.unittest.TestCase
    % MRSTExampleTests  Parameterised smoke-test runner for all MRST examples.
    %
    % Each registered MRST example is run as an individual test case. The
    % set of examples executed can be controlled in three ways:
    %
    %   1. Global skip list in getSkippedTests (backward compatible).
    %   2. Module-level ``getSkippedExamples`` function (optional per module).
    %   3. Per-file ``MRST_TEST_OPTIONS`` annotation block, e.g.::
    %
    %       %{
    %       MRST_TEST_OPTIONS
    %       interactive: true    % skip in automated runs
    %       tags:        slow    % excluded by default in CI
    %       timeout:     120     % seconds
    %       %}
    %
    % See also: mrstExampleOptions, getExampleIntegrationTestSuiteMRST

    properties (TestParameter)
        name   = getTestNames();
        module = getTestModules();
    end

    methods (Test, ParameterCombination='sequential')
        function runExample(test, name, module)
            disp(name)
            % Read per-file annotation options
            opt = mrstExampleOptions(name);
            % Skip interactive examples unconditionally in automated runs
            test.assumeFalse(opt.interactive, ...
                sprintf(['Example ''%s'' is marked interactive and ', ...
                         'must be tested manually.'], name));
            % Skip examples tagged as data-required when the env flag is set
            if getenv('MRST_SKIP_DATA_REQUIRED')
                test.assumeFalse(any(strcmpi(opt.tags, 'data-required')), ...
                    sprintf(['Example ''%s'' requires external data ', ...
                             '(tag: data-required).'], name));
            end
            [m, g, v, d, p] = clear_env();
            mrstModule('add', module);
            runScoped(name);
            restore_env(m, g, v, d, p);
        end
    end
end

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

function names = getTestNames()
    names = getTestNamesInternal();
end

function runScoped(name)
    run(name);
end

function mods = getTestModules()
    [~, mods] = getTestNamesInternal();
end

function [names, modules] = getTestNamesInternal(varargin)
    % Optional: filter by tags (cell array of tags to exclude)
    opt = struct('excludeTags', {{}});
    opt = merge_options(opt, varargin{:});

    [skip, skip_mod] = getSkippedTests();
    mods = mrstPath();
    mods = setdiff(mods, skip_mod);
    
    testNames = cell(numel(mods), 1);
    modNames = cell(numel(mods), 1);
    for i = 1:numel(mods)
        ex = mrstExamples(mods{i});
        examples = ex{1};
        % Collect module-level skip list (optional per-module function)
        modSkip = getModuleSkipList(mods{i});
        keep = true(numel(examples), 1);
        for j = 1:numel(examples)
            test_parts = strsplit(examples{j}, filesep);
            testname   = test_parts{end};
            [~, basename] = fileparts(testname);
            % Filter globally skipped tests
            toSkip = any(strcmpi(skip, testname));
            % Filter module-level skipped tests
            toSkip = toSkip || any(strcmpi(modSkip, basename)) || ...
                               any(strcmpi(modSkip, testname));
            % Filter experimental folders
            isExperimental = any(strcmpi(test_parts, 'experimental'));
            % Filter by MRST_TEST_OPTIONS annotations (interactive + tags)
            annotOpt = mrstExampleOptions(examples{j});
            isInteractive = annotOpt.interactive;
            hasExcludedTag = false;
            if ~isempty(opt.excludeTags)
                hasExcludedTag = any(ismember(lower(annotOpt.tags), ...
                                              lower(opt.excludeTags)));
            end
            keep(j) = ~(toSkip || isExperimental || isInteractive || hasExcludedTag);
            [~, examples{j}] = fileparts(examples{j});
        end
        examples = examples(keep);
        testNames{i} = examples;
        tmp = cell(size(examples));
        [tmp{:}] = deal(mods{i});
        modNames{i} = tmp;
    end
    
    names   = horzcat(testNames{:});
    modules = horzcat(modNames{:});
end

function skipList = getModuleSkipList(modname)
    % Returns per-module list of example basenames to skip.
    % Looks for a function getSkippedExamples in the module root directory.
    skipList = {};
    modPath = mrstPath('query', modname);
    if isempty(modPath), return; end
    fnpath = fullfile(modPath, 'getSkippedExamples.m');
    if exist(fnpath, 'file')
        try
            % Temporarily add module to ensure the function is callable
            oldPath = path();
            addpath(modPath);
            skipList = getSkippedExamples();
            path(oldPath);
        catch
            % Ignore errors in module-provided functions
        end
    end
end

function [names, modules] = getSkippedTests()
    names = {
        'showOptionsAMGCL', ... % Does not work due to uiwait
        'SPE10SubsetADIExample', ... % Takes too long to run, ad-fi
        'runNorneExample', ... % Takes too long to run
        'diagnosticsPostProcessorWithMRST', ... % GUI example
        'preprocessDiagnosticsEgg', ... % GUI example
        'ensembleGUIForEgg', ... % GUI example
        'trajectoryExampleEgg', ... % GUI example
        'ensemblePackedProblemsExample', ... % Launches Matlab sessions
        '', ...
        'demoPackedProblems'... % Example which launches Matlab sessions
            };
    names = cellfun(@(x) [x, '.m'], names, 'UniformOutput', false);
    modules = {'matlab_bgl', 'octave', ...
               'stokes-brinkman', 'mrst-experimental',...
               'impes', 'ad-fi'};
end

function [m, g, v, d, p] = clear_env
   close all;
   m = mrstModule;
   g = gravity;
   v = mrstVerbose;
   d = mrstDataDirectory;
   p = pause('off');

   mrstModule  clear
   gravity     reset
   mrstVerbose true
   clear       functions                                       %#ok<CLFUNC>

   mrstModule('add', m{:});

   mrstDataDirectory(d);
end

%--------------------------------------------------------------------------

function restore_env(m, g, v, d, p)
   close all;
   mrstVerbose(v);

   gravity(g)
   if norm(g) > 0
      gravity on
   end
   mrstModule('reset', m{:});
   mrstDataDirectory(d);
   pause(p);
end
