classdef MRSTExampleTests < matlab.unittest.TestCase
    properties (TestParameter)
        name = getTestNames();
        module = getTestModules();
    end
    methods
        function test = MRSTExampleTest(varargin)
            mrstModule add release
            require release
            test.examplefile = examplefile;
        end
    end
    methods (Test, ParameterCombination='sequential')
        function runExample(test, name, module)
            disp(name)
            [m, g, v, d, p] = clear_env();
            mrstModule('add', module);
            runScoped(name);
            restore_env(m, g, v, d, p);
        end
    end
end

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

function names = getTestNames()
    names = getTestNamesInternal();
end

function runScoped(name)
    run(name);
end

function mods = getTestModules()
    [~, mods] = getTestNamesInternal();
end

function [names, modules] = getTestNamesInternal()
    [skip, skip_mod] = getSkippedTests();
    mods = mrstPath();
    mods = setdiff(mods, skip_mod);
    
    testNames = cell(numel(mods), 1);
    modNames = cell(numel(mods), 1);
    for i = 1:numel(mods)
        ex = mrstExamples(mods{i});
        examples = ex{1};
        keep = true(numel(examples), 1);
        for j = 1:numel(examples)
            test_parts = strsplit(examples{j}, filesep);
            testname = test_parts{end};
            % Filter specifically skipped tests
            toSkip = any(strcmpi(skip, testname));
            % Filter experimental folders
            isExperimental = any(strcmpi(test_parts, 'experimental'));
            keep(j) = not(toSkip || isExperimental);
            [~, examples{j}] = fileparts(examples{j});
        end
        examples = examples(keep);
        testNames{i} = examples;
        tmp = cell(size(examples));
        [tmp{:}] = deal(mods{i});
        modNames{i} = tmp;
    end
    
    names = horzcat(testNames{:});
    modules = horzcat(modNames{:});
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
