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
            disp(module)
            run(name);
            restore_env(m, g, v, d, p);
        end
    end
end

function names = getTestNames()
    names = getTestNamesInternal();
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
        'showOptionsAMGCL.m', ... % Does not work due to uiwait
        'SPE10SubsetADIExample.m', ... % Takes too long to run, ad-fi
            };
    modules = {'matlab_bgl'};
end

function [m, g, v, d, p] = clear_env
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
   mrstVerbose(v);

   gravity(g)
   if norm(g) > 0
      gravity on
   end

   mrstModule('reset', m{:});

   mrstDataDirectory(d);

   pause(p);
end
