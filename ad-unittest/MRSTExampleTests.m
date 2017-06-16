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
    mods = mrstPath();
    
    testNames = cell(numel(mods), 1);
    modNames = cell(numel(mods), 1);
    for i = 1:numel(mods)
        ex = mrstExamples(mods{i});
        examples = ex{1};
        for j = 1:numel(examples)
            [~, examples{j}] = fileparts(examples{j});
        end
        testNames{i} = examples;
        tmp = cell(size(examples));
        [tmp{:}] = deal(mods{i});
        modNames{i} = tmp;
    end
    
    names = horzcat(testNames{:});
    modules = horzcat(modNames{:});
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
